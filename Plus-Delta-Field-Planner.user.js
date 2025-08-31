// ==UserScript==
// @id             iitc-plugin-plus-delta-field-planner
// @name           IITC plugin: Plus Delta Field Planner
// @category       Misc
// @version        1.0.0
// @author         gpt
// @match          https://intel.ingress.com/*
// @grant          none
// ==/UserScript==

(() => {
    'use strict';
    if (typeof window.plugin === 'undefined') window.plugin = function () {};
    // do NOT overwrite if it's already defined (function or object)

    window.plugin.deltaFieldPlanner = window.plugin.deltaFieldPlanner || {};
    const DFP = window.plugin.deltaFieldPlanner;

    DFP.openPanel = function openPanel() {
        if (DFP.dlg && DFP.dlg.dialog) {
            try { DFP.dlg.dialog('open'); } catch (_) {}
            return;
        }

        const $status = $(`
      <div id="dfp-status">
        <div><b>Triangle:</b> not scanned</div>
        <div><b>Outer score:</b> <span id="dfp-outer">-</span></div>
        <div><b>Elapsed:</b> <span id="dfp-elapsed">0s</span></div>
        <div><b>Best score:</b> <span id="dfp-best">-</span></div>
        <div><b>Faces:</b> <span id="dfp-faces">-</span></div>
      </div>
    `);


        const $controls = $('<div>').css({ marginTop: '8px', display: 'flex', gap: '6px', alignItems: 'center' });
        const $limitLbl = $('<label for="dfp-limit">Time limit (s):</label>').css({ fontSize: '11px' });
        const $limitInp = $('<input id="dfp-limit" type="number" min="1" step="1" value="120" style="width:72px;">');
        const $scanBtn  = $('<button type="button">Scan</button>').click(() => {
            try {
                const tri = DFP.readOuterTriangle(); // [{lat,lng,ref}, x3]
                const inside = DFP.getInsidePortals(tri, DFP.buildPortalList());
                // render boundary only (clear previous shapes)
                DFP.renderOuterBoundary(tri, true);

                // update panel
                $('#dfp-status').find('div').first()
                    .html(`<b>Triangle:</b> OK • Inside portals: ${inside.length}`);
                $('#dfp-outer').text(DFP.triScoreRaw(tri[0], tri[1], tri[2]).toFixed(1));
            } catch (e) { alert(e && e.message ? e.message : String(e)); }
        });
        const $runBtn   = $('<button type="button">Optimize</button>').click(() => {
            const sec = parseInt($('#dfp-limit').val(), 10) || 120;
            try {
                DFP.solveAndRender(sec);
            } catch (e) {
                alert(e && e.message ? e.message : String(e));
            }
        });

        const $stopBtn  = $('<button type="button">Stop</button>').click(() => {
            if (DFP.stopOptimization) DFP.stopOptimization();
        });

        $controls.append($limitLbl, $limitInp, $scanBtn, $runBtn, $stopBtn);

        const $body = $('<div id="dfp-panel-body">').append($status, $controls);

        DFP.dlg = dialog({
            title: 'Δ Field Planner',
            html: $body,
            id: 'dfp-panel',
            width: 360,
            closeCallback: function () { DFP.dlg = null; }
        });
    };


    // Add a simple toolbox button
    DFP.injectButton = function injectButton() {
        if (document.getElementById('dfp-open-btn')) return;
        const $btn = $('<a id="dfp-open-btn" title="Delta Field Planner">Δ-Plan</a>').click(DFP.openPanel);
        const $toolbox = $('#toolbox');
        if ($toolbox.length) $toolbox.append($btn);
    };

    // ---- utils for reading the outer triangle & counting inner portals ----

    // Tolerance when matching a drawn marker to a real portal (same as Microfield-Checker).
    DFP.TOLERANCE = 1e-5;

    DFP.key = function key(p) {
        return `${p.lat.toFixed(6)},${p.lng.toFixed(6)}`;
    };

    // Ray-casting point-in-polygon test (Microfield-Checker's pnpoly).
    DFP.pnpoly = function pnpoly(polygon, point) {
        let c = false;
        const n = polygon.length;
        for (let i = 0, j = n - 1; i < n; j = i++) {
            if (
                ((polygon[i].lat > point.lat) !== (polygon[j].lat > point.lat)) &&
                (point.lng <
                 (polygon[j].lng - polygon[i].lng) * (point.lat - polygon[i].lat) /
                 (polygon[j].lat - polygon[i].lat) +
                 polygon[i].lng)
            ) {
                c = !c;
            }
        }
        return c;
    };

    // Collect 3 drawn markers from draw-tools.
    DFP.collectDrawToolMarkers = function collectDrawToolMarkers() {
        if (!window.plugin.drawTools || !window.plugin.drawTools.drawnItems) {
            throw new Error('draw-tools is not available.');
        }
        const markers = [];
        window.plugin.drawTools.drawnItems.eachLayer(layer => {
            if (layer instanceof L.Marker) {
                const { lat, lng } = layer.getLatLng();
                markers.push({ lat, lng });
            }
        });
        return markers;
    };

    // Build a flat list of currently loaded portals (lat/lng + ref).
    DFP.buildPortalList = function buildPortalList() {
        return Object.values(window.portals).map(p => ({
            lat: p.getLatLng().lat,
            lng: p.getLatLng().lng,
            ref: p
        }));
    };

    // Find the real portal that matches a drawn marker (by tolerance).
    DFP.matchMarkerToPortal = function matchMarkerToPortal(marker, portalList) {
        const t = DFP.TOLERANCE;
        return portalList.find(p =>
                               Math.abs(p.lat - marker.lat) < t && Math.abs(p.lng - marker.lng) < t
                              );
    };

    // Read 3 anchors (A,B,C) from draw-tools and map them to real portals.
    DFP.readOuterTriangle = function readOuterTriangle() {
        const markers = DFP.collectDrawToolMarkers();
        if (markers.length !== 3) {
            throw new Error(`Expected exactly 3 markers, found ${markers.length}.`);
        }
        const plist = DFP.buildPortalList();
        const triangle = markers.map(m => {
            const match = DFP.matchMarkerToPortal(m, plist);
            if (!match) throw new Error(`No portal found at ${m.lat}, ${m.lng}`);
            return match;
        });
        return triangle; // [A,B,C], each {lat,lng,ref}
    };

    // Return portals strictly inside the triangle (exclude the 3 anchors).
    DFP.getInsidePortals = function getInsidePortals(triangle, portalList) {
        const polygon = triangle.map(p => ({ lat: p.lat, lng: p.lng }));
        const isAnchor = (p) =>
        triangle.some(v => Math.abs(v.lat - p.lat) < DFP.TOLERANCE &&
                      Math.abs(v.lng - p.lng) < DFP.TOLERANCE);
        return portalList.filter(p => !isAnchor(p) && DFP.pnpoly(polygon, p));
    };

    // Orchestrate: read outer triangle & count inside portals, then update panel.
    DFP.scanTriangleAndCount = function scanTriangleAndCount() {
        try {
            const triangle = DFP.readOuterTriangle();
            const portals  = DFP.buildPortalList();
            const inside   = DFP.getInsidePortals(triangle, portals);

            // Update panel
            const $status = $('#dfp-status');
            if ($status.length) {
                $status.html(
                    [
                        '<b>Triangle read:</b> OK',
                        'Anchors: 3',
                        `Inside portals: ${inside.length}`,
                        `Total (anchors + inside): ${inside.length + 3}`
          ].join('<br/>')
                );
            }
        } catch (err) {
            alert(err && err.message ? err.message : String(err));
        }
    };

    // ---- geometry & optimization for best triangulation (time-bounded) ----
    // ---- outer triangle helpers ----

    // score for raw points (lat/lng objects)
    DFP.triScoreRaw = function triScoreRaw(a, b, c) {
        const area = DFP.triArea(a, b, c);
        const d2 = (p, q) => {
            const dx = p.lat - q.lat, dy = p.lng - q.lng;
            return dx*dx + dy*dy;
        };
        const L2 = Math.max(d2(a,b), d2(b,c), d2(c,a));
        if (L2 <= 0 || area <= 0) return 100; // degenerate -> base score
        const eq = Math.min(1, area / ((Math.sqrt(3)/4) * L2));
        return Math.round(100 * (1 + eq));
    };

    // render ONLY the outer triangle boundary (clear layer if requested)
    DFP.renderOuterBoundary = function renderOuterBoundary(triangle, clearFirst = true) {
        if (!DFP.layerGroup) return;
        if (clearFirst) DFP.layerGroup.clearLayers();
        const A = triangle[0], B = triangle[1], C = triangle[2];

        const add = (P, Q) => {
            // distinct style to differentiate from micro edges
            const line = L.polyline([[P.lat, P.lng], [Q.lat, Q.lng]], {
                color: DFP.NEUTRAL_LINK_COLOR,
                weight: 3,
                opacity: 1.0,
                dashArray: '8,6'
            });
            DFP.layerGroup.addLayer(line);
        };
        add(A, B); add(B, C); add(C, A);
    };


    // orientation (signed double area); >0: ccw, <0: cw, =0: collinear
    DFP.orient2d = function orient2d(a, b, c) {
        return (b.lng - a.lng) * (c.lat - a.lat) - (b.lat - a.lat) * (c.lng - a.lng);
    };

    // area of triangle (positive)
    DFP.triArea = function triArea(a, b, c) {
        return Math.abs(DFP.orient2d(a, b, c)) * 0.5;
    };

    // ensure triangle indices are CCW
    DFP.ccwTri = function ccwTri(tri, pts) {
        const [i, j, k] = tri;
        return DFP.orient2d(pts[i], pts[j], pts[k]) > 0 ? tri : [i, k, j];
    };

    // max edge length squared for triangle (avoid sqrt where possible)
    DFP.maxEdgeLen2 = function maxEdgeLen2(a, b, c) {
        const d2 = (p, q) => {
            const dx = p.lat - q.lat, dy = p.lng - q.lng;
            return dx * dx + dy * dy;
        };
        const ab = d2(a, b), bc = d2(b, c), ca = d2(c, a);
        return Math.max(ab, bc, ca);
    };

    // equilaterality in [0,1]; A / ( (sqrt(3)/4) * L^2 )
    DFP.equil = function equil(a, b, c) {
        const A = DFP.triArea(a, b, c);
        const L2 = DFP.maxEdgeLen2(a, b, c);
        if (L2 <= 0 || A <= 0) return 0;
        const k = Math.sqrt(3) / 4; // ≈0.4330127019
        return Math.min(1, A / (k * L2));
    };

    // score of a single triangle for single-agent case: 100*(1+equil)
    DFP.triScore = function triScore(a, b, c) {
        return 100 * (1 + DFP.equil(a, b, c));
    };

    // compute total score and return {score, faces}
    DFP.scoreTriangulation = function scoreTriangulation(tris, pts) {
        let s = 0;
        for (let t of tris) {
            const [i, j, k] = t;
            s += DFP.triScore(pts[i], pts[j], pts[k]);
        }
        return { score: s, faces: tris.length };
    };

    // simple point-in-triangle (including boundary as inside)
    DFP.pointInTri = function pointInTri(p, a, b, c) {
        const o1 = DFP.orient2d(a, b, p);
        const o2 = DFP.orient2d(b, c, p);
        const o3 = DFP.orient2d(c, a, p);
        const hasNeg = (o1 < 0) || (o2 < 0) || (o3 < 0);
        const hasPos = (o1 > 0) || (o2 > 0) || (o3 > 0);
        return !(hasNeg && hasPos);
    };

    // ---- Microfielding DP solver (time-bounded) ----
    // Points: anchors first (indices 0,1,2), then interior portals.
    DFP.points = [];
    DFP.memo = new Map();       // key "i,j,k" (CCW) -> {score:number, choice:number|null}
    DFP.insideCache = new Map();// key "i,j,k" (CCW) -> Array<number> of interior point indices
    DFP.scoreCache = new Map(); // key "i,j,k" (sorted or CCW) -> number
    DFP.stopFlag = false;
    DFP.limitMs = 120000;
    DFP.startMs = 0;

    DFP.EPS = 1e-12;

    // neutral link color (Ingress-neutral gray)
    DFP.NEUTRAL_LINK_COLOR = '#9e9e9e';

    // orientation (signed double area); >0: ccw, <0: cw, 0: collinear
    DFP.orient2d = function orient2d(a, b, c) {
        return (b.lng - a.lng) * (c.lat - a.lat) - (b.lat - a.lat) * (c.lng - a.lng);
    };

    // enforce CCW order of a triple of indices
    DFP.ccwTriIdx = function ccwTriIdx(i, j, k) {
        const A = DFP.points[i], B = DFP.points[j], C = DFP.points[k];
        return DFP.orient2d(A, B, C) > 0 ? [i, j, k] : [i, k, j];
    };

    // canonical key for memo/cache
    DFP.keyOf = function keyOf(i, j, k) {
        const t = DFP.ccwTriIdx(i, j, k);
        return t[0] + ',' + t[1] + ',' + t[2];
    };

    // area and side helpers
    DFP.triArea = function triArea(a, b, c) {
        return Math.abs(DFP.orient2d(a, b, c)) * 0.5;
    };
    DFP.maxEdgeLen2 = function maxEdgeLen2(a, b, c) {
        const d2 = (p, q) => {
            const dx = p.lat - q.lat, dy = p.lng - q.lng;
            return dx*dx + dy*dy;
        };
        return Math.max(d2(a,b), d2(b,c), d2(c,a));
    };
    DFP.equil = function equil(a, b, c) {
        const A = DFP.triArea(a, b, c);
        const L2 = DFP.maxEdgeLen2(a, b, c);
        if (L2 <= 0 || A <= 0) return 0;
        return Math.min(1, A / ((Math.sqrt(3)/4) * L2));
    };
    DFP.triScore = function triScore(i, j, k) {
        // cache by CCW key (orientation doesn't change score because we use |area|)
        const key = DFP.keyOf(i, j, k);
        if (DFP.scoreCache.has(key)) return DFP.scoreCache.get(key);
        const a = DFP.points[i], b = DFP.points[j], c = DFP.points[k];
        const s = Math.round(100 * (1 + DFP.equil(a, b, c)));
        DFP.scoreCache.set(key, s);
        return s;
    };

    // strict point-in-triangle (exclude boundary within EPS)
    DFP.pointInTriStrict = function pointInTriStrict(p, i, j, k) {
        const A = DFP.points[i], B = DFP.points[j], C = DFP.points[k];
        const o1 = DFP.orient2d(A, B, p);
        const o2 = DFP.orient2d(B, C, p);
        const o3 = DFP.orient2d(C, A, p);
        const eps = DFP.EPS;
        const pos = (o1 > eps && o2 > eps && o3 > eps);
        const neg = (o1 < -eps && o2 < -eps && o3 < -eps);
        return pos || neg;
    };

    // lazily compute interior points for a CCW triangle key "i,j,k"
    DFP.getInsideForKey = function getInsideForKey(key) {
        const cached = DFP.insideCache.get(key);
        if (cached) return cached;
        const [i, j, k] = key.split(',').map(Number);
        const res = [];
        for (let idx = 3; idx < DFP.points.length; idx++) {
            if (idx === i || idx === j || idx === k) continue;
            if (DFP.pointInTriStrict(DFP.points[idx], i, j, k)) res.push(idx);
        }
        DFP.insideCache.set(key, res);
        return res;
    };

    // time check
    DFP.timeExceeded = function timeExceeded() {
        return DFP.stopFlag || (Date.now() - DFP.startMs) > DFP.limitMs;
    };

    // Best(i,j,k): return {score, choice} with memoization; microfielding recursion
    DFP.Best = function Best(i, j, k) {
        const key = DFP.keyOf(i, j, k);
        const hit = DFP.memo.get(key);
        if (hit) return hit;

        const cand = DFP.getInsideForKey(key);
        if (cand.length === 0) {
            const base = { score: 0, choice: null };
            DFP.memo.set(key, base);
            return base;
        }

        // if time is up, return a greedy feasible choice (only immediate gain)
        if (DFP.timeExceeded()) {
            let bestP = null, bestGain = -Infinity;
            for (const p of cand) {
                const gain = DFP.triScore(i, j, p) + DFP.triScore(j, k, p) + DFP.triScore(k, i, p);
                if (gain > bestGain) { bestGain = gain; bestP = p; }
            }
            const g = { score: bestGain, choice: bestP };
            DFP.memo.set(key, g);
            return g;
        }

        let best = -Infinity, bestP = null;
        for (const p of cand) {
            // immediate gain: three new fields at this split
            const gain = DFP.triScore(i, j, p) + DFP.triScore(j, k, p) + DFP.triScore(k, i, p);
            // recursive gains
            const s1 = DFP.Best(i, j, p).score;
            const s2 = DFP.Best(j, k, p).score;
            const s3 = DFP.Best(k, i, p).score;
            const total = gain + s1 + s2 + s3;
            if (total > best) { best = total; bestP = p; }

            // optional: opportunistic time check between iterations
            if (DFP.timeExceeded()) break;
        }

        const out = { score: best, choice: bestP };
        DFP.memo.set(key, out);
        return out;
    };

    // reconstruct edges by following memoized choices; returns {edges:Set<string>, faces:number}
    DFP.reconstructEdges = function reconstructEdges(i, j, k) {
        const edges = new Set(); // "u-v" with u<v
        let usedSplitNodes = 0;

        const addEdge = (u, v) => {
            const a = Math.min(u, v), b = Math.max(u, v);
            edges.add(a + '-' + b);
        };

        const dfs = (a, b, c) => {
            const key = DFP.keyOf(a, b, c);
            const rec = DFP.memo.get(key);
            if (!rec || rec.choice == null) return;
            const p = rec.choice;
            addEdge(a, p); addEdge(b, p); addEdge(c, p);
            usedSplitNodes += 1;
            dfs(a, b, p);
            dfs(b, c, p);
            dfs(c, a, p);
        };

        dfs(i, j, k);
        return { edges, faces: usedSplitNodes * 3 };
    };

    // solve & render best microfielding within time limit
    DFP.solveAndRender = function solveAndRender(secondsLimit) {
        // 1) collect anchors and interior portals (as before)
        const triangle = DFP.readOuterTriangle();       // 3 anchors: [{lat,lng,ref},...]
        const portals = DFP.buildPortalList();
        const inside = DFP.getInsidePortals(triangle, portals);

        // build points array: anchors first, then inside
        DFP.points = triangle.concat(inside);

        // reset caches
        DFP.memo.clear();
        DFP.insideCache.clear();
        DFP.scoreCache.clear();
        DFP.stopFlag = false;
        DFP.limitMs = Math.max(1, Math.floor((secondsLimit || 120) * 1000));
        DFP.startMs = Date.now();

        // elapsed UI ticker (note: heavy runs may only update at the end due to JS event loop)
        const $elapsed = $('#dfp-elapsed');
        const tick = setInterval(() => {
            if ($elapsed.length) {
                const sec = Math.floor((Date.now() - DFP.startMs) / 1000);
                $elapsed.text(`${sec}s`);
            }
        }, 250);

        // build points array: anchors first, then inside
        DFP.points = triangle.concat(inside);
        const outerA = DFP.points[0], outerB = DFP.points[1], outerC = DFP.points[2];

        // 2) solve from the outer triangle (0,1,2)
        const res = DFP.Best(0, 1, 2);

        // 3) draw: outer boundary first, then unique micro edges
        if (DFP.layerGroup) DFP.layerGroup.clearLayers();
        DFP.renderOuterBoundary([outerA, outerB, outerC], false);

        const { edges, faces } = DFP.reconstructEdges(0, 1, 2);
        edges.forEach(key => {
            const [u, v] = key.split('-').map(Number);
            const P = DFP.points[u], Q = DFP.points[v];
            const pl = L.polyline([[P.lat, P.lng], [Q.lat, Q.lng]], {
                color: DFP.NEUTRAL_LINK_COLOR,
                weight: 2,
                opacity: 0.9
                });
            DFP.layerGroup.addLayer(pl);
        });

        // 4) update UI
        clearInterval(tick);
        $('#dfp-outer').text(DFP.triScoreRaw(outerA, outerB, outerC));
        $('#dfp-best').text(res.score.toFixed(1));
        $('#dfp-faces').text(faces);
    };

    DFP.setup = function setup() {
        DFP.layerGroup = new L.LayerGroup();
        window.addLayerGroup('Δ Field Planner', DFP.layerGroup, true);
        DFP.injectButton();
    };

    // Hook into IITC boot
    window.bootPlugins = window.bootPlugins || [];
    window.bootPlugins.push(DFP.setup);
    if (window.iitcLoaded) DFP.setup();
})();
