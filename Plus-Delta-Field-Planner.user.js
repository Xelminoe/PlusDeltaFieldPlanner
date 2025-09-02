// ==UserScript==
// @id             iitc-plugin-plus-delta-field-planner
// @name           IITC plugin: Plus Delta Field Planner
// @category       Misc
// @version        1.0.2
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
        <div> Draw polygons to select portals, place markers to exclude them. Requires Draw-Tools. </div>
        <div><b>Elapsed:</b> <span id="dfp-elapsed">0s</span></div>
        <div><b>Best score:</b> <span id="dfp-best">-</span></div>
        <div><b>Faces:</b> <span id="dfp-faces">-</span></div>
      </div>
    `);


        const $controls = $('<div>').css({ marginTop: '8px', display: 'flex', gap: '6px', alignItems: 'center' });
        const $controls2 = $('<div>').css({ marginTop: '8px', display: 'flex', gap: '6px', alignItems: 'center' });
        const $limitLbl = $('<label for="dfp-limit">Time limit (s):</label>').css({ fontSize: '11px' });
        const $limitInp = $('<input id="dfp-limit" type="number" min="1" step="1" value="120" style="width:72px;">');
        const $polyBtn = $('<button type="button">Optimize</button>').click(() => {
            const sec = parseInt($('#dfp-limit').val(), 10) || 120;
            DFP.solveFromPolygons(sec).catch(e => alert(e && e.message ? e.message : String(e)));
        });
        const $stopBtn  = $('<button type="button">Stop</button>').click(() => {
            if (DFP.stopOptimization) DFP.stopOptimization();
        });
        const $exportBtn = $('<button type="button">Export to Draw Tools</button>').click(() => {
            try {
                DFP.exportPlanToDrawTools({
                    includeHull: true,
                    includeTriang: true,
                    includeMicro: true
                });
                alert('Exported to Draw Tools.');
            } catch (e) {
                alert(e && e.message ? e.message : String(e));
            }
        });

        $controls.append($limitLbl, $limitInp, $polyBtn, $stopBtn,);
        $controls2.append($exportBtn);

        const $body = $('<div id="dfp-panel-body">').append($status, $controls, $controls2);

        DFP.dlg = dialog({
            title: 'Δ Field Planner',
            html: $body,
            id: 'dfp-panel',
            width: 300,
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

    // --- polygons in draw-tools -> [{ outerRing:[{lat,lng}], holes:[[...]] }, ...]
    DFP.collectSearchPolygons = function collectSearchPolygons() {
        if (!window.plugin.drawTools || !window.plugin.drawTools.drawnItems)
            throw new Error('draw-tools is not available.');
        const res = [];

        window.plugin.drawTools.drawnItems.eachLayer(drawnItem => {
            if (drawnItem instanceof L.GeodesicPolygon) {
                // 参考 poly-counts2: 用 _latlngs 近似采样曲线边界  :contentReference[oaicite:3]{index=3}
                const latlngs = drawnItem._latlngs;
                if (latlngs && latlngs.length) {
                    const polys = Array.isArray(latlngs[0].lng)
                    ? latlngs.map(r => r.map(pt => ({lat:pt.lat, lng:pt.lng})))
                    : [latlngs.map(pt => ({lat:pt.lat, lng:pt.lng}))];
                    polys.forEach(poly => res.push({ type:'polygon', outerRing: poly, holes: [] }));
                }
            } else if (drawnItem instanceof L.Polygon || (typeof L.MultiPolygon==="function" && drawnItem instanceof L.MultiPolygon)) {
                // 参考 poly-counts2: toGeoJSON() -> coordinates -> [ [ring], [hole], ... ]  :contentReference[oaicite:4]{index=4}
                const g = drawnItem.toGeoJSON();
                let coords = g.geometry.coordinates;
                if (coords[0].length===2 && typeof coords[0][0]==="number") coords=[coords]; // normalize
                coords.forEach(polygonCoords => {
                    if (polygonCoords[0].length===2 && typeof polygonCoords[0][0]==="number")
                        polygonCoords=[polygonCoords];
                    const searchPoly = { type:'polygon', outerRing:[], holes:[] };
                    polygonCoords.forEach((linearRing, j) => {
                        const ring = linearRing.map(([lng,lat]) => ({lat,lng}));
                        if (j===0) searchPoly.outerRing = ring; else searchPoly.holes.push(ring);
                    });
                    res.push(searchPoly);
                });
            } else {
                // ignore markers/polyline/circle  (circle 可后续需要时再加)
            }
        });

        return res;
    };

    // --- portals ∩ polygons
    DFP.portalsInPolygons = function portalsInPolygons(portalList, polys) {
        const picked = [];
        for (const p of portalList) {
            for (const poly of polys) {
                if (poly.type==='polygon' && DFP.pointInSearchPolygonInclusive({lat:p.lat,lng:p.lng}, poly)) {
                    picked.push(p); break;
                }
            }
        }
        return picked;
    };

    DFP.collectExclusionKeys = function collectExclusionKeys(portalList) {
        const keys = new Set();
        if (!window.plugin.drawTools || !window.plugin.drawTools.drawnItems) return keys;
        window.plugin.drawTools.drawnItems.eachLayer(layer => {
            if (layer instanceof L.Marker) {
                const {lat,lng} = layer.getLatLng();
                const m = DFP.matchMarkerToPortal({lat,lng}, portalList);
                if (m) keys.add(DFP.key({lat:m.lat, lng:m.lng}));
            }
        });
        return keys;
    };

    DFP.convexHullIdx = function convexHullIdx(points) {
        // points: Array<{lat,lng,ref}>
        const n = points.length;
        if (n < 3) return [];
        const pts = points.map((p,idx)=>({x:p.lng, y:p.lat, idx}));
        pts.sort((a,b)=> a.x===b.x ? a.y-b.y : a.x-b.x);
        const cross = (o,a,b)=> (a.x-o.x)*(b.y-o.y) - (a.y-o.y)*(b.x-o.x);

        const lower=[];
        for (const p of pts) {
            while (lower.length>=2 && cross(lower[lower.length-2], lower[lower.length-1], p) <= 0) lower.pop();
            lower.push(p);
        }
        const upper=[];
        for (let i=pts.length-1;i>=0;i--) {
            const p=pts[i];
            while (upper.length>=2 && cross(upper[upper.length-2], upper[upper.length-1], p) <= 0) upper.pop();
            upper.push(p);
        }
        const hull = lower.slice(0,-1).concat(upper.slice(0,-1));
        return hull.map(h=>h.idx); // 索引相对于传入的 points 数组
    };


    // ---- geometry & optimization for best triangulation (time-bounded) ----

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

    // —— 将 "u-v" 边键转为坐标
    DFP._edgeToLatLngPair = function(edgeKey, points) {
        const [u, v] = edgeKey.split('-').map(Number);
        const P = points[u], Q = points[v];
        return [[P.lat, P.lng], [Q.lat, Q.lng]];
    };

    // —— 统一的“画一条折线到 draw-tools”
    // 优先使用 Draw Tools Plus 的 drawPolyline（与 link-prolongation 相同调用方式）
    // 否则回退到原生 draw-tools：直接往 drawnItems 加一条 L.polyline
    DFP._drawOnePolylineToDT = function(latlngs, color) {
        if (window.plugin.drawToolsPlus && typeof window.plugin.drawToolsPlus.drawPolyline === 'function') {
            // Draw Tools Plus 路径（参考 link-prolongation） :contentReference[oaicite:1]{index=1}
            window.plugin.drawToolsPlus.drawPolyline(latlngs.map(([la,ln]) => L.latLng(la,ln)), color);
        } else {
            throw new Error('Draw Tools Plus is required for exporting items.');
        }
    };

    // —— 导出当前 DFP.lastPlan 为 draw-tools 的多条 polyline
    // options: { includeHull?:true, includeTriang?:true, includeMicro?:true, colorHull?, colorTri?, colorMicro? }
    DFP.exportPlanToDrawTools = function(options = {}) {
        if (!(window.plugin.drawToolsPlus && typeof window.plugin.drawToolsPlus.drawPolyline === 'function')) {
            throw new Error('Draw Tools Plus is required. Please install/enable the Draw Tools Plus plugin.');
        }
        const plan = DFP.lastPlan;
        if (!plan) throw new Error('No plan to export. Run the planner first.');

        const {
            includeHull = true,
            includeTriang = true,
            includeMicro = true,
            colorHull = '#9e9e9e',
            colorTri = '#9e9e9e',
            colorMicro = '#9e9e9e'
        } = options;

        const points = plan.points;
        const addEdge = (edgeKey, color) => {
            const latlngs = DFP._edgeToLatLngPair(edgeKey, points);
            DFP._drawOnePolylineToDT(latlngs, color);
        };

        // 1) 凸包边
        if (includeHull && plan.hull?.length >= 3) {
            for (let i = 0; i < plan.hull.length; i++) {
                const u = plan.hull[i], v = plan.hull[(i+1) % plan.hull.length];
                addEdge(`${Math.min(u,v)}-${Math.max(u,v)}`, colorHull);
            }
        }

        // 2) 三角剖分对角线
        if (includeTriang && plan.diagonals?.length) {
            for (const key of plan.diagonals) addEdge(key, colorTri);
        }

        // 3) 内部 microfield 边
        if (includeMicro && plan.microEdges?.length) {
            for (const key of plan.microEdges) addEdge(key, colorMicro);
        }
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

    // cooperative yield: let UI timers/paint run
    DFP._lastYield = 0;
    DFP._yieldIntervalMs = 25; // try to yield every ~25ms
    DFP.yieldIfNeeded = async function yieldIfNeeded() {
        const now = Date.now();
        if (now - DFP._lastYield >= DFP._yieldIntervalMs) {
            DFP._lastYield = now;
            // yield to macrotask queue so setInterval/requestAnimationFrame can run
            await new Promise(resolve => setTimeout(resolve, 0));
        }
    };


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
    // Replace DFP.equil with a metric-safe version (geodesic edges + Heron area)
    DFP.equil = function equil(a, b, c) {
        // Use Leaflet geodesic distances (meters) to avoid lat/lng anisotropy.
        const A = L.latLng(a.lat, a.lng);
        const B = L.latLng(b.lat, b.lng);
        const C = L.latLng(c.lat, c.lng);

        const ab = A.distanceTo(B);
        const bc = B.distanceTo(C);
        const ca = C.distanceTo(A);

        // Heron area (meters^2), robust to orientation
        const s  = 0.5 * (ab + bc + ca);
        const sq = Math.max(s * (s - ab) * (s - bc) * (s - ca), 0);
        const area = Math.sqrt(sq);

        const Lmax = Math.max(ab, bc, ca);
        if (Lmax <= 0 || area <= 0) return 0;

        // E in [0,1]
        const E = (4 * area) / (Math.sqrt(3) * Lmax * Lmax);
        return Math.min(1, E);
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
        for (let idx = 0; idx < DFP.points.length; idx++) {
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
    // Best(i,j,k) -> Promise<{score:number, choice:number|null}>
    DFP.Best = async function Best(i, j, k) {
        const baseTri = DFP.triScore(i, j, k);
        const key = DFP.keyOf(i, j, k);
        const hit = DFP.memo.get(key);
        if (hit) return hit;

        const cand = DFP.getInsideForKey(key);
        if (cand.length === 0) {
            const base = { score: baseTri, choice: null };
            DFP.memo.set(key, base);
            return base;
        }

        // time cap: greedy fallback (still async so UI can update)
        if (DFP.timeExceeded()) {
            let bestP = null, bestGain = -Infinity;
            for (const p of cand) {
                const gain = DFP.triScore(i, j, p) + DFP.triScore(j, k, p) + DFP.triScore(k, i, p);
                if (gain > baseTri + bestGain) { bestGain = gain; bestP = p; }
                await DFP.yieldIfNeeded();
            }
            const g = { score: bestGain, choice: bestP };
            DFP.memo.set(key, g);
            return g;
        }

        let best = baseTri, bestP = null;
        for (const p of cand) {
            await DFP.yieldIfNeeded(); // let elapsed/UI update

            const s1 = (await DFP.Best(i, j, p)).score;
            const s2 = (await DFP.Best(j, k, p)).score;
            const s3 = (await DFP.Best(k, i, p)).score;
            const total = baseTri + s1 + s2 + s3;

            if (total > best) { best = total; bestP = p; }

            if (DFP.timeExceeded()) break; // bail early if cap hit mid-loop
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

    // 判断点P是否在有向线段AB上（带公差）
    DFP._EPS = 1e-9; // 角点/边界容差（度）
    DFP._pointOnSegment = function(A, B, P) {
        const ax=A.lat, ay=A.lng, bx=B.lat, by=B.lng, px=P.lat, py=P.lng;
        // 面积接近0（叉积）
        const cross = (bx-ax)*(py-ay) - (by-ay)*(px-ax);
        if (Math.abs(cross) > DFP._EPS) return false;
        // 投影在线段范围内（点积）
        const dot = (px-ax)*(bx-ax) + (py-ay)*(by-ay);
        if (dot < -DFP._EPS) return false;
        const len2 = (bx-ax)*(bx-ax) + (by-ay)*(by-ay);
        if (dot - len2 > DFP._EPS) return false;
        return true;
    };

    DFP._pointOnRing = function(ring, P) {
        for (let i=0, j=ring.length-1; i<ring.length; j=i++) {
            if (DFP._pointOnSegment(ring[j], ring[i], P)) return true;
        }
        return false;
    };

    // 含洞的“宽松包含”
    DFP.pointInSearchPolygonInclusive = function(point, searchPoly) {
        // 在外圈边界上 -> 视为 inside
        if (DFP._pointOnRing(searchPoly.outerRing, point)) return true;
        // 先按原 pnpoly 判外圈
        let inside = DFP.pnpoly(searchPoly.outerRing, point);
        if (!inside) return false;
        // 在任一洞边界上 -> 视为 outside
        for (const hole of searchPoly.holes) {
            if (DFP._pointOnRing(hole, point)) return false;
            if (DFP.pnpoly(hole, point)) return false;
        }
        return true;
    };


    // 权值缓存（与 memo/insideCache/scoreCache 并列）
    DFP.FCache = new Map();

    DFP.weightOfTriangle = async function weightOfTriangle(i,k,j) {
        const key = DFP.keyOf(i,k,j); // 规范化 CCW key  :contentReference[oaicite:10]{index=10}
        const hit = DFP.FCache.get(key);
        if (hit) return hit.score;
        const res = await DFP.Best(i,k,j); // 你的“节点计分”版 Best（已含本三角分数）
        DFP.FCache.set(key, res);
        return res.score;
    };

    // 凸包顶点序列 H（为 DFP.points 的索引） -> {score,tris,diagonals}
    DFP.triangulateHullMax = async function triangulateHullMax(H) {
        const m = H.length;
        if (m < 3) return { score:0, tris:[], diagonals:[] };

        const dp = Array.from({length:m}, ()=>Array(m).fill(0));
        const cut = Array.from({length:m}, ()=>Array(m).fill(-1));

        for (let len=2; len<m; len++) {
            for (let i=0; i+len<m; i++) {
                const j = i+len;
                if (j<=i+1) { dp[i][j]=0; continue; }
                let best = -Infinity, bestK = -1;
                for (let k=i+1; k<j; k++) {
                    const ii=H[i], kk=H[k], jj=H[j];
                    const w = await DFP.weightOfTriangle(ii,kk,jj);
                    const cand = dp[i][k] + dp[k][j] + w;
                    if (cand > best) { best=cand; bestK=k; }
                    await DFP.yieldIfNeeded?.();
                }
                dp[i][j]=best; cut[i][j]=bestK;
            }
        }

        // 回溯收集三角形与对角线
        const tris = [];
        const diags = new Set();
        (function rec(i,j){
            if (j<=i+1) return;
            const k = cut[i][j];
            const a=H[i], b=H[k], c=H[j];
            tris.push([a,b,c]);
            const add=(u,v)=>{ const x=Math.min(u,v), y=Math.max(u,v); diags.add(`${x}-${y}`); };
            add(a,b); add(b,c);
            rec(i,k); rec(k,j);
        })(0,m-1);

        return { score: dp[0][m-1], tris, diagonals:[...diags] };
    };

    DFP.solveFromPolygons = async function solveFromPolygons(secondsLimit) {
        const portalList = DFP.buildPortalList();              // 已有  :contentReference[oaicite:11]{index=11}
        const polys = DFP.collectSearchPolygons();             // 新
        if (!polys.length) throw new Error('No polygon found. Draw polygon(s) first.');

        // 参与集合 S
        const S = DFP.portalsInPolygons(portalList, polys);

        // 排除集合（marker 匹配真实门户）  :contentReference[oaicite:12]{index=12}
        const exKeys = DFP.collectExclusionKeys(portalList);
        const filtered = S.filter(p => !exKeys.has(DFP.key(p)));
        if (filtered.length < 3) throw new Error('Not enough portals inside polygons after exclusions.');

        // 凸包（基于 filtered 的索引）
        const hullLocalIdx = DFP.convexHullIdx(filtered);
        if (hullLocalIdx.length < 3) throw new Error('Convex hull has <3 vertices.');

        // DFP.points：将“凸包点在前，内点在后”的顺序重排，方便 DP
        const isHull = new Array(filtered.length).fill(false);
        hullLocalIdx.forEach(id => { isHull[id]=true; });
        const hullPoints = hullLocalIdx.map(id => filtered[id]);
        const innerPoints = filtered.filter((_,idx)=>!isHull[idx]);
        DFP.points = hullPoints.concat(innerPoints);      // 0..h-1 为 hull 顶点
        const H = [...Array(hullPoints.length).keys()];  // [0..h-1]

        // 计时/缓存
        DFP.memo.clear(); DFP.insideCache.clear(); DFP.scoreCache.clear(); DFP.FCache.clear();
        DFP.stopFlag = false;
        DFP.limitMs = Math.max(1, Math.floor((secondsLimit || 120) * 1000));
        DFP.startMs = Date.now(); DFP._lastYield = DFP.startMs;

        const tick = setInterval(() => {
            const el = document.getElementById('dfp-elapsed');
            if (el) el.textContent = `${Math.floor((Date.now()-DFP.startMs)/1000)}s`;
        }, 250);

        // 求解
        const plan = await DFP.triangulateHullMax(H);

        // 渲染
        if (DFP.layerGroup) DFP.layerGroup.clearLayers();

        // 画凸包边（虚线）
        for (let i=0;i<H.length;i++){
            const a=DFP.points[H[i]], b=DFP.points[H[(i+1)%H.length]];
            DFP.layerGroup.addLayer(L.polyline([[a.lat,a.lng],[b.lat,b.lng]], {
                color: DFP.NEUTRAL_LINK_COLOR, weight: 3, opacity: 1.0, dashArray:'8,6'
            }));
        }

        // 画三角剖分的内部对角线（实线，粗一点）
        for (const key of plan.diagonals) {
            const [u,v] = key.split('-').map(Number);
            const P = DFP.points[u], Q = DFP.points[v];
            DFP.layerGroup.addLayer(L.polyline([[P.lat,P.lng],[Q.lat,Q.lng]], {
                color: DFP.NEUTRAL_LINK_COLOR, weight: 3, opacity: 1.0, dashArray:'8,6'
            }));
        }

        // 三角内 microfield 边集合（合并去重）
        const edgeSet = new Set(); let faceSum = 0;
        for (const [i,k,j] of plan.tris) {
            const {edges, faces} = DFP.reconstructEdges(i,k,j);
            faceSum += faces;
            edges.forEach(e => edgeSet.add(e));
        }
        faceSum += plan.tris.length;
        edgeSet.forEach(key => {
            const [u,v] = key.split('-').map(Number);
            const P=DFP.points[u], Q=DFP.points[v];
            DFP.layerGroup.addLayer(L.polyline([[P.lat,P.lng],[Q.lat,Q.lng]], {
                color: DFP.NEUTRAL_LINK_COLOR, weight: 2, opacity: 0.9
            }));
        });

        clearInterval(tick);
        document.getElementById('dfp-best')?.replaceChildren(document.createTextNode(Math.round(plan.score)));
        document.getElementById('dfp-faces')?.replaceChildren(document.createTextNode(String(faceSum)));
        const el = document.getElementById('dfp-elapsed');
        if (el) el.textContent = `${Math.floor((Date.now()-DFP.startMs)/1000)}s`;

        // 渲染完毕后（faces/elapsed 已更新）
        // 缓存本次规划结果，供“导出到 Draw Tools”使用
        DFP.lastPlan = {
            points: DFP.points,                 // Array<{lat,lng,ref}>
            hull: H.slice(),                    // 凸包顶点索引序列（相对 points）
            tris: plan.tris.slice(),            // 三角剖分三角形索引三元组
            diagonals: plan.diagonals.slice(),  // "u-v" 字符串
            microEdges: [...edgeSet]            // "u-v" 字符串（内部微田字边，已去重）
        };

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
