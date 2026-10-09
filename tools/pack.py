# Pack OSM building footprints into the map's local metre frame (see XY() in index.html).
# Fetch raw.json from Overpass (overpass.kumi.systems if the main server is busy):
#   [out:json][timeout:240];way["building"](37.885,-122.665,38.105,-122.355);out tags geom qt;
#   [out:json][timeout:240];rel["building"](37.885,-122.665,38.105,-122.355);out geom;   (members need full `out geom`)
# Optional heights.csv (t,osm_id,h): Overture Maps building heights for OSM footprints, mostly Microsoft ML
# estimates from imagery (ODbL). Used when OSM has no height or levels tag. Export with DuckDB from
#   s3://overturemaps-us-west-2/release/<release>/theme=buildings/type=building/*  (sources[1].record_id = w123@1)
# then: python3 tools/pack.py ways.json rels.json [heights.csv] buildings.json
import json, sys, math, re
srcs, out = [a for a in sys.argv[1:-1] if not a.endswith('.csv')], sys.argv[-1]
ML = {}
for a in sys.argv[1:-1]:
    if a.endswith('.csv'):
        for line in open(a).read().splitlines()[1:]:
            t, i, h = line.split(','); ML[(t, int(i))] = float(h)
D2R = math.pi / 180; LAT0, LON0 = 37.983, -122.515
KX = 111320 * math.cos(LAT0 * D2R); KZ = 110950
def xy(lat, lon): return ((lon - LON0) * KX, -(lat - LAT0) * KZ)
KIND = {}  # 0 untagged, 1 house, 2 apartments, 3 commercial, 4 industrial, 5 minor/shed, 6 civic
for t in ['house','detached','semidetached_house','terrace','residential','bungalow','cabin','farm']: KIND[t] = 1
for t in ['apartments','dormitory','hotel']: KIND[t] = 2
for t in ['commercial','retail','office','supermarket','kiosk']: KIND[t] = 3
for t in ['industrial','warehouse','hangar','manufacture','storage_tank']: KIND[t] = 4
for t in ['garage','garages','carport','shed','roof','static_caravan','hut','greenhouse','boathouse','houseboat','shelter']: KIND[t] = 5
for t in ['school','university','college','government','hospital','church','religious','public','civic','chapel','fire_station','train_station','transportation','sports_hall']: KIND[t] = 6

def dp(pts, tol):
    if len(pts) < 3: return pts
    ax, az = pts[0]; bx, bz = pts[-1]; dx, dz = bx - ax, bz - az; L = math.hypot(dx, dz) or 1e-9
    best, bi = -1, 0
    for i in range(1, len(pts) - 1):
        px, pz = pts[i]
        d = abs(dx * (az - pz) - dz * (ax - px)) / L if L > 1e-6 else math.hypot(px - ax, pz - az)
        if d > best: best, bi = d, i
    if best <= tol: return [pts[0], pts[-1]]
    return dp(pts[:bi + 1], tol)[:-1] + dp(pts[bi:], tol)
def simplify(ring, tol):
    # split the closed ring at its farthest point pair so DP keeps the shape
    ring = ring[:-1] if ring[0] == ring[-1] else ring
    if len(ring) < 4: return ring
    far = max(range(len(ring)), key=lambda i: (ring[i][0] - ring[0][0]) ** 2 + (ring[i][1] - ring[0][1]) ** 2)
    a = dp(ring[:far + 1], tol); b = dp(ring[far:] + [ring[0]], tol)
    return a[:-1] + b[:-1]
def area(r): return 0.5 * sum(r[i][0] * r[(i + 1) % len(r)][1] - r[(i + 1) % len(r)][0] * r[i][1] for i in range(len(r)))
def num(s):
    m = re.match(r'\s*([0-9.]+)', s or '');
    try: return float(m.group(1)) if m else None
    except ValueError: return None
def stitch(ways):
    rings, segs = [], [w[:] for w in ways if len(w) > 1]
    while segs:
        cur = segs.pop()
        changed = True
        while cur[0] != cur[-1] and changed:
            changed = False
            for i, s in enumerate(segs):
                if s[0] == cur[-1]: cur += s[1:]
                elif s[-1] == cur[-1]: cur += s[::-1][1:]
                elif s[-1] == cur[0]: cur = s + cur[1:]
                elif s[0] == cur[0]: cur = s[::-1] + cur[1:]
                else: continue
                segs.pop(i); changed = True; break
        if cur[0] == cur[-1] and len(cur) >= 4: rings.append(cur)
    return rings

datas = [json.load(open(f)) for f in srcs]
outl, npts = [], 0
for e in (e for d in datas for e in d['elements']):
    t = e.get('tags', {})
    if e['type'] == 'way': rings = [[(g['lat'], g['lon']) for g in e.get('geometry', [])]]; rings = [r for r in rings if len(r) >= 4 and r[0] == r[-1]]
    else: rings = stitch([[(g['lat'], g['lon']) for g in m.get('geometry', [])] for m in e.get('members', []) if m.get('type') == 'way' and m.get('role') in ('outer', '')])
    kind = KIND.get(t.get('building'), 0)
    h = num(t.get('height')); lv = num(t.get('building:levels'))
    if not h and not lv: h = ML.get(('w' if e['type'] == 'way' else 'r', e['id']))
    hdm = int(round((h if h else lv * 3.3 if lv else 0) * 10))
    for r in rings:
        p = [xy(a, b) for a, b in r]
        p = simplify(p, 0.9)
        if len(p) < 3: continue
        ar = area(p)
        if abs(ar) < 12: continue
        if ar < 0: p = p[::-1]  # consistent winding
        q = [(int(round(x * 2)), int(round(z * 2))) for x, z in p]
        dq = [q[0][0], q[0][1]]
        for i in range(1, len(q)): dq += [q[i][0] - q[i - 1][0], q[i][1] - q[i - 1][1]]
        outl.append([kind, hdm] + dq); npts += len(q)
json.dump({ 'src': 'OpenStreetMap contributors, ODbL', 'date': datas[0]['osm3s']['timestamp_osm_base'][:10], 'q': 0.5, 'b': outl }, open(out, 'w'), separators=(',', ':'))
print(len(outl), 'buildings', npts, 'points')
