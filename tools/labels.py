# Bake labels.json: landmark, open-space and road labels from OpenStreetMap.
# Inputs: the Overpass files used for roads.py / ground.py, plus named.json:
#   [out:json];(nwr["name"~"Civic Center|Mission San Rafael|Dominican University|Marin Country Mart|
#   Corte Madera Town Center|Village at Corte Madera|Northgate|San Quentin"](37.885,-122.665,38.105,-122.355););out tags center qt;
# Usage: python3 tools/labels.py roads_raw.json named.json labels.json lc1.json [lc2.json …]
# Output: { l: [[kind, name, lat, lon], …] } with kind p (landmark), o (open space), r (road)
import sys, json, math

roads_path, named_path, out_path = sys.argv[1:4]
lc_paths = sys.argv[4:]
LONMIN, LONMAX, LATMIN, LATMAX = -122.665, -122.355, 37.885, 38.105
KX, KZ = 111320 * math.cos(37.983 * math.pi / 180), 110950
inside = lambda la, lo: LATMIN + 0.004 < la < LATMAX - 0.004 and LONMIN + 0.004 < lo < LONMAX - 0.004
out, seen = [], set()
def add(kind, name, la, lo):
    if name and name not in seen and inside(la, lo): seen.add(name); out.append([kind, name, round(la, 5), round(lo, 5)])
def ring_area_centroid(g):  # metres², (lat, lon)
    a = cx = cy = 0
    for (la1, lo1), (la2, lo2) in zip(g, g[1:]):
        x1, y1, x2, y2 = lo1 * KX, la1 * KZ, lo2 * KX, la2 * KZ; c = x1 * y2 - x2 * y1
        a += c; cx += (x1 + x2) * c; cy += (y1 + y2) * c
    if abs(a) < 1e-9: return 0, g[0]
    return abs(a) / 2, (cy / (3 * a) / KZ, cx / (3 * a) / KX)
def stitch(ways):  # join multipolygon member fragments into closed rings
    rs, segs = [], [w[:] for w in ways if len(w) > 1]
    while segs:
        cur = segs.pop(); changed = True
        while cur[0] != cur[-1] and changed:
            changed = False
            for i, q in enumerate(segs):
                if q[0] == cur[-1]: cur += q[1:]
                elif q[-1] == cur[-1]: cur += q[::-1][1:]
                elif q[-1] == cur[0]: cur = q + cur[1:]
                elif q[0] == cur[0]: cur = q[::-1] + cur[1:]
                else: continue
                segs.pop(i); changed = True; break
        if cur[0] == cur[-1] and len(cur) >= 4: rs.append(cur)
    return rs
def rings(e):
    if e['type'] == 'way': return [[(p['lat'], p['lon']) for p in e.get('geometry', [])]]
    return stitch([[(p['lat'], p['lon']) for p in m.get('geometry', [])] for m in e.get('members', []) if m.get('role') in ('outer', '')])
def area(e): return max([ring_area_centroid(r)[0] for r in rings(e) if len(r) > 2] or [0])
def where(e):
    if 'lat' in e: return e['lat'], e['lon']
    if 'center' in e: return e['center']['lat'], e['center']['lon']
    best = (0, None)
    for r in rings(e):
        if len(r) > 2: best = max(best, ring_area_centroid(r), key=lambda t: t[0])
    return best[1] if best[1] else (None, None)

# landmarks: transit, campuses, hospitals, malls, marinas, golf
feats = [e for f in lc_paths for e in json.load(open(f))['elements']]
for e in feats:
    t = e.get('tags', {}); n = t.get('name')
    if not n: continue
    k = t.get('amenity') in ('hospital', 'university', 'college', 'ferry_terminal') or t.get('railway') == 'station' or \
        t.get('shop') == 'mall' or t.get('leisure') in ('marina', 'golf_course')
    # campuses and hospitals need a real footprint; skip named wings, clinics and day-care centres
    if k and t.get('amenity') in ('hospital', 'university', 'college') and area(e) < 15000: k = False
    if k:
        la, lo = where(e)
        if la: add('p', n + ' station' if t.get('railway') == 'station' else n, la, lo)
for e in json.load(open(named_path))['elements']:
    t = e.get('tags', {}); n = t.get('name', '')
    if n in ('Mission San Rafael Arcangel', 'San Quentin State Prison', 'Marin Country Mart', 'The Village at Corte Madera'):
        la, lo = where(e)
        if la: add('p', n.replace('Arcangel', 'Arcángel').replace(' State Prison', ''), la, lo)

# open space: the largest named parks and preserves whose middle falls on the map
os_ = []
for e in feats:
    t = e.get('tags', {}); n = t.get('name')
    if not n or not (t.get('boundary') in ('protected_area', 'national_park') or t.get('leisure') in ('nature_reserve', 'park')): continue
    best = (0, None)
    for r in rings(e):
        if len(r) > 2: best = max(best, ring_area_centroid(r), key=lambda t: t[0])
    if best[1] and best[0] > 250000: os_.append((best[0], n, best[1]))
for a, n, (la, lo) in sorted(os_, reverse=True)[:22]:
    add('o', n, la, lo)

# roads: major named roads, up to three anchors each at least 3 km apart, at the middle of their longest ways
ROADCLS = ('motorway', 'trunk', 'primary', 'secondary')
byname = {}
for e in json.load(open(roads_path))['elements']:
    t = e.get('tags', {})
    if t.get('highway') not in ROADCLS: continue
    n = t.get('ref') if t.get('highway') == 'motorway' else t.get('name')
    if not n: continue
    n = n.split(';')[0].replace('US 101', 'US-101').replace('I 580', 'I-580')
    g = [(p['lat'], p['lon']) for p in e.get('geometry', [])]
    L = sum(math.hypot((b[1] - a[1]) * KX, (b[0] - a[0]) * KZ) for a, b in zip(g, g[1:]))
    byname.setdefault(n, []).append((L, g))
for n, ways in byname.items():
    total = sum(L for L, _ in ways)
    if total < 1200: continue
    anchors = []
    for L, g in sorted(ways, key=lambda w: -w[0]):
        if L < 250 or len(anchors) >= 3: break
        half, acc = L / 2, 0
        for a, b in zip(g, g[1:]):
            s = math.hypot((b[1] - a[1]) * KX, (b[0] - a[0]) * KZ)
            if acc + s >= half: t = (half - acc) / s if s else 0; mid = (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t); break
            acc += s
        if all(math.hypot((mid[1] - q[1]) * KX, (mid[0] - q[0]) * KZ) > 3000 for q in anchors) and inside(*mid):
            anchors.append(mid); out.append(['r', n, round(mid[0], 5), round(mid[1], 5)])
json.dump({ 'src': 'OpenStreetMap contributors, ODbL', 'l': out }, open(out_path, 'w'), ensure_ascii=False, separators=(',', ':'))
from collections import Counter
print(Counter(k for k, *_ in out)); print([x[1] for x in out if x[0] != 'r'])
