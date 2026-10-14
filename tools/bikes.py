# Pack San Rafael's bikeways into bikes.json in the map's local metre frame (see XY() in index.html).
# Source: the Citywide Bicycle & Pedestrian Plan update's draft-recommendations map (Alta Planning for the City),
#   https://sanrafaelrecommendations.altaplanning.cloud/  — its layers are served as one file:
#   curl -A 'Atlas2058' 'https://sanrafaelrecommendations.altaplanning.cloud/_data/features.json' -o raw_bikes.json
# Usage: python3 tools/bikes.py raw_bikes.json bikes.json
# Output:
#   e: existing bikeways, each [class, x0, z0, dx1, dz1, …] in units of q metres; class 1 path, 2 lane, 3 route, 4 separated
#   p: proposed corridors, each [class, top20, name index, x0, z0, dx1, dz1, …]; class 1–4 as above, 5 bike boulevard, 6 TBD
#   i: proposed intersection/crossing projects, each [top20, name index, x, z]
#   n: names and descriptions, each [name, description]
import json, sys, math

def dp(pts, tol):  # Douglas–Peucker
    if len(pts) < 3: return pts
    ax, az = pts[0]; bx, bz = pts[-1]; dx, dz = bx - ax, bz - az; L = math.hypot(dx, dz)
    best, bi = -1, 0
    for i in range(1, len(pts) - 1):
        px, pz = pts[i]
        d = abs(dx * (az - pz) - dz * (ax - px)) / L if L > 1e-6 else math.hypot(px - ax, pz - az)
        if d > best: best, bi = d, i
    if best <= tol: return [pts[0], pts[-1]]
    return dp(pts[:bi + 1], tol)[:-1] + dp(pts[bi:], tol)

src, out = sys.argv[1], sys.argv[2]
D2R = math.pi / 180; LAT0, LON0 = 37.983, -122.515
KX = 111320 * math.cos(LAT0 * D2R); KZ = 110950
Q = 0.5
def xy(c): return ((c[0] - LON0) * KX, -(c[1] - LAT0) * KZ)
def lines(geom):
    if geom['type'] == 'LineString': return [geom['coordinates']]
    if geom['type'] == 'MultiLineString': return geom['coordinates']
    return []
def pack(coords):
    p = dp([xy(c) for c in coords], 1.0)
    q = [(int(round(x / Q)), int(round(z / Q))) for x, z in p]
    dq = [q[0][0], q[0][1]]
    for i in range(1, len(q)): dq += [q[i][0] - q[i - 1][0], q[i][1] - q[i - 1][1]]
    return dq
def clean(s): return ' '.join((s or '').split()).strip(' ,')

layers = {}
for ft in json.load(open(src))['featureTypes']: layers[ft['name']] = ft.get('features_geojson') or []
EXISTING = { 'Bike - Existing - Class I - Shared Use Paths': 1, 'Bike - Existing - Class II - Bike Lanes': 2,
             'Bike - Existing - Class III - Bike Routes': 3, 'Bike - Existing - Class IV - Separated Bikeways': 4 }
PROPOSED = { 'Bike - Recommended - Class I - Shared Use Paths': 1, 'Bike - Recommended - Class II - Bike Lanes': 2,
             'Bike - Recommended - Class III - Bike Routes': 3, 'Bike - Recommended - Class III+ - Bike Routes': 5,
             'Bike - Recommended - Class IV - Separated Bikeways': 4, 'Bike - Recommended - TBD': 6 }
POINTS = ['Bike - Recommended - Point Projects', 'Pedestrian - Recommended - Point Projects']

names, nidx = [], {}
def name(n, d):
    k = (n, d)
    if k not in nidx: nidx[k] = len(names); names.append([n, d])
    return nidx[k]

e = []
for lay, c in EXISTING.items():
    for f in layers[lay]:
        for l in lines(f['geometry']):
            if len(l) > 1: e.append([c] + pack(l))
p = []
for lay, c in PROPOSED.items():
    for f in layers[lay]:
        pr = f['properties']
        desc = clean(pr.get('Bike Descr'))
        ped = clean(pr.get('Ped Descr'))
        if ped and ped.lower() not in desc.lower(): desc += '; walking: ' + ped[0].lower() + ped[1:]
        k = name(clean(pr.get('CORRIDOR')), desc)
        for l in lines(f['geometry']):
            if len(l) > 1: p.append([c, 1 if pr.get('Top 20') == 'Y' else 0, k] + pack(l))
pts, seen = [], set()
for lay in POINTS:
    for f in layers[lay]:
        pr = f['properties']
        if f['geometry']['type'] != 'Point': continue
        x, z = xy(f['geometry']['coordinates'])
        key = (pr.get('ID_PRO'), round(x / 5), round(z / 5))
        if key in seen: continue
        seen.add(key)
        parts = [clean(pr.get('Ped_Descr')), clean(pr.get('Bike_Descr'))]
        k = name(clean(pr.get('FullIntNam')) or clean(pr.get('Int Name')) or 'Crossing improvement', '; '.join(s for s in parts if s))
        pts.append([1 if pr.get('Top 20') == 'Y' else 0, k, int(round(x / Q)), int(round(z / Q))])

json.dump({ 'src': 'City of San Rafael Citywide Bicycle & Pedestrian Plan, draft recommendations (2026)', 'q': Q,
            'e': e, 'p': p, 'i': pts, 'n': names }, open(out, 'w'), separators=(',', ':'))
print(len(e), 'existing', len(p), 'proposed lines', len(pts), 'points', len(names), 'names')
