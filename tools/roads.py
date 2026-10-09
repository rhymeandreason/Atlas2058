# Pack OpenStreetMap streets into roads.json in the map's local metre frame (see XY() in index.html).
# Freeways (motorway, motorway_link) are left out: US-101 and I-580 stay hand-traced because the
# 101→580 connector project is built along them. Tunnels are dropped.
# Fetch roads.json from Overpass (overpass.kumi.systems if the main server is busy):
#   [out:json][timeout:240];way["highway"~"^(trunk|trunk_link|primary|primary_link|secondary|secondary_link|
#   tertiary|tertiary_link|residential|unclassified|living_street)$"](37.885,-122.665,38.105,-122.355);out tags geom qt;
# Usage: python3 tools/roads.py raw_roads.json roads.json
# Each entry is [class, bridge, x0, z0, dx1, dz1, …] in units of q metres;
# class 1 trunk/primary, 2 secondary, 3 tertiary, 4 residential/minor, 5 ramps/links.
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
CLS = { 'trunk': 1, 'primary': 1, 'secondary': 2, 'tertiary': 3, 'residential': 4, 'unclassified': 4, 'living_street': 4,
        'trunk_link': 5, 'primary_link': 5, 'secondary_link': 5, 'tertiary_link': 5 }
outl, npts = [], 0
for e in json.load(open(src))['elements']:
    t = e.get('tags', {}); c = CLS.get(t.get('highway'))
    if not c or t.get('tunnel') in ('yes', 'building_passage', 'culvert'): continue
    p = [((g['lon'] - LON0) * KX, -(g['lat'] - LAT0) * KZ) for g in e.get('geometry', [])]
    if len(p) < 2: continue
    p = dp(p, 1.5)
    q = [(int(round(x * 2)), int(round(z * 2))) for x, z in p]
    dq = [q[0][0], q[0][1]]
    for i in range(1, len(q)): dq += [q[i][0] - q[i - 1][0], q[i][1] - q[i - 1][1]]
    outl.append([c, 1 if t.get('bridge') and t.get('bridge') != 'no' else 0] + dq); npts += len(q)
json.dump({ 'src': 'OpenStreetMap contributors, ODbL', 'q': 0.5, 'r': outl }, open(out, 'w'), separators=(',', ':'))
print(len(outl), 'ways', npts, 'points')
