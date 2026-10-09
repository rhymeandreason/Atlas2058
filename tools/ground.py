# Bake ground.png, the texture draped over the terrain (W×H pixels spanning the map bbox, row 0 = north):
#   R  ground-cover class (see CLASSES; colours live in index.html as --lc-* tokens)
#   G  NOAA sea-level-rise threshold: 0 = dry even at +10 ft (or no data), k ≥ 1 = first flooded at (k-1)/2 ft
#      above today's high tide (MHHW). k = 1 is already wet at today's high tide.
#   B  255 inside parks and protected open space, else 0
#
# Inputs (fetch once; see each source's note):
#   nlcd.png    NLCD 2021 land cover, MRLC WMS GetMap layers=NLCD_2021_Land_Cover_L48, srs=EPSG:4326,
#               bbox=-122.665,37.885,-122.355,38.105, width=2600, height=2300, format=image/png
#   lc*.json    Overpass land-use, leisure, natural, amenity, protected-area and pier features. Use `out geom`
#               (not `out tags geom`, which drops relation members, i.e. every multipolygon park and preserve)
#   coast.json  Overpass natural=coastline (as for terrain.py)
#   tiles/      NOAA SLR Viewer depth tiles, zoom 14, named {level}_{y}_{x}.png, from
#               https://coast.noaa.gov/arcgis/rest/services/dc_slr/slr_{level}/MapServer/tile/14/{y}/{x}
#               for level in 0ft, 0_5ft, … 10ft over the low-lying shore
# Usage: python3 tools/ground.py nlcd.png coast.json tiles/ ground.png lc1.json [lc2.json …]
import sys, os, json, math
import numpy as np
from PIL import Image, ImageDraw
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coast import land_mask, LONMIN, LONMAX, LATMIN, LATMAX

W, H = 2600, 2300
nlcd_path, coast_path, tile_dir, out_path = sys.argv[1:5]
lc_paths = sys.argv[5:]
CLASSES = ['land', 'forest', 'scrub', 'grass', 'sand', 'wetland', 'water', 'park', 'golf', 'sports', 'campus',
           'retail', 'industrial', 'parking', 'pier', 'cemetery']
C = {n: i for i, n in enumerate(CLASSES)}

# --- base vegetation from NLCD (standard legend colours) ---
NLCD = { (70,107,159): 'water', (222,197,197): 'land', (217,146,130): 'land', (235,0,0): 'land', (171,0,0): 'land',
         (179,172,159): 'sand', (104,171,95): 'forest', (28,95,44): 'forest', (181,197,143): 'forest', (204,184,121): 'scrub',
         (223,223,194): 'grass', (220,217,57): 'grass', (171,108,40): 'grass', (184,217,235): 'wetland', (108,159,184): 'wetland' }
rgb = np.array(Image.open(nlcd_path).convert('RGB')).astype(int)
assert rgb.shape[:2] == (H, W), rgb.shape
keys = np.array(list(NLCD.keys())); vals = np.array([C[v] for v in NLCD.values()], np.uint8)
best = np.full((H, W), 1e9); cls = np.zeros((H, W), np.uint8)
for k, v in zip(keys, vals):
    d = ((rgb - k) ** 2).sum(axis=2); m = d < best; best[m] = d[m]; cls[m] = v
cls[best > 900] = C['land']  # unknown colours (edges, no data)
cls[cls == C['water']] = C['land']  # open water comes from the coastline and OSM below, not 30 m NLCD cells

# --- OSM features, painted lowest priority first ---
def px(lat, lon): return ((lon - LONMIN) / (LONMAX - LONMIN) * W, (LATMAX - lat) / (LATMAX - LATMIN) * H)
def stitch(ways):
    rings, segs = [], [w[:] for w in ways if len(w) > 1]
    while segs:
        cur = segs.pop(); changed = True
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
def geom(e):
    """(outer rings, inner rings, open lines) in pixel coords"""
    if e['type'] == 'way':
        g = [px(p['lat'], p['lon']) for p in e.get('geometry', [])]
        return ([g], [], []) if len(g) >= 4 and g[0] == g[-1] else ([], [], [g])
    if e['type'] == 'relation':
        mem = lambda role: [[px(p['lat'], p['lon']) for p in m.get('geometry', [])] for m in e.get('members', []) if m.get('type') == 'way' and m.get('role') == role]
        return stitch(mem('outer') + mem('')), stitch(mem('inner')), []
    return [], [], []
def classify(t):
    lu, le, na, am, mm = t.get('landuse'), t.get('leisure'), t.get('natural'), t.get('amenity'), t.get('man_made')
    if mm in ('pier', 'breakwater', 'groyne'): return 'pier'
    if na == 'water' or lu in ('reservoir', 'basin') or t.get('waterway') == 'riverbank' or le == 'marina': return 'water'
    if na in ('beach', 'sand', 'bare_rock'): return 'sand'
    if na == 'wetland': return 'wetland'
    if am == 'parking': return 'parking'
    if le in ('pitch', 'track', 'playground'): return 'sports'
    if le == 'golf_course': return 'golf'
    if lu in ('retail', 'commercial') or t.get('shop') == 'mall': return 'retail'
    if lu == 'industrial': return 'industrial'
    if am in ('school', 'university', 'college', 'hospital'): return 'campus'
    if lu == 'cemetery' or am == 'grave_yard': return 'cemetery'
    if na in ('wood',) or lu == 'forest': return 'forest'
    if na in ('scrub', 'heath'): return 'scrub'
    if na == 'grassland' or lu in ('grass', 'meadow', 'village_green', 'farmland', 'vineyard', 'orchard'): return 'grass'
    if le in ('park', 'garden', 'recreation_ground', 'nature_reserve') or lu == 'recreation_ground' or t.get('boundary') in ('protected_area', 'national_park'): return 'park'
    return None
ORDER = ['grass', 'scrub', 'forest', 'cemetery', 'park', 'campus', 'industrial', 'retail', 'golf', 'sports', 'parking', 'wetland', 'sand', 'water', 'pier']
masks = {n: Image.new('L', (W, H), 0) for n in ORDER}
park = Image.new('L', (W, H), 0)
for f in lc_paths:
    for e in json.load(open(f))['elements']:
        t = e.get('tags', {}); c = classify(t)
        if not c: continue
        outer, inner, lines = geom(e)
        d = ImageDraw.Draw(masks[c])
        for r in outer: d.polygon(r, fill=255)
        for r in inner: d.polygon(r, fill=0)
        for l in lines:
            if c == 'pier' and len(l) >= 2: d.line(l, fill=255, width=2)
        if c == 'park':  # parks also tint whatever vegetation they hold
            dp = ImageDraw.Draw(park)
            for r in outer: dp.polygon(r, fill=255)
            for r in inner: dp.polygon(r, fill=0)
for n in ORDER:
    if n == 'park': continue  # parks keep their vegetation; they only set the B flag
    cls[np.array(masks[n]) > 0] = C[n]

# --- bay from the OSM coastline (piers stay on top) ---
lon = LONMIN + (LONMAX - LONMIN) * (np.arange(W) + 0.5) / W
lat = LATMAX - (LATMAX - LATMIN) * (np.arange(H) + 0.5) / H
LON, LAT = np.meshgrid(lon, lat)
land, _ = land_mask(coast_path, LON, LAT)
pier = np.array(masks['pier']) > 0
cls[~land & ~pier] = C['water']

# --- NOAA connected inundation: first level at which each pixel floods ---
LEVELS = ['0ft', '0_5ft', '1ft', '1_5ft', '2ft', '2_5ft', '3ft', '3_5ft', '4ft', '4_5ft', '5ft', '5_5ft', '6ft', '6_5ft', '7ft', '7_5ft', '8ft', '8_5ft', '9ft', '9_5ft', '10ft']
Z = 14; n = 2 ** Z
gx = (LON + 180) / 360 * n * 256
gy = (1 - np.log(np.tan(np.radians(LAT)) + 1 / np.cos(np.radians(LAT))) / math.pi) / 2 * n * 256
tx, ty = (gx // 256).astype(int), (gy // 256).astype(int)
ix, iy = (gx % 256).astype(int), (gy % 256).astype(int)
flood = np.zeros((H, W), np.uint8)
tiles = sorted(set(zip(tx.ravel(), ty.ravel())))
for li, L in enumerate(LEVELS):
    for (x, y) in tiles:
        p = os.path.join(tile_dir, f'{L}_{y}_{x}.png')
        if not os.path.exists(p) or os.path.getsize(p) == 0: continue
        try: a = np.array(Image.open(p).convert('RGBA'))[:, :, 3]
        except Exception: continue
        m = (tx == x) & (ty == y)
        wet = a[iy[m], ix[m]] > 0
        cur = flood[m]; cur[(cur == 0) & wet] = li + 1; flood[m] = cur
flood[~land] = 0  # the bay itself isn't a flood

parkf = (np.array(park) > 0) & (cls != C['water'])  # refuges reaching into the bay don't tint the water
out = np.stack([cls, flood, np.where(parkf, 255, 0).astype(np.uint8)], axis=2)
Image.fromarray(out, 'RGB').save(out_path, optimize=True)
print('classes', {CLASSES[i]: int((cls == i).sum()) for i in range(len(CLASSES))})
print('flood pixels by level', np.bincount(flood.ravel(), minlength=22)[1:].tolist())
