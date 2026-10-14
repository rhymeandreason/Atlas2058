# Bake trees.bin: individual trees from Meta & WRI's 1 m global canopy height map (CC BY 4.0), in the map's local
# metre frame (see XY() in index.html).
# Input: the canopy-height tile covering the map, a uint8 GeoTIFF in EPSG:3857 (metres of canopy height):
#   https://dataforgood-fb-data.s3.amazonaws.com/forests/v1/alsgedi_global_v6_float/chm/023010203.tif  (~250 MB)
# Needs rasterio and scipy. Usage: python3 tools/trees.py 023010203.tif buildings.json trees.bin
#
# Treetops are local maxima of the lightly smoothed canopy at ~2 m resolution, at least 4 m tall and outside OSM
# building footprints (the model sometimes reads roofs as canopy). Every treetop within ~60 m of a building is kept,
# so streets and yards show their real trees; on the wooded hills, which the ground layer already shades as forest,
# only the tallest treetop per ~25 m cell is kept to give the slopes texture without millions of trees.
# Output (little-endian): uint32 n, then int16 x[n], int16 z[n] (metres), uint8 h[n] (quarter metres).
# Also prints tree cover (share of ground under canopy ≥ 3 m) around each neighbourhood label.
import sys, math, json, struct
import numpy as np, rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import maximum_filter, gaussian_filter
from PIL import Image, ImageDraw

src, bpath, out = sys.argv[1:4]
D2R = math.pi / 180; LAT0, LON0 = 37.983, -122.515
KX = 111320 * math.cos(LAT0 * D2R); KZ = 110950
R = 6378137.0
def merc(lon, lat): return (lon * D2R * R, R * math.log(math.tan(math.pi / 4 + lat * D2R / 2)))
def unmerc(mx, my): return (mx / R / D2R, (2 * math.atan(math.exp(my / R)) - math.pi / 2) / D2R)
MX0, MY0 = merc(-122.665, 37.885); MX1, MY1 = merc(-122.355, 38.105)

with rasterio.open(src) as d:
    w = from_bounds(MX0, MY0, MX1, MY1, d.transform).round_offsets().round_lengths()
    full = d.read(1, window=w)
    T = d.window_transform(w)
left, top, res = T.c, T.f, T.a
H, W = full.shape[0] // 2, full.shape[1] // 2
b = full[:H * 2, :W * 2].reshape(H, 2, W, 2).max(axis=(1, 3))   # 2×2 max: ~2.4 m Mercator ≈ 1.9 m on the ground
r2 = res * 2
def px(lat, lon):  # lat/lon → column, row in the 2× grid
    mx, my = merc(lon, lat); return ((mx - left) / r2, (top - my) / r2)

# --- building footprints, rasterised on the 2× grid ---
bd = json.load(open(bpath)); q = bd['q']
img = Image.new('L', (W, H), 0); dr = ImageDraw.Draw(img)
for r in bd['b']:
    x = z = 0; ring = []
    for k in range(2, len(r), 2):
        x += r[k] * q; z += r[k + 1] * q
        ring.append(px(LAT0 - z / KZ, LON0 + x / KX))
    if len(ring) > 2: dr.polygon(ring, fill=1)
bmask = np.asarray(img, dtype=bool)
# "in town": within ~60 m of a building, on 16-cell (~30 m) blocks
G = 16
cb = bmask[:H // G * G, :W // G * G].reshape(H // G, G, W // G, G).any(axis=(1, 3))
town = np.repeat(np.repeat(maximum_filter(cb, size=5), G, axis=0), G, axis=1)
town = np.pad(town, ((0, H - town.shape[0]), (0, W - town.shape[1])))

# --- treetops ---
rng = np.random.default_rng(2058)
s = gaussian_filter(b.astype(np.float32), 1.0) + rng.random(b.shape, dtype=np.float32) * 0.01  # noise breaks plateau ties
peak = (s == maximum_filter(s, size=7)) & (b >= 4) & ~bmask
rows, cols = np.nonzero(peak)
keep_town = town[rows, cols]
# hills: tallest treetop per 13-cell (~25 m) block
C = 13
hr, hc = rows[~keep_town], cols[~keep_town]
hh = b[hr, hc]
key = (hr // C) * (W // C + 1) + hc // C
order = np.lexsort((-hh, key))
first = np.ones(len(order), dtype=bool); first[1:] = key[order][1:] != key[order][:-1]
sel = order[first]
R_ = np.concatenate([rows[keep_town], hr[sel]]); C_ = np.concatenate([cols[keep_town], hc[sel]])
print('treetops', len(rows), 'in town', int(keep_town.sum()), 'hills kept', len(sel))

xs, zs, hs = [], [], []
for rr, cc in zip(R_, C_):
    lon, lat = unmerc(left + (cc + 0.5) * r2, top - (rr + 0.5) * r2)
    xs.append(round((lon - LON0) * KX)); zs.append(round(-(lat - LAT0) * KZ)); hs.append(min(255, int(b[rr, cc]) * 4))
n = len(xs)
with open(out, 'wb') as f:
    f.write(struct.pack('<I', n))
    f.write(np.array(xs, dtype='<i2').tobytes()); f.write(np.array(zs, dtype='<i2').tobytes()); f.write(np.array(hs, dtype=np.uint8).tobytes())
print(n, 'trees written')

# --- tree cover around neighbourhood labels (full-resolution pixels, 450 m radius) ---
HOODS = [['Downtown', 37.9748, -122.5300], ['Gerstle Park', 37.9665, -122.5370], ['Sun Valley', 37.9860, -122.5480],
         ['Dominican', 37.9830, -122.5165], ['Terra Linda', 38.0070, -122.5520], ['Lucas Valley', 38.0290, -122.5620],
         ['Santa Venetia', 38.0020, -122.5150], ['Peacock Gap', 37.9930, -122.4700], ['Glenwood', 37.9850, -122.4790],
         ['Loch Lomond', 37.9715, -122.4850], ['Canal', 37.9630, -122.5000], ['Bret Harte', 37.9555, -122.5215],
         ['Northgate', 38.0060, -122.5420], ['Francisco Blvd E', 37.9590, -122.5030]]
for name, lat, lon in HOODS:
    mx, my = merc(lon, lat); c0, r0 = (mx - left) / res, (top - my) / res
    rad = 450 / math.cos(lat * D2R) / res
    rr, cc = np.ogrid[int(r0 - rad):int(r0 + rad), int(c0 - rad):int(c0 + rad)]
    m = (rr - r0) ** 2 + (cc - c0) ** 2 <= rad * rad
    sub = full[int(r0 - rad):int(r0 + rad), int(c0 - rad):int(c0 + rad)]
    print(f'{name:18s} {100 * (sub[m] >= 3).mean():5.1f}% tree cover')
