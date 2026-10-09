# Bake terrain.bin: USGS 3DEP elevations resampled onto the map's terrain grid, with land and bay
# split along the OpenStreetMap coastline. Output is Int16 decimetres, (NX+1)*(NZ+1) values, row j
# (north to south) then column i (west to east), matching HGT in index.html.
#
# Inputs (fetch once):
#   dem.tif    3DEP exportImage with square pixels (the server reshapes non-square requests):
#              https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage
#              ?bbox=-122.67,37.88,-122.35,38.11&bboxSR=4326&imageSR=4326&size=1600,1150&format=tiff
#              &pixelType=F32&interpolation=RSP_BilinearInterpolation&f=image
#   coast.json Overpass: [out:json];way["natural"="coastline"](37.885,-122.665,38.105,-122.355);out geom;
# Usage: python3 tools/terrain.py dem.tif coast.json terrain.bin
import sys, json, math
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree

NX, NZ = 620, 554
LONMIN, LONMAX, LATMIN, LATMAX = -122.665, -122.355, 37.885, 38.105
LAT0 = 37.983
KX, KZ = 111320 * math.cos(LAT0 * math.pi / 180), 110950
DATUM = 0.9      # NAVD88 → approximate local mean sea level (m)
LAND_MIN = 1.7   # the model keeps dry land above the bay plane; flooding is handled by the zone model

dem_path, coast_path, out_path = sys.argv[1:4]
im = Image.open(dem_path)
dem = np.array(im).astype('f8')
sx, sy, _ = im.tag_v2[33550]; _, _, _, lon0, lat0, _ = im.tag_v2[33922]
dem[~np.isfinite(dem) | (np.abs(dem) > 2000)] = np.nan  # open-water nodata comes back as garbage

lon = LONMIN + (LONMAX - LONMIN) * np.arange(NX + 1) / NX
lat = LATMAX - (LATMAX - LATMIN) * np.arange(NZ + 1) / NZ
LON, LAT = np.meshgrid(lon, lat)
fx = (LON - lon0) / sx - 0.5; fy = (lat0 - LAT) / sy - 0.5  # pixel-centre coordinates
ix = np.clip(np.floor(fx).astype(int), 0, dem.shape[1] - 2); iy = np.clip(np.floor(fy).astype(int), 0, dem.shape[0] - 2)
tx = np.clip(fx - ix, 0, 1); ty = np.clip(fy - iy, 0, 1)
a, b, c, d = dem[iy, ix], dem[iy, ix + 1], dem[iy + 1, ix], dem[iy + 1, ix + 1]
elev = (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty
elev = np.where(np.isnan(elev), np.nanmax(np.stack([a, b, c, d]), axis=0), elev)  # one missing corner → use the rest

# coastline: OSM ways run with land on the left. Rasterise them at ~8 m so they cut the map into regions,
# then call each region land or bay by majority vote of the side test on pixels right beside the line
# (a nearest-segment side test alone fails at sharp corners and far from the shore).
from scipy import ndimage
pts, side_a, side_b = [], [], []
for w in json.load(open(coast_path))['elements']:
    g = [(p['lon'] * KX, p['lat'] * KZ) for p in w.get('geometry', [])]
    for (ax, ay), (bx, by) in zip(g, g[1:]):
        n = max(1, int(math.hypot(bx - ax, by - ay) / 2))
        for k in range(n + 1):
            t = k / n; pts.append((ax + (bx - ax) * t, ay + (by - ay) * t)); side_a.append((ax, ay)); side_b.append((bx, by))
pts, side_a, side_b = np.array(pts), np.array(side_a), np.array(side_b)
tree = cKDTree(pts)
def side(q):
    dist, nn = tree.query(q); A, B = side_a[nn], side_b[nn]
    return dist, (B[:, 0] - A[:, 0]) * (q[:, 1] - A[:, 1]) - (B[:, 1] - A[:, 1]) * (q[:, 0] - A[:, 0]) > 0
PX = 8.0
# the raster must stay inside the Overpass bbox: past it the coastline has gaps and regions leak together
mx0, my1 = LONMIN * KX, LATMAX * KZ
W, H = int((LONMAX - LONMIN) * KX / PX), int((LATMAX - LATMIN) * KZ / PX)
wall = np.zeros((H, W), bool)
ci = ((pts[:, 0] - mx0) / PX).astype(int); cj = ((my1 - pts[:, 1]) / PX).astype(int)
ok = (ci >= 0) & (ci < W) & (cj >= 0) & (cj < H); wall[cj[ok], ci[ok]] = True
lab, nlab = ndimage.label(~wall, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
near = ndimage.binary_dilation(wall, iterations=3) & ~wall
jj, ii = np.where(near)
_, s_near = side(np.stack([mx0 + (ii + 0.5) * PX, my1 - (jj + 0.5) * PX], axis=1))
votes = np.zeros(nlab + 1); cnt = np.zeros(nlab + 1)
np.add.at(votes, lab[jj, ii], s_near.astype(float)); np.add.at(cnt, lab[jj, ii], 1)
is_land = votes / np.maximum(cnt, 1) > 0.5
q = np.stack([LON.ravel() * KX, LAT.ravel() * KZ], axis=1)
dist, s_q = side(q)
gi = np.clip(((q[:, 0] - mx0) / PX).astype(int), 0, W - 1); gj = np.clip(((my1 - q[:, 1]) / PX).astype(int), 0, H - 1)
L = lab[gj, gi]
land = np.where(L > 0, is_land[L], s_q).reshape(elev.shape)
dist = dist.reshape(elev.shape)

e = np.where(land, np.maximum(np.nan_to_num(elev - DATUM, nan=LAND_MIN), LAND_MIN), -1.5 - np.minimum(dist, 3000) / 600)
out = np.round(e * 10).astype('<i2')
out.tofile(out_path)
print('grid', out.shape, 'land', int(land.sum()), 'max m', round(float(e.max()), 1))
