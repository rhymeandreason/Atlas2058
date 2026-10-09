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
import sys, os, math
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coast import land_mask

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

land, dist = land_mask(coast_path, LON, LAT)

e = np.where(land, np.maximum(np.nan_to_num(elev - DATUM, nan=LAND_MIN), LAND_MIN), -1.5 - np.minimum(dist, 3000) / 600)
out = np.round(e * 10).astype('<i2')
out.tofile(out_path)
print('grid', out.shape, 'land', int(land.sum()), 'max m', round(float(e.max()), 1))
