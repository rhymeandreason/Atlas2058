# Land/bay split along the OpenStreetMap coastline, shared by terrain.py and ground.py.
# OSM coastline ways run with land on the left. They are rasterised at ~8 m so they cut the map into regions,
# then each region is called land or bay by majority vote of the side test on pixels right beside the line
# (a nearest-segment side test alone fails at sharp corners and far from the shore).
import json, math
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

LONMIN, LONMAX, LATMIN, LATMAX = -122.665, -122.355, 37.885, 38.105
LAT0 = 37.983
KX, KZ = 111320 * math.cos(LAT0 * math.pi / 180), 110950

def land_mask(coast_path, LON, LAT):
    """LON/LAT: arrays of query points (degrees). Returns (land bool, distance to coast in m), same shape."""
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
    land = np.where(L > 0, is_land[L], s_q)
    return land.reshape(LON.shape), dist.reshape(LON.shape)
