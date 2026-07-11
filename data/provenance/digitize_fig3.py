import numpy as np
from PIL import Image, ImageDraw

img = Image.open('/mnt/user-data/uploads/FullSizeRender.jpeg').convert('RGB')
W, H = img.size
print("image size:", W, H)
a = np.asarray(img).astype(int)
sx, sy = W/878.0, H/815.0   # my visual reads were on 878x815 render

# --- axis calibration: 5 y-tick label centers (visual, refined by LSQ) ---
ticks_y_vis = np.array([148, 240, 352, 460, 570])   # 1e2,1e0,1e-2,1e-4,1e-6
ticks_log   = np.array([2, 0, -2, -4, -6])
ys = ticks_y_vis * sy
A = np.vstack([ys, np.ones_like(ys)]).T
(b, c), res, *_ = np.linalg.lstsq(A, ticks_log, rcond=None)
pred = A @ np.array([b, c])
print("tick fit residuals (decades):", np.round(ticks_log - pred, 3))
def y2log(y): return b*y + c

# --- masks ---
R, G, B = a[...,0], a[...,1], a[...,2]
red   = (R>140) & (G<110) & (B<110)
black = (R<70) & (G<70) & (B<70)
# plot interior; exclude legend box and axis frame
X, Y = np.meshgrid(np.arange(W), np.arange(H))
inplot = (X > 162*sx) & (X < 745*sx) & (Y > 118*sy) & (Y < 565*sy)
legend = (X > 340*sx) & (Y < 225*sy)
red   &= inplot & ~legend
black &= inplot & ~legend

stations_vis = [170, 355, 540, 718]
half = 9*sx
out = {}
for name, mask in (("high(black)", black), ("poor(red)", red)):
    vals = []
    for i, xv in enumerate(stations_vis):
        x0 = xv*sx
        m = mask & (np.abs(X-x0) < half)
        n = m.sum()
        if n < 30:
            vals.append((i, None, n)); continue
        yc = Y[m].mean()
        vals.append((i, 10**y2log(yc), n))
    out[name] = vals
    print(name, [(d, f"{v:.2e}" if v else None, n) for d, v, n in vals])

# --- overlay for visual verification ---
ov = img.copy(); dr = ImageDraw.Draw(ov)
for lg in ticks_log:                       # horizontal decade lines from fit
    yv = (lg - c)/b
    dr.line([(150*sx, yv), (760*sx, yv)], fill=(0,160,255), width=2)
for name, vals, col in (("h", out["high(black)"], (0,200,0)), ("p", out["poor(red)"], (255,0,255))):
    for d, v, n in vals:
        if v is None: continue
        x0 = stations_vis[d]*sx; yv = (np.log10(v) - c)/b
        dr.line([(x0-12, yv), (x0+12, yv)], fill=col, width=3)
        dr.line([(x0, yv-12), (x0, yv+12)], fill=col, width=3)
ov.save('/home/claude/fig3_overlay.png')
print("overlay saved")
