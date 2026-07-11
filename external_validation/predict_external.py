# -*- coding: utf-8 -*-
"""predict_external.py — disciplined external-validation protocol.

Reads external_validation/datasets_template.csv and, for every dataset with
complete numbers, runs:
  MODE A (defect transfer): use the Wu-calibrated f(d) as-is. This tests
    whether ANOTHER lab's defect population matches Wu's — expected to fail
    at the level of process variance; reported as an applicability-domain
    probe, NOT as model falsification.
  MODE B (one-point re-anchor): fit only the absolute scale f_res to the
    dataset's designated 'anchor' row (same protocol as the paper's Stage B),
    then predict every 'predict' row with zero further fitting. This is the
    transferable-structure test (tortuosity + ladder + lag), and is the
    metric that should enter the External Validation section.
Outputs per-row log-decade errors and per-dataset MAE/RMSE (decades).
Rows with extraction_status != OK are skipped and listed as TODO.
"""
import csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RG, TREF = 8.314, 298.15
K2 = lambda W: W/1e3/86400; G2 = lambda J: J*1e3*86400

# calibrated baseline (paper Table 1); organic properties per system below
CAL = dict(f0=0.399, dc=1.07, f_res=3.99e-8, fc0=1.62e-6, dcrit=43.9, r=50e-9)
ORG = {  # organic interlayer transport properties per material system
  "parylene": dict(D=5e-13, Ea=40e3, S=1.5),
  "pV3D3":    dict(D=None, Ea=40e3, S=None),   # TODO: primary-source values
}

def arr(x, Ea, T): return x*math.exp(-Ea/RG*(1/T-1/TREF))

def f_of_d(d, f_res):
    c = CAL
    return c["f0"]*math.exp(-d/c["dc"]) + f_res + c["fc0"]*max(0.,(d-c["dcrit"])/c["dcrit"])**2

def steady_wvtr(n, d_in, d_org, T, RH, f_res, org="parylene", R_sub=0.0):
    o = ORG[org]
    if o["D"] is None: raise ValueError(f"organic properties missing for {org}")
    P_org = arr(o["D"], o["Ea"], T)*o["S"]
    P_lat = arr(1e-21, 60e3, T)*0.10
    fv = f_of_d(d_in, f_res)
    P_in = P_lat + fv*P_org
    s = CAL["r"]*math.sqrt(math.pi/max(fv,1e-30))
    t2 = 1 + s**2*math.log(max(s/CAL["r"], math.e))/(2*math.pi*(d_org*1e-9)**2) if d_org>0 else 1.0
    R = R_sub + n*d_in*1e-9/P_in + max(n-1,0)*(d_org*1e-9)*t2/P_org + (d_org*1e-9)/P_org
    return G2(RH/100.0 / R)

def reanchor_f_res(row):
    """MODE B: solve f_res so the anchor row is matched exactly."""
    st = json.loads(row["structure_json"]); T = float(row["T_C"])+273.15
    target = float(row["measured_wvtr_gm2day"])
    lo, hi = 1e-11, 1e-4
    for _ in range(80):
        mid = math.sqrt(lo*hi)
        w = steady_wvtr(st["n"], st["d_in_nm"], st.get("d_org_nm",0), T,
                        float(row["RH_pct"]), mid, row.get("org","parylene"))
        (lo, hi) = (mid, hi) if w < target else (lo, mid)
    return math.sqrt(lo*hi)

def main():
    path = os.path.join(HERE, "datasets_template.csv")
    rows = list(csv.DictReader(open(path)))
    by_src = {}
    for r in rows: by_src.setdefault(r["source_id"], []).append(r)
    for src, rs in by_src.items():
        ok = [r for r in rs if r["extraction_status"] == "OK"
              and r["measured_wvtr_gm2day"] not in ("", "null")]
        todo = [r["row_id"] for r in rs if r not in ok]
        print(f"\n== {src}: {len(ok)} usable rows; TODO extraction: {todo}")
        if not ok: continue
        anchors = [r for r in ok if r["role"] == "anchor"]
        f_res = reanchor_f_res(anchors[0]) if anchors else CAL["f_res"]
        mode = "B (re-anchored)" if anchors else "A (defect transfer)"
        errs = []
        for r in ok:
            if r["role"] == "anchor": continue
            st = json.loads(r["structure_json"]); T = float(r["T_C"])+273.15
            pred = steady_wvtr(st["n"], st["d_in_nm"], st.get("d_org_nm",0), T,
                               float(r["RH_pct"]), f_res, r.get("org","parylene"))
            e = math.log10(pred/float(r["measured_wvtr_gm2day"]))
            errs.append(e)
            cflag = "" if r["censor"] == "=" else f" [censor {r['censor']}]"
            print(f"  {r['row_id']} mode {mode}: pred {pred:.2e} vs meas "
                  f"{float(r['measured_wvtr_gm2day']):.2e} -> {e:+.2f} dec{cflag}")
        if errs:
            print(f"  MAE {np.mean(np.abs(errs)):.2f} dec | RMSE {np.sqrt(np.mean(np.square(errs))):.2f} dec")

if __name__ == "__main__":
    main()
