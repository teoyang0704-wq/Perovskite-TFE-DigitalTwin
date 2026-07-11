"""Quickstart: evaluate one TFE design and print its six objectives (~1 s).
Run from repository root:  python examples/quickstart.py"""
import sys; sys.path.insert(0, "optimization"); sys.path.insert(0, "src")
from step6_optimize import evaluate
F, kpi = evaluate((100.0, 30.0, 3))   # (d_org_nm, d_inorg_nm, n_pairs)
print("Design 3 x (Al2O3 30 nm / organic 100 nm)")
for k, v in kpi.items(): print(f"  {k:12s}: {v}")
print("\nNote: T80_h uses the figure-baseline M_crit=0.01 g/m2 (DB row E009);")
print("the device-anchored floor is DB row E015 (0.97 g/m2): multiply the dose term x97; paper Sec. 3.4.")
