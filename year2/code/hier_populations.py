# -*- coding: utf-8 -*-
"""Two populations: every published film, versus films intended for use.

The hierarchical fit of `hierarchical.py` treats four films as draws from one
population.  One of them is not like the others: Wu's poor-quality film was
deposited deliberately at high base pressure to demonstrate that particulates
destroy a barrier.  Pooling it with films meant for service inflates the
between-film spread and makes the predictive distribution for a new laboratory
pessimistic.

This script repeats the inference on the three service-quality films and
compares the two populations, which is the honest way to report the spread:
one number for "any film in the literature", another for "a film you would
actually deploy".
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import hierarchical as H

CLOSURE_EPS = H.CLOSURE_EPS


def run(exclude=(), nw=52, nsteps=6000, burn=2000, seed=11):
    films = [f for f in ["Wu high-q", "Wu poor", "Lee", "Carcia"] if f not in exclude]
    H.FILMS, H.NF = films, len(films)
    H.SINGLE = [s for s in H.SINGLE_ALL if s[0] in films]
    H.LADDER = [l for l in H.LADDER_ALL if l[0] in films]
    ndim = 9 + H.NF
    rng = np.random.default_rng(seed)
    start = np.array([np.log10(0.372), np.log10(1.07), np.log10(40.0), 44.0,
                      -7.4, 0.4, *([-7.4] * H.NF), np.log10(2e4),
                      np.log10(25.0), np.log10(5e-3)])
    p0 = start + 0.02 * rng.standard_normal((nw, ndim))
    p0[:, 5] = np.abs(p0[:, 5])
    chain, arate = H.stretch_sample(H.log_posterior, p0, nsteps, rng)
    flat = chain[burn:].reshape(-1, ndim)
    return films, flat, arate


def report(tag, films, flat):
    print(f"\n=== {tag} ===")
    lo, md, hi = np.percentile(flat[:, 5], [5, 50, 95])
    print(f"  between-film spread   {md:.2f}  [{lo:.2f}, {hi:.2f}] decades")
    for i, f in enumerate(films):
        v = 10 ** np.percentile(flat[:, 6 + i], [5, 50, 95])
        print(f"  f_res {f:12s} {v[1]:.2e}  [{v[0]:.1e}, {v[2]:.1e}]")
    rn = np.random.default_rng(3).standard_normal(len(flat))
    lf_new = flat[:, 4] + flat[:, 5] * rn
    v = 10 ** np.percentile(lf_new, [5, 50, 95])
    print(f"  new film (predictive) {v[1]:.2e}  [{v[0]:.1e}, {v[2]:.1e}]")
    d_close = 10 ** flat[:, 1] * np.log(10 ** flat[:, 0] / (CLOSURE_EPS * 10 ** lf_new))
    w = flat[:, 3] - d_close
    q = np.percentile(d_close, [5, 50, 95])
    print(f"  window lower edge     {q[1]:.1f}  [{q[0]:.1f}, {q[2]:.1f}] nm")
    q2 = np.percentile(flat[:, 3], [5, 50, 95])
    print(f"  window upper edge     {q2[1]:.1f}  [{q2[0]:.1f}, {q2[2]:.1f}] nm")
    q3 = np.percentile(w, [5, 50, 95])
    print(f"  window width          {q3[1]:.1f}  [{q3[0]:.1f}, {q3[2]:.1f}] nm"
          f"   (empty in {100*np.mean(w <= 0):.1f} % of draws)")
    return lf_new, d_close


if __name__ == "__main__":
    H.SINGLE_ALL, H.LADDER_ALL = list(H.SINGLE), list(H.LADDER)

    films_a, flat_a, ra = run()
    lf_a, dc_a = report("all published films", films_a, flat_a)
    films_b, flat_b, rb = run(exclude=("Wu poor",))
    lf_b, dc_b = report("service-quality films only (Wu poor excluded)", films_b, flat_b)
    print(f"\n  acceptance rates {ra:.2f} / {rb:.2f}")

    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8), dpi=200)
    bins = np.linspace(-11, -3.5, 90)
    ax[0].hist(lf_a, bins=bins, density=True, alpha=.45, color="#8A5AA8",
               label="new film | all published")
    ax[0].hist(lf_b, bins=bins, density=True, alpha=.55, color="#1B3A5C",
               label="new film | service quality")
    for i, f in enumerate(films_a):
        ax[0].axvline(np.median(flat_a[:, 6 + i]), lw=1.1, ls=":", color="#D55E00")
        ax[0].text(np.median(flat_a[:, 6 + i]), ax[0].get_ylim()[1] * (0.92 - .12 * i),
                   f, fontsize=6, rotation=90, va="top", ha="right", color="#D55E00")
    ax[0].set(xlabel="log$_{10}$ defect area fraction $f_{res}$", ylabel="density")
    ax[0].legend(fontsize=7, frameon=False)
    ax[0].set_title("How much do films differ?", fontsize=9)

    ax[1].hist(dc_a, bins=np.linspace(8, 40, 80), density=True, alpha=.45,
               color="#8A5AA8", label="all published")
    ax[1].hist(dc_b, bins=np.linspace(8, 40, 80), density=True, alpha=.55,
               color="#1B3A5C", label="service quality")
    ax[1].axvline(22.5, color="k", ls="--", lw=1.4)
    ax[1].text(22.8, ax[1].get_ylim()[1] * .9, "Year-1 point value", fontsize=6.5)
    ax[1].set(xlabel="lower edge of the design window (nm)", ylabel="density")
    ax[1].legend(fontsize=7, frameon=False)
    ax[1].set_title("What a new laboratory should assume", fontsize=9)
    fig.tight_layout(); fig.savefig("fig_hierarchical.png")
    print("\nsaved fig_hierarchical.png")
