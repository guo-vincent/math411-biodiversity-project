import numpy as np
import matplotlib.pyplot as plt


def simpson_index(counts) -> float:
    """1 - D for an arbitrary vector of species counts."""
    counts = np.asarray(counts, dtype=np.float64)
    N = counts.sum()
    if N < 2:
        return np.nan
    D = np.sum(counts * (counts - 1.0)) / (N * (N - 1.0))
    return 1.0 - D


def simpson_even(S, n):
    """
    1 - D for a perfectly even community of S species with n individuals each.
    Closed form, vectorized over S and n:   D = (n - 1) / (S*n - 1)
    """
    S = np.asarray(S, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    N = S * n
    with np.errstate(invalid="ignore", divide="ignore"):
        out = 1.0 - (n - 1.0) / (N - 1.0)
    return np.where(N >= 2, out, np.nan)


# ------------------------ First plot (combined) ------------------------

for _S, _n in [(1, 10), (3, 7), (25, 4), (100, 100)]:
    assert np.isclose(simpson_even(_S, _n), simpson_index(np.full(_S, _n)))

S_vals = np.unique(np.logspace(0, 3, 150).round().astype(int))  # 1 .. 1000
n_vals = np.unique(np.logspace(0, 3, 150).round().astype(int))  # 1 .. 1000
SS, NN = np.meshgrid(S_vals, n_vals, indexing="ij")
Z = simpson_even(SS, NN)

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(projection="3d")
surf = ax.plot_surface(
    np.log10(SS), np.log10(NN), Z,
    cmap="viridis", vmin=0, vmax=1,
    rstride=1, cstride=1, linewidth=0, antialiased=True,
)

# ---- Red N = Sn line ----
N_target = 1000

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(projection="3d", computed_zorder=False)
surf = ax.plot_surface(
    np.log10(SS), np.log10(NN), Z,
    cmap="viridis", vmin=0, vmax=1,
    rstride=1, cstride=1, linewidth=0, antialiased=True,
    zorder=1,
)

logN = np.log10(N_target)
lo, hi = max(0.0, logN - 3.0), min(3.0, logN)
if lo >= hi:
    raise ValueError(f"N = {N_target:,} doesn't cross the plotted range (N must be between 1 and 1,000,000).")

logS_line = np.linspace(lo, hi, 400)
S_curve = 10 ** logS_line
n_curve = N_target / S_curve
Z_curve = simpson_even(S_curve, n_curve)

ax.plot(
    logS_line, np.log10(n_curve), Z_curve,
    color="red", lw=3.5, zorder=10,
    label=f"N = S x n = {N_target:,}",
)

# ---- Text label sitting on the line ----
label_frac = 0.2 
i = int(label_frac * (len(logS_line) - 1))
ax.text(
    logS_line[i], np.log10(n_curve[i]), Z_curve[i] + 0.04,
    f"N = {N_target:,}",
    color="red", fontsize=12, fontweight="bold", ha="center", va="bottom",
    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="red", alpha=0.9),
    zorder=11,
)

for axis, setter, labeller in (
    (ax.xaxis, ax.set_xticks, ax.set_xticklabels),
    (ax.yaxis, ax.set_yticks, ax.set_yticklabels),
):
    setter([0, 1, 2, 3])
    labeller(["1", "10", "100", "1000"])

ax.set_xlabel("S  (number of species)")
ax.set_ylabel("n  (individuals per species)")
ax.set_zlabel("Simpson's Index of Diversity, 1 - D")
ax.set_zlim(0, 1)
ax.set_title("1 - D for an even community (N = S x n)")
ax.view_init(elev=22, azim=-125)
fig.colorbar(surf, shrink=0.6, pad=0.12, label="1 - D")
fig.subplots_adjust(left=0.02, right=0.92, bottom=0.05, top=0.95)

# ------------------------ Second plot (Selected S values) ------------------------

fig2, ax2 = plt.subplots(figsize=(8, 5))
n_line = np.logspace(0, 4, 400)
for S in (2, 3, 5, 10, 50):
    ax2.plot(n_line, simpson_even(S, n_line), label=f"S = {S}")
    ax2.axhline(1 - 1 / S, color="grey", lw=0.6, ls=":")

ax2.set_xscale("log")
ax2.set_xlabel("n  (individuals per species)")
ax2.set_ylabel("SID (1 - D)")
ax2.set_ylim(0, 1.02)
ax2.set_title("Fixed S: Each curve saturates at its floor 1 - 1/S (evenness)")
ax2.legend()
fig2.tight_layout()

# ------------------------ Third plot (Selected N values) ------------------------
fig3, ax3 = plt.subplots(figsize=(8, 5))

# integer S values, log-spaced, 1 .. 10,000
S_line = np.unique(np.logspace(0, 4, 400).round().astype(int))

for n in (2, 3, 5, 10, 50):
    ax3.plot(S_line, simpson_even(S_line, n), label=f"n = {n}")

ax3.axhline(1.0, color="grey", lw=0.8, ls=":")   # the limit as S -> infinity

ax3.set_xscale("log")
ax3.set_xlabel("S  (number of species)")
ax3.set_ylabel("SID (1 - D)")
ax3.set_ylim(0, 1.02)
ax3.set_title("Fixed n: SID approaches 1 as S grows (richness)")
ax3.legend()
fig3.tight_layout()

plt.show()