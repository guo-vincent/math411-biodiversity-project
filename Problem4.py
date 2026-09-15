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

fig2, ax2 = plt.subplots(figsize=(8, 5))
n_line = np.logspace(0, 4, 400)
for S in (2, 3, 5, 10, 50):
    ax2.plot(n_line, simpson_even(S, n_line), label=f"S = {S}")
    ax2.axhline(1 - 1 / S, color="grey", lw=0.6, ls=":")

ax2.set_xscale("log")
ax2.set_xlabel("n  (individuals per species)")
ax2.set_ylabel("1 - D")
ax2.set_ylim(0, 1.02)
ax2.set_title("Each curve saturates at its ceiling 1 - 1/S (dotted)")
ax2.legend()
fig2.tight_layout()

plt.show()