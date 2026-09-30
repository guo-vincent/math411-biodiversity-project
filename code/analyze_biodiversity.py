"""Reproduce plots/data: python3 code/analyze_biodiversity.py.
MATH 411 Group 5; source details are in data/source_manifest.json.
Adds integer domains, feasible fixed-total points, exact checks, and exports.
"""

import csv, json, os, platform
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROWS = [
    ("Spotted Dusky Salamander", 7, 18), ("Imitator Salamander", 6, 3),
    ("Seal Salamander", 5, 15), ("Black-bellied Salamander", 7, 11),
    ("Desmognathus spp. Salamander", 4, 17),
    ("Blue Ridged Two-line Salamander", 1, 31), ("Spring Salamander", 2, 1),
    ("Northern Slimy Salamander", 0, 1), ("Santeetlah Salamander", 1, 0),
    ("Southern Red-backed Salamander", 2, 0),
]
def sid_even(S, n):
	"""Equal abundances; positive integer inputs; (1, 1) is undefined."""
	S, n = np.broadcast_arrays(np.asarray(S), np.asarray(n))
	for x in (S, n):
		if x.dtype.kind not in "iuf" or np.any(
			~np.isfinite(x) | (x < 1) | (x != np.floor(x))
		):
			raise ValueError("S and n must be positive integers")
	S, n = S.astype(float), n.astype(float)
	with np.errstate(divide="ignore", invalid="ignore"):
		return np.where(S * n >= 2, 1 - (n - 1) / (S * n - 1), np.nan)
def exact_summary(counts):
	"""Check the formula by counting unordered different-species pairs."""
	N = sum(counts)
	numerator = sum(x * (x - 1) for x in counts)
	D = Fraction(numerator, N * (N - 1))
	pair_sid = Fraction(sum(a * b for a, b in combinations(counts, 2)), comb(N, 2))
	assert 1 - D == pair_sid
	return dict(N=N, richness=sum(x > 0 for x in counts),
	same_pair_numerator=numerator, D=str(D), SID=str(1 - D),
	SID_decimal=float(1 - D), SID_rounded_4dp=f"{float(1-D):.4f}")
def main():
	for directory in ("plots", "data"):
		(ROOT / directory).mkdir(parents=True, exist_ok=True)
	streams = {name: exact_summary([row[col] for row in ROWS])
	for name, col in (("Lower Dorsey Stream", 1), ("Pig Pen Stream", 2))}
	for S in range(1, 8):
		for n in range(1, 8):
			if S * n > 1:
				expected = float(Fraction(n*(S-1), S*n-1))
				assert abs(float(sid_even(S, n)) - expected) < 1e-12
				assert exact_summary([n]*S)["SID"] == str(Fraction(n*(S-1), S*n-1))
	bad_inputs = (0, -1, 1.5, True, np.nan, np.inf, "2")
	for value in bad_inputs:
		try:
			sid_even(value, 2)
		except ValueError:
			continue
		raise AssertionError(f"Invalid input accepted: {value!r}")
	assert np.isnan(sid_even(1, 1))
	with (ROOT / "data/stream_counts.csv").open("w", newline="") as f:
		writer = csv.writer(f)
		writer.writerow(("category", "lower_dorsey", "pig_pen")); writer.writerows(ROWS)
	divisors = np.array([S for S in range(1, 1001) if 1000 % S == 0])
	abundances = 1000 // divisors
	for S, n in zip(divisors, abundances):
		expected = Fraction(1000-int(n), 999)
		assert S*n == 1000 and exact_summary([int(n)]*int(S))["SID"] == str(expected)
		assert abs(float(sid_even(S, n)) - float(expected)) < 1e-12
	with (ROOT / "data/fixed_total_points.csv").open("w", newline="") as f:
		writer = csv.writer(f); writer.writerow(("S", "n", "N", "SID_fraction", "SID"))
		writer.writerows((int(S), int(n), 1000, str(Fraction(1000-int(n), 999)),
		float(sid_even(S, n))) for S, n in zip(divisors, abundances))
	result = dict(streams=streams, four_species_100=exact_summary([100]*4),
	provenance={"source_manifest": "data/source_manifest.json",
	"commit": "fded9fa72c9b5b009eba0688004aeb02de8c90b9",
	"adapted_from": "Problem4.py", "counts": "solvedbetter.tex"},
	versions={"Python": platform.python_version(), "numpy": np.__version__,
	"matplotlib": matplotlib.__version__},
	verification={"exact_pair_checks": "passed", "equal_abundance_cases": 48,
	"invalid_input_checks": 7, "undefined_1_1": "NaN",
	"fixed_total_exact_checks": len(divisors), "randomness": "none"})
	(ROOT / "data/analysis.json").write_text(json.dumps(result, indent=2) + "\n")
	plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
	"axes.spines.right": False, "savefig.dpi": 320})
	x = np.arange(1, 1001)
	colors = ("#245d83", "#b57926", "#477454")
	fig, axes = plt.subplots(1, 2, figsize=(7.6, 2.95), sharey=True, layout="constrained")
	for S, color, marker in zip((2, 4, 10), colors, ("o", "s", "^")):
		axes[0].plot(x, sid_even(S, x), marker+"-", ms=1.7, lw=1, color=color, label=f"S = {S}")
		axes[0].axhline(1-1/S, color=color, ls="--", lw=0.9)
	for n, color, marker in zip((2, 5, 20), colors, ("o", "s", "^")):
		axes[1].plot(x, sid_even(x, n), marker+"-", ms=1.7, lw=1, color=color, label=f"n = {n}")
		axes[1].axhline(1, color="0.35", ls="--", lw=0.9)
	for ax, title, xlabel in zip(axes, ("Fixed species count S", "Fixed abundance n"),
	("n (individuals per species)", "S (species present)")):
		ax.set(xscale="log", xlim=(1, 1000), ylim=(0, 1.025), title=title, xlabel=xlabel)
		ax.grid(alpha=0.2); ax.legend(loc="lower right", frameon=False, fontsize=8, markerscale=2)
	axes[0].set_ylabel("SID (probability of different species)")
	fig.savefig(ROOT / "plots/sid_limits.png"); plt.close(fig)
	grid = np.unique(np.geomspace(1, 1000, 150).round().astype(int))
	SS, NN = np.meshgrid(grid, grid, indexing="ij")
	fig = plt.figure(figsize=(6.7, 4.6))
	ax = fig.add_subplot(projection="3d", computed_zorder=False)
	surface = ax.plot_surface(np.log10(SS), np.log10(NN), sid_even(SS, NN),
	cmap="cividis", vmin=0, vmax=1, alpha=0.82, zorder=1, rstride=1, cstride=1)
	ax.plot(np.log10(divisors), np.log10(abundances), sid_even(divisors, abundances),
	"o-", color="#b62828", ms=4, lw=1.5, zorder=10, label="N = 1000: integer pairs")
	ax.set(xlabel="S (species present)", ylabel="n (individuals per species)",
	zlim=(0, 1), title="Equal-abundance diversity")
	ax.text2D(-0.14, 0.56, "SID", transform=ax.transAxes, rotation=90, va="center")
	for setter, labeller in ((ax.set_xticks, ax.set_xticklabels),
	(ax.set_yticks, ax.set_yticklabels)):
		setter([0, 1, 2, 3]); labeller(["1", "10", "100", "1000"])
	ax.view_init(elev=23, azim=-125); ax.legend(loc="upper left", frameon=False, fontsize=8)
	fig.colorbar(surface, ax=ax, shrink=0.64, pad=0.1, label="SID")
	fig.subplots_adjust(left=0.05, right=0.92, bottom=0.09, top=0.94)
	fig.savefig(ROOT / "plots/sid_surface.png", bbox_inches="tight", pad_inches=0.1)
	plt.close(fig)

if __name__ == "__main__":
	main()