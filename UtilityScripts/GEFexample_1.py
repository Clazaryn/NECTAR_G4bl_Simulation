""" EXAMPLE TO TEST THE GEF lmd python reader """

from gef_lmd_reader import GEFLmdReader
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

# ---------- STYLE ----------
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = ["Computer Modern"]
mpl.rcParams["mathtext.fontset"] = "cm"
mpl.rcParams["font.size"] = 14
mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["xtick.labelsize"] = 13
mpl.rcParams["ytick.labelsize"] = 13
mpl.rcParams["figure.figsize"] = (8, 5)
mpl.rcParams["text.usetex"] = True

# ---------- READ DATA ----------

r = GEFLmdReader(92, 239, 2.0, 1, base_dir=".")
r.open()

KE_light = []
N = r.n_events()

for i in range(1, N + 1):
    Tl, Th = r.GetKEpost_pair(i)
    KE_light.append(Tl)

    # --- progress bar ---
    if i % 1000 == 0 or i == N:
        percent = 100.0 * i / N
        bar_len = 30
        filled = int(bar_len * i / N)
        bar = "=" * filled + "-" * (bar_len - filled)
        sys.stdout.write("\r[{0}] {1:.1f}% ({2}/{3})".format(bar, percent, i, N))
        sys.stdout.flush()

print("\nDone.")
r.close()


# ---------- PLOT ----------
plt.figure()

plt.hist(
    KE_light,
    bins=100,
    histtype="step",     # <-- line only
    linewidth=1.8,
    color="#1f77b4"
)

plt.xlabel(r"$T_{\mathrm{light}}^{\mathrm{post}}$ (MeV)")
plt.ylabel("Counts")
plt.title("$T_{\mathrm{light}}^{\mathrm{post}}$ distribution")

plt.grid(alpha=0.3)
plt.tight_layout()

plt.show()
