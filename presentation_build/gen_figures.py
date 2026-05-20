"""Generate KU Leuven-styled schematic figures for the MTOP presentation."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, FancyArrowPatch

TEAL = "#1D8DB0"
SLATE = "#2F4D5D"
SKY = "#52BDEC"
GOLD = "#E7B037"
GRAY = "#AEBFC9"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
OUT = "presentation_build"


def grid(ax, x0, y0, w, h, nx, ny, color, lw, z=2):
    for i in range(nx + 1):
        x = x0 + i * w / nx
        ax.plot([x, x], [y0, y0 + h], color=color, lw=lw,
                solid_capstyle="round", zorder=z)
    for j in range(ny + 1):
        y = y0 + j * h / ny
        ax.plot([x0, x0 + w], [y, y], color=color, lw=lw,
                solid_capstyle="round", zorder=z)


# ----------------------------------------------------------------------
# Figure 1: concrete MTOP discretisation (classical vs MTOP + 5x5 zoom)
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(10.2, 4.8), dpi=200)
fig.patch.set_facecolor("white")

axc = fig.add_axes([0.045, 0.580, 0.605, 0.300])
axm = fig.add_axes([0.045, 0.120, 0.605, 0.345])
axz = fig.add_axes([0.720, 0.135, 0.255, 0.345])
for ax in (axc, axm, axz):
    ax.set_aspect("equal")
    ax.axis("off")

# classical panel
axc.set_xlim(-0.06, 3.06)
axc.set_ylim(-0.30, 1.18)
grid(axc, 0, 0, 3, 1, 24, 8, SLATE, 0.8)
axc.add_patch(Rectangle((0, 0), 3, 1, fill=False, ec=SLATE, lw=2.4))
axc.text(1.5, 1.30, "Classical: analysis mesh = density mesh",
         ha="center", va="bottom", fontsize=12.5, fontweight="bold", color=SLATE)
axc.text(1.5, -0.16, "600 × 200 elements", ha="center", va="top",
         fontsize=9.5, color=SLATE)

# mtop panel
axm.set_xlim(-0.06, 3.06)
axm.set_ylim(-0.30, 1.18)
grid(axm, 0, 0, 3, 1, 30, 10, GRAY, 0.5, z=1)
ex0, ey0 = 2.5, 0.5
axm.add_patch(Rectangle((ex0, ey0), 0.5, 0.5, facecolor=SKY, alpha=0.40,
                        edgecolor="none", zorder=1.5))
grid(axm, 0, 0, 3, 1, 6, 2, TEAL, 2.3, z=2)
axm.add_patch(Rectangle((0, 0), 3, 1, fill=False, ec=SLATE, lw=2.4))
axm.text(1.5, 1.30, "MTOP: coarse analysis mesh + fine density mesh",
         ha="center", va="bottom", fontsize=12.5, fontweight="bold", color=TEAL)
axm.text(1.5, -0.16,
         "analysis 120 × 40          density 600 × 200",
         ha="center", va="top", fontsize=9.5, color=SLATE)

# zoom panel
axz.set_xlim(-0.13, 1.13)
axz.set_ylim(-0.34, 1.20)
axz.add_patch(Rectangle((0, 0), 1, 1, facecolor=SKY, alpha=0.16, edgecolor="none"))
grid(axz, 0, 0, 1, 1, 5, 5, TEAL, 1.4)
axz.add_patch(Rectangle((0, 0), 1, 1, fill=False, ec=TEAL, lw=2.7))
for i in range(5):
    for j in range(5):
        axz.plot((i + 0.5) / 5, (j + 0.5) / 5, "o", ms=4.6, color=GOLD,
                 mec=SLATE, mew=0.6, zorder=3)
axz.text(0.5, 1.18, "one analysis element", ha="center", va="bottom",
         fontsize=10.5, fontweight="bold", color=SLATE)
axz.text(0.5, -0.18, "5 × 5 = 25 density cells\nmidpoint integration points",
         ha="center", va="top", fontsize=9, color=SLATE)

# connectors element -> zoom
for (cx, cy), (zx, zy) in [((3.0, 1.0), (0, 1)), ((3.0, 0.5), (0, 0))]:
    fig.add_artist(ConnectionPatch(
        xyA=(cx, cy), coordsA=axm.transData,
        xyB=(zx, zy), coordsB=axz.transData,
        color=SLATE, lw=1.0, ls=(0, (4, 2)), zorder=5))

fig.savefig(f"{OUT}/gen_concept_meshes.png", dpi=200, facecolor="white")
plt.close(fig)

# ----------------------------------------------------------------------
# Figure 2: abstract idea teaser (fine design mesh vs coarse analysis mesh)
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(6.7, 5.6), dpi=200)
fig.patch.set_facecolor("white")

axt = fig.add_axes([0.07, 0.585, 0.86, 0.30])
axb = fig.add_axes([0.07, 0.115, 0.86, 0.30])
for ax in (axt, axb):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-0.06, 3.06)
    ax.set_ylim(-0.34, 1.30)

# top: fine design mesh
grid(axt, 0, 0, 3, 1, 30, 10, GRAY, 0.6)
axt.add_patch(Rectangle((0, 0), 3, 1, fill=False, ec=SLATE, lw=2.4))
axt.text(1.5, 1.34, "DESIGN  ·  fine density mesh", ha="center", va="bottom",
         fontsize=12.5, fontweight="bold", color=SLATE)
axt.text(1.5, -0.18, "600 × 200 density cells  →  full design resolution",
         ha="center", va="top", fontsize=9.5, color=SLATE)

# bottom: coarse analysis mesh
grid(axb, 0, 0, 3, 1, 6, 2, TEAL, 2.3)
axb.add_patch(Rectangle((0, 0), 3, 1, fill=False, ec=SLATE, lw=2.4))
axb.text(1.5, 1.34, "ANALYSIS  ·  coarse finite element mesh",
         ha="center", va="bottom", fontsize=12.5, fontweight="bold", color=TEAL)
axb.text(1.5, -0.18, "120 × 40 elements  →  smaller system  K u = f",
         ha="center", va="top", fontsize=9.5, color=SLATE)

# arrow between
fig.add_artist(FancyArrowPatch(
    (0.5, 0.545), (0.5, 0.470), transform=fig.transFigure,
    arrowstyle="-|>", mutation_scale=26, lw=2.6, color=GOLD))
fig.text(0.545, 0.508, "MTOP decouples\nthe two meshes", ha="left", va="center",
         fontsize=10.5, fontweight="bold", color=SLATE)

fig.savefig(f"{OUT}/gen_idea.png", dpi=200, facecolor="white")
plt.close(fig)

print("wrote gen_concept_meshes.png and gen_idea.png")
