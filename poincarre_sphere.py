import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Create figure
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection="3d")

# Define axes limits (from -1 to 1 for all Stokes parameters)
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

# Set axis labels with LaTeX formatting and ensure S_3 is properly positioned
ax.set_xlabel(r"$\mathbf{S_2}$", fontsize=18, fontweight="bold", labelpad=15)
ax.set_ylabel(r"$\mathbf{S_1}$", fontsize=18, fontweight="bold", labelpad=15)
ax.set_zlabel(r"$\mathbf{S_3}$", fontsize=18, fontweight="bold", labelpad=15)

# Generate smooth sphere surface
u = np.linspace(0, 2 * np.pi, 150)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))

# Apply soft aesthetic color with transparency
ax.plot_surface(x, y, z, color="lightpink", alpha=0.25, edgecolor="none")

# Draw thick coordinate axes with labeled ticks
axis_linewidth = 4
ax.plot([-1, 1], [0, 0], [0, 0], "k-", lw=axis_linewidth)  # S1 axis
ax.plot([0, 0], [-1, 1], [0, 0], "k-", lw=axis_linewidth)  # S2 axis
ax.plot([0, 0], [0, 0], [-1, 1], "k-", lw=axis_linewidth)  # S3 axis

# Set tick values for better readability
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
ax.set_yticks([-1, -0.5, 0, 0.5, 1])
ax.set_zticks([-1, -0.5, 0, 0.5, 1])

# Ensure S_3 is properly visible by adjusting view angle
ax.view_init(elev=20, azim=30)

# Adjust aspect ratio to make the sphere appear correctly
ax.set_box_aspect([1, 1, 1])

# Basis states at extreme points, with more spacing for better readability
label_offset = 0.3
points = {
    r"$\mathbf{|H\rangle}$": (1 + label_offset, 0, 0),
    r"$\mathbf{|V\rangle}$": (-1 - label_offset, 0, 0),
    r"$\mathbf{|D\rangle}$": (0, 1 + label_offset, 0),
    r"$\mathbf{|A\rangle}$": (0, -1 - label_offset, 0),
    r"$\mathbf{|R\rangle}$": (0, 0, 1 + label_offset),
    r"$\mathbf{|L\rangle}$": (0, 0, -1 - label_offset),
}

# Plot labels with high contrast and increased spacing
for label, (px, py, pz) in points.items():
    ax.text(
        px,
        py,
        pz,
        label,
        fontsize=20,
        fontweight="bold",
        color="black",
        ha="center",
        va="center",
    )

# Three manuscript trajectories on the Poincare sphere.
# theta is the polarization-state angle used in the paper.
num_points = 300
theta = np.linspace(0, np.pi, num_points)

# HDVA: rotation around S3, path in the S1-S2 plane.
hdva_x = np.cos(2 * theta)
hdva_y = np.sin(2 * theta)
hdva_z = np.zeros_like(theta)
ax.plot(hdva_x, hdva_y, hdva_z, color="red", linewidth=3, label="HDVA")

# HRVL: rotation around S2, path in the S1-S3 plane.
hrvl_x = np.cos(2 * theta)
hrvl_y = np.zeros_like(theta)
hrvl_z = np.sin(2 * theta)
ax.plot(hrvl_x, hrvl_y, hrvl_z, color="blue", linewidth=3, label="HRVL")

# DRAL: rotation around S1, path in the S2-S3 plane.
dral_x = np.zeros_like(theta)
dral_y = np.sin(2 * theta)
dral_z = -np.cos(2 * theta)
ax.plot(dral_x, dral_y, dral_z, color="green", linewidth=3, label="DRAL")

# Add arrows indicating increasing theta.
arrow_idx = num_points // 4
for path_x, path_y, path_z, color in [
    (hdva_x, hdva_y, hdva_z, "red"),
    (hrvl_x, hrvl_y, hrvl_z, "blue"),
    (dral_x, dral_y, dral_z, "green"),
]:
    ax.quiver(
        path_x[arrow_idx],
        path_y[arrow_idx],
        path_z[arrow_idx],
        path_x[arrow_idx + 1] - path_x[arrow_idx],
        path_y[arrow_idx + 1] - path_y[arrow_idx],
        path_z[arrow_idx + 1] - path_z[arrow_idx],
        color=color,
        linewidth=3,
        arrow_length_ratio=0.2,
    )

ax.legend(loc="upper left", fontsize=16)

# Save figure with 300 dpi resolution
plt.savefig("poincare_sphere_HDVA_HRVL_DRAL.png", dpi=300)

# Show figure
plt.show()
