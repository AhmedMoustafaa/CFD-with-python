import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
nx, ny = 20, 20
lx, ly = 1.0, 1.0
dx, dy = lx / (nx - 1), ly / (ny - 1)
Re = 100
nu = 1.0 / Re
dt = 0.001  # Time step
nt = 100  # number of time steps

# Initialize fields
u = np.zeros((ny, nx))  # x-velocity
v = np.zeros((ny, nx))  # y-velocity
p = np.zeros((ny, nx))  # Pressure
b = np.zeros((ny, nx))  # RHS of the Poisson equation


# Helper functions
def build_rhs(u, v, rho=1.0):
    """Build the right-hand side of the Poisson equation."""
    b[1:-1, 1:-1] = (
            rho * (1 / dt) * ((u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dx)
                              + (v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dy))
            - ((u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dx)) ** 2
            - 2 * ((u[2:, 1:-1] - u[:-2, 1:-1]) / (2 * dy) *
                   (v[1:-1, 2:] - v[1:-1, :-2]) / (2 * dx))
            - ((v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dy)) ** 2
    )


def pressure_poisson(p, b):
    """Solve the pressure Poisson equation."""
    for _ in range(20):  # Fewer iterations
        pn = p.copy()
        p[1:-1, 1:-1] = ((pn[1:-1, 2:] + pn[1:-1, :-2]) * dy ** 2 +
                         (pn[2:, 1:-1] + pn[:-2, 1:-1]) * dx ** 2) / (2 * (dx ** 2 + dy ** 2)) \
                        - b[1:-1, 1:-1] * dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2))
        p[:, 0] = p[:, 1]  # dp/dx = 0 at left wall
        p[:, -1] = p[:, -2]  # dp/dx = 0 at right wall
        p[0, :] = p[1, :]  # dp/dy = 0 at bottom wall
        p[-1, :] = 0  # Reference pressure at top wall


def velocity_update(u, v, p):
    """Update velocity fields."""
    un = u.copy()
    vn = v.copy()
    u[1:-1, 1:-1] = (un[1:-1, 1:-1]
                     - un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, :-2])
                     - vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[:-2, 1:-1])
                     - dt / (2 * dx) * (p[1:-1, 2:] - p[1:-1, :-2])
                     + nu * dt * ((un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, :-2]) / dx ** 2
                                  + (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[:-2, 1:-1]) / dy ** 2))
    v[1:-1, 1:-1] = (vn[1:-1, 1:-1]
                     - un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, :-2])
                     - vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[:-2, 1:-1])
                     - dt / (2 * dy) * (p[2:, 1:-1] - p[:-2, 1:-1])
                     + nu * dt * ((vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, :-2]) / dx ** 2
                                  + (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[:-2, 1:-1]) / dy ** 2))
    # Boundary conditions
    u[0, :] = 0
    u[:, 0] = 0
    u[:, -1] = 0
    u[-1, :] = 1  # Top wall moving with velocity 1
    v[0, :] = 0
    v[-1, :] = 0
    v[:, 0] = 0
    v[:, -1] = 0


# Main loop
fig, ax = plt.subplots()
im = ax.imshow(np.zeros_like(u), origin='lower', cmap='viridis', extent=[0, lx, 0, ly])


def update(frame):
    global u, v, p
    build_rhs(u, v)
    pressure_poisson(p, b)
    velocity_update(u, v, p)
    im.set_array(np.sqrt(u ** 2 + v ** 2))  # Plot velocity magnitude

    # Streamlines
    ax.streamplot(np.linspace(0, lx, nx), np.linspace(0, ly, ny), u, v, color='white', linewidth=0.5, density=1.5)

    return [im]


if __name__ == "__main__":
    ani = FuncAnimation(fig, update, frames=nt, interval=5, blit=True)
    plt.colorbar(im)
    plt.title("Lid-Driven Cavity Flow with Streamlines")

    # Save the animation as a GIF
    ani.save('flow_animation.gif', writer='pillow', fps=30)

    plt.show()
