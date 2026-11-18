
"""
Lorenz System Simulation
Author: Jordan Cruze

Packages:
- NumPy
- Matplotlib
- SciPy

Implementation:
Numerically solves the Lorenz system and generates
a 3D trajectory and individual X, Y, and Z plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Lorenz system parameters
SIGMA = 10
BETA = 1.5

# Initial conditions from the professor's example
INITIAL_CONDITIONS = [11.8, 4.4, 2.4]

# Simulation settings
TIME_SPAN = (0, 100)
TIME_POINTS = np.linspace(0, 100, 10000)


# Lorenz differential equations
def lorenz_system(t, state, r):
    x, y, z = state

    dx = SIGMA * (y - x)
    dy = x * (r - z) - y
    dz = x * y - BETA * z

    return [dx, dy, dz]


# Run numerical simulation
def simulate(r):
    solution = solve_ivp(
        lorenz_system,
        TIME_SPAN,
        INITIAL_CONDITIONS,
        args=(r,),
        t_eval=TIME_POINTS,
        method="RK45",
        rtol=1e-9,
        atol=1e-11,
        max_step=0.01
    )

    if not solution.success:
        raise RuntimeError(solution.message)

    return (
        solution.t,
        solution.y[0],
        solution.y[1],
        solution.y[2]
    )


# Display plots individually
def plot_results(r, t, x, y, z):

    # 3D trajectory
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(x, y, z, linewidth=0.6)

    ax.set_title(f"Lorenz - r: {r}")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")

    plt.tight_layout()
    plt.show()
    plt.close()

    # X versus time
    plt.figure(figsize=(9, 5))
    plt.plot(t, x)
    plt.title(f"X(t) - r: {r}")
    plt.xlabel("t - Time")
    plt.ylabel("X")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.close()

    # Y versus time
    plt.figure(figsize=(9, 5))
    plt.plot(t, y)
    plt.title(f"Y(t) - r: {r}")
    plt.xlabel("t - Time")
    plt.ylabel("Y")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.close()

    # Z versus time
    plt.figure(figsize=(9, 5))
    plt.plot(t, z)
    plt.title(f"Z(t) - r: {r}")
    plt.xlabel("t - Time")
    plt.ylabel("Z")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.close()


# Main program
def main():
    print("Lorenz System Simulation")
    print("------------------------")
    print("Enter any positive value for r.")
    print("Enter -1 to quit.")

    while True:
        try:
            r = float(input("\nEnter a value for r: "))

            if r == -1:
                break

            if r <= 0:
                print("Please enter a positive value.")
                continue

            print(f"\nRunning simulation for r = {r}...")

            t, x, y, z = simulate(r)

            plot_results(r, t, x, y, z)

            print("All four plots completed.")

        except ValueError:
            print("Please enter a valid number.")

        except RuntimeError as error:
            print("Simulation error:", error)


if __name__ == "__main__":
    main()
