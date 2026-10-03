"""
CST-305 Project 2 - Runge-Kutta-Fehlberg (RKF45) for ODEs
=========================================================
Programmer(s): <your name(s) here>
Packages used: math, time, numpy, matplotlib, scipy (optional, reference solver only)

Assigned problem (the defaults):  y' = y / (e^x - 1),  x0 = 1,  y0 = 5,  h = 0.02
Exact solution for that ODE (separable):  y(x) = y0 * (1 - e^-x) / (1 - e^-x0)
The user may instead enter any ODE y' = f(x, y), x0, y0, h and number of steps.
For other ODEs no closed form is known to the program, so a high-accuracy SciPy
reference solution (DOP853, rtol=1e-13) is used for the error comparison.

Approach to implementation
--------------------------
1. Read the ODE y' = f(x, y), x0, y0, h and #steps from the user (with defaults).
2. Advance the solution with the Runge-Kutta-Fehlberg 4(5) embedded pair using a
   fixed step h. Each step computes six slopes (k1..k6) and two estimates:
   a 4th-order value and a 5th-order value. Their difference |y5 - y4| is the
   local error estimate; the 5th-order value is carried forward.
3. Repeat for 1000 and 2000 steps and for two variations of the problem,
   counting steps and function evaluations and timing each run.
4. Solve the same problem with a "traditional" method (Euler) and with the exact
   formula, and compare against RKF45 to show the error.
5. At the end, plot (matplotlib): (1) the Runge-Kutta solution, (2) the solved
   (exact) ODE, and (3) both overlaid.
"""

import math
import time

import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Fehlberg 4(5) coefficients
# --------------------------------------------------------------------------
A = [0, 1/4, 3/8, 12/13, 1, 1/2]                       # nodes (x offsets)
B = [[],                                               # stage coefficients
     [1/4],
     [3/32, 9/32],
     [1932/2197, -7200/2197, 7296/2197],
     [439/216, -8, 3680/513, -845/4104],
     [-8/27, 2, -3544/2565, 1859/4104, -11/40]]
C4 = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0]         # 4th-order weights
C5 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]  # 5th-order weights


# --------------------------------------------------------------------------
# Core numerical routines
# --------------------------------------------------------------------------
def make_f(expr):
    """Convert a string such as 'y/(exp(x)-1)' into a function f(x, y)."""
    names = {k: getattr(math, k) for k in dir(math) if not k.startswith("_")}
    code = compile(expr, "<ode>", "eval")
    return lambda x, y: eval(code, {"__builtins__": {}}, {**names, "x": x, "y": y})


def rkf45_step(f, x, y, h):
    """One RKF45 step: returns (4th-order value, 5th-order value, slopes)."""
    k = []
    for i in range(6):
        yi = y + h * sum(B[i][j] * k[j] for j in range(i))
        k.append(f(x + A[i] * h, yi))
    y4 = y + h * sum(c * s for c, s in zip(C4, k))
    y5 = y + h * sum(c * s for c, s in zip(C5, k))
    return y4, y5, k


def rkf45_fixed(f, x0, y0, h, n):
    """Take n fixed RKF45 steps. Returns xs, ys, error estimates, #f-evaluations."""
    xs, ys, errs = [x0], [y0], [0.0]
    x, y = x0, y0
    for _ in range(n):
        y4, y5, _ = rkf45_step(f, x, y, h)
        errs.append(abs(y5 - y4))          # local error estimate
        x, y = x + h, y5                   # propagate the 5th-order value
        xs.append(x)
        ys.append(y)
    return xs, ys, errs, 6 * n


def euler(f, x0, y0, h, n):
    """Traditional first-order Euler method (for comparison)."""
    xs, ys = [x0], [y0]
    x, y = x0, y0
    for _ in range(n):
        y += h * f(x, y)
        x += h
        xs.append(x)
        ys.append(y)
    return xs, ys


def exact(x, x0, y0):
    """Exact solution of y' = y/(e^x - 1) with y(x0) = y0."""
    return y0 * (1 - math.exp(-x)) / (1 - math.exp(-x0))


# --------------------------------------------------------------------------
# Input helpers
# --------------------------------------------------------------------------
DEFAULT_EXPR = "y/(exp(x)-1)"


def ask_float(prompt, default, positive=False):
    """Ask for a number; Enter accepts the default."""
    while True:
        text = input(f"{prompt} [default {default}]: ").strip()
        try:
            val = float(text) if text else float(default)
            if positive and val <= 0:
                print("  Please enter a value greater than 0.")
                continue
            return val
        except ValueError:
            print("  Not a valid number, try again.")


def ask_steps(default=1000, lo=1000, hi=2000):
    """Number of solutions: the assignment asks for 1000 to 2000."""
    while True:
        text = input(f"Number of steps ({lo}-{hi}) [default {default}]: ").strip()
        try:
            n = int(text) if text else default
            if lo <= n <= hi:
                return n
            print(f"  Please enter a whole number from {lo} to {hi}.")
        except ValueError:
            print("  Not a valid whole number, try again.")


def ask_ode(x0_hint=1.0, y0_hint=5.0):
    """Ask for f(x,y) until it parses and can be evaluated."""
    while True:
        expr = input(f"Enter f(x,y) for y' = f(x,y) [default {DEFAULT_EXPR}]: ").strip() or DEFAULT_EXPR
        try:
            make_f(expr)(x0_hint, y0_hint)          # syntax / name check
            return expr
        except ZeroDivisionError:
            return expr                              # may be fine at other x0, checked later
        except Exception as e:
            print(f"  Could not use that expression ({e}). Use x, y and math functions "
                  f"such as exp, sin, log, sqrt, e.g. -2*y + x")


def reference_values(expr, f, x0, y0, xs):
    """Values to compare against. Returns (list or None, label)."""
    if expr.replace(" ", "") == DEFAULT_EXPR:
        return [exact(x, x0, y0) for x in xs], "Exact solution"
    try:
        from scipy.integrate import solve_ivp
        sol = solve_ivp(lambda x, y: f(x, y[0]), (xs[0], xs[-1]), [y0], method="DOP853",
                        t_eval=xs, rtol=1e-13, atol=1e-14)
        if sol.success and len(sol.t) == len(xs):
            return list(sol.y[0]), "Reference solution (SciPy DOP853)"
    except ImportError:
        pass
    return None, "Reference solution unavailable"


# --------------------------------------------------------------------------
# Reporting helpers
# --------------------------------------------------------------------------
def print_table(f, expr, x0, y0, h, n=5):
    """Print the first n RKF45 steps next to the reference values."""
    xs, ys, errs, _ = rkf45_fixed(f, x0, y0, h, n)
    ref, label = reference_values(expr, f, x0, y0, xs)
    print(f"\nFirst {n} RKF45 steps (h = {h}); comparison: {label}")
    print(f"{'n':>2} {'x_n':>8} {'y_n (RKF45)':>16} {'y_n (ref)':>16} "
          f"{'|error|':>11} {'|y5-y4|':>11}")
    for i in range(n + 1):
        if ref:
            r, e = f"{ref[i]:16.10f}", f"{abs(ys[i] - ref[i]):11.3e}"
        else:
            r, e = f"{'n/a':>16}", f"{'n/a':>11}"
        print(f"{i:>2} {xs[i]:8.4f} {ys[i]:16.10f} {r} {e} {errs[i]:11.3e}")
    _, _, k = rkf45_step(f, x0, y0, h)
    print("\nStage slopes for step 1:  " +
          "  ".join(f"k{i}={v:.10f}" for i, v in enumerate(k, 1)))


def run_case(label, expr, x0, y0, h, n):
    """Solve with RKF45 and Euler, time both, compare with the reference."""
    f = make_f(expr)

    t0 = time.perf_counter()
    xs, ys, errs, nfev = rkf45_fixed(f, x0, y0, h, n)
    t_rk = time.perf_counter() - t0

    t0 = time.perf_counter()
    _, ye = euler(f, x0, y0, h, n)
    t_eu = time.perf_counter() - t0

    ref, ref_label = reference_values(expr, f, x0, y0, xs)

    print(f"\n=== {label}: y' = {expr}, x0={x0}, y0={y0}, h={h}, steps={n} ===")
    print(f"RKF45 : steps={n}, f-evals={nfev}, time={t_rk * 1e3:.3f} ms, "
          f"final (x, y) = ({xs[-1]:.4f}, {ys[-1]:.10f})")
    print(f"Euler : steps={n}, f-evals={n}, time={t_eu * 1e3:.3f} ms, "
          f"final y = {ye[-1]:.10f}")
    if ref:
        err_rk = max(abs(a - b) for a, b in zip(ys, ref))
        err_eu = max(abs(a - b) for a, b in zip(ye, ref))
        print(f"{ref_label}: final y = {ref[-1]:.10f}")
        print(f"Max |error| over all steps:  RKF45 = {err_rk:.3e}   Euler = {err_eu:.3e}")
    else:
        print("No reference solution available (install SciPy) - error not computed.")
    print(f"Max local error estimate |y5-y4| = {max(errs):.3e}")
    return xs, ys, ref, ref_label


# --------------------------------------------------------------------------
# Plotting (called at the very end)
# --------------------------------------------------------------------------
def show_plots(xs, ys, ref, ref_label, title):
    """Three figures: Runge-Kutta plot, solved ODE, and both overlaid."""
    xs, ys = np.array(xs), np.array(ys)

    # Plot 1: Runge-Kutta (RKF45) numerical solution
    plt.figure("Runge-Kutta Plot", figsize=(8, 5))
    plt.plot(xs, ys, color="tab:blue", lw=2, label="RKF45 numerical solution")
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.title("Runge-Kutta (RKF45) Plot\n" + title, fontsize=10)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()

    if ref is not None:
        ref = np.array(ref)

        # Plot 2: Solved ODE (exact or reference solution)
        plt.figure("Solved ODE", figsize=(8, 5))
        plt.plot(xs, ref, color="tab:green", lw=2, label=ref_label)
        plt.xlabel("x")
        plt.ylabel("y(x)")
        plt.title("Solved ODE\n" + title, fontsize=10)
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()

        # Plot 3: both overlaid (thick reference line, thinner dashed RKF45 on top)
        plt.figure("Overlay", figsize=(8, 5))
        plt.plot(xs, ref, color="tab:green", lw=5, alpha=0.6, label=ref_label)
        plt.plot(xs, ys, color="tab:blue", lw=1.5, ls="--", label="RKF45 solution")
        plt.xlabel("x")
        plt.ylabel("y(x)")
        plt.title("RKF45 vs Solved ODE (Overlay)\n" + title, fontsize=10)
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
    else:
        print("Plots 2 and 3 skipped: no reference solution available.")

    plt.show()   # opens all figure windows


# --------------------------------------------------------------------------
# Main program
# --------------------------------------------------------------------------
if __name__ == "__main__":
    # 1. Read the ODE and initial conditions (Enter accepts the assigned defaults)
    expr = ask_ode()
    x0 = ask_float("Initial x0", 1.0)
    y0 = ask_float("Initial y0", 5.0)
    h = ask_float("Step size h", 0.02, positive=True)
    n_user = ask_steps(1000)
    f = make_f(expr)

    try:
        # 2. Table for the first five steps (compare with the manual calculation)
        print_table(f, expr, x0, y0, h, 5)

        # 3. Required runs: the chosen number of steps, then 2000 steps
        xs_main, ys_main, ref_main, ref_label = run_case(
            f"Original, {n_user} steps", expr, x0, y0, h, n_user)
        if n_user != 2000:
            run_case("Original, 2000 steps", expr, x0, y0, h, 2000)

        # 4. Two variations of the problem
        run_case(f"Variation 1 (smaller step h = {h / 4:g})", expr, x0, y0, h / 4, 2000)
        run_case(f"Variation 2 (shifted start x0 = {x0 + 1:g}, y0 = {0.6 * y0:g})",
                 expr, x0 + 1, 0.6 * y0, h, 2000)

        # 5. Reference: SciPy's adaptive RK45 (optional)
        try:
            from scipy.integrate import solve_ivp
            t0 = time.perf_counter()
            sol = solve_ivp(lambda x, y: f(x, y[0]), (x0, x0 + n_user * h), [y0],
                            method="RK45", rtol=1e-10, atol=1e-12)
            t = time.perf_counter() - t0
            print(f"\nSciPy RK45 (adaptive): steps={len(sol.t) - 1}, f-evals={sol.nfev}, "
                  f"time={t * 1e3:.3f} ms")
        except ImportError:
            print("\n(SciPy not installed - skipping reference solver)")

    except (ZeroDivisionError, OverflowError, ValueError) as e:
        raise SystemExit(f"\nThe ODE cannot be evaluated over this interval ({e}). "
                         "Check for a singularity (e.g. x = 0 for y/(exp(x)-1)) or overflow.")

    # 6. Plots at the end: Runge-Kutta, solved ODE, and overlay
    show_plots(xs_main, ys_main, ref_main, ref_label,
               f"y' = {expr},  y({x0:g}) = {y0:g},  h = {h:g},  {n_user} steps")