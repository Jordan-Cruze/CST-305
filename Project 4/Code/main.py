"""
CST-305 - Systems of Ordinary Differential Equations

Programmer:
    Jordan Cruze

Packages:
    NumPy
    Matplotlib

Description:
    This program models data and I/O flow between processors using
    systems of first-order linear ordinary differential equations.

    Part 1:
        Three processors (A, B, C) with data flowing between them
        and from Processor C to a network.

    Part 2:
        Two processors (A, B) with data flowing between them.

    The processor transfer rates are represented as proportions
    of the current processor contents.

    Part 1:
        A -> B = 2% per second
        B -> A = 6% per second
        B -> C = 5% per second
        C -> B = 1% per second
        C -> Network = 4% per second

    Part 2:
        A -> B = 2% per second
        B -> A = 3% per second

    Approach:
    Each processor's rate of change is (data received) - (data sent),
    with every rate divided by the 100 MB memory. This builds x' = A x.
    The program computes the eigenvalues of A, solves x(t) = e^(At) x(0)
    with the matrix exponential, checks it with odeint, and plots the result.
"""

import argparse

import numpy as np
from scipy.integrate import odeint
from scipy.linalg import expm
import matplotlib.pyplot as plt

UNITS, TIME = "MB", "s"

# Defaults from the assignment. rates[(i, j)] = MB/s sent from processor i to j.
PRESETS = {
    "1": dict(name="Part 1 - Three processors with network drain", n=3, scale=100,
              rates={(0, 1): 6, (1, 0): 2, (1, 2): 5, (2, 1): 1},
              drain=[0, 0, 4], x0=[100, 100, 100], t_start=0, t_end=300),
    "2": dict(name="Part 2 - Closed two-processor system", n=2, scale=100,
              rates={(0, 1): 2, (1, 0): 3},
              drain=[0, 0], x0=[1, -1], t_start=-100, t_end=100),
    "3": dict(name="Custom system", n=3, scale=100, rates={},
              drain=None, x0=None, t_start=0, t_end=300),
}


def ask(prompt, default, cast=float):
    """Prompt the user; Enter accepts the default."""
    while True:
        raw = input("%s [%s]: " % (prompt, default)).strip()
        if raw == "":
            return default
        try:
            return cast(raw)
        except ValueError:
            print("  Please enter a valid number.")


def read_model():
    """Interactively read the model; returns a dict describing the ODE system."""
    print("Models: 1 = Part 1 (3 processors + network), 2 = Part 2 (closed, 2 processors),"
          " 3 = custom")
    choice = str(ask("Choose model", 1, int))
    p = PRESETS[choice]
    n = ask("Number of processors", p["n"], int) if p["name"] == "Custom system" else p["n"]
    labels = [chr(ord("A") + i) for i in range(n)]
    scale = ask("Memory per processor / scaling denominator (MB)", p["scale"])

    print("\nTransfer rates (%s/%s) from one processor to another:" % (UNITS, TIME))
    rates = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                rates[i, j] = ask("  %s -> %s" % (labels[i], labels[j]), p["rates"].get((i, j), 0))

    print("\nRate sent from each processor out to the network (%s/%s):" % (UNITS, TIME))
    pd = p["drain"] or [0] * n
    drain = [ask("  %s -> network" % labels[i], pd[i]) for i in range(n)]

    print("\nDegradation (bit flips, charge leakage, media decay), applied to every processor:")
    d = ask("  Degradation rate d (1/%s)" % TIME, 0.0)

    print("\nInitial data + I/O in each processor at t = 0 (%s):" % UNITS)
    px = p["x0"] or [100] * n
    x0 = [ask("  x_%s(0)" % labels[i], px[i]) for i in range(n)]
    t_start = ask("\nPlot start time (%s)" % TIME, p["t_start"])
    t_end = ask("Plot end time (%s)" % TIME, p["t_end"])
    return dict(name=p["name"], labels=labels, scale=scale, rates=rates,
                drain=drain, d=d, x0=np.array(x0, float), t_start=t_start, t_end=t_end,
                two_figures=(choice == "2"))   # instructor shows both graphs only for Part 2


def build_matrix(m):
    """A from: rate of change = input rate - output rate (each scaled by memory)."""
    n, s = len(m["labels"]), m["scale"]
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                A[j, i] += m["rates"][i, j] / s    # input to j from i
                A[i, i] -= m["rates"][i, j] / s    # output from i to j
        A[i, i] -= m["drain"][i] / s               # output to the network
        A[i, i] -= m["d"]                          # degradation loss
    return A


def print_equations(labels, A):
    print("\nSystem of ODEs  x' = A x  (f = 0):")
    for i, li in enumerate(labels):
        terms = " ".join("%+.4f x_%s" % (A[i, j], lj) for j, lj in enumerate(labels) if A[i, j] != 0)
        print("  dx_%s/dt = %s" % (li, terms or "0"))


def solve_exact(A, x0, t):
    """x(t) = e^{At} x0 (valid because f = 0)."""
    return np.array([expm(A * ti) @ x0 for ti in t]).T


def solve_numeric(A, x0, t):
    """Numerical check with odeint. x0 is the state at t = 0, so integrate forward
    over t >= 0 and backward over t <= 0 (odeint expects rhs(x, t))."""
    rhs = lambda x, _t: A @ x
    out = np.zeros((len(x0), len(t)))
    for mask, order in ((t >= 0, 1), (t <= 0, -1)):
        tt = t[mask][::order]                      # start at t = 0 and move away from it
        if len(tt):
            out[:, mask] = odeint(rhs, x0, tt, rtol=1e-10, atol=1e-10).T[:, ::order]
    return out


def closed_form_labels(A, x0, labels):
    """If x(0) is an eigenvector of A (A x0 = lam x0) then x_i(t) = x0_i e^(lam t);
    return labels such as 'e^-0.05t' / '-e^-0.05t'. Otherwise use the processor names."""
    lam = (x0 @ (A @ x0)) / (x0 @ x0)
    if lam != 0 and np.allclose(A @ x0, lam * x0, atol=1e-9):
        out = []
        for c in x0:
            coef = "" if c == 1 else "-" if c == -1 else "%g " % c
            out.append("%se^%gt" % (coef, round(lam, 6)))
        return out
    return ["x%d(t): Processor %s" % (i + 1, l) for i, l in enumerate(labels)]


def plot_solution(m, A, t, X):
    """Part 2 shows two figures, styled like the instructor's examples (Part 1 shows only Figure 1):
       Figure 1 (left example)  - titled, gridded, curves labeled by their formula
       Figure 2 (right example) - 'x(t) = e^(At)*c', Solution 1 / Solution 2."""
    # Figure 1
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    for i, lab in enumerate(closed_form_labels(A, m["x0"], m["labels"])):
        ax1.plot(t, X[i], label=lab)
    ax1.set_title("Degradation of data integrity - ODE Solutions by Matrix Method")
    ax1.set_xlabel("t")
    ax1.set_ylabel("x")
    if m["two_figures"]:
        ax1.set_ylim(-40, 40)      # same y-axis range as the instructor's example
    ax1.grid(True)
    ax1.legend(loc="upper right")
    fig1.tight_layout()

    if not m["two_figures"]:                       # Part 1 / custom: one graph only
        return

    # Figure 2
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    for i in range(len(X)):
        ax2.plot(t, X[i], label="Solution %d" % (i + 1))
    ax2.set_title("x(t) = e^(At)*c")
    ax2.set_xlabel("t")
    ax2.set_ylabel("x(t)")
    ax2.legend(loc="upper right")
    fig2.tight_layout()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-show", action="store_true", help="do not open the plot window (testing)")
    args = ap.parse_args()

    m = read_model()
    A = build_matrix(m)
    print_equations(m["labels"], A)
    print("\nA =\n", np.round(A, 5))

    # Eigenvalues -> stability and time constants
    lam = np.linalg.eigvals(A)
    print("Eigenvalues (1/%s):" % TIME, np.round(np.sort_complex(lam)[::-1], 5))
    nz = np.abs(lam.real)[np.abs(lam.real) > 1e-12]
    if nz.size:
        print("Slowest time constant: %.1f %s" % (1 / nz.min(), TIME))
    print("Stability:", "asymptotically stable (all data eventually lost)"
          if np.all(lam.real < -1e-12) else "stable, data conserved (a zero eigenvalue)")

    # Solve exactly, verify numerically
    t = np.unique(np.append(np.linspace(m["t_start"], m["t_end"], 600), 0.0))
    t = t[(t >= m["t_start"]) & (t <= m["t_end"])]
    X = solve_exact(A, m["x0"], t)
    err = np.max(np.abs(X - solve_numeric(A, m["x0"], t)))
    print("Max |exact - odeint| = %.2e" % err)
    print("Final state (%s):" % UNITS, np.round(X[:, -1], 3), "| total:", round(X[:, -1].sum(), 3))

    plot_solution(m, A, t, X)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()