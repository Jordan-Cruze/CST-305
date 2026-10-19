"""
main.py
Programmers : <Team member names here>
Course      : <Course / Assignment name>
Packages    : sympy (symbolic Green's function, integration, ODE solving),
              numpy (numerics), scipy (numerical check with solve_ivp),
              matplotlib (static plots + animation)

Approach
--------
For a linear constant-coefficient IVP
        a*y'' + b*y' + c*y = f(t),   y(0) = 0,  y'(0) = 0,   t >= 0
1. The ODE is read as text (e.g. "y''+y=4") and parsed with sympy.
2. The homogeneous equation is solved symbolically (dsolve).
3. Green's function G(t,tau) is the solution of the homogeneous equation
   (in t) with G(tau,tau)=0 and dG/dt(tau,tau)=1/a.
4. The solution is  y(t) = Integral_0^t G(t,tau) f(tau) d tau.
5. The result is verified (substituted back into the ODE and compared with a
   numerical solve_ivp solution), plotted, and animated.

Usage
-----
  python green_ode.py                      # prompts you for the ODE(s), t_max, etc.
  python green_ode.py "y''+y=4"            # one ODE given as text
  python green_ode.py "y''+4*y=t" --tmax 15 --no-show   # prints results only, no plots
Units: t in seconds (s), y in meters (m) -- a generic forced oscillator
       (y'' + w^2 y = F(t)); change UNITS below to fit your context.
"""
import argparse
import warnings
import re
import numpy as np
import sympy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp

warnings.filterwarnings("ignore", message="Animation was deleted")  # only relevant when no window is shown
t, tau = sp.symbols("t tau", real=True)
y = sp.Function("y")
T_UNIT, Y_UNIT = "s", "m"          # terminology / units for the plots

# ----------------------------------------------------------------------------
def parse_ode(text):
    """Parse e.g. "y''+4*y=x" -> (a, b, c, f(t)) for a*y''+b*y'+c*y=f(t)."""
    if text.count("=") != 1:
        raise ValueError("the ODE must contain exactly one '=' sign")
    s = text.replace("^", "**").replace("y''", "@2@").replace("y'", "@1@")
    s = re.sub(r"\by\b", "@0@", s)                     # remaining plain y
    s = s.replace("@2@", "Dd2").replace("@1@", "Dd1").replace("@0@", "Dd0")
    lhs, rhs = s.split("=")
    D2, D1, D0 = sp.symbols("Dd2 Dd1 Dd0")
    # 'x' is accepted as the independent variable (as in the assignment)
    loc = {"Dd2": D2, "Dd1": D1, "Dd0": D0, "t": t, "x": t}
    expr = sp.expand(sp.sympify(lhs, locals=loc) - sp.sympify(rhs, locals=loc))
    a, b, c = (expr.coeff(D, 1) for D in (D2, D1, D0))
    f = -expr.subs({D2: 0, D1: 0, D0: 0})
    return sp.nsimplify(a), sp.nsimplify(b), sp.nsimplify(c), sp.simplify(f)

def homogeneous_solution(a, b, c, y0=0, v0=0):
    """General homogeneous solution and the one that meets given ICs."""
    ode = sp.Eq(a*y(t).diff(t, 2) + b*y(t).diff(t) + c*y(t), 0)
    general = sp.dsolve(ode, y(t))
    ics = sp.dsolve(ode, y(t), ics={y(0): y0, y(t).diff(t).subs(t, 0): v0})
    return general.rhs, ics.rhs

def greens_function(a, b, c):
    """G(t,tau): homogeneous solution with G(tau,tau)=0, G_t(tau,tau)=1/a."""
    ode = sp.Eq(a*y(t).diff(t, 2) + b*y(t).diff(t) + c*y(t), 0)
    sol = sp.dsolve(ode, y(t)).rhs
    C = sorted(sol.free_symbols - {t}, key=str)
    eqs = [sol.subs(t, tau), sp.diff(sol, t).subs(t, tau) - 1/a]
    cs = sp.solve(eqs, C, dict=True)[0]
    return sp.simplify(sol.subs(cs))

def particular_and_total(a, b, c, f):
    G = greens_function(a, b, c)
    y_total = sp.simplify(sp.integrate(sp.expand(G*f.subs(t, tau)), (tau, 0, t)))
    return G, y_total

def check(a, b, c, f, ysol):
    res = sp.simplify(a*ysol.diff(t, 2) + b*ysol.diff(t) + c*ysol - f)
    return res == 0, ysol.subs(t, 0), ysol.diff(t).subs(t, 0)

# ----------------------------------------------------------------------------
def run(ode_text, tmax, show, tag):
    a, b, c, f = parse_ode(ode_text)
    print("="*70)
    print(f"ODE: {ode_text}   ->  {a}*y'' + {b}*y' + {c}*y = {f}")
    G, ysol = particular_and_total(a, b, c, f)
    ok, y0, v0 = check(a, b, c, f, ysol)
    gen, hom_ic = homogeneous_solution(a, b, c)
    print("Green's function  G(t,tau) =", G)
    print("Homogeneous (general)      =", gen)
    print("Green's solution  y(t)     =", ysol)
    print(f"Verification: ODE residual = 0 -> {ok};  y(0)={y0}, y'(0)={v0}")

    # Split y = y_h + y_p. y_p = polynomial/forced part, y_h = homogeneous part
    # fixed by the ICs (shown in the plots as the 'homogeneous part').
    ys = sp.expand(ysol)
    yp = sum((term for term in sp.Add.make_args(ys)
              if not term.has(sp.sin, sp.cos, sp.exp)), sp.Integer(0))
    yh = sp.simplify(ys - yp)
    print("Particular part y_p        =", yp)
    print("Homogeneous part y_h       =", yh)

    fy = sp.lambdify(t, ysol, "numpy")
    fyp = sp.lambdify(t, yp, "numpy")
    fyh = sp.lambdify(t, yh, "numpy")
    ff = sp.lambdify(t, f, "numpy")
    fG = sp.lambdify((t, tau), G, "numpy")
    tt = np.linspace(0, tmax, 600)
    ones = np.ones_like(tt)

    # numerical check
    fn = sp.lambdify(t, f, "numpy")
    num = solve_ivp(lambda s, u: [u[1], (float(fn(s)) - float(b)*u[1] - float(c)*u[0])/float(a)],
                    [0, tmax], [0, 0], t_eval=tt, rtol=1e-9, atol=1e-11)
    err = np.max(np.abs(num.y[0] - fy(tt)*ones))
    print(f"Max |Green - numerical (solve_ivp)| = {err:.2e} {Y_UNIT}")

    # ---- Figure 1: homogeneous part (+ general homogeneous basis) ----------
    fig1, ax = plt.subplots(1, 2, figsize=(11, 4))
    C = sorted(gen.free_symbols - {t}, key=str)
    for k, Ck in enumerate(C):
        basis = sp.lambdify(t, gen.subs({Cj: (1 if Cj == Ck else 0) for Cj in C}), "numpy")
        ax[0].plot(tt, basis(tt)*ones, label=f"basis {k+1}: {gen.subs({Cj: (1 if Cj == Ck else 0) for Cj in C})}")
    ax[0].set_title("Homogeneous equation: basis solutions")
    ax[0].set_xlabel(f"t ({T_UNIT})"); ax[0].set_ylabel(f"y ({Y_UNIT})")
    ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
    ax[1].plot(tt, fyh(tt)*ones, "r")
    ax[1].set_title(f"Homogeneous part of the IVP solution\n$y_h$ = {sp.latex(yh) and yh}")
    ax[1].set_xlabel(f"t ({T_UNIT})"); ax[1].set_ylabel(f"y_h ({Y_UNIT})"); ax[1].grid(alpha=.3)
    fig1.suptitle(f"{ode_text},  y(0)=y'(0)=0")
    fig1.tight_layout()

    # ---- Figure 2: Green's function solution, particular, homogeneous ------
    fig2, ax2 = plt.subplots(figsize=(8, 4.5))
    ax2.plot(tt, fy(tt)*ones, "b", lw=2.2, label=f"Green's solution y(t) = {ysol}")
    ax2.plot(tt, fyp(tt)*ones, "g--", label=f"particular y_p = {yp}")
    ax2.plot(tt, fyh(tt)*ones, "r:", label=f"homogeneous y_h = {yh}")
    ax2.plot(tt[::25], num.y[0][::25], "ko", ms=3, label="numerical check (solve_ivp)")
    ax2.set_xlabel(f"time t ({T_UNIT})"); ax2.set_ylabel(f"displacement y ({Y_UNIT})")
    ax2.set_title(f"IVP: {ode_text}, y(0)=y'(0)=0"); ax2.grid(alpha=.3); ax2.legend(fontsize=8)
    fig2.tight_layout()

    # ---- Animation: Green's integral built up as tau sweeps 0 -> t ---------
    fig3, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7))
    a1.set_xlim(0, tmax); a1.set_xlabel(f"time t ({T_UNIT})"); a1.set_ylabel(f"y ({Y_UNIT})")
    ymin, ymax = np.min(fy(tt)*ones), np.max(fy(tt)*ones); pad = 0.1*(ymax-ymin or 1)
    a1.set_ylim(ymin-pad, ymax+pad); a1.grid(alpha=.3)
    a1.plot(tt, fy(tt)*ones, color="lightgray", label="full y(t)")
    line, = a1.plot([], [], "b", lw=2, label="y(t) built so far")
    dot, = a1.plot([], [], "ro")
    a1.legend(loc="upper left"); a1.set_title("Solution y(t) = ∫₀ᵗ G(t,τ) f(τ) dτ")
    a2.set_xlim(0, tmax); a2.set_xlabel(f"τ ({T_UNIT})"); a2.grid(alpha=.3)
    a2.set_ylabel("integrand G(t,τ)·f(τ)")
    integ, = a2.plot([], [], "g"); fill = [None]
    frames = np.linspace(0, tmax, 120)
    def update(t0):
        i = np.searchsorted(tt, t0)
        line.set_data(tt[:i+1], (fy(tt)*ones)[:i+1]); dot.set_data([t0], [fy(t0)])
        ts = np.linspace(0, t0, 200)
        vals = fG(t0, ts)*ones[:200]*ff(ts)
        integ.set_data(ts, vals)
        if fill[0]: fill[0].remove()
        fill[0] = a2.fill_between(ts, 0, vals, color="g", alpha=.25)
        m = max(1e-9, np.max(np.abs(vals))); a2.set_ylim(-1.2*m, 1.2*m)
        a2.set_title(f"t = {t0:.2f} {T_UNIT}: area under curve = y(t) = {fy(t0):.3f} {Y_UNIT}")
        return line, dot, integ
    ani = FuncAnimation(fig3, update, frames=frames, interval=50, blit=False)
    fig3.tight_layout()
    if show:
        plt.show()      # figures and animation are displayed only; nothing is saved to disk
    plt.close("all")

def ask(prompt, default):
    """Print a prompt, show the default, and return the user's answer."""
    ans = input(f"{prompt} [default: {default}]: ").strip()
    return ans if ans else default

def interactive_inputs():
    """Tell the user exactly what to type, then read it."""
    print("=" * 70)
    print(" Green's Function IVP Solver:  a*y'' + b*y' + c*y = f(t),  y(0)=y'(0)=0")
    print("=" * 70)
    print("Enter the ODE as text, using y'' for y double prime and y' for y prime.")
    print("The independent variable can be written as t or x.")
    print("Use * for multiplication and ** (or ^) for powers. Examples:")
    print("    y''+y=4")
    print("    y''+4*y=x")
    print("    y''+2*y'+5*y=sin(t)")
    print("Initial conditions are fixed at y(0)=0 and y'(0)=0, for t >= 0.")
    print("Press Enter with nothing typed to use the two assigned ODEs.\n")
    odes = []
    while True:
        o = input(f"ODE #{len(odes)+1} (blank to {'finish' if odes else 'use the assigned ODEs'}): ").strip()
        if not o:
            break
        try:
            parse_ode(o)                      # validate right away
            odes.append(o)
        except Exception as e:
            print(f"  Could not read that ODE ({e}). Please try again, e.g. y''+y=4")
    if not odes:
        odes = ["y''+y=4", "y''+4*y=x"]
        print("Using the assigned ODEs:", odes)
    while True:
        try:
            tmax = float(ask("Final time t_max in seconds (a positive number)", f"{4*np.pi:.4f}"))
            if tmax > 0:
                break
        except ValueError:
            pass
        print("  Please enter a positive number, e.g. 12.5")
    show = ask("Open plot windows? (y/n)", "y").lower().startswith("y")
    return odes, tmax, show

def main():
    p = argparse.ArgumentParser(description="Solve y'' + ... = f(t), y(0)=y'(0)=0 via Green's function")
    p.add_argument("ode", nargs="*", help="ODE text, e.g. \"y''+y=4\" (omit to be prompted)")
    p.add_argument("--tmax", type=float, default=None, help="final time (s)")
    p.add_argument("--no-show", action="store_true", help="no plot windows (print results only)")
    args = p.parse_args()
    if args.ode:                                   # command-line mode
        odes, tmax, show = args.ode, args.tmax or 4*np.pi, not args.no_show
    else:                                          # interactive mode
        odes, tmax, show = interactive_inputs()
        if args.tmax: tmax = args.tmax
        if args.no_show: show = False
    if not show:
        matplotlib.use("Agg")
    for i, o in enumerate(odes, 1):
        run(o, tmax, show, f"ode{i}")

if __name__ == "__main__":
    main()