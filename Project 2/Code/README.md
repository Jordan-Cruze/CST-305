# CST-305 Project 2 – Runge-Kutta-Fehlberg (RKF45) for ODEs

## Overview

This project uses the Runge-Kutta-Fehlberg 4(5) method (RKF45) to solve an ordinary differential equation (ODE) and to measure the performance of the computer running it.

The assigned problem is:

$$
y' = \frac{y}{e^{x}-1}, \qquad x_0 = 1,\quad y_0 = 5,\quad h = 0.02
$$

The program implements RKF45 from scratch with a fixed step size, solves the ODE for 1,000 to 2,000 steps, and reports the number of computational steps, the number of function evaluations, and the actual computing time. It compares the RKF45 result with the exact solution and with a traditional method (Euler's method) to show the numerical error, and then plots the results.

The user may enter any ODE of the form `y' = f(x, y)` along with the initial conditions, step size, and number of steps. Pressing Enter at any prompt uses the assigned problem.

## Mathematical Model

Each RKF45 step computes six slopes from the current point `(xₙ, yₙ)`:

$$
\begin{aligned}
k_1 &= f(x_n,\ y_n) \\
k_2 &= f\left(x_n+\tfrac{h}{4},\ y_n+h\tfrac{k_1}{4}\right) \\
k_3 &= f\left(x_n+\tfrac{3h}{8},\ y_n+h\left(\tfrac{3k_1}{32}+\tfrac{9k_2}{32}\right)\right) \\
k_4 &= f\left(x_n+\tfrac{12h}{13},\ y_n+h\,\tfrac{1932k_1-7200k_2+7296k_3}{2197}\right) \\
k_5 &= f\left(x_n+h,\ y_n+h\left(\tfrac{439k_1}{216}-8k_2+\tfrac{3680k_3}{513}-\tfrac{845k_4}{4104}\right)\right) \\
k_6 &= f\left(x_n+\tfrac{h}{2},\ y_n+h\left(-\tfrac{8k_1}{27}+2k_2-\tfrac{3544k_3}{2565}+\tfrac{1859k_4}{4104}-\tfrac{11k_5}{40}\right)\right)
\end{aligned}
$$

and combines them into a 4th-order and a 5th-order estimate:

$$
y_{n+1}^{(4)} = y_n + h\left(\tfrac{25}{216}k_1+\tfrac{1408}{2565}k_3+\tfrac{2197}{4104}k_4-\tfrac{1}{5}k_5\right)
$$

$$
y_{n+1}^{(5)} = y_n + h\left(\tfrac{16}{135}k_1+\tfrac{6656}{12825}k_3+\tfrac{28561}{56430}k_4-\tfrac{9}{50}k_5+\tfrac{2}{55}k_6\right)
$$

The 5th-order value is carried forward to the next step. The local error estimate is `|y⁽⁵⁾ − y⁽⁴⁾|`.

For the assigned ODE, the equation is separable and has the exact solution:

$$
y(x) = \frac{y_0\left(1-e^{-x}\right)}{1-e^{-x_0}} = \frac{5\left(1-e^{-x}\right)}{1-e^{-1}}
$$

Example input:
```text
y/(exp(x)-1)
```

## Requirements

* Python 3
* NumPy
* Matplotlib
* SciPy (used for the reference solution and the adaptive RK45 comparison)

## Installation

Create and activate a Python virtual environment:

```text
python -m venv .venv
```

Activate it on Windows PowerShell:

```text
.venv\Scripts\Activate.ps1
```

Activate it on macOS or Linux:

```text
source .venv/bin/activate
```

Install the required packages:

```text
pip install -r requirements.txt
```

## Running the Program

Run:

```text
python main.py
```

The program will ask for:

1. The ODE `f(x, y)` (default `y/(exp(x)-1)`).
2. The initial value `x0` (default `1`).
3. The initial value `y0` (default `5`).
4. The step size `h` (default `0.02`).
5. The number of steps, from 1000 to 2000 (default `1000`).

Press Enter at any prompt to accept the default, or enter your own values.

Example using a different ODE:

```text
Enter f(x,y) for y' = f(x,y): -2*y + x

Initial x0: 0

Initial y0: 1

Step size h: 0.01

Number of steps (1000-2000): 1500
```

Expressions may use `x`, `y`, and the functions from Python's `math` module, such as `exp`, `sin`, `cos`, `log`, and `sqrt`.

The program then:

1. Prints a table of the first five RKF45 steps next to the exact values, along with the six stage slopes for the first step (used to check the manual calculation).
2. Solves the ODE for the chosen number of steps and for 2000 steps.
3. Solves two additional variations of the problem:
   * Variation 1: a step size four times smaller (`h/4`).
   * Variation 2: a shifted starting point (`x0 + 1`, `0.6 * y0`), which is `(2, 3)` for the assigned problem.
4. Reports, for each run, the number of steps, the number of function evaluations, the computing time in milliseconds, the final value, and the maximum error of RKF45 and Euler's method.
5. Runs SciPy's adaptive RK45 for comparison.
6. Opens three Matplotlib windows at the end:
   * **Runge-Kutta Plot** – the RKF45 numerical solution.
   * **Solved ODE** – the exact (or reference) solution.
   * **Overlay** – both solutions on the same axes.

The graphs are displayed only and are not saved to disk. Close the plot windows to end the program.

## Packages

### NumPy

Used for array handling when plotting.

### Matplotlib

Used to create the Runge-Kutta, solved ODE, and overlay plots.

### SciPy

Used for a high-accuracy reference solution (`solve_ivp` with the DOP853 method) when the entered ODE is not the assigned one, and for the adaptive RK45 comparison.

### math and time (standard library)

`math` is used to evaluate the entered function and the exact solution, and `time` (`perf_counter`) is used to measure computing time.

## Error Estimate

The program reports error in three ways:

* **Local error estimate:** the difference `|y⁽⁵⁾ − y⁽⁴⁾|` computed at every step.
* **Comparison with the exact solution:** for the assigned ODE, the maximum difference between RKF45 and the exact formula over all steps. For any other ODE, the exact formula is not available, so the program compares against a very tight SciPy reference solution (DOP853, relative tolerance 1e-13) and labels it as such.
* **Comparison with Euler's method:** the traditional first-order method, shown for contrast with RKF45.

For the assigned problem, the maximum RKF45 error is on the order of 1e-11 over 1,000 to 2,000 steps, while Euler's error is on the order of 1e-2.

## Accuracy and Computation Time

Each RKF45 step costs 6 function evaluations, compared with 1 for Euler's method. Reducing the step size `h` improves accuracy (the global error of RKF45 is O(h⁵)) but requires more steps and more computing time. Below roughly 1e-13, floating-point round-off dominates and a smaller `h` no longer improves accuracy.

## Project Structure

```text
Project2/
├── rkf45_project2.py
├── README.md
├── requirements.txt
├── screenshots/
└── .venv/
```