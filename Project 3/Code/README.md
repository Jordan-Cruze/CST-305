# CST-305 Project 3 – Green's Function ODE Solver

## Overview

This project solves second-order ordinary differential equations (ODEs) using Green's functions and compares the resulting symbolic solution with a numerical solution.

The program accepts a second-order linear ODE of the form

**a y'' + b y' + c y = f(t)**

with the initial conditions

**y(0) = 0**  
**y'(0) = 0**

The program uses SymPy to parse and solve the differential equation symbolically, construct the Green's function, evaluate the Green's-function integral, and verify the resulting solution. SciPy is used to generate a numerical solution for comparison, NumPy is used for numerical calculations, and Matplotlib is used to generate graphs and an animation showing the construction of the Green's-function solution.

## Mathematical Approach

For a second-order linear ODE,

**a y'' + b y' + c y = f(t)**

the program first solves the associated homogeneous equation:

**a y'' + b y' + c y = 0**

Two linearly independent homogeneous solutions, y₁(t) and y₂(t), are used to construct the Green's function.

The Green's function is given by

**G(t, τ) = [y₁(τ)y₂(t) − y₂(τ)y₁(t)] / W(τ)**

where **W(τ)** is the Wronskian of the two homogeneous solutions.

For the zero initial conditions used by this project, the solution is calculated using

**y(t) = ∫₀ᵗ G(t, τ)f(τ)dτ**

The program then simplifies and verifies the resulting symbolic solution.

## Assigned Equations

The program was tested using the following assigned equations.

### Equation 1

**y'' + y = 4**

The Green's function is

**G(t, τ) = sin(t − τ)**

and the resulting solution is

**y(t) = 4 − 4cos(t)**

### Equation 2

**y'' + 4y = t**

The Green's function is

**G(t, τ) = ½sin(2t − 2τ)**

and the resulting solution is

**y(t) = t/4 − sin(2t)/8**

Both solutions satisfy the corresponding differential equation and the initial conditions y(0) = 0 and y'(0) = 0.

## Requirements

- Python 3
- NumPy
- SciPy
- SymPy
- Matplotlib

## Installation

Create and activate a Python virtual environment:

```text
python -m venv .venv
```

Activate it on Windows PowerShell:

```text
.venv\Scripts\Activate.ps1
```

Install the required packages:

```text
pip install -r requirements.txt
```

## Running the Program

Run the program with:

```text
python main.py
```

The program provides the option to enter an ODE manually or use the assigned equations.

Example input:

```text
y''+y=4
```

or

```text
y''+4*y=x
```

The program calculates the symbolic homogeneous solution, Green's function, Green's-function solution, and numerical solution. It also performs symbolic verification and calculates the maximum difference between the symbolic and numerical solutions.

The program generates three visualizations:

1. The homogeneous solution and homogeneous component.
2. The Green's-function solution compared with the numerical solution.
3. An animation showing the Green's-function integral being built over time.

## Program Output

For each equation, the program displays:

- The parsed differential equation.
- The homogeneous solution.
- The Green's function.
- The Green's-function solution.
- Verification of the differential equation.
- Verification of the initial conditions.
- The particular and homogeneous components.
- The maximum difference between the symbolic and numerical solutions.

For the assigned equations, the numerical comparison produced the following maximum differences:

| Equation | Green's-Function Solution | Maximum Difference |
|---|---|---:|
| y'' + y = 4 | y(t) = 4 − 4cos(t) | 9.15 × 10⁻⁹ m |
| y'' + 4y = t | y(t) = t/4 − sin(2t)/8 | 1.83 × 10⁻⁹ m |

These small differences demonstrate that the numerical solutions closely agree with the symbolic Green's-function solutions.

## Packages

### NumPy

NumPy is used for numerical calculations, creating arrays of time values, and processing the numerical solution.

### SciPy

SciPy is used to generate the numerical solution with `solve_ivp()` and compare it with the symbolic Green's-function solution.

### SymPy

SymPy is used to parse the differential equation, solve the homogeneous equation, construct the Green's function, evaluate the integral, simplify expressions, and verify the resulting solution.

### Matplotlib

Matplotlib is used to create the solution graphs and the animation showing the Green's-function integral being constructed.

## Verification

The program performs two types of verification.

### Symbolic Verification

The calculated solution is substituted back into the original differential equation. The program calculates the residual

**a y'' + b y' + c y − f(t)**

and verifies that it simplifies to zero.

The initial conditions are also checked to verify that

**y(0) = 0**

and

**y'(0) = 0**

### Numerical Verification

The symbolic Green's-function solution is compared with a numerical solution generated using SciPy's `solve_ivp()` function.

The maximum absolute difference between the two solutions is calculated. The resulting errors for the assigned equations are on the order of 10⁻⁹ meters, demonstrating close agreement between the two methods.

## Command-Line Usage

The program can also be run with an ODE provided directly as a command-line argument:

```text
python main.py "y''+y=4"
```

The maximum simulation time can be specified with:

```text
python main.py "y''+y=4" --tmax 20
```

The graphical display can be disabled with:

```text
python main.py "y''+y=4" --no-show
```

## Project Structure

```text
Project1/
├── main.py
├── README.md
├── requirements.txt
├── screenshots/
└── .venv/
```

## Project Limitations

The current implementation is designed for second-order linear ODEs with constant coefficients and zero initial conditions.

The input parser supports equations in the expected second-order format, such as:

```text
y''+y=4
```

and

```text
y''+4*y=x
```

The program is intended primarily for the assigned equations and similar linear ODEs rather than arbitrary differential equations.

## References

The project uses the official documentation for SymPy, NumPy, SciPy, and Matplotlib. These references are listed in the project documentation.