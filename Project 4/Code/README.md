# CST-305 – Project 4: Degradation of Data Integrity

**Grand Canyon University**  
**Course:** CST-305 – Ordinary Differential Equations  
**Programmer:** Jordan Cruze  
**Language:** Python

## Project Description

This project uses systems of ordinary differential equations (ODEs) to model data movement and degradation in digital storage systems.

The program models two different processor configurations:

1. A three-processor system where data moves between processors A, B, and C. Processor C also transfers data to an external network.
2. A closed two-processor system where data moves between processors A and B without entering or leaving the overall system.

The program uses matrix calculations, eigenvalue analysis, matrix exponentials, and numerical integration to examine how the systems behave over time.

## Features

- Three-processor mathematical model.
- Closed two-processor mathematical model.
- Custom processor model option.
- Automatic coefficient matrix construction.
- Differential equation generation.
- Eigenvalue calculations.
- Matrix exponential calculations.
- Numerical ODE integration using SciPy.
- Comparison of analytical and numerical solutions.
- Graph generation to visualize system behavior.
- Optional degradation rate for custom models.

## Technologies and Libraries

The project uses Python and the following libraries:

| Library | Purpose |
|---|---|
| NumPy | Matrix creation, array operations, and mathematical calculations. |
| SciPy | Matrix exponential calculations and numerical ODE integration. |
| Matplotlib | Generating graphs of the differential equation solutions. |
| argparse | Supporting command-line program options. |

## Requirements

Install Python 3 and the required libraries before running the program.

Install dependencies using:

```bash
pip install numpy scipy matplotlib
```

The `argparse` module is included with Python and does not require a separate installation.

## How to Run the Program

1. Download or clone the repository.
2. Open a terminal in the project directory.
3. Install the required Python libraries.
4. Run the Python program.

Example:

```bash
python main.py
```

Replace `main.py` with the actual Python filename if it is different.

Follow the prompts to select a mathematical model and enter the requested parameters.

## Mathematical Models

### Part 1: Three-Processor System

The first model consists of three interconnected processors, each with a capacity of 100 MB.

The assumed transfer rates are:

| Transfer | Rate |
|---|---:|
| A → B | 6 MB/s |
| B → A | 2 MB/s |
| B → C | 5 MB/s |
| C → B | 1 MB/s |
| C → Network | 4 MB/s |

The system is represented by:

\[
x'(t)=Ax(t)
\]

Where:

\[
A=
\begin{bmatrix}
-0.06&0.02&0\\
0.06&-0.07&0.01\\
0&0.05&-0.05
\end{bmatrix}
\]

The initial conditions are:

\[
x(0)=
\begin{bmatrix}
100\\100\\100
\end{bmatrix}
\]

The approximate eigenvalues are:

- -0.1054
- -0.05315
- -0.02141

All three eigenvalues are negative, indicating that the modeled system is asymptotically stable.

### Part 2: Closed Two-Processor System

The second model consists of two processors, A and B, with data moving between them.

The assumed transfer rates are:

| Transfer | Rate |
|---|---:|
| A → B | 2 MB/s |
| B → A | 3 MB/s |

The coefficient matrix is:

\[
A=
\begin{bmatrix}
-0.02&0.03\\
0.02&-0.03
\end{bmatrix}
\]

The initial state vector is:

\[
x(0)=
\begin{bmatrix}
1\\-1
\end{bmatrix}
\]

The eigenvalues are:

- 0
- -0.05

The analytical solution is:

\[
x(t)=
\begin{bmatrix}
e^{-0.05t}\\
-e^{-0.05t}
\end{bmatrix}
\]

The zero eigenvalue represents the conserved total state of the closed mathematical system.

**Note:** The initial state vector follows the mathematical initial conditions used in the assignment and program. The negative second component is not a physically meaningful quantity of stored data.

## Program Implementation

The program is organized into functions that handle different parts of the mathematical and computational process.

| Function | Description |
|---|---|
| `read_model()` | Collects model information and user input. |
| `build_matrix()` | Constructs the coefficient matrix from processor transfer rates. |
| `print_equations()` | Displays the resulting differential equations. |
| `solve_exact()` | Calculates the analytical solution using the matrix exponential. |
| `solve_numeric()` | Calculates the numerical solution using SciPy's ODE integration. |
| `plot_solution()` | Generates graphs displaying the system behavior over time. |

## Numerical Verification

The program compares two methods of solving the differential equation systems:

1. The analytical matrix exponential method using `scipy.linalg.expm()`.
2. The numerical integration method using `scipy.integrate.odeint()`.

The maximum absolute difference between the two solutions is used to evaluate numerical accuracy.

The recorded maximum errors were:

| Model | Maximum Absolute Error |
|---|---:|
| Three-processor system | \(1.07 \times 10^{-8}\) |
| Two-processor system | \(1.31 \times 10^{-7}\) |

These small differences indicate that the numerical solutions closely match the analytical solutions.

## Repository Contents

The repository contains the Python implementation and supporting project documentation.

- Python source code.
- Project documentation.
- README.md.
- Generated graphs and screenshots, if included.