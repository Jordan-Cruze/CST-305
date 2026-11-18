# CST-305 – Project 5: Lorenz System Simulation

**Grand Canyon University**  
**Course:** CST-305 – Ordinary Differential Equations  
**Programmer:** Jordan Cruze  
**Language:** Python

## Project Description

This project uses the Lorenz system of differential equations to simulate deterministic chaotic behavior. The Lorenz system models the behavior of a simplified atmospheric convection system and demonstrates how small changes in initial conditions and system parameters can produce significantly different outcomes.

The program numerically solves the Lorenz equations and generates visualizations of the system's behavior in three-dimensional space and across time. Users can enter different values for the parameter \(r\) to observe how the system changes under different conditions.

## Features

- Simulates the Lorenz system using numerical integration.
- Allows users to enter different values for the parameter \(r\).
- Generates a three-dimensional trajectory of the system.
- Generates individual X, Y, and Z position-versus-time graphs.
- Uses SciPy's `solve_ivp()` function with the RK45 numerical integration method.
- Displays plots one at a time for easier viewing.
- Allows users to run multiple simulations without restarting the program.

## Technologies and Libraries

| Technology | Purpose |
|---|---|
| Python | Primary programming language |
| NumPy | Numerical calculations and time arrays |
| SciPy | Numerical integration of differential equations |
| Matplotlib | Graphing and visualization |

## Requirements

- Python 3.10 or newer
- pip (Python package installer)

The required Python libraries are listed in `requirements.txt`:

```text
numpy
scipy
matplotlib
```

## How to Run the Program

### 1. Download the Repository

Clone or download the project repository and open a terminal in the project directory.

### 2. Create a Virtual Environment

Create a virtual environment to keep the project's Python dependencies separate from other Python projects.

```bash
python -m venv venv
```

### 3. Activate the Virtual Environment

Choose the command that matches your operating system and terminal.

**Windows – Command Prompt:**

```bash
venv\Scripts\activate
```

**Windows – PowerShell:**

```powershell
.\venv\Scripts\Activate.ps1
```

**macOS/Linux:**

```bash
source venv/bin/activate
```

Once activated, the terminal should display `(venv)` before the command prompt.

### 4. Install Dependencies

Install all required libraries using the provided requirements file:

```bash
pip install -r requirements.txt
```

### 5. Run the Program

Start the simulation:

```bash
python main.py
```

### 6. Enter Simulation Parameters

When prompted, enter a positive numerical value for the Lorenz parameter \(r\).

For example:

```text
Lorenz System Simulation
------------------------

Enter a value for r (or -1 to quit): 10
```

The program will generate four plots:

1. Three-dimensional Lorenz system trajectory.
2. X position over time.
3. Y position over time.
4. Z position over time.

Each plot appears in a separate window. Close the current plot window to display the next one.

After all four plots have been displayed, the program prompts the user for another value of \(r\). Enter `-1` to exit.

### 7. Deactivate the Virtual Environment

When finished running the program, deactivate the virtual environment:

```bash
deactivate
```

**Note:** The virtual environment only needs to be created once. When returning to the project, activate the existing environment before running the program. The `venv` directory should not be included in the GitHub repository.

## Mathematical Model

The Lorenz system consists of three coupled nonlinear differential equations:

$$
\frac{dx}{dt} = \sigma(y-x)
$$

$$
\frac{dy}{dt} = x(r-z)-y
$$

$$
\frac{dz}{dt} = xy-\beta z
$$

Where:

- \(x\) represents the system's first state variable.
- \(y\) represents the system's second state variable.
- \(z\) represents the system's third state variable.
- \(σ\) represents the Prandtl number.
- \(r\) represents the Rayleigh number.
- \(β\) represents a geometric parameter.

### Simulation Parameters

The program uses the following default parameters:

| Parameter | Value |
|---|---:|
| σ | 10 |
| β | 1.5 |
| Initial X | 11.8 |
| Initial Y | 4.4 |
| Initial Z | 2.4 |
| Start Time | 0 |
| End Time | 100 |
| Evaluation Points | 10,000 |

The value of \(r\) is entered by the user when the program runs.

## Program Implementation

The program is organized into several functions:

- **`lorenz_system()`** – Defines the three Lorenz differential equations and calculates the derivatives of the system's state variables.
- **`simulate()`** – Uses SciPy's `solve_ivp()` function to numerically solve the differential equations and return the simulation results.
- **`plot_results()`** – Generates the three-dimensional trajectory and individual X, Y, and Z graphs.
- **`main()`** – Handles user input, runs simulations, displays results, and allows repeated executions.

### Numerical Integration

The program uses the Runge-Kutta method (RK45) through SciPy's `solve_ivp()` function.

The numerical solver uses the following settings:

| Setting | Value |
|---|---:|
| Method | RK45 |
| Relative Tolerance | \(10^{-9}\) |
| Absolute Tolerance | \(10^{-11}\) |
| Maximum Step Size | 0.01 |

These settings provide tighter numerical tolerances and smaller integration steps to improve the accuracy of the simulation.

## Numerical Visualization

The program produces four visualizations for each selected value of \(r\):

1. **3D Trajectory:** Displays the system's movement through three-dimensional state space.
2. **X vs. Time:** Shows how the X variable changes throughout the simulation.
3. **Y vs. Time:** Shows how the Y variable changes throughout the simulation.
4. **Z vs. Time:** Shows how the Z variable changes throughout the simulation.

Different values of \(r\) demonstrate how the Lorenz system can exhibit different behaviors, including convergence toward an equilibrium and chaotic motion.