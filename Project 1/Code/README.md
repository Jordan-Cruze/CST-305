# CST-305 Project 1 – Visualize ODE With SciPy

## Overview

This project models CPU temperature after a heavy computer workload using an ordinary differential equation (ODE).

The program is designed specifically for a CPU-temperature model based on Newton's Law of Cooling. The user enters the equation used by the model, along with the initial CPU temperature and simulation time.

The program uses SymPy to interpret the equation, SciPy to numerically solve the ODE, NumPy for numerical calculations, and Matplotlib to visualize the result.

## Mathematical Model

The CPU-temperature model is:

$$
\frac{dT}{dt}=-k(T-T_a)
$$

where:

* `T` is CPU temperature in °C
* `t` is time in seconds
* `k` is the cooling coefficient
* `Ta` is ambient temperature

Example input:

```text
-0.08 * (T - 25)
```

## Requirements

* Python 3
* NumPy
* SciPy
* SymPy
* Matplotlib

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

Run:

```text
python main.py
```

The program will ask for:

1. The CPU-temperature ODE.
2. The initial CPU temperature.
3. The starting simulation time.
4. The ending simulation time.

Example:

```text
ODE: -0.08 * (T - 25)

Initial CPU temperature T(0) in °C: 85

Starting time in seconds: 0

Ending time in seconds: 60
```

The program then calculates the numerical solution, displays the results and estimated numerical error, and generates a graph of CPU temperature versus time.

## Packages

### NumPy

Used for numerical calculations and generating simulation time points.

### SciPy

Used to numerically solve the ODE with `solve_ivp()`.

### SymPy

Used to interpret the mathematical equation entered by the user.

### Matplotlib

Used to create the CPU-temperature visualization.

## Error Estimate

The program solves the ODE using two different accuracy settings and compares the results. The largest difference between the solutions is reported as an estimated maximum numerical error in degrees Celsius.

## Project Structure

```text
Project1/
├── main.py
├── README.md
├── requirements.txt
├── screenshots/
└── .venv/
```