# CST-305 Project 4 – Digital Storage Data Degradation ODE Model

## Overview

This project uses a system of ordinary differential equations (ODEs) to model the degradation of data stored in a digital storage system over time.

The model represents five factors that can affect digital data integrity:

1. Bit-flip errors
2. Electrical and storage errors
3. Write and update degradation
4. Physical media degradation
5. Backup and recovery errors

Each factor is represented by a state variable in the ODE system. The individual degradation factors are combined using weighted contributions to calculate the overall degradation of the stored data. The remaining data quality is then calculated as a percentage.

The project uses Python to parse the ODE system, solve the equations numerically, calculate the overall data degradation, and visualize the results. Matplotlib is used to display the mathematical results, while Pillow is used to create progressively degraded versions of a sample image.

The image degradation is a visual representation of the mathematical model and is not intended to simulate the exact physical behavior of a specific storage device.

## Mathematical Approach

The project uses a system of five first-order linear ODEs.

The general form of the system is:

```text
x'(t) = A x(t) + f(t)
```

where `x(t)` is the vector containing the five degradation factors, `A` is the coefficient matrix, and `f(t)` is the forcing vector.

The five state variables are defined as:

```text
x₁(t) = bit-flip degradation
x₂(t) = electrical/storage degradation
x₃(t) = write/update degradation
x₄(t) = physical media degradation
x₅(t) = backup/recovery degradation
```

The default system used by the program is:

```text
x₁' = -0.010x₁ + 0.010
x₂' = -0.015x₂ + 0.015
x₃' = -0.020x₃ + 0.020
x₄' = -0.008x₄ + 0.008
x₅' = -0.012x₅ + 0.012
```

The matrix form is:

```text
        [ -0.010   0       0       0       0   ]
        [  0      -0.015   0       0       0   ]
x'  =   [  0       0      -0.020   0       0   ] x
        [  0       0       0      -0.008   0   ]
        [  0       0       0       0      -0.012]
```

with the forcing vector:

```text
        [0.010]
        [0.015]
f(t) =  [0.020]
        [0.008]
        [0.012]
```

The initial conditions are:

```text
x₁(0) = x₂(0) = x₃(0) = x₄(0) = x₅(0) = 0
```

Each equation follows the form:

```text
xᵢ' = rᵢ(1 - xᵢ)
```

where `rᵢ` represents the degradation rate for the corresponding factor.

## Overall Data Degradation

The individual degradation factors are combined using the following weights:

```text
Bit-flip errors:          0.30
Electrical/storage:       0.15
Write/update:             0.20
Physical media:           0.25
Backup/recovery:          0.10
```

The weights sum to 1.0.

The overall degradation is calculated using:

```text
D(t) = 0.30x₁ + 0.15x₂ + 0.20x₃ + 0.25x₄ + 0.10x₅
```

The remaining data quality is calculated as:

```text
Q(t) = 1 - D(t)
```

The program displays data quality as a percentage:

```text
Data Quality (%) = 100(1 - D(t))
```

A higher degradation value therefore corresponds to a lower remaining data quality.

## Simulation

The default simulation runs for 100 years using 500 time points.

The different degradation rates cause the five factors to increase at different speeds. The overall degradation is determined by both the individual rates and their assigned weights.

At the end of the default 100-year simulation, the model produces approximately:

```text
Overall degradation: 68.9%
Remaining data quality: 31.1%
```

These values are outputs of the mathematical model and should not be interpreted as a prediction for a specific physical storage medium.

## Requirements

- Python 3
- NumPy
- SciPy
- SymPy
- Matplotlib
- Pillow

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

Run the program from the project directory with:

```text
python main.py
```

The program first allows the user to choose whether to use the default ODE system or enter a custom system.

For the default model, the program uses the five equations described in the Mathematical Approach section.

The program then solves the ODE system and displays the results.

## Program Output

The program produces several outputs.

### ODE Results

The program displays the coefficient matrix and forcing vector used by the ODE system.

It then numerically solves the system using SciPy's `solve_ivp()` function.

### Mathematical Visualization

The first visualization contains two graphs:

1. Individual Data Degradation Factors
2. Overall Data Integrity Over Time

The first graph shows the change in each of the five degradation factors.

The second graph shows both overall degradation and remaining data quality over the simulation period.

### Image Visualization

After the mathematical graphs are closed, the program displays five versions of the sample image:

```text
Original
25 Years
50 Years
75 Years
100 Years
```

The amount of visual noise added to the image is based on the calculated degradation value at each point in time.

The image progression provides an intuitive visual representation of how increasing degradation can affect stored digital data.

The program expects the image file to be named:

```text
sample_image.jpg
```

and placed in the same directory as `main.py`.

### Final Numerical Results

The program also displays the final degradation values for each factor, the overall degradation, and the remaining data quality.

The program provides an option to save the generated degraded images.

## Packages

### NumPy

NumPy is used for numerical calculations, arrays, matrix operations, time values, and processing the ODE results.

### SciPy

SciPy is used to numerically solve the system of ordinary differential equations using the `solve_ivp()` function.

The default numerical method is the RK45 method.

### SymPy

SymPy is used to parse the ODE equations and convert the equations into the matrix and forcing-vector representation used by the numerical solver.

### Matplotlib

Matplotlib is used to create the graphs showing the individual degradation factors, overall degradation, and remaining data quality.

### Pillow

Pillow is used to load the sample image and create progressively degraded versions based on the calculated degradation values.

## Implementation

The program follows these general steps:

1. Start the program.
2. Load the default ODE system or accept a custom system.
3. Parse the equations using SymPy.
4. Convert the equations into matrix form.
5. Define the initial conditions.
6. Define the simulation time.
7. Solve the ODE system using SciPy.
8. Calculate the weighted overall degradation.
9. Calculate the remaining data quality.
10. Display the mathematical graphs.
11. Load `sample_image.jpg`.
12. Generate progressively degraded image versions.
13. Display the image progression.
14. Display the final numerical results.
15. Optionally save the degraded images.
16. End the program.

## Project Structure

```text
DigitalStorage/
├── main.py
├── README.md
├── requirements.txt
├── sample_image.jpg
└── .venv/
```

The `.venv` directory contains the Python virtual environment and does not need to be included in the GitHub repository.

## Project Limitations

The model is a simplified mathematical representation of digital data degradation.

The five degradation factors and their rates are assumptions used to demonstrate how a system of ODEs can represent multiple sources of data degradation. Actual degradation depends on the specific storage technology, environmental conditions, storage duration, write cycles, error-correction methods, and other physical factors.

The image visualization is also a simplified representation. The program uses artificial image noise to visually demonstrate increasing degradation rather than simulating actual bit-level corruption of a JPEG file.

The current implementation is primarily designed around the five-factor model used for this project. Custom equations may be entered, but they must follow the expected first-order linear system format supported by the program.

## References

The project uses the official documentation for NumPy, SciPy, SymPy, Matplotlib, and Pillow. Additional references concerning numerical ODE methods and digital storage are included in the project documentation.