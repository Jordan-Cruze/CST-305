"""
CST-305 Project 1 - Visualize ODE With SciPy

Programmers:
    Jordan Cruze

Packages:
    NumPy
    SciPy
    SymPy
    Matplotlib

Approach:
    This program solves one specific ordinary differential equation
    (ODE) related to computer system performance: CPU temperature
    after a heavy workload.

    The CPU temperature model is based on Newton's Law of Cooling:

        dT/dt = -k(T - T_ambient)

    The program asks the user to enter the equation used for the
    CPU temperature model, along with the initial temperature and
    simulation time. SymPy is used to interpret the equation,
    SciPy's solve_ivp function is used to numerically solve the ODE,
    and Matplotlib is used to visualize the CPU temperature over time.

    The program is designed specifically for the CPU temperature
    model and is not intended to solve arbitrary types of ODEs.

    Example equation input:
        -0.08 * (T - 25)
"""


# NumPy is used for numerical calculations and generating time points.
import numpy as np

# Matplotlib is used to create the CPU temperature graph.
import matplotlib.pyplot as plt

# SciPy's solve_ivp function numerically solves the ODE.
from scipy.integrate import solve_ivp

# SymPy functions are used to interpret the equation entered by the user.
from sympy import sympify, symbols, lambdify


# ------------------------------------------------------------
# Display project information
# ------------------------------------------------------------
# Display the project title and explain what the program models.
print("=" * 60)
print("CST-305 Project 1 - ODE Visualization")
print("=" * 60)

print("\nComputer System Application:")
print("CPU temperature cooling after a heavy workload.")

# Explain the general form of the equation used by the program.
print("\nThe program solves an ODE of the form:")
print("dT/dt = f(t, T)")

# Give the user an example of an equation they can enter.
print("\nExample:")
print("dT/dt = -0.08 * (T - 25)")


# ------------------------------------------------------------
# Create mathematical variables
# ------------------------------------------------------------
# Create the symbolic variables t and T for time and temperature.
# These variables allow SymPy to interpret the user's equation.
t, T = symbols("t T")


# ------------------------------------------------------------
# Get the ODE from the user
# ------------------------------------------------------------
# Tell the user how the equation should be entered.
# Only the right side of dT/dt is entered.
print("\nEnter the right side of the ODE.")
print("Use 't' for time and 'T' for CPU temperature.")
print("Do not enter 'dT/dt ='.")

# Read the ODE entered by the user as text.
ode_input = input("\nODE: ")


# ------------------------------------------------------------
# Convert the user's input into a SymPy expression
# ------------------------------------------------------------
# SymPy converts the text entered by the user into a mathematical
# expression that can be used by the rest of the program.
try:
    ode_expression = sympify(ode_input)

except Exception:
    # If SymPy cannot interpret the user's input, display an error
    # and stop the program.
    print("\nError: The ODE could not be interpreted.")
    print("Example input:")
    print("-0.08 * (T - 25)")
    raise SystemExit


# ------------------------------------------------------------
# Convert the SymPy expression into a numerical function
# ------------------------------------------------------------
# Convert the symbolic expression into a Python function that
# accepts numerical values for time and temperature.
try:
    ode_function = lambdify(
        (t, T),
        ode_expression,
        modules="numpy"
    )

except Exception:
    # If the expression cannot be converted into a numerical
    # function, display an error and stop the program.
    print("\nError: The ODE could not be converted into a numerical function.")
    raise SystemExit


# ------------------------------------------------------------
# Get initial condition and simulation settings
# ------------------------------------------------------------
# Ask the user for the initial CPU temperature and the amount
# of time that should be simulated.
try:
    # The initial condition specifies the CPU temperature
    # at the beginning of the simulation.
    initial_temperature = float(
        input("\nEnter the initial CPU temperature T(0) in °C: ")
    )

    # Get the starting time for the simulation.
    start_time = float(
        input("Enter the starting time in seconds: ")
    )

    # Get the ending time for the simulation.
    end_time = float(
        input("Enter the ending time in seconds: ")
    )

except ValueError:
    # If the user enters something that is not a number,
    # display an error and stop the program.
    print("\nError: Please enter numerical values.")
    raise SystemExit


# ------------------------------------------------------------
# Validate simulation time
# ------------------------------------------------------------
# The ending time must be greater than the starting time
# for the simulation to be valid.
if end_time <= start_time:
    print("\nError: Ending time must be greater than starting time.")
    raise SystemExit


# ------------------------------------------------------------
# Generate time values for the solution
# ------------------------------------------------------------
# Create 300 evenly spaced time points between the starting
# and ending times. These points are used by solve_ivp to
# provide values for the graph.
time_points = np.linspace(
    start_time,
    end_time,
    300
)


# ------------------------------------------------------------
# Define the function used by SciPy's solve_ivp
# ------------------------------------------------------------
def ode_system(time, temperature):
    """
    Return the rate of change of CPU temperature.

    SciPy passes the current time and current temperature
    into this function.
    """

    # The user's equation calculates dT/dt.
    # temperature[0] contains the current CPU temperature.
    #
    # The result is placed inside a list because solve_ivp
    # expects the function to return an array-like value.
    return [ode_function(time, temperature[0])]


# ------------------------------------------------------------
# Solve the ODE using SciPy
# ------------------------------------------------------------
# solve_ivp numerically integrates the ODE from the starting
# time to the ending time.
#
# [initial_temperature] provides the initial condition.
# t_eval specifies the time points where the solution should
# be returned.
# rtol and atol control the numerical accuracy of the solver.
solution = solve_ivp(
    ode_system,
    (start_time, end_time),
    [initial_temperature],
    t_eval=time_points,
    rtol=1e-6,
    atol=1e-8
)


# ------------------------------------------------------------
# Check whether SciPy successfully solved the ODE
# ------------------------------------------------------------
# solve_ivp returns a success value indicating whether the
# numerical integration completed successfully.
if not solution.success:
    print("\nThe ODE solver failed.")
    print(solution.message)
    raise SystemExit


# ------------------------------------------------------------
# Display numerical results
# ------------------------------------------------------------
# Get the final temperature from the last value in the solution.
final_temperature = solution.y[0, -1]

print("\n" + "=" * 60)
print("ODE SOLUTION")
print("=" * 60)

# Display the equation and simulation information.
print(f"ODE entered: dT/dt = {ode_input}")
print(f"Initial temperature: {initial_temperature:.2f} °C")
print(f"Starting time: {start_time:.2f} seconds")
print(f"Ending time: {end_time:.2f} seconds")
print(f"Final temperature: {final_temperature:.2f} °C")

# Display the numerical solver and its accuracy settings.
print("\nSolver settings:")
print("Method: SciPy solve_ivp")
print("Relative tolerance: 1e-6")
print("Absolute tolerance: 1e-8")


# ------------------------------------------------------------
# Calculate an error estimate
# ------------------------------------------------------------
#
# A second solution is calculated using tighter tolerances.
# The difference between the normal and higher-accuracy
# solutions provides an estimate of the numerical error.
#
# This provides an estimate of the numerical error without
# requiring a known analytical solution for every possible ODE.
# ------------------------------------------------------------

# Solve the same ODE again with stricter accuracy settings.
high_accuracy_solution = solve_ivp(
    ode_system,
    (start_time, end_time),
    [initial_temperature],
    t_eval=time_points,
    rtol=1e-10,
    atol=1e-12
)

# Only calculate the error if the higher-accuracy solution
# was successfully generated.
if high_accuracy_solution.success:

    # Calculate the absolute difference between the two
    # numerical solutions at every time point.
    error = np.abs(
        solution.y[0] - high_accuracy_solution.y[0]
    )

    # Find the largest difference between the two solutions.
    maximum_error = np.max(error)

    # Display the estimated maximum numerical difference.
    print(f"Estimated maximum numerical error: "
          f"{maximum_error:.8f} °C")


# ------------------------------------------------------------
# Create the visualization
# ------------------------------------------------------------
# Create a figure that is 10 inches wide and 6 inches tall.
plt.figure(figsize=(10, 6))

# Plot CPU temperature against time.
# solution.t contains the time values.
# solution.y[0] contains the calculated temperatures.
plt.plot(
    solution.t,
    solution.y[0],
    label="CPU Temperature"
)

# Label the horizontal and vertical axes with their units.
plt.xlabel("Time (seconds)")
plt.ylabel("CPU Temperature (°C)")

# Create a graph title that also displays the equation
# entered by the user.
plt.title(
    "CPU Temperature Over Time\n"
    f"dT/dt = {ode_input}"
)

# Add grid lines to make the graph easier to read.
plt.grid(True)

# Display the CPU Temperature label on the graph.
plt.legend()

# Automatically adjust the spacing around the graph.
plt.tight_layout()

# Display the completed graph to the user.
plt.show()