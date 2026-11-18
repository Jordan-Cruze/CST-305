"""
CST-305 - Digital Storage Data Degradation Model

Project:
    Modeling Data Degradation in a Digital Storage System Using a System of ODEs

Description:
    This program models the degradation of a digital image stored in a
    hypothetical digital storage system. Five degradation factors are modeled
    using a linear system of ordinary differential equations (ODEs).

    The five factors are:
        1. Bit-flip errors
        2. Electrical/storage errors
        3. Write/update degradation
        4. Physical media degradation
        5. Backup/recovery errors

    Each degradation factor is represented by an ODE of the form:

        x_i' = -r_i*x_i + r_i

    where:
        x_i = degradation level of factor i
        r_i = degradation rate

    The individual degradation factors are combined using weighted values
    to calculate the overall degradation of the stored data.

    Data quality is calculated as:

        Q(t) = 1 - D(t)

    where D(t) is the overall degradation.

    The program also uses a sample image to visually demonstrate the effect
    of increasing data degradation over time. The image visualization is
    based on the mathematical degradation model and is intended as a
    representation of data corruption rather than a physical storage model.

Packages:
    NumPy
    SciPy
    SymPy
    Matplotlib
    Pillow

Implementation:
    1. Define or enter the ODE system.
    2. Convert the equations into matrix form.
    3. Enter initial conditions and simulation duration.
    4. Numerically solve the ODE system using SciPy's solve_ivp.
    5. Calculate overall weighted data degradation.
    6. Calculate remaining data quality.
    7. Display degradation factors and overall quality charts.
    8. Apply increasing visual noise to sample_image.jpg.
    9. Display progressively degraded versions of the image.

Programmers:
    Jordan Cruze
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from PIL import Image


# ============================================================
# CONSTANTS
# ============================================================

NUM_VARIABLES = 5

VARIABLE_NAMES = [
    "Bit-Flip Errors",
    "Electrical/Storage Errors",
    "Write/Update Degradation",
    "Physical Media Degradation",
    "Backup/Recovery Errors"
]

# The weights determine how strongly each degradation factor
# contributes to the overall degradation of the stored data.
DEGRADATION_WEIGHTS = np.array([
    0.30,   # Bit-flip errors
    0.15,   # Electrical/storage errors
    0.20,   # Write/update degradation
    0.25,   # Physical media degradation
    0.10    # Backup/recovery errors
])

# Image used for the visual degradation demonstration.
IMAGE_FILENAME = "sample_image.jpg"


# ============================================================
# PROGRAM HEADER
# ============================================================

def print_header():
    """Display the program title and project information."""

    print("=" * 70)
    print("CST-305 - DIGITAL STORAGE DATA DEGRADATION MODEL")
    print("=" * 70)
    print()
    print("This program models the degradation of a digital image")
    print("using a system of ordinary differential equations (ODEs).")
    print()
    print("Five degradation factors are modeled:")
    for i, name in enumerate(VARIABLE_NAMES, start=1):
        print(f"  x{i}: {name}")
    print()
    print("=" * 70)
    print()


# ============================================================
# EQUATION PARSING
# ============================================================

def parse_equation(equation):
    """
    Convert a user-entered equation into a SymPy expression.

    The user may enter either:

        x1' = -0.010*x1 + 0.010

    or simply:

        -0.010*x1 + 0.010

    The derivative notation is converted to a standard expression
    before SymPy parses the equation.
    """

    equation = equation.strip()

    # Remove spaces to simplify processing.
    equation = equation.replace(" ", "")

    # If the user enters x1' = ..., keep only the right side.
    if "=" in equation:
        equation = equation.split("=")[1]

    # Define the model variables.
    x1, x2, x3, x4, x5 = sp.symbols("x1 x2 x3 x4 x5")

    variables = {
        "x1": x1,
        "x2": x2,
        "x3": x3,
        "x4": x4,
        "x5": x5
    }

    try:
        expression = sp.sympify(equation, locals=variables)
        return expression

    except Exception as error:
        raise ValueError(
            f"Unable to parse equation: {equation}\n"
            f"Error: {error}"
        )


# ============================================================
# DEFAULT ODE SYSTEM
# ============================================================

def get_default_equations():
    """
    Return the default linear ODE system.

    Each equation follows:

        x_i' = -r_i*x_i + r_i

    The different rates represent different levels of vulnerability
    to degradation.
    """

    equations = [
        "-0.010*x1 + 0.010",
        "-0.015*x2 + 0.015",
        "-0.020*x3 + 0.020",
        "-0.008*x4 + 0.008",
        "-0.012*x5 + 0.012"
    ]

    return [parse_equation(eq) for eq in equations]


# ============================================================
# GET ODE EQUATIONS
# ============================================================

def get_ode_equations():
    """
    Ask the user whether to use the default ODE system or enter
    a custom system.
    """

    print("ODE SYSTEM")
    print("-" * 70)

    print("The program has a default system of five equations.")
    print()
    print("Default equations:")

    default_equations = [
        "x1' = -0.010*x1 + 0.010",
        "x2' = -0.015*x2 + 0.015",
        "x3' = -0.020*x3 + 0.020",
        "x4' = -0.008*x4 + 0.008",
        "x5' = -0.012*x5 + 0.012"
    ]

    for equation in default_equations:
        print(f"  {equation}")

    print()
    choice = input(
        "Use the default equations? [Y/n]: "
    ).strip().lower()

    if choice == "" or choice == "y" or choice == "yes":
        print()
        print("Using default ODE system.")
        print()
        return get_default_equations()

    # Custom equations
    print()
    print("Enter five equations.")
    print("Example: -0.010*x1 + 0.010")
    print()

    equations = []

    for i in range(NUM_VARIABLES):
        while True:
            try:
                equation = input(
                    f"Equation for x{i + 1}': "
                )

                parsed = parse_equation(equation)
                equations.append(parsed)
                break

            except ValueError as error:
                print(error)
                print("Please enter the equation again.")

    print()

    return equations


# ============================================================
# EXTRACT MATRIX FORM
# ============================================================

def extract_linear_system(equations):
    """
    Convert the symbolic equations into matrix form:

        x'(t) = A*x(t) + f

    where:
        A = coefficient matrix
        f = forcing/input vector
    """

    x1, x2, x3, x4, x5 = sp.symbols("x1 x2 x3 x4 x5")

    variables = [x1, x2, x3, x4, x5]

    A = np.zeros((NUM_VARIABLES, NUM_VARIABLES))
    f = np.zeros(NUM_VARIABLES)

    for i, equation in enumerate(equations):

        for j, variable in enumerate(variables):
            coefficient = equation.coeff(variable)

            # Require numerical coefficients for the linear system.
            try:
                A[i, j] = float(coefficient)
            except (TypeError, ValueError):
                raise ValueError(
                    "All coefficients must be numerical constants."
                )

        # Remove all variable terms to find the forcing term.
        forcing = equation

        for variable in variables:
            forcing = forcing.subs(variable, 0)

        try:
            f[i] = float(forcing)
        except (TypeError, ValueError):
            raise ValueError(
                "Forcing terms must be numerical constants."
            )

    return A, f


# ============================================================
# DISPLAY MATRIX
# ============================================================

def display_matrix(A, f):
    """Display the ODE system in matrix form."""

    print("MATRIX FORM")
    print("-" * 70)

    print("The system is represented as:")
    print()
    print("       x'(t) = A x(t) + f")
    print()

    print("Coefficient Matrix A:")
    print(A)

    print()
    print("Forcing Vector f:")
    print(f)

    print()


# ============================================================
# INITIAL CONDITIONS
# ============================================================

def get_initial_conditions():
    """
    Get initial conditions for all five degradation variables.

    The default is x_i(0) = 0, meaning no degradation at the
    beginning of the simulation.
    """

    print("INITIAL CONDITIONS")
    print("-" * 70)

    print(
        "Enter the starting degradation level for each factor."
    )

    print(
        "Values should be between 0 and 1."
    )

    print(
        "Press Enter to use 0."
    )

    print()

    initial_conditions = []

    for i, name in enumerate(VARIABLE_NAMES):

        while True:
            value = input(
                f"x{i + 1}(0) - {name} [0]: "
            ).strip()

            if value == "":
                value = 0.0
            else:
                try:
                    value = float(value)
                except ValueError:
                    print("Please enter a numerical value.")
                    continue

            if 0 <= value <= 1:
                initial_conditions.append(value)
                break

            print("Value must be between 0 and 1.")

    print()

    return np.array(initial_conditions)


# ============================================================
# SIMULATION TIME
# ============================================================

def get_simulation_time():
    """Get the length of the simulation in years."""

    print("SIMULATION TIME")
    print("-" * 70)

    while True:

        value = input(
            "Simulation duration in years [100]: "
        ).strip()

        if value == "":
            value = 100.0
        else:
            try:
                value = float(value)
            except ValueError:
                print("Please enter a numerical value.")
                continue

        if value > 0:
            print()
            return value

        print("Simulation time must be greater than zero.")


# ============================================================
# SOLVE ODE SYSTEM
# ============================================================

def solve_ode_system(A, f, initial_conditions, simulation_time):
    """
    Numerically solve the system:

        x'(t) = A*x(t) + f

    using SciPy's solve_ivp with the RK45 method.
    """

    def ode_function(t, x):
        return A @ x + f

    # Generate time values for the simulation.
    time = np.linspace(
        0,
        simulation_time,
        500
    )

    # Solve the system numerically.
    solution = solve_ivp(
        ode_function,
        [0, simulation_time],
        initial_conditions,
        t_eval=time,
        method="RK45"
    )

    if not solution.success:
        raise RuntimeError(
            f"ODE solver failed: {solution.message}"
        )

    return solution.t, solution.y


# ============================================================
# CALCULATE OVERALL DEGRADATION
# ============================================================

def calculate_overall_degradation(solution):
    """
    Calculate overall data degradation using weighted factors.

        D(t) = w1*x1 + w2*x2 + ... + w5*x5
    """

    degradation = (
        DEGRADATION_WEIGHTS @ solution
    )

    # Keep values within the valid range [0,1].
    degradation = np.clip(
        degradation,
        0,
        1
    )

    return degradation


# ============================================================
# CALCULATE DATA QUALITY
# ============================================================

def calculate_quality(degradation):
    """
    Convert degradation into remaining data quality.

        Q(t) = 1 - D(t)

    Returned as a percentage.
    """

    quality = (
        1 - degradation
    ) * 100

    return quality


# ============================================================
# DISPLAY MATHEMATICAL RESULTS
# ============================================================

def plot_simulation_results(
    time,
    solution,
    degradation,
    quality
):
    """
    Display both mathematical charts in the same figure.

    Chart 1:
        Individual degradation factors.

    Chart 2:
        Overall degradation and remaining data quality.
    """

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16, 6)
    )

    # --------------------------------------------------------
    # Chart 1 - Individual Degradation Factors
    # --------------------------------------------------------

    for i, name in enumerate(VARIABLE_NAMES):

        axes[0].plot(
            time,
            solution[i] * 100,
            linewidth=2,
            label=name
        )

    axes[0].set_title(
        "Individual Data Degradation Factors"
    )

    axes[0].set_xlabel(
        "Time (years)"
    )

    axes[0].set_ylabel(
        "Degradation (%)"
    )

    axes[0].set_ylim(
        0,
        100
    )

    axes[0].grid(
        True,
        alpha=0.3
    )

    axes[0].legend(
        fontsize=8
    )

    # --------------------------------------------------------
    # Chart 2 - Overall Degradation and Data Quality
    # --------------------------------------------------------

    axes[1].plot(
        time,
        degradation * 100,
        linewidth=3,
        label="Overall Degradation"
    )

    axes[1].plot(
        time,
        quality,
        linewidth=3,
        label="Remaining Data Quality"
    )

    axes[1].set_title(
        "Overall Data Integrity Over Time"
    )

    axes[1].set_xlabel(
        "Time (years)"
    )

    axes[1].set_ylabel(
        "Percentage (%)"
    )

    axes[1].set_ylim(
        0,
        100
    )

    axes[1].grid(
        True,
        alpha=0.3
    )

    axes[1].legend()

    fig.suptitle(
        "Digital Storage Data Degradation Model",
        fontsize=16,
        fontweight="bold"
    )

    plt.tight_layout()

    # One show() call displays both charts simultaneously.
    plt.show()


# ============================================================
# LOAD SAMPLE IMAGE
# ============================================================

def load_image():
    """
    Load the project's sample image.

    The image file should be stored in the same folder as main.py.
    """

    print("IMAGE VISUALIZATION")
    print("-" * 70)

    print(
        f"Loading image: {IMAGE_FILENAME}"
    )

    try:
        image = Image.open(
            IMAGE_FILENAME
        ).convert("RGB")

        print(
            f"Image loaded successfully: "
            f"{image.width} x {image.height} pixels"
        )

        print()

        return image

    except FileNotFoundError:

        print()
        print(
            f"ERROR: Could not find {IMAGE_FILENAME}."
        )

        print(
            "Make sure the image is in the same folder as main.py."
        )

        print()

        return None

    except Exception as error:

        print()
        print(
            f"ERROR: Unable to load image: {error}"
        )

        print()

        return None


# ============================================================
# CREATE DEGRADED IMAGE
# ============================================================

def create_degraded_image(
    image,
    degradation,
    seed=42
):
    """
    Create a visually degraded version of the image.

    Increasing Gaussian noise is added according to the
    degradation level calculated by the ODE model.

    This is a visualization of the mathematical model rather
    than a physical simulation of how storage hardware corrupts
    image files.
    """

    # Convert image into a NumPy array.
    image_array = np.asarray(
        image,
        dtype=np.float32
    )

    # Use a consistent seed so the visualization is repeatable.
    rng = np.random.default_rng(seed)

    # Generate random Gaussian noise.
    noise = rng.normal(
        loc=0,
        scale=degradation * 255,
        size=image_array.shape
    )

    # Apply the noise to the original image.
    degraded_array = image_array + noise

    # Keep RGB values between 0 and 255.
    degraded_array = np.clip(
        degraded_array,
        0,
        255
    ).astype(np.uint8)

    return Image.fromarray(
        degraded_array
    )


# ============================================================
# DISPLAY IMAGE PROGRESSION
# ============================================================

def display_image_progression(
    image,
    time,
    degradation,
    quality
):
    """
    Display the original image and four progressively degraded
    versions.

    Five points are selected from the simulation:
        0%
        25%
        50%
        75%
        100%
    """

    # Select five evenly spaced simulation points.
    indices = np.linspace(
        0,
        len(time) - 1,
        5
    ).astype(int)

    fig, axes = plt.subplots(
        1,
        5,
        figsize=(20, 5)
    )

    for position, index in enumerate(indices):

        current_degradation = degradation[index]
        current_quality = quality[index]

        if position == 0:
            # The first image is the original image.
            display_image = image
            title = (
                f"Original\n"
                f"Year {time[index]:.0f}\n"
                f"Quality: {current_quality:.1f}%"
            )

        else:
            # Later images use the calculated degradation.
            display_image = create_degraded_image(
                image,
                current_degradation,
                seed=42 + position
            )

            title = (
                f"Degraded\n"
                f"Year {time[index]:.0f}\n"
                f"Quality: {current_quality:.1f}%"
            )

        axes[position].imshow(
            display_image
        )

        axes[position].set_title(
            title,
            fontsize=10
        )

        axes[position].axis(
            "off"
        )

    fig.suptitle(
        "Visual Demonstration of Digital Data Degradation",
        fontsize=16,
        fontweight="bold"
    )

    plt.tight_layout()

    plt.show()


# ============================================================
# SAVE DEGRADED IMAGES
# ============================================================

def save_degraded_images(
    image,
    time,
    degradation,
    quality
):
    """
    Ask whether the user wants to save the generated degraded
    image stages to disk.
    """

    choice = input(
        "Save degraded image stages? [y/N]: "
    ).strip().lower()

    if choice not in ["y", "yes"]:
        return

    indices = np.linspace(
        0,
        len(time) - 1,
        5
    ).astype(int)

    for position, index in enumerate(indices):

        if position == 0:
            output_image = image

        else:
            output_image = create_degraded_image(
                image,
                degradation[index],
                seed=42 + position
            )

        filename = (
            f"degraded_stage_{position + 1}.jpg"
        )

        output_image.save(
            filename,
            quality=95
        )

        print(
            f"Saved: {filename}"
        )

    print()


# ============================================================
# PRINT FINAL RESULTS
# ============================================================

def print_final_results(
    time,
    solution,
    degradation,
    quality
):
    """
    Display the final numerical results from the simulation.
    """

    print("=" * 70)
    print("FINAL SIMULATION RESULTS")
    print("=" * 70)

    final_index = -1

    print()
    print(
        f"Simulation time: {time[final_index]:.2f} years"
    )

    print()

    print("Individual degradation factors:")

    for i, name in enumerate(VARIABLE_NAMES):

        print(
            f"  {name}: "
            f"{solution[i, final_index] * 100:.2f}%"
        )

    print()

    print(
        f"Overall data degradation: "
        f"{degradation[final_index] * 100:.2f}%"
    )

    print(
        f"Remaining data quality: "
        f"{quality[final_index]:.2f}%"
    )

    print()

    print("=" * 70)
    print("Simulation complete.")
    print("=" * 70)


# ============================================================
# MAIN PROGRAM
# ============================================================

def main():

    # --------------------------------------------------------
    # Step 1 - Display program information
    # --------------------------------------------------------

    print_header()

    # --------------------------------------------------------
    # Step 2 - Get ODE system
    # --------------------------------------------------------

    equations = get_ode_equations()

    # --------------------------------------------------------
    # Step 3 - Convert equations to matrix form
    # --------------------------------------------------------

    try:
        A, f = extract_linear_system(
            equations
        )

    except ValueError as error:

        print(
            f"ERROR: {error}"
        )

        return

    display_matrix(
        A,
        f
    )

    # --------------------------------------------------------
    # Step 4 - Get initial conditions
    # --------------------------------------------------------

    initial_conditions = (
        get_initial_conditions()
    )

    # --------------------------------------------------------
    # Step 5 - Get simulation duration
    # --------------------------------------------------------

    simulation_time = (
        get_simulation_time()
    )

    # --------------------------------------------------------
    # Step 6 - Solve ODE system
    # --------------------------------------------------------

    print("SOLVING ODE SYSTEM")
    print("-" * 70)

    try:

        time, solution = solve_ode_system(
            A,
            f,
            initial_conditions,
            simulation_time
        )

    except RuntimeError as error:

        print(
            f"ERROR: {error}"
        )

        return

    print(
        "ODE system solved successfully."
    )

    print()

    # --------------------------------------------------------
    # Step 7 - Calculate degradation
    # --------------------------------------------------------

    degradation = (
        calculate_overall_degradation(
            solution
        )
    )

    # --------------------------------------------------------
    # Step 8 - Calculate remaining quality
    # --------------------------------------------------------

    quality = (
        calculate_quality(
            degradation
        )
    )

    # --------------------------------------------------------
    # Step 9 - Display mathematical charts
    # --------------------------------------------------------

    print(
        "Displaying mathematical results..."
    )

    print(
        "Both charts will appear in the same window."
    )

    print()

    plot_simulation_results(
        time,
        solution,
        degradation,
        quality
    )

    # --------------------------------------------------------
    # Step 10 - Load sample image
    # --------------------------------------------------------

    image = load_image()

    # --------------------------------------------------------
    # Step 11 - Display image degradation
    # --------------------------------------------------------

    if image is not None:

        print(
            "Displaying visual degradation progression..."
        )

        print()

        display_image_progression(
            image,
            time,
            degradation,
            quality
        )

        # ----------------------------------------------------
        # Step 12 - Optionally save image stages
        # ----------------------------------------------------

        save_degraded_images(
            image,
            time,
            degradation,
            quality
        )

    else:

        print(
            "Image visualization skipped because the sample "
            "image could not be loaded."
        )

        print()

    # --------------------------------------------------------
    # Step 13 - Print final results
    # --------------------------------------------------------

    print_final_results(
        time,
        solution,
        degradation,
        quality
    )


# ============================================================
# PROGRAM ENTRY POINT
# ============================================================

if __name__ == "__main__":
    main()