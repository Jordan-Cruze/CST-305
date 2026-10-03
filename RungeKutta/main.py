
#* Runge-Kutta Solver *#

print("Runge-Kutta Solver\n")
print("Enter the initial conditions and parameters for the Runge-Kutta solver:")

print("Enter the function f(x, y) for the differential equation:")
function_input = input()

def f(x, y):
    return eval(function_input)

print("Initial condition x0:")
x0 = float(input())

print("Initial condition y0:")
y0 = float(input())

print("Step size h:")
h = float(input())

print("Number of steps N:")
N = int(input())

print("Solving with Runge-Kutta Method")
print("-------------------------------")

def runge_kutta(f, x0, y0, h, N):
    x = x0
    y = y0
    for i in range(N):
        if x >= 10:
            break

        k1 = f(x, y)
        print("k1 = ", k1)
        
        k2 = f(x + (0.5 * h), y + (0.5 * k1 * h))
        print("k2 = ", k2)

        k3 = f(x + (0.5 * h), y + (0.5 * k2 * h))
        print("k3 = ", k3)

        k4 = f(x + h, y + (k3 * h))
        print("k4 = ", k4)

        y = y + (h / 6) * (k1 + (2 * k2) + (2 * k3) + k4)
        x = x + h

        print()
        print("x", i, " =  ", x, "y ", i, " = ", y)
        print("-------------------------------")

    return x, y

x_final, y_final = runge_kutta(f, x0, y0, h, N)

'''
from scipy.integrate import odeint
import numpy as np

print("ODE Solver\n")
print("Enter the differential equation:")

function_input = input()

def f(y, x):
    return eval(function_input)

print("Initial condition x0:")
x0 = float(input())

print("Initial condition y0:")
y0 = float(input())

print("Step size h:")
h = float(input())

print("Number of steps N:")
N = int(input())

# Create the x values
x_values = np.arange(x0, x0 + (N + 1) * h, h)

# Solve the ODE
y_values = odeint(f, y0, x_values)

print("\nODE Solution")
print("-------------------------------")

for i in range(len(x_values)):
    print("x", i, "=", x_values[i], "y", i, "=", y_values[i][0])

print("-------------------------------")

print("Final x =", x_values[-1])
print("Final y =", y_values[-1][0])
'''