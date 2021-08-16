import matplotlib.pyplot as plt
from differential_eq_py import DifferentialEquationSystem

PARAM_P = 5

def finit(x):
    if -PARAM_P < x < PARAM_P:
        return 5
    return 0

init_funcs = [
    lambda t: t,
    lambda t: t + 10
]

des = DifferentialEquationSystem(init_funcs=init_funcs, 
                                 h=0.1,
                                 T=1,
                                 L=100,
                                 relations=[[0, .1], [.1, 0]],
                                 finit_func=finit)
des.Solve(100)
solution = des.GetSolutionData()

plt.title("DifferentialEquationSystem example solution")
plt.xlabel("t")
plt.ylabel("x(t)")

for t, x in solution:
    plt.plot(t, x)

plt.show()