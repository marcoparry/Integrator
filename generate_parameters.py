import numpy as np
import itertools
import os

initial_parameters = {
    "a1": 1.0,
    "a2": 25.0,
    "b1": 0.6,
    "b2": 3.0,
    "c1": 0.646,
    "c2": 0.6,
    "c34": 0.2,
    "m1": 0.46,
    "Q": 5.2e-3,
    "q1": 0.08,
    "w1": 0.24,
    "w2": 0.02,
    "q54": 0.6
}

q1_range = np.linspace(0, 1, 200)

res = ["0\t" + "\t".join(map(str, list(initial_parameters.values()))) + "\n"]

i = 0
# for m1, q1, w1 in itertools.product(m1_range, q1_range, w1_range):
for q1 in q1_range:
    params = initial_parameters.copy()
    params["q1"] = q1
    params_values = list(params.values())
    res.append("\t".join(map(str, [i + 1] + params_values)) + "\n")
    i += 1

current_location = os.path.dirname(os.path.abspath(__file__))
with open(f"{current_location}/data/parameters_different_q1.dat", "w+") as fp:
    fp.writelines(res)



