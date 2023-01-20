import os
import shutil
import json

import numpy as np
import matplotlib.pyplot as plt

from gsemo import GSEMO

def main():
    budget = 100_000
    dim = 100
    lamb = 10
    pid = 1
    nreps = 100
    logged = True
    path = os.path.realpath("./data")
    shutil.rmtree(path)

    for gtype in ("static", "TwoRate", "varctrl", "logNormal"):
        alg = GSEMO(budget, True, 1/dim, lamb, algorithm_name=gtype)
        result = alg.evaluate_problem(pid, dim, logged, nreps, path)

    for f in os.listdir(path):
        runs = read_data(os.path.join(path, f))
        ert = ert_front_onemm(runs, dim, budget)
        plt.semilogy(ert, label=f)

    plt.grid()
    plt.legend()
    plt.show()


def read_data(path):
    for f in os.listdir(path):
        if f.endswith(".json"):
            with open(os.path.join(path, f), "r") as h:
                data = json.loads(h.read())
                for scen in data['scenarios']:
                    data_file = os.path.join(path, scen['path'])
                    with open(data_file) as d:
                        runs = []
                        run = None
                        for line in d:
                            if line.startswith("evaluations"):
                                if run is not None:
                                    runs.append(np.array(run))
                                run = []
                                continue
                            run.append(list(map(float, line.strip().split())))
                        runs.append(np.array(run))
    return runs

def ert_front_onemm(runs, dim, budget):
    ggrid = np.zeros(dim + 1)
    success = np.zeros(dim + 1)
    for run in runs:
        grid = np.ones(dim + 1) * -1
        for time, x, y in run:
            if grid[int(x)] == -1:
                grid[int(x)] = time
                success[int(x)] += 1
        ggrid += grid

    ggrid += (len(runs) - success) * budget
    ert = ggrid / len(runs)
    return ert

if __name__ == "__main__":
    main()