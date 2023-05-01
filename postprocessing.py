# coding=UTF-8

import numpy as np
import os, json, sys
import pandas as pd

class solution:
    def __init__(self,y):
        self.y = y
    def __eq__(self, another):
        r = True
        for i in range(len(self.y)):
            r &= (self.y[i] == another.y[i])
        return r
    def __str__(self):
        r = ''
        for i in range(len(self.y)):
            r += str(self.y[i]) + ' '
        return r

def pareto_y(pname, dimension):
    y = []
    if pname == "OneMax":
        for i in range(dimension):
            y.append(solution([i,dimension - i]))
        y.append(solution([dimension, 0]))
    if pname == "LeadingOnes":
        y.append(solution([0, dimension]))
        for i in range(dimension):
            for j in range(dimension - i + 1):
                #if j != dimension - i:
                y.append(solution([i,j]))
        y.append(solution([dimension, 0]))
    if pname == "onejumpzerojump":
        k = 2
        for i in range(1, dimension - k + 1):
            if i < k:
                y.append(solution([i+k,i]))
            else:
                y.append(solution([i+k,dimension - i + k]))
        for i in range(dimension - k + 1, dimension):
            y.append(solution([dimension - i, dimension - i + k]))
        y.append(solution([k, dimension + k]))
        y.append(solution([dimension + k, k]))
    if pname == "COCZ":
        for i in range(dimension+1):
            if i < dimension / 2:
                for j in range(int(dimension / 2 + i + 1)):
                    #if j != dimension - i:
                    y.append(solution([i,j]))
            else:
                for j in range(int(dimension / 2 + (dimension - i) + 1)):
                    #if j != dimension - i:
                    y.append(solution([i,j]))
    return y


def read_table(jsondata, runs, pname, dimension):
    eval_table = pd.DataFrame()
    eval_table['pareto'] = pareto_y(pname,dimension)
    eval_table['algorithm'] = jsondata['algorithm']['name']
    eval_table['func'] = jsondata['function_name'] 
    eval_table['dimension'] = dimension

    

    for index in range(len(runs)):
        t = pd.DataFrame(runs[index])
        t.columns =['Evaluation', 'x1', 'x2']
        found_tag = 'found'+str(index)
        hit_tag = 'First_hit' + str(index)
        eval_table[found_tag] = False
        eval_table[hit_tag] = -1

        for index, row in eval_table.iterrows():
            if len(t[(t.x1 == eval_table.loc[index,'pareto'].y[0]) & (t.x2==eval_table.loc[index,'pareto'].y[1])]['Evaluation']) > 0:
                eval_table.loc[index,hit_tag] = min(t[(t.x1 == eval_table.loc[index,'pareto'].y[0]) & (t.x2==eval_table.loc[index,'pareto'].y[1])]['Evaluation'])
                eval_table.loc[index,found_tag] = True
            t = t.drop(t[(t.x1 == eval_table.loc[index,'pareto'].y[0]) & (t.x2==eval_table.loc[index,'pareto'].y[1])].index)
    return eval_table

def read_data(path):
    eval_table = []
    for f in os.listdir(path):
        if f.endswith(".json"):
            with open(os.path.join(path, f), "r") as h:
                # print(f,flush==True)
                data = json.loads(h.read())
                for scen in data['scenarios']:
                    data_file = os.path.join(path, scen['path'])
                    dimension = int(scen['dimension'])
                    pname = data['function_name']
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
                    eval_table = read_table(data,runs, pname, dimension)
                    
                    file_name = data['algorithm']['name'] + data['function_name'] + 'D' + str(dimension) + '.csv'
                    eval_table.to_csv(file_name,index=False)
    return eval_table

def main():
    pathname = sys.argv[1]
    foldername = sys.argv[2]
    print(os.path.join(pathname,foldername))
    read_data(os.path.join(pathname,foldername))

if __name__ == "__main__":
    main()

