import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from IPython import display
import itertools
import imageio
import os,sys


open_file = sys.argv[1]
print(open_file)
gif_name = sys.argv[2]
run_id = int(sys.argv[3])

file = open(open_file)
run_info = []
run = 1
for line in file:
    if line[0] == 'e':
        run += 1
        if run == run_id + 1:
            break
    else:
        if run == run_id:
            info = line.strip().split(' ')
            run_info.append(info)

t = pd.DataFrame(run_info)
t.columns =['evaluation', 'x1', 'x2']
t['x1'] = t['x1'].astype(float)
t['x2'] = t['x2'].astype(float)
t['x1'] = t['x1'].astype(int)
t['x2'] = t['x2'].astype(int)
t['evaluation'] = t['evaluation'].astype(float)
t['evaluation'] = t['evaluation'].astype(int)

x = list(range(0,101))
dt = []
for i in x:
    for j in range(0,101-i):
        dt.append([i,j,-1])
dt = pd.DataFrame(dt)
dt.columns = ['x1','x2','evaluation']

for i in range(len(t)):
    tt = t.loc[i]
    if dt[(dt['x1'] == tt['x1']) & (dt['x2'] == tt['x2'])]['evaluation'].item() == -1:
        dt[(dt['x1'] == tt['x1']) & (dt['x2'] == tt['x2'])] = dt[(dt['x1'] == tt['x1']) & (dt['x2'] == tt['x2'])].replace(-1,tt['evaluation'])
        if (len(dt[dt['evaluation']==-1]) == 0):
            break

tmp = dt
tmp = tmp[tmp['evaluation']!= -1]
tmp = tmp.sort_values(by=['evaluation'])
tmp

run_info = tmp
filenames = []
flag = True
index = 0
tmp_folder = gif_name+'-'+str(run_id)
os.mkdir(tmp_folder)
marker_size = 25
while index < len(run_info):
        # 清除原有图像
    # get current and next coordinates
    x = run_info['x1'][:index]
    y = run_info['x2'][:index]
    c = run_info['evaluation'][:index]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw = dict(aspect="equal"))
    plt.scatter(x, y, c=c,s=8,marker = "s")
    plt.title("tmp_folder")
    plt.xlim(0,100)
    plt.ylim(0,100)
        # remove spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # grid
    ax.set_axisbelow(True)
    # ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.7)
    # ax.xaxis.grid(color='gray', linestyle='dashed', alpha=0.7)
        # build file name and append to list of file names
    filename = f'{tmp_folder}/frame_{index}.png'
    filenames.append(filename)
       
    
    # save img
    plt.savefig(filename, dpi=96)
    plt.close()
    index += 15
    if (index > len(run_info)) & flag:
        index = len(run_info)-1
        flag = False
# Build GIF
print('creating gif\n')
with imageio.get_writer(f'{gif_name}.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
print('gif complete\n')
print('Removing Images\n')
# Remove files
for filename in set(filenames):
    os.remove(filename)
print('done')