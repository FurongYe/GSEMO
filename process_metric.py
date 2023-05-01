import numpy as np
import os, json
import pandas as pd
import copy
# import seaborn as sns
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', 110) #replace n with the number of columns you want to see completely
pd.set_option('display.max_rows', 110) #replace n with the number of rows you want to see completely
import warnings
warnings.filterwarnings('ignore')
import sys

def get_PF(pf,y1,y2):
    for i in range(len(pf)):
        if (pf['y1'][i] >= y1) & (pf['y2'][i] >= y2) :
            return pf
        if (pf['y1'][i] < y1) & (pf['y2'][i] < y2) :
            pf.iloc[i,pf.columns.get_loc('y1')] = y1
            pf.iloc[i,pf.columns.get_loc('y2')] = y2
    pf = pf.append(pd.DataFrame([[y1,y2]],columns=['y1','y2']))
    pf = pf.drop_duplicates()
    pf = pf.sort_values(by=['y1','y2']).reset_index(drop=True)
    return pf

def HV(dt):
    y1_ref = -1
    y2_ref = -1
    result = 0
    for i in range(len(dt)):
        result += ((dt.iloc[i]['y1'] - y1_ref) * (dt.iloc[i]['y2'] - y2_ref))
        y1_ref = dt.iloc[i]['y1']
    return result

def NUM_PF(dt,id):
    if (id == 1):
        tmp = dt[(dt['y1'] + dt['y2']) == 100]
    if (id == 2):
        tmp = dt[(dt['y1'] + dt['y2']) == 100]
    if (id == 4):
        tmp = dt[(dt['y1'] >= 50) & (dt['y1'] + dt['y2'] == 150)]
    return len(tmp)

def pareto_txt_1(txt):
    y = txt.split(' ')
    y[0] = int(y[0])
    y[1] = int(y[1])
    return y[0]

def pareto_txt_2(txt):
    y = txt.split(' ')
    y[0] = int(y[0])
    y[1] = int(y[1])
    return y[1]

def get_per_run_metric(file,run,pid):
    rawdata = pd.read_csv(file)
    rawdata['y1'] = rawdata['pareto'].apply(lambda x: pareto_txt_1(x))
    rawdata['y2'] = rawdata['pareto'].apply(lambda x: pareto_txt_2(x))
    tmp_data = rawdata[["First_hit"+str(run),'y1','y2']]
    tmp_data = tmp_data[tmp_data["First_hit"+str(run)] != -1]

    tmp_data['y1_m'] = 0
    tmp_data['y2_m'] = 0
    tmp_data['hv_m'] = 0
    tmp_data['num_m'] = 0
    tmp_data = tmp_data.sort_values(by=['First_hit'+str(run)]).reset_index(drop=True)
    pf = pd.DataFrame([],columns=['y1','y2'])
    for i in range(len(tmp_data)):
        # print(i)
        if i == 0:
            tmp_data.iloc[i,tmp_data.columns.get_loc('y1_m')] = tmp_data['y1'][i]
            tmp_data.iloc[i,tmp_data.columns.get_loc('y2_m')] = tmp_data['y2'][i]
        else:
            tmp_data.iloc[i,tmp_data.columns.get_loc('y1_m')] = max(tmp_data['y1'][i],tmp_data['y1_m'][i-1])
            tmp_data.iloc[i,tmp_data.columns.get_loc('y2_m')] = max(tmp_data['y2'][i],tmp_data['y2_m'][i-1])
        pf = get_PF(pf,tmp_data['y1'][i],tmp_data['y2'][i])
        tmp_data.iloc[i,tmp_data.columns.get_loc('hv_m')] = HV(pf)
        tmp_data.iloc[i,tmp_data.columns.get_loc('num_m')] = NUM_PF(pf,pid)

    return tmp_data


def get_metric(file,pid):
    maxFE = -1
    rawdata = pd.read_csv(file)
    for i in range(100):
        maxFE = max(maxFE,max(rawdata['First_hit'+str(i)]))
    metric_pd = pd.DataFrame(list(range(1,maxFE+1)),columns=['FE'])
    metric_pd['y1_m'] = 0
    metric_pd['y2_m'] = 0
    metric_pd['hv_m'] = 0
    metric_pd['num_m'] = 0
    nrun = 100
    for i in range(nrun):
        metric_tmp = get_per_run_metric(file,i,pid)
        print(i)
        last_FE = 0
        for j in range(len(metric_tmp)):
            
            metric_pd.loc[last_FE: metric_tmp['First_hit'+str(i)][j],'y1_m'] += metric_tmp['y1_m'][j]
            metric_pd.loc[last_FE: metric_tmp['First_hit'+str(i)][j],'y2_m'] += metric_tmp['y2_m'][j]
            metric_pd.loc[last_FE: metric_tmp['First_hit'+str(i)][j],'hv_m'] += metric_tmp['hv_m'][j]
            metric_pd.loc[last_FE: metric_tmp['First_hit'+str(i)][j],'num_m'] += metric_tmp['num_m'][j]
            last_FE = metric_tmp['First_hit'+str(i)][j] + 1
        
        
        metric_pd.loc[last_FE: maxFE,'y1_m'] += metric_tmp['y1_m'][len(metric_tmp)-1]
        metric_pd.loc[last_FE: maxFE,'y2_m'] += metric_tmp['y2_m'][len(metric_tmp)-1]
        metric_pd.loc[last_FE: maxFE,'hv_m'] += metric_tmp['hv_m'][len(metric_tmp)-1]
        metric_pd.loc[last_FE: maxFE,'num_m'] += metric_tmp['num_m'][len(metric_tmp)-1]
        
    metric_pd['y1_m'] = metric_pd['y1_m'] / nrun
    metric_pd['y2_m'] = metric_pd['y2_m'] / nrun
    metric_pd['hv_m'] = metric_pd['hv_m'] / nrun
    metric_pd['num_m'] = metric_pd['num_m'] / nrun

    metric_pd.to_csv(file[:-4]+"-metric.csv",index=False)

def main():
    pathname = sys.argv[1]
    foldername = sys.argv[2]
    pid = (int)(sys.argv[3])
    print(os.path.join(pathname,foldername))
    get_metric(os.path.join(pathname,foldername),pid)

if __name__ == "__main__":
    main()

