#https://www.youtube.com/watch?v=2AFGPdNn4FM

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import *
import pandas as pd
import csv
import plotly.plotly as py
import plotly.tools as tls
import seaborn as sns


#df = pd.read_csv('/Users/maryam.evari/Desktop/CBIC_qap_functional_temporal.csv')
#subject = pd.read_csv('/Users/maryam.evari/Desktop/qap_functional_temporal.csv')
#sub_id = subject.iloc[0].Participant  #

df = pd.read_csv('/Users/maryam.evari/Desktop/qap_functional_temporal.csv') #input of the program is QAP csv file
output_dir = '/Users/maryam.evari/Desktop/FD/'
sublist = '/Users/maryam.evari/Desktop/sublist.csv' #input of the program. it is list of subjects
sub_df = pd.read_csv(sublist)

    #not all subject have all functionals
for row in range(len(sub_df)):
    sub_id = sub_df.iloc[row,0]
    os.makedirs(os.path.join(output_dir,sub_id))
    sub_dir = output_dir + sub_id + '/'
    a,b = sub_id.split('-')
    sub_id = a+'_'+b
    subject =df[df.Participant == sub_id ]
    print(sub_dir)
    funcs = []
    for j in range(len(subject)):
        funcs.append(subject.iloc[j].Series)
        #print(funcs[j],sub_id,j)
        #funcs = ['func_peer_run_1','func_peer_run_2','func_peer_run_3','func_movieDM','func_movieTP','func_rest_run_1','func_rest_run_2']
    for i in range(len(funcs)):
            plt.clf
            print( sub_id,funcs[i])
            df1 = df[(df.Series == funcs[i]) & (df['RMSD (Mean)'] <= 1.5 )]
            x1 = df1['RMSD (Mean)'].drop(df1.index[0])
            ax = sns.distplot(x1, hist=True,rug=False, kde=True,kde_kws={'clip': (0.0, 1.5),'color': 'green'},bins= 54, color = 'lightblue', hist_kws={'edgecolor':'lightblue'})
            n = ax.get_yticks()
            m1 = int(max(n))
            fd_val1 = round(subject[(subject.Participant == sub_id)& (subject.Series == funcs[i])].xs('RMSD (Mean)',axis=1).iloc[0],4)
            sub_x1 =[]
            p=0
            for p in range(m1):
                sub_x1.append(fd_val1)
            
            sub_y1 =[]
            k=0
            for k in range(m1):
                sub_y1.append(k)
        
            l=plt.plot(sub_x1,sub_y1,color='r',linestyle='dashed',lw=2)
            plt.title(funcs[i])
            plt.xlabel('Mean FD')
            
            #drow point
            sub_x_dot1=[fd_val1]
            sub_y_dot1=[fd_val1]
            
            plt.scatter(sub_x_dot1,sub_y_dot1)
            plt.annotate(sub_x_dot1[0],(sub_x_dot1[0],sub_y_dot1[0]))
            #print(output_dir +sub_id +'_MeanFD_'+funcs[i]+'.png')
            plt.savefig(sub_dir+'MeanFD_'+funcs[i]+'.png')
            #   plt.show()
            plt.close()
