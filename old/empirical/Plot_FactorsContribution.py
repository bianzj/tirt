import matplotlib.pyplot as plt
import pandas  as pd
import math

# 1.画组内1-10站RF的贡献率图
importances_path = r"D:\data\lstSimulate\Result\RF\indicator\1\importance.xlsx"
data = pd.read_excel(importances_path,sheet_name='1',header=0,index_col=0)
labels = ['M1', 'M2', 'M3', 'M4', 'M5','M6','M7','M8','M9','M10','all']
features =  ["MS","LAI","Ta","WS","RH","DR","DLR"]
MS = data.iloc[0]
LAI = data.iloc[1]
Ta = data.iloc[2]
WS = data.iloc[3]
RH = data.iloc[4]
DS = data.iloc[5]
DL = data.iloc[6]
width = 0.75    # the width of the bars: can also be len(x) sequence
fig, ax = plt.subplots(figsize=(10,6))
ax.bar(labels,Ta, width, label='Ta')
ax.bar(labels,DS, width, label='DS',bottom = Ta)
ax.bar(labels,LAI, width, label='LAI',bottom = Ta +DS)
ax.bar(labels,RH, width, label='RH',bottom = Ta +DS + LAI)
ax.bar(labels,DL, width, label='DL',bottom = Ta +DS + LAI+RH)
ax.bar(labels,MS, width, label='MS',bottom = Ta +DS + LAI+RH+DL)
ax.bar(labels,WS, width, label='WS',bottom = Ta +DS + LAI+RH+DL+MS)
ax.set_ylabel('Factors Contribution')
ax.set_xlabel('Stations')
plt.ylim(0,1.1)
ax.legend(loc = 9,ncol=7)
plt.show()

# 3.画不同年不同气候类型RF的贡献率图
importances_path = r"D:\data\lstSimulate\Result\RF\indicator\3\importance.xlsx"
data = pd.read_excel(importances_path,sheet_name='1',header=0,index_col=0)
print(type(data))
labels = ['2013\nDM','2014\nDM','2015\nDM','2016\nDM','2017\nDM','all\nDM','2013\nAR', '2014\nAR','2015\nAR','2016\nAR','2017\nAR', 'all\nAR','2013\nSDQ','2014\nSDQ','2015\nSDQ','2016\nSDQ','2017\nSDQ','all\nSDQ']
features =  ["MS","LAI","Ta","WS","RH","DR","DLR"]
MS = data.iloc[0]
LAI = data.iloc[1]
Ta = data.iloc[2]
WS = data.iloc[3]
RH = data.iloc[4]
DS = data.iloc[5]
DL = data.iloc[6]
width = 0.75      # the width of the bars: can also be len(x) sequence
fig, ax = plt.subplots(figsize=(10,6))
ax.bar(labels,Ta, width, label='Ta')
ax.bar(labels,DS, width, label='DS',bottom = Ta)
ax.bar(labels,LAI, width, label='LAI',bottom = Ta +DS)
ax.bar(labels,RH, width, label='RH',bottom = Ta +DS + LAI)
ax.bar(labels,DL, width, label='DL',bottom = Ta +DS + LAI+RH)
ax.bar(labels,MS, width, label='MS',bottom = Ta +DS + LAI+RH+DL)
ax.bar(labels,WS, width, label='WS',bottom = Ta +DS + LAI+RH+DL+MS)
ax.set_ylabel('Factors Contribution')
ax.set_xlabel('Stations')
plt.ylim(0,1.1)
ax.legend(loc = 9,ncol=7)
plt.show()