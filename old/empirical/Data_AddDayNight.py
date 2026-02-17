import pandas as pd
import numpy as np

trainyear = [13,14,15,16,17]
stationlist = ['AR','DM','SDQ']
# 添加白天和夜晚标识
for year in trainyear:
	for station in stationlist:
		xfilepath = r"D:\data\lstSimulate\inputX_20" + "%02d" % year + "_" + station + ".csv"
		data = pd.read_csv(xfilepath)
		# data.columns=['1','2','3','4','5','6','7','8','9']
		# print(data)
		# 下行短波辐射
		ds = data.values[:,5]
		# 白天夜晚标识
		dns = []
		for i in range(len(ds)):
			if ds[i]>10:
				dns.append(1)
			else:
				dns.append(0)
		dns = pd.DataFrame(dns)
		dataframe = data.join(dns)
		dataframe.to_csv(xfilepath,mode = 'w',index =None,sep=',')
		print(xfilepath)

