import numpy as np
import pandas as pd
from scipy import interpolate

LAI_file_path = r"D:\data\lstSimulate\\SDQ_Lai_500m_2013_19.xlsx"
file_path = r"D:\data\lstSimulate\\inputX_20"
out_path = r"D:\data\lstSimulate\SDQ\\inputX_20"
HHL_15_path = r"D:\data\lstSimulate\2015\AWS_HHL_2015.csv"
startline = 48*150

years = [13,14,15,16,17,18,19]
for year in years:
	data = pd.read_excel(LAI_file_path,sheet_name=str(year))
	y = data["mean"]
	# 插值为全年52560个数据
	x = np.linspace(1, 17520, len(y))
	xnew = np.linspace(1, 17520, 17520)
	f = interpolate.interp1d(x, y, kind='linear')
	ynew = f(xnew)
	meteo_path = file_path + str(year) + "_SDQ.csv"
	meteo = pd.read_csv(meteo_path,header=None)
	LAI = ynew[startline:startline+len(meteo)]
	# meteo.loc[:,1] = LAI

	if year == 15 or year == 13:
		hhl = pd.read_csv(HHL_15_path, low_memory=False)
		ms = hhl["A.Ms_2(%)"]
		index = []
		for i in range(len(ms)):
			if (i + 1) % 6 != 1 and (i + 1) % 6 != 4:
				index.append(i)
		xdf = ms.drop(index=index)
		xdf = xdf.reset_index(drop=True)
		ms = xdf[startline:startline+len(meteo)]
		meteo.loc[:,0] = np.array(ms)


	out_file_path = out_path + str(year) + "_SDQ.csv"
	meteo.to_csv(out_file_path,header=None,index= None,mode="w",sep = ",")
	print(out_file_path)

