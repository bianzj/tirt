import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

years = [14,16,17,18]
file_path = r"D:\data\lstSimulate\SDQ\inputX_20"

# 将nan值转化为列平均值
def DeleteNan(t1):
	for i in range(t1.shape[0]):  # 将过小的异常值转化为nan值
		for j in range(t1.shape[1]):
			if t1[i, j] < -100:
				t1[i, j] = np.nan
			elif t1[i, j] <0:
				t1[i, j] =0
			# if t1[i, 10] < 200: # 排除极端温度数据情况
			# 	t1[i, 10] = np.nan
	for i in range(t1.shape[1]):  # 遍历每一列（每一列中的nan替换成该列的均值）
		temp_col = t1[:, i]  # 当前的一列
		nan_num = np.count_nonzero(temp_col != temp_col)
		if nan_num != 0:  # 不为0，说明当前这一列中有nan
			temp_not_nan_col = temp_col[temp_col == temp_col]  # 去掉nan的ndarray
			# 选中当前为nan的位置，把值赋值为不为nan的均值
			temp_col[np.isnan(temp_col)] = temp_not_nan_col.mean()  # mean()表示求均值。
	return t1


for i in range(9):
	x = range(0,10*48)
	y = []
	title = ["MS","LAI","Ta","WS","RH","DR","DLR","UR","ULR"]
	for year in years:
		csv_file_path = file_path + str(year) + "_SDQ.csv"
		xdata = pd.read_csv(csv_file_path,header=None)
		xdata = DeleteNan(np.array(xdata))
		y = xdata[0:10*48,i]
		plt.plot(x,y,label = year)
		plt.legend(loc='best')
		plt.title(title[i])
	plt.show()




