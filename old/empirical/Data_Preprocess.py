import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale

# 数据预处理
# 1.读入数据，2.消除数据异常值，3.对数据进行标准化，4.准备训练集

def Preprocess(batchsize):
	# 1.读入数据 2\3\5\6\7\10\12\13\14\15
	# 读取自变量为xlist,土壤水，叶面积指数，空气温度，风速，空气湿度，下行短波辐射，下行长波辐射，上行短波辐射，上行长波辐射，白天夜晚标识
	# 因变量为地表温度
	# global dataset
	for i in [2,3,4,5,6,7,10,12,13,14,15]:
		xfilepath = r"D:\data\lstSimulate\inputX_2012_" + "%02d" % i + ".csv"
		yfilepath = r"D:\data\lstSimulate\inputY_2012_" + "%02d" % i + ".csv"
		xtmp = pd.read_csv(xfilepath, header=None).values
		ytmp = pd.read_csv(yfilepath, header=None).values
		if i == 2:
			dataset = np.append(xtmp, ytmp, axis=1)
			continue
		datatmp = np.append(xtmp, ytmp, axis=1)
		dataset = np.append(dataset, datatmp, axis=0)

	# 2.消除异常值
	# 删除所有y值为nan的数据行
	datalist = dataset[~np.isnan(dataset).any(axis=1)]

	# 3.对数据进行标准化
	# 两种最常用方法，min-max标准化，z-score标准化
	m1 = 0
	if m1:
		for i in range(0, 10):
			# min-max标准化
			max_value = np.max(datalist[:, i])
			min_value = np.min(datalist[:, i])
			scalar = max_value - min_value
			datalist[:, i] = list(map(lambda x: (x - min_value) / scalar, datalist[:, i]))
	else: # z-score方法
		datalist = scale(X=datalist,with_mean=True,with_std=True,copy=True)

	# 4.准备训练集
	# 选取训练集方法：1.随机选取，2.按顺序选取
	m2 = 1
	test_rate = 0.15
	x = datalist[:, 0:10]
	y = datalist[:, 10]

	if m2 == 1:
		# 随机选取
		x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_rate, random_state=1)
	elif m2 ==2:
		# 顺序选取
		datalist = datalist[:(len(datalist) // (batchsize )) * batchsize -1, :]
		data_size = len(datalist)
		x_train = datalist[:-int(test_rate * data_size), 0:10]
		y_train = datalist[:-int(test_rate * data_size), 10]
		x_test = datalist[-int(test_rate * data_size):, 0:10]
		y_test = datalist[-int(test_rate * data_size):, 10]
	elif m2 ==3:
		# 顺序选取，全部用来训练，全部用来验证
		datalist = datalist[:(len(datalist) // (batchsize )) * batchsize , :]
		x_train = datalist[:, 0:10]
		y_train = datalist[:, 10]
		x_test = datalist[:, 0:10]
		y_test = datalist[:, 10]
	x_train = x_train[:(len(x_train) // (batchsize)) * batchsize , :]
	y_train = y_train[:(len(y_train) // (batchsize)) * batchsize]
	x_test = x_test[:(len(x_test) // (batchsize)) * batchsize , :]
	y_test = y_test[:(len(y_test) // (batchsize)) * batchsize]
	return x_train,y_train,x_test,y_test