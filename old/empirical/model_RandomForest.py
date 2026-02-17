from Data_Preprocess import *
import matplotlib.pyplot as plt

# 2013-2017年
trainyear = [13, 14, 15, 16]
testyear = [17]
stationlist = ["AR", "DM", "SDQ"]

# 将nan值转化为列平均值
def DeleteNan(t1):
	for i in range(t1.shape[0]):  # 将过小的异常值转化为nan值
		for j in range(t1.shape[1]):
			if t1[i, j] < -100:
				t1[i, j] = np.nan
			if t1[i, 10] < 200: # 排除极端温度数据情况
				print(1)
				t1[i, 10] = np.nan
	for i in range(t1.shape[1]):  # 遍历每一列（每一列中的nan替换成该列的均值）
		temp_col = t1[:, i]  # 当前的一列
		nan_num = np.count_nonzero(temp_col != temp_col)
		if nan_num != 0:  # 不为0，说明当前这一列中有nan
			temp_not_nan_col = temp_col[temp_col == temp_col]  # 去掉nan的ndarray
			# 选中当前为nan的位置，把值赋值为不为nan的均值
			temp_col[np.isnan(temp_col)] = temp_not_nan_col.mean()  # mean()表示求均值。
	return t1

# 处理时序数据,使每一天数据在一起，48行数据为一天
def ReformData(dataset):
	index = []
	for i in range(len(dataset) - 1):
		if dataset[i + 1, 9] == 1 and dataset[i, 9] == 2:  # 记录日夜变化的标签
			index.append(i)
	dataset = np.delete(dataset, range(index[-1] + 1, len(dataset)), axis=0)  # 删除尾部多余数据
	dataset = np.delete(dataset, range(0, index[0] + 1), axis=0)  # 输出头部多余数据
	return dataset


# 数据标准化
def Normalization(datalist, m1):
	# 两种最常用方法，min-max标准化，z-score标准化
	if m1:
		global max_value, min_value, scalar
		max_value = []
		min_value = []
		scalar = []
		for i in range(0, 11):
			# min-max标准化
			max_value.append(np.max(datalist[:, i]))
			min_value.append(np.min(datalist[:, i]))
			scalar.append(max_value[i] - min_value[i])
			datalist[:, i] = list(map(lambda x: (x - min_value[i]+0.00001) / scalar[i] + 0.00001, datalist[:, i]))
	else:  # z-score方法
		datalist = scale(X=datalist, with_mean=True, with_std=True, copy=True)
	return datalist

# 训练集数据预处理
def Preprocess2(trainstation):
	# 读入训练集数据
	for i in trainstation:
		xfilepath = r"D:\data\lstSimulate\inputX_2012_" + "%02d" % i + ".csv"
		yfilepath = r"D:\data\lstSimulate\inputY_2012_" + "%02d" % i + ".csv"
		xtmp = pd.read_csv(xfilepath, header=None).values
		ytmp = pd.read_csv(yfilepath, header=None).values
		if i == trainstation[0]:
			dataset = np.append(xtmp, ytmp, axis=1)
			dataset = DeleteNan(dataset)
			# dataset = ReformData(dataset)
			continue
		datatmp = np.append(xtmp, ytmp, axis=1)
		datatmp = DeleteNan(datatmp)
		dataset = np.append(dataset, datatmp, axis=0)
	# 标准化
	datalist = Normalization(dataset, 1)

	# 训练集
	# x_train = datalist[:, 0:10]
	# y_train = datalist[:, 10]
	train_dataset = datalist
	return train_dataset

# 测试集预处理
def PreprocessTest(teststation):
	xfilepath = r"D:\data\lstSimulate\inputX_2012_" + "%02d" % teststation + ".csv"
	yfilepath = r"D:\data\lstSimulate\inputY_2012_" + "%02d" % teststation + ".csv"
	xtmp = pd.read_csv(xfilepath, header=None).values
	ytmp = pd.read_csv(yfilepath, header=None).values

	dataset = np.append(xtmp, ytmp, axis=1)
	dataset = DeleteNan(dataset)
	# dataset = ReformData(dataset)

	testdatalist = dataset
	for i in range(0, 11):
		testdatalist[:, i] = (testdatalist[:, i] - min_value[i]) / scalar[i] + 0.00001
	# x_test = datalist[:, 0:10]
	# y_test = datalist[:, 10]
	test_dataset = testdatalist
	return test_dataset

# 导入训练集和验证集
train_dataset = Preprocess2(trainstation)
x_train = train_dataset[:,0:10]
y_train = train_dataset[:,10]
# 训练随机森林解决回归问题
from sklearn.ensemble import RandomForestRegressor

regressor = RandomForestRegressor(n_estimators=200, random_state=0)
regressor.fit(x_train, y_train)

for station in teststation:
	test_dataset = PreprocessTest(station)
	x_test = test_dataset[:, 0:10]
	y_test = test_dataset[:, 10]
	y_pred = regressor.predict(x_test)

	# 评估回回归性能
	from sklearn import metrics

	# print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
	# print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
	print(station,'Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test * scalar[10]+min_value[10], y_pred * scalar[10]+min_value[10])))
	# 画图检验
	plt.plot(y_pred* scalar[10]+ min_value[10], 'b', label='prediction', linewidth=2)
	plt.plot(y_test* scalar[10]+ min_value[10], 'r', label='real', linewidth=1)
	plt.legend(loc='best')
	plt.ylim((250, 350))
	plt.show()

