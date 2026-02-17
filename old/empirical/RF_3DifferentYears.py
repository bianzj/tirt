from Data_Preprocess import *
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame

# 2.13-17 年 AR、DM、SDQ三个超级站
# 14,16,17,18,19
# DM,AR,SDQ
trainyear = [13,14,15,16,17]
trainstation = ["SDQ"]
testyear = [13,14,15,16,17]
teststation = ["SDQ"]

if np.array(trainyear).shape[0] == 1:
	yearname  = "20" + str(trainyear[0])
else:
	yearname = "2013-17"

# 输出结果excel保存位置
save_excel_path = r"D:\data\lstSimulate\Result\RF\3\\"

# RMSE\R2\BIAS 保存位置
save_RMSE_path = r"D:\data\lstSimulate\Result\RF\indicator\3\3_" +str(trainstation[0]) + r"\\RMSE.xlsx"
save_r2_path = r"D:\data\lstSimulate\Result\RF\indicator\3\3_" +str(trainstation[0]) + r"\\r2.xlsx"
save_MAE_path = r"D:\data\lstSimulate\Result\RF\indicator\3\3_" +str(trainstation[0]) + r"\\MAE.xlsx"
save_bias_path = r"D:\data\lstSimulate\Result\RF\indicator\3\3_" +str(trainstation[0]) + r"\\bias.xlsx"
save_importance_path = r"D:\data\lstSimulate\Result\RF\indicator\3\\importance.xlsx"
bias_save = []
RMSE_save = []
r2_save = []
MAE_save = []
importance_save = []

x_file_path =  r"D:\data\lstSimulate\inputX_20"
y_file_path =  r"D:\data\lstSimulate\inputY_20"

# 将nan值转化为列平均值
def DeleteNan(t1):
	for i in range(t1.shape[0]):  # 将过小的异常值转化为nan值
		for j in range(t1.shape[1]):
			if t1[i, j] < -100:
				t1[i, j] = np.nan
			elif t1[i, j] <0:
				t1[i, j] =0
			if t1[i, 7] < 200: # 排除极端温度数据情况
				t1[i, 7] = np.nan
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
		for i in range(0, 8):
			# min-max标准化
			max_value.append(np.max(datalist[:, i]))
			min_value.append(np.min(datalist[:, i]))
			scalar.append(max_value[i] - min_value[i])
			if scalar[i] == 0:
				datalist[:, i] = 0
				continue
			datalist[:, i] = list(map(lambda x: (x - min_value[i]) / (scalar[i]+0.00001) + 0.00001, datalist[:, i]))
	else:  # z-score方法
		datalist = scale(X=datalist, with_mean=True, with_std=True, copy=True)
	return datalist

# 训练集数据预处理
def Preprocess2(trainstation,trainyear):
	# 读入训练集数据
	for i in trainyear:
		for station in trainstation:
			xfilepath = x_file_path + "%02d" % i + "_" + station + ".csv"
			yfilepath = y_file_path + "%02d" % i + "_" + station + ".csv"
			xtmp = pd.read_csv(xfilepath, header=None).values[:,0:7]
			ytmp = pd.read_csv(yfilepath, header=None).values
			if station == trainstation[0] and i == trainyear[0]:
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
def PreprocessTest(station,testyear):
	xfilepath = x_file_path + "%02d" % testyear + "_" + station + ".csv"
	yfilepath = y_file_path + "%02d" % testyear + "_" + station + ".csv"
	xtmp = pd.read_csv(xfilepath, header=None).values[:,0:7]
	ytmp = pd.read_csv(yfilepath, header=None).values

	dataset = np.append(xtmp, ytmp, axis=1)
	dataset = DeleteNan(dataset)
	# dataset = ReformData(dataset)

	testdatalist = dataset
	for i in range(0, 8):
		if scalar[i] == 0:
			testdatalist[:, i] = 0
			continue
		testdatalist[:, i] = (testdatalist[:, i] - min_value[i]) / scalar[i] + 0.00001
	# x_test = datalist[:, 0:10]
	# y_test = datalist[:, 10]
	test_dataset = testdatalist
	return test_dataset

# 测试集预处理，全部站点
def PreprocessTestAll(teststation,testyear):
	for i in testyear:
		for station in teststation:
			xfilepath = x_file_path + "%02d" % i + "_" + station + ".csv"
			yfilepath = y_file_path + "%02d" % i + "_" + station + ".csv"
			xtmp = pd.read_csv(xfilepath, header=None).values[:,0:7]
			ytmp = pd.read_csv(yfilepath, header=None).values
			if i == testyear[0] and station == teststation[0]:
				dataset = np.append(xtmp, ytmp, axis=1)
				dataset = DeleteNan(dataset)
				# dataset = ReformData(dataset)
				continue
			datatmp = np.append(xtmp, ytmp, axis=1)
			datatmp = DeleteNan(datatmp)
			dataset = np.append(dataset, datatmp, axis=0)
	testdatalist = dataset
	for i in range(0, 8):
		if scalar[i] == 0:
			testdatalist[:, i] = 0
			continue
		testdatalist[:, i] = (testdatalist[:, i] - min_value[i]) / scalar[i] + 0.00001
	# x_test = datalist[:, 0:10]
	# y_test = datalist[:, 10]
	test_dataset = testdatalist
	return test_dataset


# 导入训练集和验证集
train_dataset = Preprocess2(trainstation,trainyear)
x_train = train_dataset[:,0:7]
y_train = train_dataset[:,7]
# 训练随机森林解决回归问题
from sklearn.ensemble import RandomForestRegressor

regressor = RandomForestRegressor(n_estimators = 200, random_state=0,n_jobs=-1,oob_score=True,max_features="auto",max_depth=80,min_samples_split=2,min_samples_leaf=1)
regressor.fit(x_train, y_train)

importances =regressor.feature_importances_
features = ["MS","LAI","Ta","WS","RH","DR","DLR","UR","ULR" , "DNS"]
indices = np.argsort(importances)[::-1]
for f in range(x_train.shape[1]):
	print(features[f]+ ':'+str(importances[f]))
# # 保存模型
# import joblib
# if np.array(trainyear).shape[0] == 1:
# 	savefilepath = r"D:\lyf\LSTResult\RF\2\outputY_20" + "%02d" % trainyear[0]  +  trainstation[0]+'train_model.m'
# else:
# 	savefilepath = r"D:\lyf\LSTResult\RF\2\outputY_" + "all" + 'train_model.m'
# print(savefilepath)
# joblib.dump(regressor, savefilepath) # 存储
# # clf = joblib.load("train_model.m") # 调用

for i in testyear:
	for station in teststation:
		test_dataset = PreprocessTest(station, i)
		x_test = test_dataset[:, 0:7]
		y_test = test_dataset[:, 7]
		y_pred = regressor.predict(x_test)

		# 评估回回归性能
		from sklearn import metrics
		y_test = y_test * scalar[7] + min_value[7]
		y_pred = y_pred * scalar[7] + min_value[7]
		rmse = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
		RMSE_save.append(rmse)
		r2 = metrics.r2_score(y_test, y_pred)
		r2_save.append(r2)
		mae = metrics.mean_absolute_error(y_test, y_pred)
		MAE_save.append(mae)
		bias = np.mean(y_pred - y_test)
		bias_save.append(bias)
		# print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
		# print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
		print(trainyear, trainstation, '---', i,station, 'Root Mean Squared Error:',rmse)
		# # 画图检验
		# plt.plot(y_pred * scalar[10] + min_value[10], 'b', label='prediction', linewidth=2)
		# plt.plot(y_test * scalar[10] + min_value[10], 'r', label='real', linewidth=1)
		# plt.legend(loc='best')
		# plt.ylim((250, 350))
		# plt.show()
		savefilepath = save_excel_path + str(yearname) + "_" + trainstation[0] + "-" + "20" + str(i) + "_" + \
					   teststation[0] + ".csv"
		save = []
		save = np.append(np.reshape(y_test, (len(y_test), 1)), np.reshape(y_pred, (len(y_test), 1)), axis=1)
		save = pd.DataFrame(np.array(save).reshape(len(y_test), 2), columns=["measured", "predicted"])
		save.to_csv(savefilepath, index=False)


# 全部检验
test_dataset = PreprocessTestAll(teststation,testyear)
x_test = test_dataset[:, 0:7]
y_test = test_dataset[:, 7]
y_pred = regressor.predict(x_test)

y_test = y_test * scalar[7] + min_value[7]
y_pred = y_pred * scalar[7] + min_value[7]
rmse = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
RMSE_save.append(rmse)
r2 = metrics.r2_score(y_test, y_pred)
r2_save.append(r2)
mae = metrics.mean_absolute_error(y_test, y_pred)
MAE_save.append(mae)
bias = np.mean(y_pred - y_test)
bias_save.append(bias)

# print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
# print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
print(trainyear, trainstation, '---', 'all','Root Mean Squared Error:', rmse)
# # 画图检验
# plt.plot(y_pred* scalar[10]+ min_value[10], 'b', label='prediction', linewidth=2)
# plt.plot(y_test* scalar[10]+ min_value[10], 'r', label='real', linewidth=1)
# plt.legend(loc='best')
# plt.ylim((250, 350))
# plt.show()
savefilepath = save_excel_path + str(yearname) + "_" +trainstation[0] + "-" + "all_"+teststation[0] + ".csv"
save = []
save = np.append(np.reshape(y_test, (len(y_test), 1)), np.reshape(y_pred, (len(y_test), 1)), axis=1)
save = pd.DataFrame(np.array(save).reshape(len(y_test), 2), columns=["measured", "predicted"])
save.to_csv(savefilepath, index=False)

ind = yearname
rmse_xlsx = pd.read_excel(save_RMSE_path,sheet_name="1")
RMSE_save = pd.Series(RMSE_save,name =ind)
rmse_xlsx[ind] = RMSE_save
DataFrame(rmse_xlsx).to_excel(save_RMSE_path,sheet_name="1",index= None ,header=True)

r2_xlsx = pd.read_excel(save_r2_path,sheet_name="1")
r2_save = pd.Series(r2_save,name =ind)
r2_xlsx [ind] = r2_save
DataFrame(r2_xlsx ).to_excel(save_r2_path,sheet_name="1",index= None ,header=True)

mae_xlsx = pd.read_excel(save_MAE_path,sheet_name="1")
MAE_save = pd.Series(MAE_save,name =ind)
mae_xlsx[ind] = MAE_save
DataFrame(mae_xlsx).to_excel(save_MAE_path,sheet_name="1",index= None ,header=True)

bias_xlsx = pd.read_excel(save_bias_path,sheet_name="1")
bias_save = pd.Series(bias_save,name =ind)
bias_xlsx[ind] = bias_save
DataFrame(bias_xlsx).to_excel(save_bias_path,sheet_name="1",index= None ,header=True)

importance_xlsx = pd.read_excel(save_importance_path,sheet_name="1")
importance_save = pd.Series(importances,name =str(ind+trainstation[0]))
importance_xlsx[str(ind+trainstation[0])] = importance_save
DataFrame(importance_xlsx).to_excel(save_importance_path,sheet_name="1",index= None ,header=True)
print(ind+trainstation[0])

print("Save rmse,r2,mae successful")
