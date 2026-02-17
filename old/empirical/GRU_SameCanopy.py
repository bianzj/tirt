import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
from Data_Preprocess import *
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import matplotlib.pyplot as plt

# 修改参数处
# 学习效率
learning_rate = 0.02
# 特征个数
in_dim = 10
# 线性层隐藏层神经元个数
hidden_dim = 128
# GRU隐藏层层数
n_layers = 2
# 输出维度
output_dim = 1
# 训练次数*scalar
times = 100000
batchsize = 10
trainstation = [3]
teststation = [6]

# 组内10组用在自身
# 采用2012,2、3、5、6、7、10、12、13、14、15，共5个站的数据进行训练，并在组内进行验证
# 训练模型的参数为200-64-2

# 消除nan值
def Deletenan(dataset):
	datalist = dataset[~np.isnan(dataset).any(axis=1)]
	return datalist

# 数据标准化
def Normalization(datalist,m1):
	# 两种最常用方法，min-max标准化，z-score标准化
	if m1:
		global max_value, min_value, scalar
		max_value = []
		min_value = []
		scalar = []
		for i in range(0, 11):
			# min-max标准化
			max_value.append( np.max(datalist[:, i]))
			min_value.append(np.min(datalist[:, i]))
			scalar.append(max_value[i] - min_value[i])
			datalist[:, i] = list(map(lambda x: (x - min_value[i]) / scalar[i] + 0.00001, datalist[:, i]))
	else:  # z-score方法
		datalist = scale(X=datalist, with_mean=True, with_std=True, copy=True)
	return datalist


def Preprocess(bathsize):
	# 1.读入数据 2\3\5\6\7\10\12\13\14\15
	# 读取自变量为xlist,土壤水，叶面积指数，空气温度，风速，空气湿度，下行短波辐射，下行长波辐射，上行短波辐射，上行长波辐射，白天夜晚标识
	# 因变量为地表温度
	# global dataset
	flag = 0
	for i in trainstation:
		xfilepath = r"D:\data\lstSimulate\inputX_2012_" + "%02d" % i + ".csv"
		yfilepath = r"D:\data\lstSimulate\inputY_2012_" + "%02d" % i + ".csv"
		xtmp = pd.read_csv(xfilepath, header=None).values
		ytmp = pd.read_csv(yfilepath, header=None).values
		if flag == 0:
			dataset = np.append(xtmp, ytmp, axis=1)
			flag = flag + 1
			continue
		datatmp = np.append(xtmp, ytmp, axis=1)
		dataset = np.append(dataset, datatmp, axis=0)
		flag = flag+1

	# 2消除nan
	datalist = Deletenan(dataset)

	# 3标准化,min-max标准化或Z score
	datalist = Normalization(datalist, 1)

	# 4.准备训练集
	# 选取训练集方法：1.随机选取，2.按顺序选取 3.全部
	m2 = 2
	test_rate = 0.20
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
		# datalist = datalist[:(len(datalist) // (batchsize )) * batchsize , :]
		x_train = datalist[:, 0:10]
		y_train = datalist[:, 10]
		x_test = datalist[:, 0:10]
		y_test = datalist[:, 10]
	x_train = x_train[:(len(x_train) // (batchsize)) * batchsize , :]
	y_train = y_train[:(len(y_train) // (batchsize)) * batchsize]
	x_test = x_test[:(len(x_test) // (batchsize)) * batchsize , :]
	y_test = y_test[:(len(y_test) // (batchsize)) * batchsize]
	return x_train,y_train,x_test,y_test


# Preprocess函数用来选取站来训练和检验
def Preprocess2(batchsize,trainstation,teststation):
	# 读入训练集数据
	for i in trainstation:
		xfilepath = r"D:\data\lstSimulate\inputX_2012_" + "%02d" % i + ".csv"
		yfilepath = r"D:\data\lstSimulate\inputY_2012_" + "%02d" % i + ".csv"
		xtmp = pd.read_csv(xfilepath, header=None).values
		ytmp = pd.read_csv(yfilepath, header=None).values
		if i == trainstation[0]:
			dataset = np.append(xtmp, ytmp, axis=1)
			continue
		datatmp = np.append(xtmp, ytmp, axis=1)
		dataset = np.append(dataset, datatmp, axis=0)

	# 消除nan
	datalist = Deletenan(dataset)
	# 标准化,min-max标准化或Z score
	global max_value, min_value, scalar
	max_value = []
	min_value = []
	scalar = []
	for i in range(0, 11):
		# min-max标准化
		max_value.append(np.max(datalist[:, i]))
		min_value.append( np.min(datalist[:, i]))
		scalar.append(max_value[i] - min_value[i])
		datalist[:, i] = list(map(lambda x: (x - min_value[i]) / scalar[i] + 0.00001, datalist[:, i]))
	x_train = datalist[:, 0:10]
	y_train = datalist[:, 10]
	# 消除不整除batchsize的余项部分
	x_train = x_train[:(len(x_train) // (batchsize)) * batchsize, :]
	y_train = y_train[:(len(y_train) // (batchsize)) * batchsize]

	# 读入测试数据
	flag = 0
	for i in teststation:
		xfilepath = r"D:\data\lstSimulate\inputX_2012_" + "%02d" % i + ".csv"
		yfilepath = r"D:\data\lstSimulate\inputY_2012_" + "%02d" % i + ".csv"
		xtmp = pd.read_csv(xfilepath, header=None).values
		ytmp = pd.read_csv(yfilepath, header=None).values
		if flag == 0:
			dataset = np.append(xtmp, ytmp, axis=1)
			continue
		datatmp = np.append(xtmp, ytmp, axis=1)
		dataset = np.append(dataset, datatmp, axis=0)
		flag = flag + 1

	# 消除nan
	testdatalist = Deletenan(dataset)
	# 标准化,min-max标准化或Z score
	for i in range(0, 11):
		testdatalist[:, i] = testdatalist[:, i] - min_value[i] / scalar[i] + 0.00001
	x_test = datalist[:, 0:10]
	y_test = datalist[:, 10]
	x_test = x_test[:(len(x_test) // (batchsize)) * batchsize, :]
	y_test = y_test[:(len(y_test) // (batchsize)) * batchsize]
	return x_train,y_train,x_test,y_test


# 使用GPU进行及计算
device = torch.device('cuda' if torch.cuda.is_available else 'cpu')
# 导入训练集和验证集
x_train ,y_train,x_test,y_test = Preprocess(batchsize)
# x_train ,y_train,x_test,y_test = Preprocess2(batchsize,trainstation,teststation)
# 将其转化为张量
x_train = torch.tensor(x_train.reshape(-1, 10, 10), device=device).to(torch.float32)
y_train = torch.tensor(y_train.reshape(-1, 10, 1), device=device).to(torch.float32)
x_test = torch.tensor(x_test.reshape(-1, 1, 10), device=device).to(torch.float32)
# torch.Size([54961, 10])

# 构建神经网络，第一层为线性层，第二层为GRU曾，第三层为dropout层，用于避免过拟合
# 第四层为线性层，最后一层为线性分类层
class Gru(nn.Module):
	def __init__(self, in_dim, hidden_dim, n_layers, output_dim):
		super(Gru, self).__init__()
		self.hidden_dim = hidden_dim
		self.ready = nn.Linear(in_dim, hidden_dim)
		self.gru = nn.GRU(hidden_dim, 2 * hidden_dim, n_layers, batch_first=True)
		self.outer = nn.Linear(2 * hidden_dim, hidden_dim)
		self.num_directions = 1
		self.linear = nn.Linear(hidden_dim * self.num_directions, output_dim)
		self.drop = nn.Dropout(p=0.5)

	def forward(self, x):
		out = self.ready(x)
		out, _ = self.gru(out)
		out = self.outer(out)
		out = self.linear(out)
		# out = self.drop(out)
		return out

# 把参数传入神经网络
model = Gru(in_dim, hidden_dim, n_layers, output_dim)
model.to(device)
# model = model.cuda()
# 定义loss和optimizer，使用交叉熵衡量模型损失，使用SGD优化器优化网络内部参数
criterion = nn.CrossEntropyLoss()
optimizer = optim.SGD(model.parameters(), lr=learning_rate)
#optimizer = torch.optim.Adam(model.parameters(), lr=0.02)
# 进行训练

loss_func = nn.MSELoss()

# # 训练次数
# for i in range(times):
# 	out = model(x_train)
# 	loss = loss_func(out, y_train)
# 	optimizer.zero_grad()
# 	loss.backward()
# 	optimizer.step()
# 	if (i+1)%2==0:
# 		print('Epoch:{}, Loss:{:.5f}'.format(i+1, loss.item()))
#
# #  对模型进行测试
# pred = model(x_test)
# pred_test = pred.view(-1).data.cpu().numpy()
# from sklearn import metrics
# print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, pred_test))
# print('Mean Squared Error:', metrics.mean_squared_error(y_test, pred_test))
# print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, pred_test)))

for i in range(times):
	out = model(x_train)
	loss = loss_func(out, y_train)
	optimizer.zero_grad()
	loss.backward()
	optimizer.step()
	if (i + 1) % 100 == 0:
		print('Epoch:{}, Loss:{:.5f}'.format(i + 1, loss.item()))

#  对模型进行测试
pred = model(x_test)
pred_test = pred.view(-1).data.cpu().numpy()
from sklearn import metrics

# print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, pred_test))
# print('Mean Squared Error:', metrics.mean_squared_error(y_test, pred_test))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test*scalar[10], pred_test*scalar[10])))


# 画图检验
plt.plot(y_test*scalar[10]+min_value[10], 'b', label='real', linewidth=1)
plt.plot(pred_test*scalar[10]+min_value[10], 'r', label='prediction',linewidth=1)
#print(s,miny)
plt.legend(loc='best')
plt.ylim((250, 350))
plt.show()