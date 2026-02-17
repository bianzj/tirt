from Data_Preprocess import *
import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import torch.nn.functional as F
import torch.optim as optim
import matplotlib.pyplot as plt





import torch.optim as optim
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
import torchvision.transforms as transforms
import torch.utils.data as Data

device = torch.device('cuda' if torch.cuda.is_available else 'cpu')
# 导入训练集和验证集
x_train, y_train, x_test, y_test = Preprocess(10)

# 将其转化为张量
x_train = torch.tensor(x_train.reshape(-1, 10, 10), device=device).to(torch.float32)
y_train = torch.tensor(y_train.reshape(-1, 10, 1), device=device).to(torch.float32)
x_test = torch.tensor(x_test.reshape(-1, 1, 10), device=device).to(torch.float32)
# torch.Size([54961, 10])

# 学习效率
learning_rate = 0.02
# 特征个数
in_dim = 10
# 线性层隐藏层神经元个数
hidden_dim = 100
# GRU隐藏层层数
n_layers = 2
# 输出维度
output_dim = 1
# 训练次数*scalar
times = 100
# 训练集和测试集选取密度（数值越大，密度越小）
k = 2


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
	# def __init__(self, in_dim, hidden_dim, output_dim=1, n_layers=1, batch_size=1):
	# 	super(Rnn, self).__init__()
	#
	# 	self.batch_size = batch_size
	# 	self.hidden_size = hidden_dim
	# 	self.n_layers = n_layers
	# 	self.out_size = output_dim
	#
	# 	# 这里指定了BATCH FIRST,所以输入时BATCH应该在第一维度
	# 	self.gru = torch.nn.GRU(in_dim, hidden_dim, n_layers, batch_first=True, bidirectional=True)
	#
	# 	# 加了一个线性层，全连接
	# 	self.fc1 = torch.nn.Linear(hidden_dim * 2, 300)
	# 	# 加入了第二个全连接层
	# 	self.fc2 = torch.nn.Linear(300, output_dim)
	#
	# def forward(self, word_inputs, hidden):
	# 	# hidden 就是上下文输出，output 就是 RNN 输出
	# 	output, hidden = self.gru(word_inputs, hidden)
	# 	# output是所有隐藏层的状态，hidden是最后一层隐藏层的状态
	# 	output = self.fc1(output)
	# 	output = self.fc2(output)
	#
	# 	# 仅仅获取 time seq 维度中的最后一个向量
	# 	# the last of time_seq
	# 	output = output[:, -1, :]
	#
	# 	return output, hidden


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

# 训练次数
for i in range(times):
	out = model(x_train)
	loss = loss_func(out, y_train)
	optimizer.zero_grad()
	loss.backward()
	optimizer.step()
	if (i + 1) % 1 == 0:
		print('Epoch:{}, Loss:{:.5f}'.format(i + 1, loss.item()))

#  对模型进行测试
pred = model(x_test)
pred_test = pred.view(-1).data.cpu().numpy()
from sklearn import metrics

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, pred_test))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, pred_test))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, pred_test)))
print(pred_test)
# 画图检验
plt.plot(y_test, 'b', label='real', linewidth=1)
plt.plot(pred.view(-1).data.cpu().numpy(), 'r', label='prediction', linewidth=1)

plt.legend(loc='best')
# plt.ylim((250, 350))
plt.show()
