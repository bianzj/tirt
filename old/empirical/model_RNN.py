import torch
from torch import nn
import numpy as np
from Data_Preprocess import *
import matplotlib.pyplot as plt

# torch.manual_seed(1)    # reproducible

# Hyper Parameters
TIME_STEP = 10  # rnn time step
INPUT_SIZE = 10  # rnn input size
LR = 0.02  # learning rate

# show data
# steps = np.linspace(0, np.pi * 2, 100, dtype=np.float32)  # float32 for converting torch FloatTensor
# x_np = np.sin(steps)
# y_np = np.cos(steps)
# plt.plot(steps, y_np, 'r-', label='target (cos)')
# plt.plot(steps, x_np, 'b-', label='input (sin)')
# plt.legend(loc='best')
# plt.show()
# 导入训练集和验证集
x_train ,y_train,x_test,y_test = Preprocess()
# 将其转化为张量
x_train = torch.tensor(x_train.reshape(-1, 10, 10)).to(torch.float32)
y_train = torch.tensor(y_train.reshape(-1, 10, 1)).to(torch.float32)
x_test = torch.tensor(x_test.reshape(-1, 1, 10)).to(torch.float32)
# torch.Size([54961, 10])


class RNN(nn.Module):
	def __init__(self):
		super(RNN, self).__init__()

		self.rnn = nn.RNN(
			input_size=INPUT_SIZE,
			hidden_size=32,  # rnn hidden unit
			num_layers=1,  # number of rnn layer
			batch_first=True,  # input & output will has batch size as 1s dimension. e.g. (batch, time_step, input_size)
		)
		self.out = nn.Linear(32, 1)

	def forward(self, x, h_state):
		# x (batch, time_step, input_size)
		# h_state (n_layers, batch, hidden_size)
		# r_out (batch, time_step, hidden_size)
		r_out, h_state = self.rnn(x, h_state)

		outs = []  # save all predictions
		for time_step in range(r_out.size(1)):  # calculate output for each time step
			outs.append(self.out(r_out[:, time_step, :]))
		return torch.stack(outs, dim=1), h_state

	# instead, for simplicity, you can replace above codes by follows
	# r_out = r_out.view(-1, 32)
	# outs = self.out(r_out)
	# outs = outs.view(-1, TIME_STEP, 1)
	# return outs, h_state

	# or even simpler, since nn.Linear can accept inputs of any dimension
	# and returns outputs with same dimension except for the last
	# outs = self.out(r_out)
	# return outs


rnn = RNN()
print(rnn)

optimizer = torch.optim.Adam(rnn.parameters(), lr=LR)  # optimize all cnn parameters
loss_func = nn.MSELoss()

h_state = None  # for initial hidden state

plt.figure(1, figsize=(12, 5))
plt.ion()  # continuously plot

# for step in range(100):
# 	start, end = step * np.pi, (step + 1) * np.pi  # time range
# 	# use sin predicts cos
# 	steps = np.linspace(start, end, TIME_STEP, dtype=np.float32,
# 						endpoint=False)  # float32 for converting torch FloatTensor
# 	x_np = np.sin(steps)
# 	y_np = np.cos(steps)
#
# 	# x = torch.from_numpy(x_np[np.newaxis, :, np.newaxis])  # shape (batch, time_step, input_size)
# 	# y = torch.from_numpy(y_np[np.newaxis, :, np.newaxis])
#
# 	prediction, h_state = rnn(x_test, h_state)  # rnn output
# 	# !! next step is important !!
# 	h_state = h_state.data  # repack the hidden state, break the connection from last iteration
#
# 	loss = loss_func(prediction, y_test)  # calculate loss
# 	optimizer.zero_grad()  # clear gradients for this training step
# 	loss.backward()  # backpropagation, compute gradients
# 	optimizer.step()  # apply gradients
#
# 	# plotting
# 	plt.plot(steps, y_np.flatten(), 'r-')
# 	plt.plot(steps, prediction.data.numpy().flatten(), 'b-')
# 	plt.draw();
# 	plt.pause(0.05)

# 训练次数
for i in range(100):
	out,h_state= rnn(x_train,h_state)
	h_state = h_state.data
	loss = loss_func(out, y_train)
	optimizer.zero_grad()
	loss.backward()
	optimizer.step()
	if (i+1)%2==0:
		print('Epoch:{}, Loss:{:.5f}'.format(i+1, loss.item()))

#  对模型进行测试
pred = rnn(x_test,h_state)
pred_test = pred.view(-1).data.numpy()
from sklearn import metrics
print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, pred_test))
print('Mean Squared Error:', metrics.mean_squared_error(y_test, pred_test))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, pred_test)))
# 画图检验
plt.plot(y_test, 'b', label='real', linewidth=1)
plt.plot(pred.view(-1).data.numpy(), 'r', label='prediction',linewidth=1)

plt.legend(loc='best')
#plt.ylim((250, 350))
plt.show()