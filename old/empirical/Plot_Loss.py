import matplotlib.pyplot as plt
import pandas  as pd
import math

### 此文件用于生成Loss曲线图，判断模型训练程度

loss = []
val_loss = []

txt_file_path = r"D:\data\lstSimulate\model\loss\\2012_2_loss.txt"


with open(txt_file_path,"r") as f:
	data = f.readlines()
# 通过循环消除回车，并以’，‘作为分隔
for line in data:
	line = line.rstrip()
	line = line.split(',')
	loss.append(float(line[0]))
	val_loss.append(float(line[1]))


# 画图
epoch = len(loss)
plt.figure()
plt.plot(loss, 'b', label='Loss', linewidth=2)
plt.plot(val_loss, 'r', label='val_loss', linewidth=1)
plt.legend(loc='best')
# plt.ylim((270, 320))
plt.xlabel('Epochs')
plt.ylim(-0.010,0.150)
plt.show()
