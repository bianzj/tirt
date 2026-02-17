import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import torch.nn.functional as F
import torch.optim as optim
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
import torchvision.transforms as transforms
import torch.utils.data as Data

# 修改参数处
times = 10
hiddensize = 128
nlayer = 2
inputsize = 10
outputsize = 1
seq_len = 48
batch_size = 64
# trainstation 2, 3, 5,6, 7, 10, 12, 13, 14, 15
trainstation = [2]
teststation = [2, 3,  5,6, 7, 10, 12, 13, 14, 15,1,4,7]

# 组内10组用在自身
# 采用2012,2、3、5、6、7、10、12、13、14、15，共10个站的数据进行训练，并在组内进行验证


# 将nan值转化为列平均值
def DeleteNan(t1):
    for i in range(t1.shape[0]):  # 将过小的异常值转化为nan值
        for j in range(t1.shape[1]):
            if t1[i, j] < -100:
                t1[i, j] = np.nan
            if t1[i, 10] < 200: # 排除极端温度数据情况
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

# 设置dataset数据格式
class DemoDatasetLSTM(Data.Dataset):
    """
        Support class for the loading and batching of sequences of samples
        Args:
            dataset (Tensor): Tensor containing all the samples
            sequence_length (int): length of the analyzed sequence by the LSTM
            transforms (object torchvision.transform): Pytorch's transforms used to process the data
    """

    ##  Constructor
    def __init__(self, dataset, sequence_length=1, transforms=None):
        self.dataset = dataset
        self.seq_len = sequence_length
        self.transforms = transforms

    ##  Override total dataset's length getter
    def __len__(self):
        return self.dataset.__len__()

    ##  Override single items' getter
    def __getitem__(self, idx):
        if idx + self.seq_len > self.__len__():
            if self.transforms is not None:
                item = torch.zeros(self.seq_len, self.dataset[0].__len__(),device =device)
                item[:self.__len__() - idx] = self.transforms(self.dataset[idx:])
                return item, item
            else:
                item = []
                item[:self.__len__() - idx] = self.dataset[idx:]
                return item, item
        else:
            if self.transforms is not None:
                return self.transforms(self.dataset[idx:idx + self.seq_len]), self.transforms(
                    self.dataset[idx:idx + self.seq_len])
            else:
                return self.dataset[idx:idx + self.seq_len], self.dataset[idx:idx + self.seq_len]


# list 转化为 Tensor
def listToTensor(list):
    tensor = torch.empty(list.__len__(), list[0].__len__(),device = device)
    for i in range(list.__len__()):
        tensor[i, :] = torch.FloatTensor(list[i])
    return tensor

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

# 测试集预处理，全部站点
def PreprocessTestAll(teststation):
    for i in teststation:
        xfilepath = r"D:\data\lstSimulate\inputX_2012_" + "%02d" % i + ".csv"
        yfilepath = r"D:\data\lstSimulate\inputY_2012_" + "%02d" % i + ".csv"
        xtmp = pd.read_csv(xfilepath, header=None).values
        ytmp = pd.read_csv(yfilepath, header=None).values
        if i == teststation[0]:
            dataset = np.append(xtmp, ytmp, axis=1)
            dataset = DeleteNan(dataset)
            # dataset = ReformData(dataset)
            continue
        datatmp = np.append(xtmp, ytmp, axis=1)
        datatmp = DeleteNan(datatmp)
        dataset = np.append(dataset, datatmp, axis=0)
        # print(i,dataset.shape)

    testdatalist = dataset
    for i in range(0, 11):
        testdatalist[:, i] = (testdatalist[:, i] - min_value[i]) / scalar[i] + 0.00001
    # x_test = datalist[:, 0:10]
    # y_test = datalist[:, 10]
    test_dataset = testdatalist
    return test_dataset

# 测试集预处理
def PreprocessTestDifYear(station,testyear):
    xfilepath = r"D:\data\lstSimulate\inputX_20" + "%02d" % testyear + "_" + station + ".csv"
    yfilepath = r"D:\data\lstSimulate\inputY_20" + "%02d" % testyear + "_" + station + ".csv"
    xtmp = pd.read_csv(xfilepath, header=None).values
    ytmp = pd.read_csv(yfilepath, header=None).values

    dataset = np.append(xtmp, ytmp, axis=1)
    dataset = DeleteNan(dataset)
    # dataset = ReformData(dataset)

    testdatalist = dataset
    for i in range(0, 11):
        if scalar[i] == 0:
            testdatalist[:, i] = 0
            continue
        testdatalist[:, i] = (testdatalist[:, i] - min_value[i]) / scalar[i] + 0.00001
    # x_test = datalist[:, 0:10]
    # y_test = datalist[:, 10]
    test_dataset = testdatalist
    return test_dataset

# 测试集预处理，全部站点
def PreprocessTestAllDifYear(teststation,testyear):
    for i in testyear:
        for station in teststation:
            xfilepath = r"D:\data\lstSimulate\inputX_20" + "%02d" % i + "_" + station + ".csv"
            yfilepath = r"D:\data\lstSimulate\inputY_20" + "%02d" % i + "_" + station + ".csv"
            xtmp = pd.read_csv(xfilepath, header=None).values
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
    for i in range(0, 11):
        if scalar[i] == 0:
            testdatalist[:, i] = 0
            continue
        testdatalist[:, i] = (testdatalist[:, i] - min_value[i]) / scalar[i] + 0.00001
    # x_test = datalist[:, 0:10]
    # y_test = datalist[:, 10]
    test_dataset = testdatalist
    return test_dataset

# 使用GPU进行及计算
device = torch.device('cuda' if torch.cuda.is_available else 'cpu')
# 导入训练集和验证集
# x_train ,y_train,x_test,y_test = Preprocess(trainstation)
train_dataset = Preprocess2(trainstation)

y_tr = train_dataset[:,10]
# 将其转化为张量
data_transform = transforms.Lambda(lambda x: listToTensor(x))

train_dataset = DemoDatasetLSTM(train_dataset, seq_len, transforms=data_transform)
train_data_loader = Data.DataLoader(train_dataset, batch_size, shuffle=True,drop_last=False)


# 构建一个简单的1D CNN模型
class CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1d = nn.Conv1d(in_channels=10, out_channels=1, kernel_size=3)
        self.relu = nn.ReLU(inplace=True)
        self.fc1 = nn.Linear(64 * 64, 50)
        self.fc2 = nn.Linear(50, 1)

    def forward(self, x):
        # 该模型的网络结构为 一维卷积层 -> Relu层 -> Flatten -> 全连接层1 -> 全连接层2
        x = self.conv1d(x)
        x = self.relu(x)
        x = x.view(-1)
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        return x


model = CNN()
model.to(device)
optimizer = torch.optim.Adam(model.parameters(),lr = 0.01)
loss_func = nn.MSELoss()
plt.ion()
if __name__=='__main__':
    for i in range(times):
        for j, data_ in enumerate(train_data_loader):
            x, _ = data_
            # print(i,j,x)
            a, b, c = x.shape
            if x[-1, -1, 0] == 0 or a * b * c != batch_size * seq_len * (inputsize + 1):
                continue
            x_train = x[:, :, 0:10].reshape([batch_size, seq_len, inputsize])
            y_train = x[:, :, 10].reshape([batch_size, seq_len, 1])
            # y_tr = x[:, :, 10].view(-1).data.cpu().numpy()
            out = model(x_train)
            loss = loss_func(out, y_train)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            if (i + 1) % 1 == 0 and j == 0:
                print('Epoch:{}, Loss:{:.5f}'.format(i + 1, loss.item()))
                plt.cla()
                # plt.ylim((250, 350))
                plt.plot(y_tr, 'b', label='real', linewidth=1)
                plt.plot(out.view(-1).data.cpu().numpy(), 'r', label='prediction', linewidth=1)
                # plt.ylim((250, 350))
                plt.legend(loc='best')
                plt.pause(0.0001)


# 检验十个普通站
for station in teststation:
    test_dataset = PreprocessTest(station)
    x_test = test_dataset[:, 0:10]
    y_test = test_dataset[:, 10]
    x_test = x_test[:(len(x_test) // (batch_size)) * batch_size, :]
    y_test = y_test[:(len(y_test) // (batch_size)) * batch_size]

    x_test = torch.tensor(x_test.reshape(-1, batch_size, inputsize), device=device).to(torch.float32)
    #  对模型进行测试
    pred = model(x_test)
    pred_test = pred.view(-1).data.cpu().numpy()
    from sklearn import metrics
    print(station,'Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test * scalar[10], pred_test * scalar[10])))

    # 画图检验
    plt.plot(y_test * scalar[10] + min_value[10], 'b', label='real', linewidth=1)
    plt.plot(pred_test * scalar[10] + min_value[10], 'r', label='prediction', linewidth=1)
    # print(s,miny)
    plt.legend(loc='best')
    plt.ylim((250, 350))
    plt.show()