import numpy as np
import pandas as pd
from osgeo import gdal
from osgeo import osr
import os
from scipy import interpolate
import matplotlib.pyplot as plt

scalar = 0.1
coordinates = [101.1374,42.0012]
file_list_path = r"D:\data\\lstSimulate\\2018LAI2\\25_4"
# csv_file_path = r"D:\data\LSTSIM\\2019HHLmeteo.xlsx"
excel_file_path = r"D:/data/LSTSIM/2018SDQmeteo.xlsx"
x_file_path = r"D:\data\lstSimulate\inputX_2018_SDQ.csv"
y_file_path = r"D:\data\lstSimulate\inputY_2018_SDQ.csv"
LAI_csv_path = r"D:\data\\lstSimulate\\2018LAI2\\statistics_Lai_500m.csv"

def get_file_info(in_file_path):
	'''
	根据图像文件路径，只读打开
	:param in_file_path: tif文件输入路径
	:return: 数据集、地理坐标系、投影坐标系、栅格影像大小
	'''
	pcs = None
	gcs = None
	shape = None

	if in_file_path.endswith(".tif") or in_file_path.endswith(".TIF"):
		dataset = gdal.Open(in_file_path)
		pcs = osr.SpatialReference()
		pcs.ImportFromWkt(dataset.GetProjection())
		gcs = pcs.CloneGeogCS()
		extend = dataset.GetGeoTransform()
		shape = (dataset.RasterXSize,dataset.RasterYSize)
	else:
		raise("Unsupported file formet!")
	return dataset,gcs,pcs,extend,shape

def lonlat_to_xy(gcs,pcs,lon,lat):
	"""
	经纬度坐标转为投影坐标
	:param gcs: 地理空间信息
	:param pcs: 投影信息
	:param lon: 经度
	:param lat: 纬度
	:return: 投影坐标
	"""
	ct = osr.CoordinateTransformation(gcs,pcs)
	coordinates = ct.TransformPoint(lon,lat)

	return coordinates[0],coordinates[1],coordinates[2]

def xy_to_rowcol(extend,x,y):
	'''
	根据六参数模型将投影坐标转为图上坐标
	:param extend: 图像空间范围
	:param x: 投影坐标x
	:param y: 投影坐标y
	:return: 行row，列col
	'''
	a = np.array([[extend[1],extend[2]],[extend[4],extend[5]]])
	b = np.array([x-extend[0],y-extend[3]])
	row_col = np.linalg.solve(a,b)
	row = int(np.floor(row_col[1]))
	col = int(np.floor(row_col[0]))
	return row , col

def get_value_by_coordinates(file_path,coordinates,coordiantes_type = 'rowcol'):
	'''
	根据图像坐标，返回对应像元的值

	:param file_path:图像文件的路径
	:param coordinates: 坐标
	:param coordiantes_type:坐标类型，rowcol，xy，lonlat
	:return: 像元值
	'''
	dataset,gcs,pcs,extend,shape =get_file_info(file_path)
	img = dataset.GetRasterBand(1).ReadAsArray()
	value = None
	if coordiantes_type == "rowcol":
		value = img[coordinates[0],coordinates[1]]
	elif coordiantes_type =='lonlat':
		x,y,_ = lonlat_to_xy(gcs,pcs,coordinates[0],coordinates[1])
		row,col = xy_to_rowcol(extend,x,y)
		value = int(img[row+1,col+1])+int(img[row,col+1])+int(img[row+1,col])+int(img[row,col])+int(img[row+1,col-1])+int(img[row-1,col+1])+int(img[row-1,col])+int(img[row,col-1])+int(img[row-1,col-1])
		print(int(img[row+1,col+1]),int(img[row,col+1]),int(img[row+1,col]),int(img[row,col]),int(img[row+1,col-1]),int(img[row-1,col+1]),int(img[row-1,col]),int(img[row,col-1]),int(img[row-1,col-1]))
		value = value/9
		# value = img[row-1,col]
	elif coordiantes_type =='xy':
		row,col = xy_to_rowcol(extend,coordinates[0],coordinates[1])
		value = img[row, col]
	else:
		raise("Coordinates_type:Wrong para input")
	return value


# # 获取输入文件夹中所有的文件名
# NameList = os.listdir(file_list_path)
# # 遍历文件名列表中所有文件
# # print("经纬度坐标：",coordinates[0],',',coordinates[1])
# y = []
# for i in range(len(NameList)):
# 	# 获取文件名后缀
# 	dirname, basename = os.path.split(NameList[i])
# 	filename, txt = os.path.splitext(basename)
# 	# 判断文件后缀是否为 .hdf
# 	if txt == '.tif' or txt == '.TIF':
# 		# dataset = get_file_info(NameList[i])
# 		in_file_path = file_list_path +'/' + NameList[i]
# 		pix_value = get_value_by_coordinates(in_file_path,coordinates,'lonlat') * scalar
# 		y.append(pix_value)

# plt.plot(y)
# plt.ylim(2.5,4.0)
# plt.show()

LAI = pd.read_csv(LAI_csv_path)
y = LAI["mean"]

x = np.linspace(1,52560,len(y))
xnew = np.linspace(1,52560,52560)

f=interpolate.interp1d(x,y,kind='linear')
LAI = f(xnew)


# 读取自变量为xlist,土壤水，叶面积指数，空气温度，风速，空气湿度，下行短波辐射，下行长波辐射，上行短波辐射，上行长波辐射，白天夜晚标识
meteodata = pd.read_excel(excel_file_path)
col1 = pd.Series(meteodata["Ms_2cm"].values,name = "Ms_2cm")
col2 = pd.Series(LAI,name ="LAI")
col3 = pd.Series(meteodata["Ta_5m"].values,name ="Ta_5m")
col4 = pd.Series(meteodata["WS_5m"].values,name ="WS_5m")
col5 = pd.Series(meteodata["RH_5m"].values,name ="RH_5m")
col6 = pd.Series(meteodata["DR"].values,name ="DR")
col7 = pd.Series(meteodata["DLR_Cor"].values,name ="DLR_Cor")
col8 = pd.Series(meteodata["UR"].values,name ="UR")
col9 = pd.Series(meteodata["ULR_Cor"].values,name ="ULR_Cor")
DR = meteodata["DR"]
DNS = []
for i in range(len(DR)):
	if DR[i] > 10:
		DNS.append(1)
	else:
		DNS.append(0)
col10 = pd.Series(DNS,name ="DNS")

df = pd.DataFrame({col1.name:col1, col2.name:col2,col3.name:col3,col4.name:col4,col5.name:col5,col6.name:col6,col7.name:col7,col8.name:col8,col9.name:col9,col10.name:col10})

index = []
for i in range(len(df)):
	if (i+1)%6 != 1 and (i+1)%6!=4:
		index.append(i)
xdf = df.drop(index = index)
xdf = xdf.reset_index(drop=True)
outxdf = xdf[150*48:150*48+4800]
outxdf.to_csv(x_file_path,mode = 'w',header = None,index =None,sep=',')
print(x_file_path)

IRT_1 = meteodata["IRT_1"]
IRT_2 = meteodata["IRT_2"]
IRT = []
for i in range(len(IRT_1)):
	IRT.append ((IRT_1[i] + IRT_2[i])/2+273.15)
IRT = pd.DataFrame(IRT)
IRT = IRT.drop(index = index)
IRT = IRT.reset_index(drop=True)
outIRT = IRT[150*48:150*48+4800]
outIRT.to_csv(y_file_path,mode='w',header=None,index = None,sep = ",")
print(y_file_path)





