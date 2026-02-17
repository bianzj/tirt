from osgeo import gdal, osr
from osgeo import gdal_array
import os
import time
import numpy as np

start = time.time()
'''
脚本目的:批量读取tif文件并进行几何校正,定标和投影转换,裁剪
'''
'''
第一个函数:批量读取文件夹中hdf格式的modis产品数据并进行几何校正(WGS84)
'''


def readHdfWithGeo(hdfFloder, saveFloder):
	# 获取输入文件夹中所有的文件名
	hdfNameList = os.listdir(hdfFloder)
	# 遍历文件名列表中所有文件
	for i in range(len(hdfNameList)):
		# 获取文件名后缀
		dirname, basename = os.path.split(hdfNameList[i])
		filename, txt = os.path.splitext(basename)
		# 判断文件后缀是否为 .hdf
		if txt == '.hdf':
			hdfPath = hdfFloder + os.sep + hdfNameList[i]
			# 打开hdf文件
			datasets = gdal.Open(hdfPath)
			# 打开子数据集
			dsSubDatasets = datasets.GetSubDatasets()
			# 打开ndvi数据
			laiRaster = gdal.Open(dsSubDatasets[1][0])
			# 获取元数据
			metaData = datasets.GetMetadata()
			for key, value in metaData.items():
				print('{key}:{value}'.format(key = key, value = value))
			# 获取数据时间
			time = metaData['RANGEBEGINNINGDATE']
			# 命名输出完整路径文件名
			outName = saveFloder + os.sep + time + '.tif'
			# 进行几何校正
			geoData = gdal.Warp(outName, laiRaster,
								dstSRS='EPSG:4326', format='GTiff',
								resampleAlg=gdal.GRA_Bilinear)
			del geoData
			print('{outname} deal end'.format(outname=outName))


start = time.time()

readHdfWithGeo(r'D:\data\lstSimulate\2019LAI2',
			   r'D:\data\lstSimulate\2019LAI2\25_4')
end = time.time()
print('deal spend: {s} s'.format(s=end - start))

'''
第二个函数:批量读取tif文件(适用于多波段tif文件)存入数组,并进行定标
'''


def readTifAsArray(tifPath):
	dataset = gdal.Open(tifPath)
	if dataset == None:
		print(tifPath + "文件错误")
		return tifPath

	image_datatype = dataset.GetRasterBand(1).DataType
	row = dataset.RasterYSize
	col = dataset.RasterXSize
	nb = dataset.RasterCount
	proj = dataset.GetProjection()
	gt = dataset.GetGeoTransform()

	if nb != 1:
		array = np.zeros((row, col, nb),
						 dtype=gdal_array.GDALTypeCodeToNumericTypeCode(
							 image_datatype))
		for b in range(nb):
			band = dataset.GetRasterBand(b + 1)
			nan = band.GetNoDataValue()
			array[:, :, b] = band.ReadAsArray()
	else:
		array = np.zeros((row, col),
						 dtype=gdal_array.GDALTypeCodeToNumericTypeCode(
							 image_datatype))
		band = dataset.GetRasterBand(1)
		nan = band.GetNoDataValue()
		array = band.ReadAsArray()
	return array, nan, gt, proj


'''  
第三个函数:写出tif文件
'''


def writeTiff(im_data, nan, im_geotrans, im_proj, path):
	if 'int8' in im_data.dtype.name:
		datatype = gdal.GDT_Byte
	elif 'int16' in im_data.dtype.name:
		datatype = gdal.GDT_UInt16
	else:
		datatype = gdal.GDT_Float32
	if len(im_data.shape) == 3:
		im_bands, im_height, im_width = im_data.shape
	elif len(im_data.shape) == 2:
		im_data = np.array([im_data])
		im_bands, im_height, im_width = im_data.shape

	driver = gdal.GetDriverByName("GTiff")
	dataset = driver.Create(path, int(im_width), int(im_height), int(im_bands), datatype)

	if (dataset != None):
		dataset.SetGeoTransform(im_geotrans)  # 写入仿射变换参数
		dataset.SetProjection(im_proj)  # 写入投影
	for i in range(im_bands):
		dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
		outBand = dataset.GetRasterBand(i + 1)
		outBand.FlushCache()
		outBand.SetNoDataValue(nan)
	del dataset


'''
第四个函数:批量进行定标计算
'''


def cal(tifFloder, saveFloder, scale):
	tifNameList = os.listdir(tifFloder)

	for i in range(len(tifNameList)):

		filename, txt = os.path.splitext(tifNameList[i])

		if txt == '.tif':
			tifPath = tifFloder + os.sep + tifNameList[i]
			array = readTifAsArray(tifPath)
			outName = saveFloder + os.sep + filename + '_cal' + txt
			calData = writeTiff(array[0] / scale, array[1] / scale,
								array[2], array[3], outName)
			print('{filename} deal end '.format(filename=outName))
			del calData


# start = time.time()
# cal(r'F:\GDAL learning\MODIS\MOD13A3\NDVItif_withGeo',
# 	r'F:\GDAL learning\MODIS\MOD13A3\NDVItif_cal', 10000)
# end = time.time()
# print('deal spend: {s} s'.format(s=end - start))

'''
第五个函数:批量进行重投影
'''


def reproject(tifFloder, saveFloder, proj4):
	tifNameList = os.listdir(tifFloder)

	srs = osr.SpatialReference()
	srs.ImportFromProj4(proj4)

	for i in range(len(tifNameList)):

		filename, txt = os.path.splitext(tifNameList[i])

		if txt == '.tif':
			tifPath = tifFloder + os.sep + tifNameList[i]
			outName = saveFloder + os.sep + filename + '_rep' + txt
			reproData = gdal.Warp(outName, tifPath, dstSRS=srs,
								  xRes=1000, yRes=1000,
								  resampleAlg=gdal.GRA_Bilinear,
								  outputType=gdal.GDT_Float32)
			print('{filename} deal end '.format(filename=outName))
			del reproData


# start = time.time()
# reproject(r'F:\GDAL learning\MODIS\MOD13A3\NDVItif_cal',
# 		  r'F:\GDAL learning\MODIS\MOD13A3\NDVItif_albers',
# 		  '+proj=aea +lat_0=0 +lon_0=105 +lat_1=25 +lat_2=47 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
# end = time.time()
# print('deal spend: {s} s'.format(s=end - start))

'''
第六个函数:批量进行裁剪
'''


def cut(tifFloder, saveFloder, shpFile):
	tifNameList = os.listdir(tifFloder)

	for i in range(len(tifNameList)):

		filename, txt = os.path.splitext(tifNameList[i])

		if txt == '.tif':
			tifPath = tifFloder + os.sep + tifNameList[i]
			outName = saveFloder + os.sep + filename + '_cut' + txt
			cut_ds = gdal.Warp(outName, tifPath,
							   cutlineDSName=shpFile,
							   cropToCutline=True)
			print('{filename} deal end '.format(filename=outName))
			del cut_ds


# start = time.time()
# cut(r'F:\GDAL learning\MODIS\MOD13A3\NDVItif_albers',
# 	r'F:\GDAL learning\MODIS\MOD13A3\NDVItif_wrap',
# 	r'F:\GDAL learning\shp\BeiJingBei.shp')
# end = time.time()
# print('deal spend: {s} s'.format(s=end - start))
