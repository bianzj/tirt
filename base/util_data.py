import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from base.util_packages import *


'''util_data'''

def resize_data(preArray, nl, ns, method = 'cv2.INTER_NEAREST'):
    '''
    数组的缩放，用到了cv2库，方法主要有：
    ‘最邻近‘cv2.INTER_NEAREST，‘双线性插值’cv2.INTER_LINEAR，‘三次立方’cv2.INTER_CUBIC，’等面积’cv2.INTER_AREA
    :param preArray: 原有数组
    :param nl: 目标行数
    :param ns: 目标列数
    :param method:  缩放方法,默认是双线性插值
    :return: 返回缩放后的数组
    '''
    ns = np.int_(ns)
    nl = np.int_(nl)
    data = cv2.resize(preArray,(ns,nl),interpolation=0)
    return data

def resize_data_ratio(preArray, ratio = 0.5, method = 'cv2.INTER_LINEAR'):
    '''
    并不是给定行列号，而是行列号的比例确定新数据的行列号
    :param preArray: 原有数组
    :param ratio: 转换比例
    :param method:  缩放方法，默认是双线性插值
    :return: 新数组
    '''
    [pre_nl,pre_ns] = np.shape(preArray)
    ns = np.int_(pre_ns*ratio)
    nl = np.int_(pre_nl*ratio)
    data = cv2.resize(preArray,(ns,nl),interpolation=method)
    return data

def getPointfromImage(data, imagex_, imagey_, dist= 0, dn = 3,  minThreshold = -100, maxThreshold = 100):
    '''
    从图像中找点位对应的值，如有必要进行简单的统计
    :param data: 图像数据
    :param imagex_: 行列号，x坐标
    :param imagey_: 行列号，y坐标
    :param dist: 距离，默认0
    :param dn: 绝对值大小，默认3倍
    :param minThreshold: 绝对最小阈值
    :param maxThreshold: 绝对最大阈值
    :return: 图像值或统计结果数组
    '''
    size = np.size(imagex_)
    result = np.zeros(size)
    for k in range(size):
        imagex = imagex_[k]
        imagey = imagey_[k]
        x1 = imagex - dist
        x2 = imagex + dist + 1
        y1 = imagey - dist
        y2 = imagey + dist + 1
        temp = data[x1:x2, y1:y2]
        ### 绝对大小阈值判断
        ind = (temp > minThreshold) * (temp < maxThreshold)
        if np.sum(ind) < 1: continue
        ### 相对大小阈值判断
        std = np.std(temp[ind])
        ave = np.average(temp[ind])
        indd = (temp > minThreshold) * (temp < maxThreshold) *(temp < ave+dn*std) *(temp > ave-dn*std)
        if np.sum(indd) < 1: continue
        result[k] = np.average(temp[indd])
    return result


'''
util_map
'''
def getSRSPair(dataset):
    '''
    获得给定数据的投影参考系和地理参考系
    :param dataset: GDAL地理数据
    :return: 投影参考系和地理参考系
    '''
    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    return prosrs, geosrs

def geo2lonlat(dataset, x, y):
        '''
        将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
        :param dataset: GDAL地理数据
        :param x: 投影坐标x
        :param y: 投影坐标y
        :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
        '''
        prosrs, geosrs = getSRSPair(dataset)
        ct = osr.CoordinateTransformation(prosrs, geosrs)
        x = np.reshape(x, [-1])
        y = np.reshape(y, [-1])
        temp = np.asarray([x, y])
        temp = np.transpose(temp)
        coords = np.asarray(ct.TransformPoints(temp))
        return coords[:,0],coords[:,1]

def lonlat2geo(dataset,lat,lon):
    '''
        将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定），尤其注意，不同版本经度和纬度是反的
        :param dataset: GDAL地理数据
        :param lon: 地理坐标lon经度
        :param lat: 地理坐标lat纬度
        :return: 经纬度坐标(lon, lat)对应的投影坐标
    '''
    # dataset = gdal.Open(fileName, gdal.GA_ReadOnly)
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    lon = np.reshape(lon,[-1])
    lat = np.reshape(lat,[-1])
    temp = np.asarray([lon,lat])
    temp = np.transpose(temp)
    # temp = np.asarray([lat[0:2],lon[0:2]])
    coords = np.asarray(ct.TransformPoints(temp))

    return coords[:,0],coords[:,1]

def geo2imagexy(dataset, x, y):
    '''
    根据GDAL的六 参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
    :param dataset: GDAL地理数据
    :param x: 投影或地理坐标x
    :param y: 投影或地理坐标y
    :return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col) nl ns
    '''
    trans = dataset.GetGeoTransform()
    a = np.array([[trans[2], trans[1]], [trans[5], trans[4]]])
    b = np.array([x - trans[0], y - trans[3]])
    return np.linalg.solve(a, b)  # 使用numpy的linalg.solve进行二元一次方程的求解

def imagexy2geo(dataset, row,col):
    '''
    根据GDAL的六参数模型将影像图上坐标（行列号）转为投影坐标或地理坐标（根据具体数据的坐标系统转换）
    :param dataset: GDAL地理数据
    :param row: 像素的行号
    :param col: 像素的列号
    :return: 行列号(row, col)对应的投影坐标或地理坐标(x, y)
    '''
    trans = dataset.GetGeoTransform()
    px = trans[0] + col * trans[1] + row * trans[2]
    py = trans[3] + col * trans[4] + row * trans[5]
    return px, py

def reproj_image_gdal(dataset_destination, dataset_source):
    '''
    源数据集重投影，提取与目标数据集对应的数据，数据来自源数据，结果的坐标与投影与目标数据一致
    :param dataset_destination: 目标数据集
    :param dataset_source: 源数据集
    :return: 数据集，其数据来自源，坐标与投影与目标一致
    '''
    ns1 = dataset_destination.RasterXSize
    nl1 = dataset_destination.RasterYSize
    nb1 = dataset_destination.RasterCount
    # data1,ns1,nl1,nb1,trans1,proj1,dataset1 = read_image_gdal_dataset(infile1)
    xi = np.zeros([nl1,ns1])
    yi = np.zeros([nl1,ns1])
    temp = np.linspace(0,nl1-1,nl1)
    for k in range(ns1):
        xi[:,k] = temp
        yi[:,k] = k
    temp11,temp12 = imagexy2geo(dataset_destination, xi, yi)
    lon,lat = geo2lonlat(dataset_destination, temp11, temp12)
    ns2 = dataset_source.RasterXSize
    nl2 = dataset_source.RasterYSize
    nb2 = dataset_source.RasterCount
    data2 = dataset_source.ReadAsArray(0, 0, ns2, nl2)
    temp21, temp22 = lonlat2geo(dataset_source, lat, lon)
    imagex, imagey = geo2imagexy(dataset_source, temp21, temp22)
    imagex = np.asarray(imagex+0.5, np.int_)
    imagey = np.asarray(imagey+0.5, np.int_)
    imagex[imagex < 0] = 0
    imagey[imagey < 0] = 0
    imagex[imagex >= nl2-1] = nl2-1
    imagey[imagey >= ns2-1] = ns2-1
    if nb2 <=1:
        data = np.zeros([nl1, ns1])
        data[:] = np.reshape(data2[imagex,imagey],[nl1,ns1])
    else:
        data = np.zeros([nb2,nl1, ns1])
        for k in range(nb2):
            temp = data2[k,:,:]*1.0
            data[k,:] =  np.reshape(temp[imagex,imagey],[nl1,ns1])
    return data

def reproj_image_ref_gdal(dataset_destination, dataset_reff, dataset_source):
    '''
    数据重投影到新区域，原数据并不直接到目标数据，原数据先数值变换到中间参考，然后通过中间参考投影信息到目标投影信息
    :param dataset_destination: 目标数据
    :param dataset_reff: 中间参考
    :param dataset_source: 原始数据
    :return: 新图像，其坐标和投影信息来自目标，数据信息来自源数据
    '''
    ns1 = dataset_destination.RasterXSize
    nl1 = dataset_destination.RasterYSize
    nb1 = dataset_destination.RasterCount
    # data1,ns1,nl1,nb1,trans1,proj1,dataset1 = read_image_gdal_dataset(infile1)
    xi = np.zeros([nl1,ns1])
    yi = np.zeros([nl1,ns1])
    temp = np.linspace(0,nl1-1,nl1)
    for k in range(ns1):
        xi[:,k] = temp
        yi[:,k] = k
    temp11,temp12 = imagexy2geo(dataset_destination, xi, yi)
    lon,lat = geo2lonlat(dataset_destination, temp11, temp12)
    ns2 = dataset_reff.RasterXSize
    nl2 = dataset_reff.RasterYSize
    nb2 = dataset_reff.RasterCount
    ns3 = dataset_source.RasterXSize
    nl3 = dataset_source.RasterYSize
    nb3 = dataset_source.RasterCount
    temp = dataset_source.ReadAsArray(0, 0, ns3, nl3)
    data2 = np.zeros([nb3,nl2, ns2])
    for k in range(nb3):
        data2[k,:,:] = resize_data(temp[k,:,:],nl2,ns2)
    temp21, temp22 = lonlat2geo(dataset_reff, lat, lon)
    imagex, imagey = geo2imagexy(dataset_reff, temp21, temp22)
    imagex = np.asarray(imagex+0.5, np.int_)
    imagey = np.asarray(imagey+0.5, np.int_)

    imagex[imagex < 0] = 0
    imagey[imagey < 0] = 0
    imagex[imagex >= nl2-1] = nl2-1
    imagey[imagey >= ns2-1] = ns2-1
    if nb3 <=1:
        data = np.zeros([nl1, ns1])
        data[:] = np.reshape(data2[imagex,imagey],[nl1,ns1])
    else:
        data = np.zeros([nb3,nl1, ns1])
        for k in range(nb3):
            temp = data2[k,:,:]*1.0
            data[k,:] =  np.reshape(temp[imagex,imagey],[nl1,ns1])
    return data

def loc2map(lat, lon, dataset):
    '''
    经纬度转图上行列号
    :param lat: 纬度
    :param lon: 经度
    :param dataset: 数据集，提供了地理和投影信息
    :return: 行列号
    '''
    ns2 = dataset.RasterXSize
    nl2 = dataset.RasterYSize
    nb2 = dataset.RasterCount
    data2 = dataset.ReadAsArray(0, 0, ns2, nl2)
    temp21, temp22 = lonlat2geo(dataset, lat, lon)
    imagex, imagey = geo2imagexy(dataset, temp21, temp22)
    # imagex = np.asarray(imagex+0.5, np.int_)
    # imagey = np.asarray(imagey+0.5, np.int_)
    imagex = np.asarray(imagex, np.int_)
    imagey = np.asarray(imagey, np.int_)
    return imagex,imagey

def calc_azimuth(lat1, lon1, lat2, lon2):
    '''
    计算地球上两点间的相对方位角
    :param lat1: 点1的纬度
    :param lon1: 点1的经度
    :param lat2: 点2的纬度
    :param lon2: 点2的经度
    :return: 相对方位角
    '''
    lat1_rad = lat1 * np.pi / 180
    lon1_rad = lon1 * np.pi / 180
    lat2_rad = lat2 * np.pi / 180
    lon2_rad = lon2 * np.pi / 180
    y = np.sin(lon2_rad - lon1_rad) * np.cos(lat2_rad)
    x = np.cos(lat1_rad) * np.sin(lat2_rad) - \
        np.sin(lat1_rad) * np.cos(lat2_rad) * np.cos(lon2_rad - lon1_rad)
    brng = np.arctan2(y, x) * 180 / np.pi
    return float((brng + 360.0) % 360.0)

def write_vrt(vrtfile, datafile, xfile, yfile, xscale, yscale, noDataValue = 65535, bandnumber = 1):
    '''
    导出用于gdal vrt的配置文件
    # 与函数配合使用，outfile为转成的文件路径，resampleAlg为重采样的
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_NearestNeighbour)
    :param vrtfile:  该配置文件名称
    :param datafile:  所要操作的文件名称
    :param xfile:  经度坐标文件 lon
    :param yfile:  纬度坐标文件 lat
    :param xscale:  列数
    :param yscale:  行数
    :param noDataValue: 无值填补
    :param bandnumber:  波段数据
    :return: 无
    '''
    f = open(vrtfile, 'w')
    f.write(r'<VRTDataset rasterXSize="%d"'%xscale + 'rasterYSize="%d">'%yscale)
    f.write('\n')
    f.write(r'  <Metadata domain="GEOLOCATION">')
    f.write('\n')
    f.write(r'    <MDI key="SRS">GEOGCS["WGS 84(DD)",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AXIS["Long",EAST],AXIS["Lat",NORTH]]</MDI>')
    # f.write(r'    <MDI key="SRS">GEOGCS["WGS 84(DD)",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG", "9122"],AUTHORITY["EPSG", "4490"]</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="X_DATASET">' + xfile + r'</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="X_BAND">1</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="PIXEL_OFFSET">0</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="PIXEL_STEP">1</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="Y_DATASET">' + yfile + r'</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="Y_BAND">1</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="LINE_OFFSET">0</MDI>')
    f.write('\n')
    f.write(r'    <MDI key="LINE_STEP">1</MDI>')
    f.write('\n')
    f.write(r'  </Metadata>')
    f.write('\n')
    f.write(r'  <VRTRasterBand dataType = "Float32" band = "%d">'%bandnumber)
    f.write('\n')
    f.write(r'    <ColorInterp>Gray</ColorInterp >')
    f.write('\n')
    f.write(r'    <NoDataValue>%d</NoDataValue >'%noDataValue)
    f.write('\n')
    f.write(r'    <SimpleSource>')
    f.write('\n')
    f.write(r'      <SourceFilename relativeToVRT = "1" >' + datafile + r'</SourceFilename>')
    f.write('\n')
    f.write(r'      <SourceBand>1</SourceBand>')
    f.write('\n')
    f.write(r'    </SimpleSource>')
    f.write('\n')
    f.write(r'  </VRTRasterBand>')
    f.write('\n')
    f.write('</VRTDataset>')
    f.close()
    return 1

def eliminate_edge(data,edge=1):
    '''
    消除重采样后边界的异常值
    :param data: 数据
    :param edge: 要剔除数据的步长
    :return:  剔除边界的新数据
    '''
    dataa = data*1.0
    [nl,ns] = np.shape(data)
    data1 = np.ones([nl,ns])
    data2 = np.ones([nl,ns])
    data1[:,:ns-edge] = data[:,edge:ns]
    data2[:,edge:ns] = data[:,:ns-edge]
    temp = (data1)*data*data2
    dataa[temp ==0] = 0
    return dataa


'''util_time'''
def date2DOY(year,month,day):
    '''
    日期到Day of Year 的转换
    :param year: 年
    :param month:  月
    :param day: 日
    :return: doy
    '''
    days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    monthsum = np.zeros(12)
    year = np.int_(year)
    month = np.int_(month)
    day = np.int_(day)
    if isinstance(day * 1.0, float):
        total = 0
        for index in range(month - 1):
            total += days_of_month[index]
        temp = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
        if month > 2 and temp:
            total += 1
        return total + day
    else:
        for index in range(1,12):
            monthsum[index] = monthsum[index-1]+days_of_month[index-1]
        month = np.asarray(month,dtype=np.int_)
        DOY = monthsum[month-1] + day
        ind = ((year % 4 == 0) * (year % 100 != 0)) * (month>2)
        DOY[ind] = DOY[ind] + 1
        ind = ((year % 400 == 0))
        DOY[ind] = DOY[ind] + 1
    return DOY

def doy2date(year, doy):
    '''
    doy 到日期的转换
    :param year: 年
    :param doy: day of year
    :return:  日期（月,日）
    '''
    year = np.int_(year)
    doy = np.int_(doy)
    month_leapyear = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_notleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        for i in range(0, 12):
            if doy > month_leapyear[i]:
                doy -= month_leapyear[i]
                continue
            if doy <= month_leapyear[i]:
                month = i + 1
                day = doy
                break
    else:
        for i in range(0, 12):
            if doy > month_notleap[i]:
                doy -= month_notleap[i]
                continue
            if doy <= month_notleap[i]:
                month = i + 1
                day = doy
                break
    return month, day