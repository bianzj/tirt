import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from base.util_packages import *




'''util_file'''
def search(dir):
    '''
    查找路径下所有的文件和文件夹
    :param dir: 路径
    :return: 文件和文件夹名字
    '''
    results = os.listdir(dir)
    return results

def search_file(dir,specstr):
    '''
    查找某个路径下，含有特定标识的文件，需要list把标识括起来，即使只有1个
    :param dir: 路径
    :param specstr: 特定的标识
    :return: 文件名
    '''
    results = []
    num = np.size(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        specstr0 = specstr[0]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) and
                    specstr0 in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x]
    return results

def search_dir(dir,specstr):
    '''
    查找某个路径下，含有特定标识的文件夹，需要list把标识括起来，即使只有1个
    :param dir: 路径
    :param specstr:
    :return:
    '''
    results = []
    num = len(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) and
                    specstr in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x]
    return results

def search_file_rej(dir,specstr,rejstr):
    '''
    查找某个路径下，含有特定标识,又不含有某个标识的文件，需要list把标识括起来，即使只有1个，可以有多个“有”标识，但是只有1个“没有”标识
    :param dir: 路径
    :param specstr: “有”标识
    :param rejstr: “没有”标识
    :return: 文件名数组
    '''
    results = []
    num = len(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        specstr0 = specstr[0]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x))
                    and specstr0 in x
                    and rejstr not in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and rejstr not in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x
                    and rejstr not in x]
    return results

def search_dir_rej(dir,specstr,rejstr):
    '''
    查找某个路径下，含有特定标识,又不含有某个标识的文件夹，需要list把标识括起来，即使只有1个，可以有多个“有”标识，但是只有1个“没有”标识
    :param dir: 路径
    :param specstr: “有”标识
    :param rejstr: “没有”标识
    :return: 文件夹名数组
    '''
    results = []
    num = len(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        specstr0 = specstr[0]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x))
                    and specstr0 in x
                    and rejstr not in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and rejstr not in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x
                    and rejstr not in x]
    return results

def remove_file(dir, specstr):
    '''
    删除某个路径下，含有某个标识的所有文件
    :param dir: 路径
    :param specstr:标识
    :return: 没有
    '''
    for x in os.listdir(dir):
        fp = os.path.join(dir, x)
        # 如果文件存在，返回true
        if re.search(specstr, x) is not None:
            print(fp)
            os.remove(fp)

def rename_file(dir, specstr):
    '''
    对路径下，含有某个标识的文件批量修改名称
    :param dir: 路径
    :param specstr: 标识
    :return: 没有
    '''
    for x in os.listdir(dir):
        fp = os.path.join(dir, x)
        # print(x)
        # 如果文件存在，返回true
        if re.search(specstr, x) is not None:
            [filename, hz] = os.path.splitext(x)
            outfile = dir + filename + '_test.tif'
            print(fp)
            print(outfile)
            if os.path.exists(outfile) == 1:
                os.remove(outfile)
            os.rename(fp, outfile)

def move_file(infile,outfile):
    '''
    修改名称，其实是修改路径
    :param infile: 旧名称
    :param outfile: 新名称
    :return: 没有
    '''
    os.rename(infile,outfile)

def read_txt_float(filename):
    '''
    读取txt文件，并将其按照float存储，行列号作为索引的矩阵
    :param filename: 文件路径
    :return: float数组
    '''
    mydata = []
    with open(filename) as f:
        lines = f.readline()
        while lines:
            line = lines.split()
            mydata.append(line)
            lines = f.readline()
    mydata = np.asarray(mydata,dtype=float)
    return mydata

def read_txt_str(filename, spt):
    '''
    读取文件数据，然后按照特定的约束进行解析，随后按照行列进行存储
    :param filename: 文件名称
    :param spt: 解析表示
    :return: 解析后数组
    '''
    mydata = []
    with open(filename) as f:
        lines = f.readline()
        while lines:
            line = lines.split(spt)
            mydata.append(line)
            lines = f.readline()
    return mydata

def read_txt_array(filename, num_pass, num_col):
    '''
    限定读取txt文本数据，考虑了要跳过的行，但是必须要给定列数
    :param filename: 文件名
    :param num_pass: 要跳过的行数
    :param num_col: 数据的列数
    :return: 返回数组
    '''
    f = open(filename, 'r')
    temp = f.readlines()
    temp = np.asarray(temp)
    temp = temp[num_pass:]
    num = len(temp)
    lut = np.zeros([num, num_col])
    for k in range(num):
        if (temp[k] == ""): continue
        tempp = (re.split(r'\s+', temp[k].strip()))
        lut[k, :] = np.asarray(tempp)
    return lut

def read_excel_sheet(filename,sheetname='Sheet1'):
    '''
    读取excel数据，返回册数据集，这里用了xlrd
    :param filename:文件名
    :param sheetname: 册名
    :return: 册数据集
    '''
    ExcelFile = xlrd.open_workbook(filename)
    ExcelFile.sheet_names()
    sheet = ExcelFile.sheet_by_name(sheetname)
    return sheet

def read_excel_sheet_col(filename,col,sheetname='Sheet1'):
    '''
    读取excel数据，返回某个册的某列数据，这里用了pandas
    :param filename: 文件名
    :param col: 某列
    :param sheetname:册名
    :return: 列数组
    '''
    df = pd.read_excel(filename,sheet_name=sheetname)
    colvalue = df.ix[:,col]
    return colvalue

def read_excel_sheet_row(filename,row,sheetname='Sheet1'):
    '''
    读取excel数据，然后某册某行数据，这里用了pandas
    :param filename: 文件名
    :param row: 某行
    :param sheetname:某册
    :return: 行数组
    '''
    ExcelFile = xlrd.open_workbook(filename)
    ExcelFile.sheet_names()
    sheet = ExcelFile.sheet_by_name(sheetname)
    rowvalue = sheet.row_values(row)
    return rowvalue

def read_binary(filename, type = np.int_):
    '''
    读取二进制文件，并转成特定格式，默认是整型
    :param filename:文件名
    :param type: 数据格式
    :return: 数组
    '''
    fin = open(filename, 'rb')
    temp =[]
    while True:
        fileContent = fin.read(4)
        num = len(fileContent)
        if num !=4:
            break
        tp = struct.unpack('l',fileContent)
        temp.extend(tp)
    fin.close()
    temp = np.array(temp,dtype=type)
    return temp

def read_image_gdal(filename):
    '''
    使用gdal,读取图像文件，并格式输出
    :param filename:文件名
    :return: 文件，列数，行数，波段数，地理信息，投影信息
    '''
    dataset = gdal.Open(filename)
    if dataset == None:
        print(filename + "文件无法打开")
        return
    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize
    im_bands = dataset.RasterCount
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
    im_geotrans = dataset.GetGeoTransform()
    im_proj = dataset.GetProjection()
    return im_data,im_width,im_height,im_bands,im_geotrans,im_proj

def read_image_dataset_gdal(filename):
    '''
    读取图像文件，特别的，会返回数据集
    :param filename: 文件名
    :return: 文件，列数，行数，波段数，地理信息，投影信息，数据集
    '''
    dataset = gdal.Open(filename)
    if dataset == None:
        print(filename + "文件无法打开")
        return
    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize
    im_bands = dataset.RasterCount
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
    im_geotrans = dataset.GetGeoTransform()
    im_proj = dataset.GetProjection()
    return im_data,im_width,im_height,im_bands,im_geotrans,im_proj,dataset

def read_image_raw(filename, ns, nl, nb, type = np.int_):
    '''
    在给定行列号的情况下，读取二进制文件，需要指定数据类型，确定每次读取的步长
    :param filename:文件名
    :param ns:列数
    :param nl: 行数
    :param nb:波段数
    :param type:数据类型
    :return: 数组
    '''
    fb = open(filename, 'rb')
    mydata = np.zeros((ns,nl,nb))
    for kb in range(nb):
        for ks in range(ns):
            for kl in range(nl):
                if type == np.int_:
                    arr = fb.read(2)
                else:
                    arr = fb.read(4)
                elem = struct.unpack('h',arr)[0]
                mydata[ks][kl][kb] = elem
    return mydata

def read_image_Nc_group(fileName, groupName, objectName,ifscale=0):
    '''
    读取NC格式的数据，需要指定是否需要缩放转换
    :param fileName:文件名
    :param groupName: 组名
    :param objectName: 目标名
    :param ifscale: 是否缩放
    :return: 数组
    '''
    dataset = netCDF4.Dataset(fileName)
    if ifscale ==0:
        dataset.groups[groupName].variables[objectName].set_auto_maskandscale(False)
    predata = np.asarray(dataset.groups[groupName].variables[objectName][:])
    return predata

def read_image_Nc(fileName, objectName,ifscale):
    '''
    读取NC格式的数据，没有组名，也需要指定是否需要缩放
    :param fileName:  文件名
    :param objectName:  目标名
    :param ifscale:  是否缩放
    :return:  数组
    '''
    dataset = netCDF4.Dataset(fileName)
    if ifscale ==0:
        dataset.variables[objectName].set_auto_maskandscale(False)
    predata = np.asarray(dataset.variables[objectName][:])
    return predata

def getDatafromNc(fileName, objectName):
    '''
    读取NC文件，文件名和项目名
    :param fileName: 文件名
    :param objectName: 项目名
    :return: 值
    '''
    try:
        dataset = netCDF4.Dataset(fileName)
    except (OSError, RuntimeError):
        print(fileName + ' File Error')
        predata = np.asarray(None)
    else:
        try:
            predata = np.asarray(dataset.variables[objectName][:])
        except RuntimeError:
            print(fileName + ' Attribute Error')
            predata = np.asarray(None)
        finally:
            dataset.close();
    return predata

def write_image_gdal(im_data, im_width, im_height, im_bands, im_trans, im_proj, path, imageType = 'GTiff'):
    '''
    保存数据，输出成tif格式
    :param im_data: 数组
    :param im_width:  列数
    :param im_height:  行数
    :param im_bands:  波段数
    :param im_trans:  地理坐标
    :param im_proj:  投影坐标
    :param path:  路径
    :return: 无
    '''
    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'uint16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_Int16
    else:
        datatype = gdal.GDT_Float32
    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    elif len(im_data.shape) == 2:
        im_data = np.array([im_data])
    else:
        im_bands, (im_height, im_width) = 1, im_data.shape
    driver = gdal.GetDriverByName(imageType)
    dataset = driver.Create(path, im_width, im_height, im_bands, datatype)
    if (dataset != None and im_trans != '' and im_proj != ''):
        dataset.SetGeoTransform(im_trans)
        dataset.SetProjection(im_proj)
    for i in range(im_bands):
        dataset.GetRasterBand(i+1).WriteArray(im_data[i])
    del dataset

def read_dataset_gdal(fileName):
    '''
    读取数据，返回数据集，多用于基于数据集的数据转换
    :param fileName: 数据名
    :return: 数据集
    '''
    dataset = gdal.Open(fileName)
    return dataset

def mkdir_if_not_exist(dir_name, is_delete=False):
    """
    创建文件夹
    :param dir_name: 文件夹
    :param is_delete: 是否删除
    :return: 是否成功
    """
    try:
        if is_delete:
            if os.path.exists(dir_name):
                shutil.rmtree(dir_name)
                print("Path exists, remove and recreate")

        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            print("Created path:" + dir_name)
        return True
    except Exception as e:
        return False

def un_gz(file_name,infile,outdir):
    '''
    解压gz文件
    :param file_name:要解压文件名称
    :param infile:要解压文件位置
    :param outdir:输出位置
    '''
    mkdir_if_not_exist(outdir)
    f_name = outdir + file_name.replace(".gz", "")
    # 获取文件的名称，去掉
    g_file = gzip.GzipFile(infile + file_name)
    # 创建gzip对象
    open(f_name, "wb+").write(g_file.read())
    # gzip对象用read()打开后，写入open()建立的文件里。
    g_file.close()  # 关闭gzip对象
    print("Unzip " + file_name + " finished!")

