import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from base.util_packages import *


'''util_plot'''

'''
plot with plt and seaborn
plt 的用法，有多少列数据，就是多少列
seaborn 的用法，有很多列数据，但是只有1列表现出来，用标签标识不同的列
'''
def plt_scatter(data1,data2,title='',dif=10,min1 = 225,min2 = 225,max1 = 335,max2 = 335):
    '''
    画散点图，并进行数据1和数据2 的统计信息
    :param data1: 数据1
    :param data2: 数据2
    :param title: 名称
    :param dif:  偏差阈值
    :param min1:  数据1的最小值
    :param min2:  数据2的最小值
    :param max1:  数据1的最大值
    :param max2:  数据2的最大值
    :return:  无
    '''
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['figure.figsize'] = (4.0, 3.2)
    ind = (data1>min1)*(data1<max1)*(data2>min2)*(data2<max2)*(np.abs(data1-data2)<dif)
    ind = np.where(ind > 0)
    dif = data2[ind]-data1[ind]
    rmse = np.sqrt(np.mean(dif*dif))
    std = np.std(dif)
    bias = np.mean(dif)
    r = np.corrcoef(data1[ind],data2[ind])
    r2 = r[0,1]*r[0,1]
    plt.figure(figsize=[4.5,4])
    plt.plot(data1[ind],data2[ind],'ko',markersize=2.5)
    plt.plot([min1,max1],[min2,max2],'k-.')
    plt.title([rmse,r2])
    plt.xlim([min1,max1])
    plt.ylim([min2,max2])
    plt.text(min1+(max1-min1)*3/5, min2+(max2-min2)*1/5, "$RMSE$ = %2.2f$\\degree$C" % rmse + "\n$Bias$ = %2.2f$\\degree$C" % bias +
             "\n$r^2$ = %2.2f" % r2 + "\n$\\sigma$ = %2.2f$\\degree$C" % std, fontsize=14)
    plt.title(title)
    plt.show()
    return 0

def plt_hist(data,bins = 30, alpha = 0.5):
    '''
    单个数据的分布图
    :param data:数据
    :param bins:大小
    :param alpha:透过率
    :return: 无
    '''
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['figure.figsize'] = (4.0, 3.2)
    kwargs = dict(histtype='stepfilled', alpha=alpha, bins=bins)
    fig, axs = plt.subplots(ncols=1, figsize=(5, 4))
    plt.hist(data, **kwargs, label='$\\Delta$', color='orange')
    plt.legend()
    plt.xlabel('Difference or Value', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.show()

def plt_hist_2col(data1,data2,bins = 30, alpha = 0.5):
    '''
    数据的分布图,有两列数据，但是单独计算
    :param data1:数据1
    :param data2:数据2
    :param bins:大小
    :param alpha:透过率
    :return: 无
    '''
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['figure.figsize'] = (4.0, 3.2)
    data = np.transpose(np.stack([data1,data2]))
    kwargs = dict(histtype='stepfilled', alpha=alpha, bins=bins)
    fig, axs = plt.subplots(ncols=1, figsize=(5, 4))
    plt.hist(data, **kwargs, label='$\\Delta$', color='orange')
    plt.legend()
    plt.xlabel('Difference or Value', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.show()

def plt_hist_scatter(data1,data2,bins = 30, alpha = 0.5):
    '''
    数据的分布图,两个数据的相关散点图的密度图
    :param data1:数据1
    :param data2:数据2
    :param bins:大小
    :param alpha:透过率
    :return: 无
    '''
    fig, axs = plt.subplots(ncols=1, figsize=(5, 4))
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['figure.figsize'] = (4.0, 3.2)
    plt.hist2d(data1, data2, norm=LogNorm(), cmap='jet',histtype='stepfilled', alpha=alpha, bins=bins)
    plt.legend()
    plt.xlabel('Difference or Value', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.plot([0, 1], [0, 1], 'k-')
    colorbar()
    plt.show()

def plt_pie(labels,sizes,colors):
    '''
    饼图的画法，这里是个例子
    :param labels: 饼图的图例
    :param sizes: 饼图的尺寸
    :param colors: 饼图的颜色
    :return: 无
    '''
    labers = ['Raytran', 'DART', 'RGM', 'FLIGHT', 'librat', 'FliES']
    sizes = [131, 200, 93, 180, 84, 72]
    # sizes2 = [3789, 5138, 2232, 5623, 2702, 1476]
    colors = ['red', 'blue', 'green', 'yellow', 'cyan', 'orange']
    explode = 0, 0, 0, 0, 0, 0
    patches, l_text, p_text = plt.pie(sizes, explode=explode, labels=labers,
           colors=colors, autopct='%1.1f%%', shadow=False, startangle=50)
    plt.axis('equal')
    plt.show()

def plt_coutourPolar(theta,rho,z,n,maxrho=50,step = 5):
    '''
    极坐标图的画法
    :param theta: 极坐标图的轴角
    :param rho: 极坐标图的轴长
    :param z: 各值的大小
    :param n: 值的区间数
    :param maxrho: 最大的轴长的阈值
    :return:无
    '''

    # ###transform data to Cartesian coordinates.
    # delta = maxrho/50
    # xx = rho*np.cos((90-theta)*(np.pi/180))
    # yy = rho*np.sin((90-theta)*(np.pi/180))
    # xi = np.linspace(-maxrho,maxrho,2*maxrho/delta)
    # yi = xi
    # [xi,yi] = np.meshgrid(xi,yi)
    # zi = griddata((xx,yy),z,(xi,yi),'cubic')
    # # fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    # plt.contourf(xi,yi,zi,n,cmap='jet')
    # plt.show()

    delta = maxrho/50
    xx = np.radians(theta)
    yy = rho
    xi = np.radians(np.arange(0,365,step))
    yi = np.arange(0,maxrho,step)
    [xi,yi] = np.meshgrid(xi,yi)
    # xi = np.radians(xi)
    zi = griddata((xx,yy),z,(xi,yi))
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    plt.autumn()
    cax = ax.contourf(xi, yi, zi, n,cmap='jet')
    plt.autumn()
    cb = fig.colorbar(cax)
    # cb.set_label("Pixel reflectance") #颜色条的说明
    plt.show()

def sns_reg(data1,data2,min1 = 250,min2 = 350,max1 = 350, max2 = 350):
    '''
    画相关关系图, 同时画线性关系的线，及其不确定性，相关关系是该图的重点,点的大小和颜色并不是显示的重点，使用了seaborn库
    :param data1: 数据1
    :param data2: 数据2
    :param min1:  最小
    :param min2:  最小
    :param max1: 最大
    :param max2: 最大
    :return:
    '''
    ####################################
    #### statistical result
    ####################################
    ind = (data2 > 0)*(data1>0)
    dif = data2[ind] - data1[ind]
    rmse = np.sqrt(np.mean(dif * dif))
    std = np.std(dif)
    bias = np.mean(dif)
    r = np.corrcoef(data1[ind], data2[ind])
    r2 = r[0, 1] * r[0, 1]
    slope, intercept, r_value3, p_value3, std_err = stats.linregress(data1[ind],data2[ind])
    ####################################
    ### core code
    ####################################
    figure(figsize=[4.5, 3.2])
    sns.set_theme(style="white")
    ax = sns.regplot(x=data1, y=data2, ci=95, label='XXX')
    #ax = sns.regplot(x=data1, y=data3, ci=95, label='XXX')
    plt.tight_layout() ### 图像紧致
    plt.legend()
    plt.xlim([min1,max1])
    plt.ylim([min2,max2])
    plt.text(min1+(max1-min1)*3/5, min2+(max2-min2)*1/5, "$RMSE$ = %2.2f$\\degree$C" % rmse + "\n$Bias$ = %2.2f$\\degree$C" % bias +
             "\n$r^2$ = %2.2f" % r2 + "\n$\\sigma$ = %2.2f$\\degree$C" % std, fontsize=14)
    plt.show()

def sns_scatter(data1,data2,lai,smc):
    '''
    画散点图，重点是点的大小和颜色的展示,这里给出了一个例子
    :param data1: 数据1
    :param data2: 数据2
    :param lai:  类别1
    :param smc: 类别2
    :return:
    '''
    sns.set_theme(style="whitegrid")
    # sns.set_theme(style="white")
    # sns.set_palette(sns.color_palette("Paired", 4))
    cmap = sns.cubehelix_palette(rot=-.1, as_cmap=True)
    data = {'data1': data1, 'data2': data2, 'lai': lai, 'smc': smc}
    df = pd.DataFrame(data)
    g = sns.relplot(x='data1',y='data2',hue='SMC',size ='LAI',palette = cmap,data = df)
    g.set(xlabel ="data1",ylabel = "data2")
    # g.set()
    plt.show()

def sns_density(data1,data2,xmin = 250,xmax = 350, ymin = 250, ymax = 350):
    '''
    画散点图，重点是展示密度图
    :param data1:
    :param data2:
    :param xmin:
    :param xmax:
    :param ymin:
    :param ymax:
    :return:
    '''
    df = pd.DataFrame({'data1': data1, 'data2': data2})
    g = sns.jointplot(x='data1', y='data2', data=df, kind='hex', size=4)
    plt.plot([xmin, xmax], [ymin, ymax], '--')
    # g.set(xlabel='SCOPE LST (K)',ylabel='Fitted LST (K)')
    plt.xlabel("SCOPE LST")
    plt.ylabel("Fitted LST")
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    # plt.colorbar()
    plt.show()

def sns_bar(data1,data2,data3, label1='data1',label2='data2',label3 ='class'):
    '''
    画直方图，通过dataframe把数据组织起来，然后通过标签进行分割分类，进行展示
    :param data1: 数据
    :param data2: 数据
    :param data3:  数据
    :param label1: 。。。
    :param label2: 。。。
    :param label3: 。。。
    :return:
    '''
    sns.set_theme(style="white")
    # sns.set_palette(sns.color_palette("Paired", 4))
    data = {label1: data1, label2:data2, label3:data3}
    df = pd.DataFrame(data)
    plt.figure(figsize=[12.0, 4.5])
    sns.barplot(data,x='data1',y='data2',hue='class', alpha=1.0)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.show()

def sns_hist(data, label='data',):
    '''
    画密度直方图，可以有颜色，分布，累计，颜色等修饰
    :param data: 数据
    :param label: 描述
    :return:
    '''
    g = sns.histplot(x=data,  bins=30, norm=LogNorm(), log_scale=(False),fill=True,
    cumulative=False, stat="density")
    plt.tight_layout()
    plt.show()

