import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from base.util_packages import *

'''
util_fitting
求解方程函数
y = ax + b
y = f(x1,x2,x3...)
'''

def fitting_unknown_linear_ls(y,A):
    '''
    线性反演未知数，y=ax1+bx2+cx3, x3=1
    :param y: 因变量
    :param x: 自变量
    :return:
    '''
    coeffs = np.asarray(lstsq(A, y))[0]
    return coeffs

def fitting_unknown_nonlinear_fmin(fun,y,A,xn,x0):
    '''
    非线性反演/回归，基于fmin函数，能给初值，不能给边界条件
    :param fun: 函数
    :param y: 因变量
    :param x: 自变量，是个多列数组
    :param x0: 初值
    :param xn: 列数
    :return: 回归系数
    '''
    res = fmin(fun, x0=x0, args=(y,A,xn), disp=0)
    return res

def fitting_unknown_nonlinear_min(fun,y,A,xn,x0,bnds,method = 'Powell'):
    '''
    非线性反演/回归，基于minimize函数，能给初值，也能给边界条件
    :param fun: 函数
    :param y: 因变量
    :param x: 自变量，是个多列数组
    :param x0: 初值
    :param xn: 列数
    :param bnds: 边界
    :param method: 方法，默认为powell，还有slsqp，L-BFGS-B，TNC，
    :return: 回归系数
    '''
    res = minimize(fun, x0,
                   args=(y,A,xn),
                   method=method,
                   bounds=bnds, )


    return res.x

def fitting_unknown_nonlinear_ls(fun,y,A,xn,x0,bnds,method = 'trf'):
    '''
    非线性反演/回归，基于least_squares函数，能给初值，也能给边界条件
    :param fun: 函数
    :param y: 因变量
    :param A: 矩阵
    :param x0: 初值
    :param xn: 列数
    :param bnds: 边界
    :param method: 方法，trf', 'dogbox', 'lm'
    :return: 回归系数
    '''
    res = least_squares(fun, x0,
                   args=(y,A,xn),
                   method=method,
                   bounds=bnds, )
    return res.x

def fitting_unknown_linear_sklearn(y,A,method='ridge'):
    '''
    线性反演未知数，
    :param y: 因变量
    :param x: 自变量
    :param method: 方法，推荐ridge，然后是lasso，因为lasso会倾向于拟合0
    :return: 回归系数
    '''
    alphas = [0.0001,0.001,0.01,0.1,1,10,100]
    coeffs = 1.0
    if method=='ridge':
        clf = RidgeCV(alphas=alphas, fit_intercept=False)
        clf.fit(A, y)
        coeffs = clf.coef_
    elif method =='lasso':
        lasso = LassoCV(alphas=alphas, fit_intercept=False)
        lasso.fit(A, y)
        coeffs = lasso.coef_
    elif method =='bayesianridge':
        clf = BayesianRidge(fit_intercept=False)
        clf.fit(A,y)
        coeffs = clf.coef_

    return coeffs

from scipy.sparse.linalg import lsqr
def fitting_unknown_linear_sparse(y,A,iter=100):
    '''
    找到大型稀疏线性方程组的least-squares 解。
    :param y:
    :param x:
    :param iter:
    :return:
    '''
    coeffs = lsqr(A,y,iter_lim=iter)[0]
    return coeffs

from scipy.sparse.linalg import lsqr

'''
from here user defiend specific functions
'''

def fitting_unknown_nonlinear_min_condition(fun, y, A, xn, x0, bnds, method='Powell', constraints=None, **kwargs):
    '''
    非线性反演/回归，基于minimize函数，能给初值，也能给边界条件
    参数:
        fun: 目标函数，格式为fun(x, *args)
        y: 因变量数组
        A: 自变量矩阵
        xn: 变量个数
        x0: 初始值
        bnds: 边界条件列表
        method: 优化方法('Powell','SLSQP'等)
        constraints: 约束条件
        **kwargs: 传递给目标函数的额外参数
    返回:
        优化结果向量
    '''
    # 合并固定参数和额外参数
    fixed_args = (y, A, xn)
    extra_args = tuple(kwargs.values())
    all_args = fixed_args + extra_args

    res = minimize(fun, x0,
                   args=all_args,
                   method=method,
                   bounds=bnds,
                   constraints=constraints)
    return res.x


