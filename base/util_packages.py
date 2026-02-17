import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy.linalg import lstsq
from osgeo import gdal
from osgeo import osr,ogr
import re
import scipy
import xlrd
import os
import cv2
import struct
import netCDF4
from sklearn.linear_model import LinearRegression
import seaborn as sns
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from pylab import *
import gzip
import shutil
import netCDF4 as nc
from scipy import stats
from scipy.optimize import minimize
from scipy.optimize import least_squares
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LassoCV
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import BayesianRidge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import BayesianRidge
import pysolar.solar as Sun

# os.environ['PROJ_LIB'] =  r"C:\ProgramData\anaconda3\envs\python37\Lib\site-packages\osgeo\data"

'''推荐使用python311版本'''