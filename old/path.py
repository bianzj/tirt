import numpy as np
import scipy.integrate as sci


'''
路径长度
角度数目*体素数目*半球数目
a1,a1,a1...a2,a2,a2....an,an,an
v1,v2,vn...v1,v2,vn....v1,v2,vn
。。。
'''

#### 垄行结构

def path_length_row_voxel(length_voxel, height_voxel, row_width, row_blank, row_height, xza, raa, vsa = 0, ifinverse=0):
	'''
	垄行场景,每个植被体元到达植被表面的路径长度
	:param length_voxel: l 方向坐标
	:param height_voxel: h 方向坐标
	:param xza: 天顶角
	:param raa: 垄行方位角
	:param vsa: 太阳-观测方位角
	:param row_width: 垄植被宽度
	:param row_blank: 垄空白宽度
	:param row_height: 垄高度
	:param ifinverse: 是否求逆
	:return: 路径长度
	'''

	rw = row_width
	rb = row_blank
	rh = row_height
	rs = rw + rb


	# number_voxel = np.size(length_voxel)
	# number_angle = np.size(vza)
	# length_voxel = np.tile(length_voxel,number_angle)
	# height_voxel = np.tile(height_voxel,number_angle)
	# vza = np.repeat(vza, number_voxel)
	# raa = np.repeat(raa, number_voxel)
	# vsa = np.repeat(vsa, number_voxel)
	### 当求拟方向的时候，判断太阳方向和观测方向是否在垄行的同侧，如果在同侧则其路径长度相同，如果不是同侧，则路径水平距离正好相反
	if(ifinverse):
		### 垄行间空白区域体素
		ind = (vsa >= 90) * (length_voxel > rw)
		if np.sum(ind) > 0:
			length_voxel[ind] = (rs - length_voxel[ind]) + rw
		### 垄行内植被区域体素
		ind = (vsa >= 90) * (length_voxel < rw)
		if np.sum(ind) > 0:
			length_voxel[ind] = (rw - length_voxel[ind])


	'''
	对于垂直的有不同的处理方式 vza = 0
	对于观测和垄行同向的有不同的处理方式 raa = 0
	'''
	ind_vza_0 = (xza == 0)
	ind_raa_0 = (raa == 0)
	ind_vza_soil = (xza == 0) * (length_voxel > rw)
	ind_raa_soil = (raa == 0) * (length_voxel > rw)
	raa[ind_raa_0] = 1
	xza[ind_vza_0] = 1

	tgthv = np.tan(np.deg2rad(xza))
	sthetx = np.sin(np.deg2rad(xza))
	sphir = np.sin(np.deg2rad(raa))
	cthetx = np.cos(np.deg2rad(xza))


	###从冠层顶到该体素的距离'''
	###体素到冠层出射点的水平距离'''
	length_voxel_up = height_voxel * tgthv * sphir
	length_up = length_voxel + length_voxel_up
	###有多少个垄行个体，小数结果'''
	l_row = (length_up / rs)
	###垄行个体的整数部分'''
	n_row = np.floor(l_row)
	###垄行个体的小数部分，在植被部分的占比'''
	temp = (l_row - n_row)/ (rw/rs)
	###占比只能是1'''
	temp[temp > 1]  = 1
	###体素的前置小数部分'''
	###如果是大于植被部分，则只能是植被部分'''
	length_res = length_voxel*1.0
	length_res[length_res > rw] = rw
	l_row0 = (n_row + temp) * rw - length_res
	###垂直垄方向距离，变换到倾斜角度方向距离'''
	pl = l_row0 / sphir / sthetx
	pl[ind_raa_0] = height_voxel[ind_raa_0] / cthetx[ind_raa_0]
	pl[ind_vza_0] = height_voxel[ind_vza_0]
	pl[ind_vza_soil] = 0
	pl[ind_raa_soil] = 0

	return pl

#### 树冠结构

def path_length_crown_inside_voxel(xi,yj,zk,za,aa,hcr,rcr):
	'''
	计算树冠内植被体元到达树冠表面的路径长度, 体元数目要与角度数目相同，如果不相同计算困难
	核心是求解一元二次方程组，aaa,bbb,ccc 分别是其系数，这一步有点慢
	from frt in rlips.f
	:param xi:  x 方向距离
	:param yj:  y 方向距离
	:param zk:  z 方向距离
	:param za:  天顶角
	:param aa:  方位角
	:param hcr:  垂直方向半径
	:param rcr:  水平方向半径
	:param eps:  最小值
	:return: 向上方向和向下方向的路径长度（向下的路径长度是负值）
	'''
	eps = 0.0001
	number_angle = np.size(za)
	number_voxel = np.size(xi)
	number = number_voxel
	if number_angle != number_voxel:
		number = number_voxel * number_angle
		za = np.repeat(za, number_voxel)
		aa = np.repeat(aa, number_voxel)
		xi = np.tile(xi, number_angle)
		yj = np.tile(yj, number_angle)
		zk = np.tile(zk, number_angle)
	uz = np.cos(np.deg2rad(za) )
	sz = np.sin(np.deg2rad(za) )
	sa = np.sin(np.deg2rad(aa) )
	ua = np.cos(np.deg2rad(aa) )
	rlout = np.zeros(number)
	rlout_ = np.zeros(number)
	a2 = rcr * rcr
	c2 = hcr * hcr
	aaa = c2 * (sz * sz) + a2 * (uz * uz)
	bbb = 2.0 * (c2 * sz * (xi * ua + yj * sa) + a2 * uz * zk)
	ccc = c2 * (xi * xi + yj * yj) + a2 * (zk * zk - c2)
	det = (bbb * bbb - 4. * aaa * ccc)
	### 等于0，一个解
	ind = np.abs(det) < eps
	if np.sum(ind) > 0:
		rlout[ind] = -bbb[ind] * 0.50 / aaa[ind]
		rlout_[ind] = rlout[ind] * 1.0
	### 大于0，两个解
	ind = det > eps
	if (np.sum(ind) > 0):
		rlout[ind] = (np.sqrt(bbb[ind] * bbb[ind] - 4.0 * aaa[ind] * ccc[ind]) - bbb[ind]) * .50 / aaa[ind]
		rlout_[ind] = (-np.sqrt(bbb[ind] * bbb[ind] - 4.0 * aaa[ind] * ccc[ind]) - bbb[ind]) * .50 / aaa[ind]
	### 小于0，没有解
	# ind = det < eps
	# if np.sum(ind) > 0:
	# 	rlout[ind] = 0.0
	# 	rlout_[ind] = 0.0
	return rlout, rlout_

def path_length_crown_outside_voxel(xi, yj, zk, hc, hcr, rcr, xza, xaa):
	'''
	树冠表面到达树冠顶的路径长度，没有方位角概念，因为是对称的
	:param zk: 体元z方向高度
	:param xza: 天顶角
	:param xaa: 方位角
	:param hc: 高度
	:param hcr: 树冠高度方向半径
	:param rcr: 树冠水平方向半径
	:return: 路径长度 和投影面积
	'''
	### 象征性判断体素问题
	number_angle = np.size(xza)
	number_voxel = np.size(zk)
	number = number_voxel
	if number_angle != number_voxel:
		number = number_voxel * number_angle
		xza = np.repeat(xza, number_voxel)
		xaa = np.repeat(xaa, number_voxel)
		zk = np.tile(zk, number_angle)

	cthetx = np.cos(np.deg2rad(xza))
	[upArea, upVol] = crosscutting_ellipsoid(zk, xza, hc, hcr, rcr)
	pl_outside = np.zeros(number)
	ind = (zk <= hc) * (upArea > 0)
	if np.sum(ind) > 0:
		### 路径长度问题，这里直接采用了FRT的结果 在 spooj.f 80 行左右
		pl_outside[ind] = upVol[ind] / upArea[ind] / cthetx[ind]
	return pl_outside,upArea

def path_length_crown_outside_voxel_corrected(xi, yj, zk, hc, hcr, rcr, xza, xaa):
	'''
	树冠表面到树冠顶的路径长度， 这里修正了树与树重合的问题，意义同上
	:param zk: 树冠高度
	:param xza: 天顶角
	:param xaa: 方位角
	:param hc: 高度
	:param hcr: 树冠半径 垂直方向
	:param rcr: 树冠半径 水平方向
	:return: 路径长度
	'''
	number_angle = np.size(xza)
	number_voxel = np.size(zk)
	number = number_voxel
	if number_angle != number_voxel:
		number = number_voxel * number_angle
		xza = np.repeat(xza, number_voxel)
		xaa = np.repeat(xaa, number_voxel)
		zk = np.tile(zk, number_angle)
	uz = np.cos(np.deg2rad(xza))
	rk = np.sqrt(xi*xi + yj * yj)
	[upArea, upVol] = oblique_ellipsoid(zk, rk, xza, hc, hcr, rcr)
	pl_outside = np.zeros(number)
	ind = (zk <= hc) * (upArea > 0)
	if np.sum(ind) > 0:
		pl_outside[ind] = upVol[ind] / upArea[ind] / uz[ind]
	return pl_outside,upArea

def oblique_ellipsoid(zk,rk, za, hc, hcr, rcr):
	'''
	高度zk 以上沿着天顶角方向的，切出的体积和投影面积
	The projection and volume of the upper part
    of an ellipsoid at the level zk in direction thx
    it is from frt model in pi11u.f
	:param zk: height level
	:param za: zenith angle
	:param hc: tree height = trunk height + crown height
	:param hcr: crown radi / cell b
	:param rcr: crown radi / cell a
	:return: projection and volume
	'''

	hcc = hc - hcr
	hcb = hc - hcr * 2
	eps = 0.0001

	number = np.size(zk)
	tv = np.tan(np.deg2rad(za))
	### default value: above the crown
	upArea = np.zeros(number)
	upVol = np.zeros(number)
	### under the crown
	ind = (zk >= 0) * (zk <= hcb)
	if np.sum(ind) > 0:
		upArea[ind] = np.sqrt(rcr * rcr + np.power(hcr * tv, 2)) * np.pi * rcr
		upVol[ind] = 4.0 * np.pi * rcr * rcr * hcr / 3.0

	###
	ind = (tv > eps)
	dk = np.zeros(number)+hcr
	if np.sum(ind) > 0:
		rd = 2*rcr - rk
		hk = rd[ind] / np.tan(np.deg2rad(za[ind]))
		hkh = hk + (zk[ind] - hcc)
		dk[ind] = hkh * np.sin(np.deg2rad(za[ind]))
	dk[dk > hcr] = hcr
	dk[dk < -hcr] = -hcr

	# zx = zk - hcc
	zx = dk * 1.0
	ind0 = (zk <= hc) * (zk > hcb)
	if (np.sum(ind0) > 0):
		### crown volume above the level zk
		upVol[ind0] = np.pi * rcr * rcr * (hcr * 2 / 3.0 -
					zx[ind0] + zx[ind0] * zx[ind0] * zx[ind0] / 3.0 / (hcr * hcr))
		### zo height of intersection from 0 to hcr
		z0 = np.zeros(number)
		ind = (tv > eps) * ind0
		if np.sum(ind):
			z0[ind] = hcr / np.sqrt(rcr / np.power(hcr * tv[ind], 2) + 1.0)
		### crown radii at height zx
		rz = np.sqrt(hcr * hcr - zx * zx) * rcr / hcr
		### top of crown above the intersection point
		ind = (abs(zx) >= z0) * (zk > hcc) * ind0
		if np.sum(ind) > 0:
			upArea[ind > 0] = np.pi * rz[ind] * rz[ind]
		### bottom of crown under the intersection point
		ind = (abs(zx) >= z0) * (zk < hcc) * ind0
		if np.sum(ind) > 0:
			upArea[ind] = np.sqrt(rcr * rcr + np.power(hcr * tv[ind], 2)) * np.pi * rcr
		### two parts for the crown
		ind = (abs(zx) < z0) * ind0
		if (np.sum(ind) > 0):
			xyz = np.sqrt(rcr * rcr + np.power((hcr * tv[ind]), 2))
			beta = np.arccos(zx[ind] * tv[ind] / xyz)
			sel1 = rcr * xyz * beta - rz[ind] * zx[ind] * tv[ind]
			sel2 = np.pi * rz[ind] * rz[ind] / 2.0
			upArea[ind] = sel1 + sel2


	return np.abs(upArea), np.abs(upVol)

def crosscutting_ellipsoid(zk, xza, hc, hcr, rcr):
	'''
	冠层表面以上部分的体积和投影面积
	The projection and volume of the upper part
    of an ellipsoid at the level zk in direction thx
    it is from frt model in pi11u.f
	:param zk: height level
	:param xza: zenith angle
	:param hc: tree height = trunk height + crown height
	:param hcr: crown radi / cell b
	:param rcr: crown radi / cell a
	:return: projection and volume
	'''
	hcc = hc - hcr
	hcb = hc - hcr * 2
	eps = 0.0001

	number = np.size(zk)
	tgthx = np.tan(np.deg2rad(xza))
	### default value: above the crown
	upArea = np.zeros(number)
	upVol = np.zeros(number)
	### under the crown
	ind = (zk >= 0) * (zk <= hcb)
	if np.sum(ind) > 0:
		upArea[ind] = np.sqrt(rcr * rcr + np.power(hcr * tgthx[ind], 2)) * np.pi * rcr
		upVol[ind] = 4.0 * np.pi * rcr * rcr * hcr / 3.0

	zx = zk - hcc
	ind0 = (zk <= hc) * (zk > hcb)*(np.abs(zx)<hcr)
	if (np.sum(ind0) > 0):
		### crown volume above the level zk
		upVol[ind0] = np.pi * rcr * rcr * (hcr * 2 / 3.0 -
					zx[ind0] + zx[ind0] * zx[ind0] * zx[ind0] / 3.0 / (hcr * hcr))
		### zo height of intersection from 0 to hcr
		z0 = np.zeros(number)
		ind = (tgthx > eps) * ind0
		if np.sum(ind):
			z0[ind] = hcr / np.sqrt(rcr / np.power(hcr * tgthx[ind], 2) + 1.0)
		### crown radii at height zx
		rz = np.zeros(number)
		rz[ind0] = np.sqrt(hcr * hcr - zx[ind0] * zx[ind0]) * rcr / hcr
		### top of crown above the intersection point
		ind = (abs(zx) >= z0) * (zk > hcc) * ind0
		if np.sum(ind) > 0:
			upArea[ind > 0] = np.pi * rz[ind] * rz[ind]
		### bottom of crown under the intersection point
		ind = (abs(zx) >= z0) * (zk < hcc) * ind0
		if np.sum(ind) > 0:
			upArea[ind] = np.sqrt(rcr * rcr + np.power(hcr * tgthx[ind], 2)) * np.pi * rcr
		### two parts for the crown
		ind = (abs(zx) < z0) * ind0
		if (np.sum(ind) > 0):
			xyz = np.sqrt(rcr * rcr + np.power((hcr * tgthx[ind]), 2))
			beta = np.arccos(zx[ind] * tgthx[ind] / xyz)
			sel1 = rcr * xyz * beta - rz[ind] * zx[ind] * tgthx[ind]
			sel2 = np.pi * rz[ind] * rz[ind] / 2.0
			upArea[ind] = sel1 + sel2

	return np.abs(upArea), np.abs(upVol)

