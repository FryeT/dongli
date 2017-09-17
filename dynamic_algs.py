# coding: utf-8
import numpy as np
from numpy.linalg import eig, inv

# 简化法求滞回圈面积及阻尼比
# loop: 一个周期滞回圈np.ndarray对象数据，[epsilon,sigma]
def area_of_loop(loop):
    """滞回圈的面积"""
    area = 0
    for i, val in enumerate(loop):
        area += np.linalg.det([loop[i-1], val])
    return -0.5*area    # 顺时针遍历是负值，所以乘以-0.5


# 椭圆拟合法求滞回圈面积 阻尼比
def fitEllipse(x, y):
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    D = np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T, D)
    C = np.zeros([6, 6])
    C[0, 2] = C[2, 0] = 2
    C[1, 1] = -1
    E, V = eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:, n]
    return a

def ellipse_center(aa):
    b, c, d, f, g, a = aa[1]/2, aa[2], aa[3]/2, aa[4]/2, aa[5], aa[0]
    num = b*b-a*c
    x0 = (c*d-b*f)/num
    y0 = (a*f-b*d)/num
    return np.array([x0, y0])

def ellipse_axis_length(aa):
    b, c, d, f, g, a = aa[1]/2, aa[2], aa[3]/2, aa[4]/2, aa[5], aa[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1 = (b*b-a*c)*((c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2 = (b*b-a*c)*((a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1 = np.sqrt(up/down1)
    res2 = np.sqrt(up/down2)
    return np.array([res1, res2])

def ellipse_angle_of_rotation2(aa):
    b, c, d, f, g, a = aa[1]/2, aa[2], aa[3]/2, aa[4]/2, aa[5], aa[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else:
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2

def ellipse_angle_of_rotation(aa):
    b, c, d, f, g, a = aa[1]/2, aa[2], aa[3]/2, aa[4]/2, aa[5], aa[0]
    phi = 0.5*np.arctan(2*b/(a-c))
    theta = phi+np.pi/2 if phi < 0 else phi-np.pi/2 if phi > np.pi/2 else phi
    return theta

def damping_ratio_byfit(loop, xscale=1, yscale=1):
    """
    椭圆拟合法求阻尼比
    返回(阻尼比:damping ratio, 拟合点：(xx, yy),(epsilon_m,sigma_m), 椭圆中心：center, 椭圆面积：area_e, 三角形面积：area_t)
    """
    x, y = loop[:, 0]*xscale, loop[:, 1]*yscale
    R = np.arange(0, 2*np.pi, 0.01)
    aa = fitEllipse(x, y)
    center = ellipse_center(aa)
    phi = ellipse_angle_of_rotation(aa)
    a, b = ellipse_axis_length(aa)    #长短半轴
    if a < b:
        a, b = b, a
    xx = (center[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi))
    yy = (center[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi))
    epsilon_m = xx.max()-center[0]    #椭圆顶点的x值，动应变顶点
    sigma_m = yy.max()-center[1]      #椭圆定点的y值，动应力顶点
    area_t = 0.5 * epsilon_m * sigma_m    #小三角形面积
    area_e = np.pi*a*b                    #椭圆面积
    damping_ratio = area_e/area_t/4/np.pi
    #print([area_e, area_t, damping_ratio])
    return (damping_ratio, (xx, yy), center, area_e/xscale/yscale, area_t/xscale/yscale)
