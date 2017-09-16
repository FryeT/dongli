# coding: utf-8
import numpy as np
from numpy.linalg import eig, inv

class Ellipse:
    def __init__(self,aa):
        """aa="""
        self.center = ellipse_center(aa)
        self.axis_length = ellipse_axis_length(aa)
        self.angle_of_rotation = ellipse_angle_of_rotation(aa)

    def scatters_to_center(self):
        R = np.arange(0, 2*np.pi, 0.01)
        phi = self.angle_of_rotation
        a, b = self.axis_length
        if a < b:
            a, b = b, a
        xx = a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
        yy = a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)
        return (xx, yy)

    def scatters(self):
        xx, yy = self.scatters_to_center()
        return (xx+self.center[0], yy+self.center[1])

    def ellipse_area(self):
        a, b = self.axis_length
        return np.pi*a*b

    def ellipse_area_scatters(self):
        x, y = self.scatters_to_center()
        loop = np.stack((x, y), axis=-1)
        return area_of_loop(loop)

    def tri_area(self):
        x, y = self.scatters_to_center()
        return 0.5*x.max()*y.max()

    def xmax(self):
        return self.scatters_to_center()[0].max()

    def ymax(self):
        rturn self.scatters_to_center()[1].max()

class DynamicLoop:
    def __init__(self, epsilon_d, sigma_d):
        """
        DataFrame type epsilon_d and sigma_d in one loop
        """
        self.loop_x = epsilon_d.values
        self.loop_y = sigma_d.values
        self.loop = np.stack((epsilon_d.values, sigma_d.values), axis=-1)

    def ellipse_byfit(self):
        aa = fitEllipse(self.loop_x, self.loop_y)
        return Ellipse(aa)

    def area_byfit(self):
        return self.ellipse_byfit().ellipse_area()

    def tri_area_byfit(self):
        return self.ellipse_byfit().tri_area()

    def damping_ratio_byfit(self):
        return self.area_byfit()/self.tri_area_byfit()/4/np.pi

    def sigma_m_byfit(self):
        return self.ellipse_byfit().ymax()

    def epsilon_m_byfit(self):
        return self.ellipse_byfit().xmax()
    
    def modulus_byfit(self):
        return self.sigma_m_byfit()/self.epsilon_m_byfit()

    def sigma_m(self):
        return self.loop_y.max()

    def epsilon_m(self):
        return self.loop_x.max()

    def modulus(self):
        return self.sigma_m()/self.epsilon_m()

    def area(self):
        return area_of_loop(self.loop)

    def tri_area(self):
        return 
    
    def damping_ratio(self):
        e_fit = self.ellipse_byfit()
        return e_fit.ellipse_area()/e_fit.tri_area()/4/np.pi

    


# 简化法求滞回圈面积及阻尼比
# loop: 一个周期滞回圈np.ndarray对象数据，[epsilon,sigma]
def area_of_loop(loop):
    """滞回圈的面积"""
    area = 0
    for i, val in enumerate(loop):
        area += np.linalg.det([loop[i-1], val])
    return -0.5*area    # 顺时针遍历是负值，所以乘以-0.5

def damping_ratio_nofit(loop):
    """
    返回(阻尼比:damping_ratio, 顶点坐标1：(x1,y1), 顶点坐标3：(x3,y3))
    """
    x1, y1 = max(loop[:, 0]), max(loop[:, 1])    #见吴世明《土动力学》p85
    x3, y3 = min(loop[:, 0]), min(loop[:, 1])
    area_t = (x1-x3)*(y1-y3)/2
    damping_ratio = area_of_loop(loop)/area_t/np.pi
    return (damping_ratio, (x1, y1), (x3, y3))


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

def modulus_nofit(loop):
    epsilon_m, sigma_m = max(loop[:, 0]), max(loop[:, 1])
    return (epsilon_m, sigma_m, sigma_m/epsilon_m)


