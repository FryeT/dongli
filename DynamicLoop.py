# coding: utf-8
import numpy as np
from dynamic_algs import *
import matplotlib.pyplot as plt

class Ellipse:
    def __init__(self, aa):
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

    def area(self):
        a, b = self.axis_length
        return np.pi*a*b

    def area_scatters(self):
        x, y = self.scatters_to_center()
        loop = np.stack((x, y), axis=-1)
        return area_of_loop(loop)

    def tri_area(self):
        x, y = self.scatters_to_center()
        return 0.5*x.max()*y.max()

    def xmax_to_center(self):
        return self.scatters_to_center()[0].max()

    def ymax_to_center(self):
        return self.scatters_to_center()[1].max()

    def xmax(self):
        return self.xmax_to_center()+self.center[0]

    def ymax(self):
        return self.ymax_to_center()+self.center[1]

    def info(self):
        print('Ellipse class')
        print('center: ({0:.2e}, {1:.2f})'.format(self.center[0], self.center[1]))
        print('axis length: ({0:.2e}, {1:.2f})'.format(self.axis_length[0], self.axis_length[1]))
        print('angle_of_rotation: {0:.2f})'.format(self.angle_of_rotation))
        print('ellipse area: {0:.2e}'.format(self.area()))
        print('ellipse area of scatters: {0:.2e}'.format(self.area_scatters()))
        print('tri area: {0:.2e}'.format(self.tri_area()))
        print('x_max, y_max: ({0:.2e}, {1:.2f})'.format(self.xmax(), self.ymax()))
        print('xc_max, yc_max: ({0:.2e}, {1:.2f})'.format(self.xmax_to_center(), self.ymax_to_center()))

class DynamicLoop:
    def __init__(self, epsilon_d, sigma_d, xscale=1, yscale=1):
        """
        DataFrame type epsilon_d and sigma_d in one loop
        """
        self.loop_x = epsilon_d.values
        self.loop_y = sigma_d.values
        self.loop = np.stack((epsilon_d.values, sigma_d.values), axis=-1)
        self.xscale = xscale
        self.yscale = yscale
    def ellipse_byfit(self):
        aa = fitEllipse(self.loop_x*self.xscale, self.loop_y*self.yscale)
        return Ellipse(aa)

    def area_byfit(self):
        return self.ellipse_byfit().area()

    def tri_area_byfit(self):
        return self.ellipse_byfit().tri_area()

    def damping_ratio_byfit(self):
        return self.area_byfit()/self.tri_area_byfit()/4/np.pi

    def sigma_m_byfit(self):
        return self.ellipse_byfit().ymax_to_center()

    def epsilon_m_byfit(self):
        return self.ellipse_byfit().xmax_to_center()
    
    def modulus_byfit(self):
        return self.sigma_m_byfit()/self.epsilon_m_byfit()

    def center(self):
        return self.ellipse_byfit().center
    
    def loop_to_center(self):
        xc, yc = self.center()
        return (self.loop_x-xc, self.loop_y-yc) 

    def sigma_m(self):
        return self.loop_to_center()[1].max()

    def epsilon_m(self):
        return self.loop_to_center()[0].max()

    def modulus(self):
        return self.sigma_m()/self.epsilon_m()

    def area(self):
        return area_of_loop(self.loop)

    def tri_area(self):
        x1, y1 = max(self.loop_x), max(self.loop_y)    #见吴世明《土动力学》p85
        x3, y3 = min(self.loop_x), min(self.loop_y)    #见吴世明《土动力学》p85
        return 0.5*(x1-x3)*(y1-y3)

    def damping_ratio(self):
        return self.area()/self.tri_area()/np.pi

    def plot_to_center(self, plot_loop=True, plot_fit=True):
        e = self.ellipse_byfit()    #ellipse
        ex, ey = e.scatters_to_center()
        lx, ly = self.loop_to_center()
        if plot_loop == True:
            plt.plot(lx, ly, label='loop')
        if plot_fit == True:
            plt.plot(ex, ey, label='fit')
        