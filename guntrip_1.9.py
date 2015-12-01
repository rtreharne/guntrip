#!/usr/bin/env python

""" -------------------------------------------------------------------------
GunTrip 1.9
Modeling superimposed profiles from sputtering guns
Profile generated from 2D gaussian function
----------------------------------------------------------------------------"""

import pygtk
import gtk
import numpy as np
from pylab import *
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
import tkFileDialog
import re
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.mlab import griddata
import scipy.interpolate as interp

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

class Profile():
    
    def generate(self, index):
    
        
        self.d = window.params[index]
        gamma = (window.params[-1]*np.pi)/180
        
        
        res = 33
        xi = np.linspace(-8, 8, res)
        yi = np.linspace(-8, 8, res)
        X, Y = np.meshgrid(xi, yi)
        Z = twoD_Gaussian((X, Y), 1134.9, -10.1, -5.2, 6.1, 6.6, 1.24, 76.1)*(self.d/332.76)
        
        theta = self.angle+gamma
        
        grid_x = X
        grid_y = Y
        z_values = Z
        
        x_irr = (X*cos(theta)-Y*sin(theta)).ravel()
        y_irr = (X*sin(theta)+Y*cos(theta)).ravel()
        
        grid_z = interp.griddata((x_irr, y_irr), Z, (X, Y))
        
        data = [X.reshape(res,res), Y.reshape(res, res), grid_z]
        
        return data
        
    def __init__(self, angle, d):
        self.angle = angle+(34*np.pi/180)
        self.d = d
        

class Wind(gtk.Window):

    def combo_change(self, widget):
        index = self.combo.get_active()
        self.adj1.set_value(self.params[index])
        self.adj1.set_lower(self.params_range[index][0])
        self.adj1.set_upper(self.params_range[index][1]+1)

    def scale_change(self, widget):
        index = self.combo.get_active()
        value = self.adj1.get_value()
        self.params[index] = value
        if (self.func == 0):
            self.plot(widget)
        if (self.func == 1):
            self.plot_comp(widget)
    
    def destroy(self, widget, data=None):
        gtk.main_quit()
        
    def plot(self, widget):
        self.func = 0
        A = [profile0.generate(0), profile1.generate(1), profile2.generate(2), profile3.generate(3), profile4.generate(4)]
        
        Z_sum = A[0][2]
        
        for i in range (1, len(A)):
            Z_sum += A[i][2]
        
        
        self.fig1.clf()
        self.ax1 = self.fig1.add_subplot(111)
        self.ax1.set_aspect('equal')
        self.ax1.set_xlim(-5, 5)
        self.ax1.set_ylim(-5, 5)
        self.ax1.plot(A[0][0], A[0][1], 'o', color='white', markersize=2)
        
        inc = 20
        
        Z_max = np.nanmax(Z_sum.ravel())/1.1
        #if (Z_max < 200):
        #   Z_max = 300*1.2
        
        Z_1 = Z_max
          
        if (Z_max < 1000):
            inc = 100
            Z_1 = 1000
              
        if (Z_max < 100):
            inc = 10
            Z_1 = 100
            
        if (Z_max < 10):
            inc = 1
            Z_1 = 10
        
        lvls = np.arange(0, Z_1, inc)
        #plot1 = self.ax1.contourf(A[0][0], A[0][1], Z_sum, zdir='z', cmap=cm.jet, alpha = 1.0, levels = lvls)
        plot1 = self.ax1.contour(A[0][0], A[0][1], Z_sum, colors='black', linewidths=1, levels=lvls)
        self.ax1.clabel(plot1,  fmt='%1.0f', fontsize = 12)
        #color bar added to figure not axes
        #cbar = self.fig1.colorbar(plot1)
        self.canvas.draw()
        
    def plot_comp(self, widget):
        self.func = 1
        self.fig1.clf()
        self.ax1 = self.fig1.add_subplot(111)
        
        A = [profile0.generate(0), profile1.generate(1), profile2.generate(2), profile3.generate(3), profile4.generate(4)]
        
        x = A[0][0].ravel()
        y = A[0][1].ravel()
        
        Z0 = A[0][2].ravel()
        Z1 = A[1][2].ravel()
        Z2 = A[2][2].ravel()
        Z_sum = Z0 + Z1 + Z2
        A = 3.58
        B = 3.95
        C = 5.61
        
        x_zoom = []
        y_zoom = []
        z_zoom =[]
        
        
        
        X = (A*Z0/(B*Z1+C*Z2+A*Z0))*100
        Y = (B*Z1/(B*Z1+C*Z2+A*Z0))*100
        
        for i in range (0,len(x)):
            if (-4 <= x[i] <=4 and -4 <= y[i] <=4):
                x_zoom.append(X[i])
                y_zoom.append(Y[i])
                z_zoom.append(Z_sum[i])
                
        print len(x_zoom)
        
        plot1 = self.ax1.scatter(x_zoom, y_zoom, s=100, c=z_zoom, alpha = 0.5)
        self.ax1.set_xlabel('% wt. MgO')
        self.ax1.set_ylabel('% wt. Al2O3')
        cbar = self.fig1.colorbar(plot1)
        cbar.set_label('thickness, $d$, nm')
        self.canvas.draw()

    def __init__(self):
        self.func = 0
        self.gamma = 0
        self.params = [profile0.d, profile1.d, profile2.d, profile3.d, profile4.d, self.gamma]
        self.params_range = [(0, 500), (0, 500), (0, 500), (0, 500), (0, 500), (-180, 180)]
    
        self.fig1 = Figure()
        self.canvas = FigureCanvas(self.fig1)

        mb = gtk.MenuBar()
  
        filemenu = gtk.Menu()
        filem = gtk.MenuItem("File")
        filem.set_submenu(filemenu)

        exit = gtk.MenuItem("Exit")
        exit.connect("activate", gtk.main_quit)
        filemenu.append(exit)

        mb.append(filem)

        self.plot_but = gtk.Button('Plot')
        self.plot_but.connect("clicked", self.plot)
        
        self.plot_comp_but = gtk.Button('Composition')
        self.plot_comp_but.connect("clicked", self.plot_comp)
        #self.plot_but.set_size_request(700, 50)

        self.adj1 = gtk.Adjustment(self.params[0], self.params_range[0][0], self.params_range[0][1], 1.0, 1.0, 1.0)
        self.adj1.connect("value_changed", self.scale_change)
        self.hscale = gtk.HScale(self.adj1)
        self.hscale.set_size_request(800, 35)

        self.combo = gtk.combo_box_entry_new_text()
        self.combo.set_size_request(200, 50)
        self.combo.append_text('GUN 1')
        self.combo.append_text('GUN 2')
        self.combo.append_text('GUN 3')
        self.combo.append_text('GUN 4')
        self.combo.append_text('GUN 5')
        self.combo.append_text('Angular Offset')
        self.combo.set_active(0)
        self.combo.connect("changed", self.combo_change)

        self.win = gtk.Window()
        self.win.set_size_request(700,650)
        self.win.set_position(gtk.WIN_POS_CENTER)
        self.win.set_title('GunTrip 1.0')    
       
                
        self.box1 = gtk.HBox()
        self.box1.set_size_request(980,500)
        self.box1.pack_start(self.canvas)

        self.box2 = gtk.HBox()
        self.box2.pack_start(self.plot_but)
        self.box2.pack_start(self.plot_comp_but)
        self.box2.pack_start(self.combo)
        
        self.box3 = gtk.HBox()
        self.box3.set_size_request(980, 50)
        self.box3.pack_start(self.hscale)
        
    
        self.box6 = gtk.VBox(False, 2)
        self.box6.pack_start(mb, False, False, 0)
          
        self.box4 = gtk.VBox()
        self.box4.pack_start(self.box1)
        self.box4.pack_start(self.box2)
        self.box4.pack_start(self.box3)
     
        self.box4.pack_start(self.box6)
                       
        self.win.add(self.box4)        
                
        self.win.show_all()
        self.win.connect("destroy", self.destroy)
        

    def main(self):
        gtk.main()

if __name__ == "__main__":
    profile0 = Profile((2*np.pi/5)*0, 250)
    profile1 = Profile((2*np.pi/5)*1, 5)
    profile2 = Profile((2*np.pi/5)*2, 250)
    profile3 = Profile((2*np.pi/5)*3, 0)
    profile4 = Profile((2*np.pi/5)*4, 0)
    window = Wind()
    window.main()
    
    
    
    
        

