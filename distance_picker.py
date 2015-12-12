#!/usr/bin/env python
"""
This script draws distance lines to a center point on a global map. The center point is chosen
interactively by clicking on the map
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
import numpy as np

#==== MAIN =====
def main():
    fig,ax = plt.subplots(1,1)

    #---- basemap setup ----
    lon_0 = 0.
    bmap = Basemap(projection='robin',lon_0=lon_0)
    bmap.fillcontinents()
    bmap.drawcoastlines()

    #---- event handler class that draws lines around a points ----
    evhandler = PickEvent(fig,ax,bmap)

    #---- add event handler to figure and show plot ----
    cid = fig.canvas.mpl_connect('button_press_event', evhandler)
    plt.show()

class PickEvent(object):
    def __init__(self, fig, ax, bmap):
        self.showlines = False
        self.ax   = ax
        self.fig  = fig
        self.bmap = bmap

    def __call__(self,event):
        if self.showlines:
            self.lines.remove()
            del self.lines
        lon,lat     = self.bmap(event.xdata,event.ydata,inverse=True)
        self.lines  = self.get_distance_lines(lat,lon)
        self.ax.add_collection(self.lines)
        self.showlines = True
        self.fig.canvas.draw()
        print("lon,lat: ", lon,lat)

    def get_distance_lines(self,clat,clon,npts90=200):
        """
        returns a set of distance line segments around the center point clat,clon
        in map coordinates.
        """
        #line setup
        thetas = np.radians(np.arange(30.,180.,30.))
    
        #get rotation matrix for clat,clon
        ctheta  = -np.radians(90.-clat)
        cphi    = np.radians(clon)
        rmatrix = np.dot(Rz(cphi),Ry(ctheta))
    
        #make lines
        segmentlist = []
        for theta in thetas:
            #get xyz coordinates of distance line
            npts   = npts90 * np.sin(theta) + 1
            phis   = np.linspace(0.,2*np.pi,npts)
            points = np.vstack( tp2xyz(theta*np.ones(npts),phis) )
    
            #rotate xyz coordinates and transform back to lat,lon
            points_rot = np.dot(rmatrix,points)
            thetas_rot,phis_rot = xyz2tp(points_rot[0],points_rot[1],points_rot[2])
            lats_rot, lons_rot = 90.-np.degrees(thetas_rot), np.degrees(phis_rot)
    
            #find longitude jumps at the map boundaries
            lons_rot[lons_rot>self.bmap.lonmax] -= 360.
            lons_rot[lons_rot<self.bmap.lonmin] += 360.
            threshold = 90.
            isplit = np.nonzero(np.abs(np.diff(lons_rot)) > threshold)[0]
    
            #split line at jumps
            lonlat_rot = np.vstack( (lons_rot,lats_rot) )
            subsegs    = np.split(lonlat_rot,isplit+1,axis=1)
    
            #add segments to the list of all lines
            for seg in subsegs:
                x = np.transpose(np.vstack(self.bmap(seg[0],seg[1])))
    
                segmentlist.append(x)
    
            segments = LineCollection(segmentlist)
        return segments

#---- spherical coordinates ----
def tp2xyz(t,p):
    x = np.sin(t)*np.cos(p)
    y = np.sin(t)*np.sin(p)
    z = np.cos(t)
    return x,y,z

def xyz2tp(x,y,z):
    t = np.arccos(z)
    p = np.arctan2(y,x)
    return t,p

#---- rotation functions ----
def Rz(angle):
    Rmatrix = np.array([[np.cos(angle),-np.sin(angle),0],
                        [np.sin(angle), np.cos(angle),0],
                        [0             , 0           ,1]])
    return Rmatrix

def Ry(angle):
    Rmatrix = np.array([[np.cos(angle),0, -np.sin(angle)],
                        [0             , 1           ,0],
                        [np.sin(angle), 0, np.cos(angle)]])
    return Rmatrix

#==== EXECUTION ====
if __name__ == "__main__":
    main()
