# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:12:59 2016

@author: tmed2
"""


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import numpy as np


def nasa_reader(date):
    """
    Reads in, and returns the all the data corresponding to a particular 'date'
    """
    
    concs = np.array([])
    lats = np.array([])
    longs = np.array([])
    
    with open('NASA/NasaSatellite/NasaData.dat', 'rb') as file:
        for row in file:
            string_row = row.decode()
            string_row = string_row.split()
            
            if string_row[0] == date:
                if float(string_row[1]) > 0:
                    concs = np.append(concs, float(string_row[1]))
                    lats = np.append(lats, float(string_row[2]))
                    longs = np.append(longs, float(string_row[3]))
                else:
                    pass
            else:
                pass
    return concs, lats, longs


c, lt, lg = nasa_reader('15/10/2004')
bmap = Basemap(projection='kav7', lon_0=0, resolution='c')
bmap.drawcoastlines()
#bmap.fillcontinents(color='coral', lake_color='aqua')
#bmap.drawmapboundary(fill_color='aqua')
x, y = bmap(lg, lt)
cs = bmap.contourf(x, y, c, tri=True)
cbar = bmap.colorbar(cs, location='bottom', pad="5%")
plt.show()

"""
#define the allowed dates here:
DATES = []
start_month = 10
for month_since in range(111):
    #probably could find a more elegant way
    d = 15
    m =  (start_month + month_since) % 12
    y = 2004 + ((start_month + month_since) // 12)
    if m == 0:
        m = 12
        y -= 1
    DATES.append(str(d) + '/' + str(m) + '/' + str(y))


map_ani = plt.figure(figsize = (10,10))

bmap = Basemap(projection='kav7', lon_0=0, resolution='c')
bmap.drawcoastlines()
clevels = [0 + 10*i for i in range(16)]
cs0 = None


def contour_find(i):
    date = DATES[i]
    c, lt, lg = nasa_reader(date)
    x, y = bmap(lg, lt)
    
    cs = bmap.contourf(x, y, c, levels = clevels, tri = True)
    
    if i == 0:
        global cs0
        cs0 = cs
    else:
        pass
    return cs.collections

contour_list = []
for i in range(111):
    contour_list.append(contour_find(i))
    
cbar = bmap.colorbar(cs0, label='Ozone / ppb', location='bottom', pad="5%")
plt.title('Worldwide Tropospheric Ozone 2004-2013', fontsize = 22)

ani = animation.ArtistAnimation(map_ani, contour_list, interval=2000, 
                                blit=False, repeat_delay=5000)
map_ani.show()
"""