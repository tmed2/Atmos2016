# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:28:40 2016

@author: tmed2
"""


from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import matplotlib.colors as mpl_colors
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

import datetime as dtm

import atmosgastools.atmosgasts as ag


START = dtm.datetime(2003, 1, 1)
END = dtm.datetime(2014, 12, 28)

H = ag.WDCGG_TS('Heimaey', 'O3', 'Iceland O3 2003-2010')
I = ag.WDCGG_TS('Izana', 'O3', 'Izana O3 2008-2014')
MH = ag.WDCGG_TS('Mace Head', 'O3', 'MaceHeadO3')
PdM = ag.WDCGG_TS('Pic du Midi', 'O3', 'Pic du Midi O3 2008-2012')
TH = ag.WDCGG_TS('Tudor Hill', 'O3', 'TudorHillO3 2003-2014')
RP = ag.WDCGG_TS('Ragged Point', 'O3', 'BarbadosO3 2006-2015')

gas_arr_full = [H, I, MH, PdM, RP, TH]
gas_arr = [H, I, PdM, RP, TH]

for gas in gas_arr_full:
    gas.find_t_rel()
    assert gas.t_rel is not None

MH.regularise(START, END, 1)

padN = 100000
interval = 60*60
for gas in gas_arr:    
    tl = gas.t_rel[0] - (np.arange(1, padN+1, 1)[::-1]) * interval
    tr = gas.t_rel[-1] + (np.arange(1, padN+1, 1)) * interval
    
    old_t_rel = gas.t_rel
    new_t_rel = np.append(tl, old_t_rel)
    new_t_rel = np.append(new_t_rel, tr)
    gas.t_rel = new_t_rel
    
    #Compared to the above, why didn't this work...
    #gas.t_rel = np.append(tl, gas.t_rel)
    #gas.t_rel = np.append(gas.t_rel, tr)

    
    gas.c = np.pad(gas.c, padN, mode = 'constant')
    
    gas.find_interpol()
    gas.set_interpol_data(MH.t)
    assert gas.t_rel != old_t_rel

for gas in gas_arr_full:
    gas.moving_mean(24*7)
ag.ts_multi_plot(gas_arr_full, join_data = True)

###############################################################################
#define map extent
lat_l = 0
lat_h = 70
lng_l = -100
lng_h = 30
lat0 = 0.5*(lat_l + lat_h)
lng0 = 0.5*(lng_l + lng_h)

#creates the map and calls animation methods
map_ani = plt.figure(figsize = (12, 9))
m = Basemap(llcrnrlon = lng_l, llcrnrlat = lat_l, urcrnrlon = lng_h,
            urcrnrlat = lat_h, resolution = 'i', projection = 'merc',
            lat_0 = lat0, lon_0 = lng0)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='white')

# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,15.))
m.drawmeridians(np.arange(-180.,181.,15.))
m.drawmapboundary(fill_color='white')
plt.title('Concentration of Surface Ozone Over\nthe North Atlantic',
          fontsize = 22)

#not needed; replaced with concentration colour bar
#creates a mapping between scatter colour and year, and a legend also
#years = np.array([2003 + i for i in range(13)])
#colours = cm.ScalarMappable().to_rgba(years)
#cymap = {}
#patches = []
#for i, y in enumerate(years):
#    c = colours[i][:3]
#    cymap[y] = c
#    patches.append(mpatches.Patch(color = tuple(c), label = y))
#plt.legend(handles=patches)

v_min = 0
v_max = 100
normalise = mpl_colors.Normalize(v_min, v_max)
colmap = cm.ScalarMappable(norm = normalise, cmap = 'plasma')
colmap.set_array(np.linspace(v_min, v_max))
cbar = map_ani.colorbar(colmap, label = 'Ozone / ppb')

lats = np.array([])
longs = np.array([])
for gas in gas_arr_full:
    lats = np.append(lats, gas.latitude)
    longs = np.append(longs, gas.longitude)

x, y = m(longs, lats)
sctr = m.scatter(x, y, s=[], alpha = 0.8, zorder=10)
time_text = plt.text(0, 0, '',  fontsize = 20)

def scatter_find(i):
    concs = np.array([])
    cols = []
    for gas in gas_arr_full:
        concs = np.append(concs, gas.c[i])
        c = colmap.to_rgba(gas.c[i])
        cols.append(c[:3])
        
    concs *= concs
    sctr.set_sizes(concs)
    #using the last object in the list...
    sctr.set_color(cols)
    time_text.set_text(str(gas.t[i].month) + '/' + str(gas.t[i].year))
    return sctr, time_text


def scatter_init():
    sctr.set_sizes([])
    sctr.set_color([])
    time_text.set_text('')
    return sctr, time_text

frms = len(MH.t)
ani = animation.FuncAnimation(map_ani, scatter_find, init_func = scatter_init,
                              frames = frms, interval=10, blit = False,
                              repeat_delay = 5000)

ani.save('NorthAtlanticOzone.gif', fps=12, writer='imagemagick')
plt.show()
###############################################################################