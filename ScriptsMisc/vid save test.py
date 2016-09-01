# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:46:57 2016

@author:
http://stackoverflow.com/questions/23074484/cannot-save-matplotlib-animation-with-ffmpeg
"""


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
plt.rcParams['animation.ffmpeg_path'] = '/home/tmed2/ffmpeg-3.1.1/ffmpeg-3.1.1/ffmpeg'

fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = np.linspace(0, 2, 1000)
    y = np.sin(2 * np.pi * (x - 0.01 * i))
    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=20, blit=True)
plt.show()
FFwriter = animation.FFMpegWriter()
anim.save('basic_animation.mp4', writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'])