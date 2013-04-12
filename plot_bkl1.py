#!/usr/local/bin/python

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from scipy import interpolate, signal
import matplotlib.font_manager as fm
import os
import urllib2
if not os.path.exists('Humor-Sans.ttf'):
    fhandle = urllib2.urlopen('http://antiyawn.com/uploads/Humor-Sans.ttf')
    open('Humor-Sans.ttf', 'wb').write(fhandle.read())
def xkcd_line(x, y, xlim=None, ylim=None, mag=1.0, f1=30, f2=0.05, f3=15):
    x = np.asarray(x)
    y = np.asarray(y)

    if xlim is None:
        xlim = (x.min(), x.max())
    if ylim is None:
        ylim = (y.min(), y.max())

    if xlim[1] == xlim[0]:
        xlim = ylim
    
    if ylim[1] == ylim[0]:
        ylim = xlim

    x_scaled = (x-xlim[0])*1./(xlim[1]-xlim[0])
    y_scaled = (y-ylim[0])*1./(ylim[1]-ylim[0])

    dx = x_scaled[1:] - x_scaled[:-1]
    dy = y_scaled[1:] - y_scaled[:-1]
    dist_tot = np.sum(np.sqrt(dx * dx + dy * dy))

    Nu = int(200*dist_tot)
    u = np.arange(-1, Nu + 1)*1./(Nu - 1)

    k = min(3, len(x) - 1)
    res = interpolate.splprep([x_scaled, y_scaled], s=0, k=k)
    x_int, y_int = interpolate.splev(u, res[0])

    dx = x_int[2:] - x_int[:-2]
    dy = y_int[2:] - y_int[:-2]
    dist = np.sqrt(dx*dx+dy*dy)

    coeffs = mag * np.random.normal(0, 0.01, len(x_int) - 2)
    b = signal.firwin(f1, f2 * dist_tot, window=('kaiser', f3))
    response = signal.lfilter(b, 1, coeffs)

    x_int[1:-1] += response * dy / dist
    y_int[1:-1] += response * dx / dist

    return x_int, y_int

def XKCDify(ax, mag=1.0, f1=50, f2=0.01, f3=15, bgcolor='w', xaxis_loc=None, yaxis_loc=None, xaxis_arrow='+', yaxis_arrow='+', ax_extend=0.1, expand_axes=False):
    ext = ax.get_window_extent().extents
    aspect = (ext[3] - ext[1]) / (ext[2] - ext[0])

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    xspan = xlim[1] - xlim[0]
    yspan = ylim[1] - xlim[0]

    xax_lim = (xlim[0] - ax_extend * xspan, xlim[1] + ax_extend * xspan)
    yax_lim = (ylim[0] - ax_extend * yspan, ylim[1] + ax_extend * yspan)

    if xaxis_loc is None:
        xaxis_loc = ylim[0]

    if yaxis_loc is None:
        yaxis_loc = xlim[0]
    xaxis = pl.Line2D([xax_lim[0], xax_lim[1]], [xaxis_loc, xaxis_loc], linestyle='-', color='k')
    yaxis = pl.Line2D([yaxis_loc, yaxis_loc], [yax_lim[0], yax_lim[1]], linestyle='-', color='k')
    ax.text(xax_lim[1], xaxis_loc - 0.02 * yspan, ax.get_xlabel(), fontsize=14, ha='right', va='top', rotation=12)
    ax.text(yaxis_loc - 0.02 * xspan, yax_lim[1], ax.get_ylabel(), fontsize=14, ha='right', va='top', rotation=78)
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.text(0.5 * (xax_lim[1] + xax_lim[0]), yax_lim[1], ax.get_title(), ha='center', va='bottom', fontsize=16)
    ax.set_title('')

    Nlines = len(ax.lines)
    lines = [xaxis, yaxis] + [ax.lines.pop(0) for i in range(Nlines)]

    for line in lines:
        x, y = line.get_data()

        x_int, y_int = xkcd_line(x, y, xlim, ylim, mag, f1, f2, f3)
        lw = line.get_linewidth()
        line.set_linewidth(2 * lw)
        line.set_data(x_int, y_int)

        if (line is not xaxis) and (line is not yaxis):
            line_bg = pl.Line2D(x_int, y_int, color=bgcolor, linewidth=8 * lw)

            ax.add_line(line_bg)
        ax.add_line(line)

    arr1 = 0.03 * np.array([-1, 0, -1])
    arr2 = 0.02 * np.array([-1, 0, 1])

    arr1[::2] += np.random.normal(0, 0.005, 2)
    arr2[::2] += np.random.normal(0, 0.005, 2)

    x, y = xaxis.get_data()
    if '+' in str(xaxis_arrow):
        ax.plot(x[-1] + arr1 * xspan * aspect, y[-1] + arr2 * yspan, color='k', lw=2)
    if '-' in str(xaxis_arrow):
        ax.plot(x[0] - arr1 * xspan * aspect, y[0] - arr2 * yspan, color='k', lw=2)

    x, y = yaxis.get_data()
    if '+' in str(yaxis_arrow):
        ax.plot(x[-1] + arr2 * xspan * aspect, y[-1] + arr1 * yspan, color='k', lw=2)
    if '-' in str(yaxis_arrow):
        ax.plot(x[0] - arr2 * xspan * aspect, y[0] - arr1 * yspan, color='k', lw=2)
    prop = fm.FontProperties(fname='Humor-Sans.ttf', size=16)
    for text in ax.texts:
        text.set_fontproperties(prop)
    leg = ax.get_legend()
    if leg is not None:
        leg.set_frame_on(False)
                                
        for child in leg.get_children():
            if isinstance(child, pl.Line2D):
                x, y = child.get_data()
                child.set_data(xkcd_line(x, y, mag=10, f1=100, f2=0.001))
                child.set_linewidth(2 * child.get_linewidth())
            if isinstance(child, pl.Text):
                child.set_fontproperties(prop)

    ax.set_xlim(xax_lim[0] - 0.1 * xspan, xax_lim[1] + 0.1 * xspan)
    ax.set_ylim(yax_lim[0] - 0.1 * yspan, yax_lim[1] + 0.1 * yspan)
    ax.set_xticks([])
    ax.set_yticks([])

    if expand_axes:
        ax.figure.set_facecolor(bgcolor)
        ax.set_axis_off()
        ax.set_position([0, 0, 1, 1])

    return ax
tlsb = []
nmb = []
origb = []
name = []

def readdata(filename, data):
    with open(filename, 'r') as f:
        line = f.readline()
        while line:
            data.append(float(line)) 
            line = f.readline()
    f.closed
    return

def readname(filename, data):
    with open(filename, 'r') as f:
        line = f.readline()
        while line:
            data.append(line)
            line = f.readline()
    f.closed
    return

readdata('tlsbkl', tlsb)
readdata('nmbkl', nmb)
readdata('origbkl', origb)

readname('list', name)

#for i in range(0, len(bkl_nm)):
#    bkl_tls[i] =2*(bkl_tls[i] - bkl_nm[i])/(bkl_tls[i]+bkl_nm[i])
np.random.seed(0)
ax = pl.axes()
x = np.arange(0, len(tlsb), 1)
ax.plot(x, tlsb, c='red', label='tls bkl')
ax.plot(x, nmb, c='c',ls='--',label='nm bkl')
ax.plot(x, origb, c='blue',ls=':', label='origin bkl')


ax.set_title('bkl distribution')
ax.set_xlabel('protein No.')
ax.set_ylabel('bkl')
ax.legend(loc='upper right')
ax.set_xlim(0, len(tlsb)-1)
ax.set_ylim(0.0, 14.0)
#ax.set_xlim(0, 0.4)
#ax.set_ylim(-0.01, 0.01)
xtickName = ax.set_xticklabels(name)
plt.setp(xtickName, rotation = 45, fontsize=11)
ax.set_xticks(x)

plt.show()

