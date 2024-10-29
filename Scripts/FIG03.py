#same script as subplots in ARTIVLE but I added the particule abundance 

import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from datetime import timedelta
import xarray as xr
import cmocean
import matplotlib.dates as mdates
from pathlib import Path

cm=cmocean.cm.balance

def contour_levels_func(min_contour_level, max_contour_level, levels):
    """Function to define contour levels for contourf"""
    distance_levels = max_contour_level / levels
    contour_levels = np.arange(min_contour_level, max_contour_level, distance_levels)
    return contour_levels


def nonlinear_colormap():
    import pylab as pyl
    #import numpy as np
    levels1 = [0, 1, 2]

    ####################################################################
    ###                      non linear colormap                     ###
    ####################################################################

    """
    nlcmap - a nonlinear cmap from specified levels

    Copyright (c) 2006-2007, Robert Hetland <hetland@tamu.edu>
    Release under MIT license.

    Some hacks added 2012 noted in code (@MRR)
    """

    from matplotlib.colors import LinearSegmentedColormap


    class nlcmap(LinearSegmentedColormap):
        """Nonlinear colormap.
           Needs the input of a linear colormap e.g. pylab.cm.jet
           and a list of levels e.g. [0,1,2,3,6]
           in order to change the interpolation, change the levels vector,
           e.g ad 20, 50 to make it more unlinear."""
        import numpy as np
        name = 'nlcmap'

        def __init__(self, cmap, levels):
            import numpy as np
            self.cmap = cmap
            # @MRR: Need to add N for backend
            self.N = cmap.N
            self.monochrome = self.cmap.monochrome
            self.levels = np.asarray(levels, dtype='float64')
            self._x = self.levels / self.levels.max()
            self._y = np.linspace(0.0, 1.0, len(self.levels))

        #@MRR Need to add **kw for 'bytes'
        def __call__(self, xi, alpha=1.0, **kw):
            import numpy as np
            """docstring for fname"""
            # @MRR: Appears broken?
            # It appears something's wrong with the
            # dimensionality of a calculation intermediate
            #yi = stineman_interp(xi, self._x, self._y)
            yi = np.interp(xi, self._x, self._y)
            return self.cmap(yi, alpha)


    cmap_nonlin = nlcmap(pyl.cm.CMRmap, levels1)
    return cmap_nonlin



cm1 = nonlinear_colormap()

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )

M158_lon=pd.read_csv("Time_float100.csv", sep=',')
M158_depth=pd.read_csv("DEPTH_float100.csv", sep=',')
M158_Cflux=pd.read_csv("FluxC_float100.csv", sep=',')
M158_Mip=pd.read_csv("Mip_float100.csv", sep=',')
M158_Map=pd.read_csv("Map_float100.csv", sep=',')

lon158=np.squeeze(np.array(M158_lon))
depth158=np.squeeze(np.array(M158_depth))

#lon158=mdates.num2date(lon158)

#mean or median of the data to try detecting anomaly 
#trying map
mean_158map=np.mean(M158_Map,1) #200 levels of depth 
mean_158mip=np.mean(M158_Mip,1) #200 levels of depth 
mean_158C=np.mean(M158_Cflux,1) #200 levels of depth 

#value-mean to conduct anomaly 

New_value_map=M158_Map.sub(mean_158map,axis=0)
New_value_mip=M158_Mip.sub(mean_158mip,axis=0)
New_value_C=M158_Cflux.sub(mean_158C,axis=0)



# Convert numeric dates to datetime objects
reference_date = datetime(1970, 1, 1)  # Choose your reference date
real_dates = [reference_date + timedelta(days=date) for date in lon158]




levels = 20
min_contour_level = 0
max_contour_level = 30

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)

jet = plt.get_cmap('BrBG')

width, height = 0.8, 0.7
fig = plt.figure(1, figsize=(3,5))
ax=fig.add_subplot(615)
p1 = plt.contourf(real_dates, depth158 ,M158_Cflux, contour_levels, cmap=cm1, alpha=1,extend='both')

ax.set_ylim([-2000, 0])


ax.get_xaxis().set_ticklabels([])
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(e)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')
ax.set_ylabel('Depth (m)', fontsize=7)

#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('C flux \n (mgC m$^{-3}$)', fontsize=5)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)

for p1 in p1.collections:
    p1.set_rasterized(True)





### anomaly 1
levels = 15
min_contour_level = -4
max_contour_level = 4

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)
ax=fig.add_subplot(616)

jet = plt.get_cmap('BrBG')
p1 = plt.contourf(real_dates, depth158 ,New_value_C, contour_levels, cmap=cm, alpha=1,extend='both')
ax.set_ylim([-2000, 0])
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(f)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')


ax.set_ylabel('Depth (m)', fontsize=7)
#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('Cflux \n anomaly \n (mgC m$^{-3}$)', fontsize=5)
cbar.ax.locator_params(nbins=5)
cbar.ax.tick_params(labelsize=7) 

# from matplotlib.dates import DateFormatter
# # Define the date format
# date_form = DateFormatter('%b %y')
# ax.xaxis.set_major_formatter(date_form)
# Customize the x-axis labels
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %y')) # %d


plt.xticks(rotation=90,fontsize=7)


    
for p1 in p1.collections:
    p1.set_rasterized(True)





levels = 100
min_contour_level = 0
max_contour_level = 70

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)
ax=fig.add_subplot(611)


p1 = plt.contourf(real_dates, depth158 ,M158_Mip, contour_levels, cmap=cm1, alpha=1,extend='both')
ax.set_ylim([-2000, 0])
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(a)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')

ax.set_ylabel('Depth (m)', fontsize=7)

#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('Mip \n (# L$^{-1}$)', fontsize=5)
cbar.ax.locator_params(nbins=5)
cbar.ax.tick_params(labelsize=7) 
ax.get_xaxis().set_ticklabels([])

for p1 in p1.collections:
    p1.set_rasterized(True)



### anomaly 1
levels = 15
min_contour_level = -4
max_contour_level = 4

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)
ax=fig.add_subplot(612)

jet = plt.get_cmap('BrBG')
p1 = plt.contourf(real_dates, depth158 ,New_value_mip, contour_levels, cmap=cm, alpha=1,extend='both')
ax.set_ylim([-2000, 0])
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(b)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')


ax.set_ylabel('Depth (m)', fontsize=7)
#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('Mip \n anomaly \n (# L$^{-1}$)', fontsize=5)
cbar.ax.locator_params(nbins=5)
cbar.ax.tick_params(labelsize=7) 
ax.get_xaxis().set_ticklabels([])

for p1 in p1.collections:
    p1.set_rasterized(True)









levels = 15
min_contour_level = 0
max_contour_level = 1.5

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)

jet = plt.get_cmap('BrBG')

width, height = 0.8, 0.7

ax=fig.add_subplot(613)
p1 = plt.contourf(real_dates, depth158 ,M158_Map, contour_levels, cmap=cm1, alpha=1,extend='both')
    
ax.set_ylim([-2000, 0])
ax.set_ylabel('Depth (m)', fontsize=7)

ax.get_xaxis().set_ticklabels([])
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(c)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')


#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('Map \n  (# L$^{-1}$)', fontsize=5)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)

for p1 in p1.collections:
    p1.set_rasterized(True)
    

### anomaly 1
levels = 15
min_contour_level = -0.1
max_contour_level = 0.1

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)
ax=fig.add_subplot(614)

jet = plt.get_cmap('BrBG')
p1 = plt.contourf(real_dates, depth158 ,New_value_map, contour_levels, cmap=cm, alpha=1,extend='both')
ax.set_ylim([-2000, 0])
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(d)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')


ax.set_ylabel('Depth (m)', fontsize=7)
#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('Map \n anomaly  \n (# L$^{-1}$)', fontsize=5)
cbar.ax.locator_params(nbins=5)
cbar.ax.tick_params(labelsize=7) 
ax.get_xaxis().set_ticklabels([])

#eqution 2, 2 sinking speed 
lon158_2=np.squeeze(np.array(M158_lon))
b=lon158_2[10]*30
line_1=lon158_2[10:40]
Y1=-30*line_1+b
plt.plot(line_1, Y1,'r')

#equation 2
b=lon158_2[22]*30
line_2=lon158_2[22:66]
Y2=-30*line_2+b
plt.plot(line_2, Y2,'r')


#first event 
lon158_2=np.squeeze(np.array(M158_lon))
b=lon158_2[77]*30
line_4=lon158_2[77:100]
Y4=-30*line_4+b
plt.plot(line_4, Y4,'b')



#equation 2
b=lon158_2[60]*30  #speed of 28 is better
line_3=lon158_2[60:87]
Y3=-30*line_3+b
plt.plot(line_3, Y3,'b')



for p1 in p1.collections:
    p1.set_rasterized(True)


    
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("Fig03" + ".png")
plt.savefig(fig_name_pdf,dpi=300,bbox_inches="tight")


mean_flux_between=np.nanmean(M158_Map,1)