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



import calendar
import time
from datetime import date,datetime


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
           
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/" )

############trajectory files#######
df_158 = pd.read_csv("traj_M158.csv", sep=',')
df_181 = pd.read_csv("traj_M181.csv", sep=',')
df_float22 = pd.read_csv("traj_float2022.csv", sep=',')



date_158=[]
date_181=[]
date_float22=[]

for k in range(0,df_158.Lon.size):
    python_date = datetime.fromordinal(int(df_158.Time[k])) + timedelta(days=df_158.Time[k]%1) - timedelta(days = 366)
    date_158.append(python_date)
    
for k in range(0,df_181.Lon.size):
    python_date = datetime.fromordinal(int(df_181.Time[k])) + timedelta(days=df_181.Time[k]%1) - timedelta(days = 366)
    date_181.append(python_date)

for k in range(0,df_float22.Lon.size):
    python_date = datetime.fromordinal(int(df_float22.Time[k])) + timedelta(days=df_float22.Time[k]%1) - timedelta(days = 366)
    date_float22.append(python_date)

###############################
######1- SST data #######
###############################

dataset = xr.open_dataset('OI_SST_TROPATL.nc')

#I extract lon/lat
LAT_SST = dataset.variables['LAT'].values
LON_SST = dataset.variables['LON'].values
SST_TROPATL = np.squeeze(dataset.variables['SST_TROPATL'].values)
SST_TROPATL_FILTERED=np.squeeze(dataset.variables['SST_TROPATL_FILTERED'].values)
TIME_SST=dataset.variables['TIME'].values


#put a loop to open each SST_TROPATL, and squeeze bcs 3D variable
#take lon, lat, time
#should maybe either conduct a spatial mean between -8 and 8 
# or more restrained 


mean_sst_tropatl_filtered=np.zeros((TIME_SST.size,LON_SST.size))
mean_sst_tropatl=np.zeros((TIME_SST.size,LON_SST.size))

#transform matlab to python datetime
date_sst=[]


for k in range(0,TIME_SST.size):
    #I filter the data based on the lat 
    #option 1 conduct the mean between -8 and 8 
   # sst_tropatl=np.squeeze(SST_TROPATL[k,:,:])
   # sst_tropatl_filtered=np.squeeze(SST_TROPATL_FILTERED[k,:,:])
    
    sst_tropatl=np.squeeze(SST_TROPATL[k,28:35,:])
    sst_tropatl_filtered=np.squeeze(SST_TROPATL_FILTERED[k,28:35,:])
    
    mean_sst_tropatl[k]=np.nanmean(sst_tropatl,axis=0) #for mean based on row= axis=1
    mean_sst_tropatl_filtered[k]=np.nanmean(sst_tropatl_filtered,axis=0) #for mean based on row= axis=1
    python_date = datetime.fromordinal(int(TIME_SST[k])) + timedelta(days=TIME_SST[k]%1) - timedelta(days = 366)
    date_sst.append(python_date)
    
    
#dates of the beginning and end of the export 
#lon158_2[10 lon158_2[22] 2021-8-8 2021-9-8 

#lon158_2[60  lon158_2[77] 2021-12-13 2022-1-26



date_list1 = [
    datetime.strptime("2021-8-8", "%Y-%m-%d"),
    datetime.strptime("2021-9-8", "%Y-%m-%d")]



date_list2 = [
    datetime.strptime("2021-12-13", "%Y-%m-%d"),
    datetime.strptime("2022-1-26", "%Y-%m-%d")]
    
##################
###CHLA
###################
    
dataset = xr.open_dataset('CHLA_TROPATL.nc')

#I extract lon/lat
LAT_SST = dataset.variables['LAT'].values
LON_CHL = dataset.variables['LON'].values
CHL_TROPATL = np.squeeze(dataset.variables['CHLA_TROPATL'].values)
CHL_TROPATL_FILTERED=np.squeeze(dataset.variables['CHLA_TROPATL_FILTERED'].values)
TIME_SST=dataset.variables['TIME'].values

#put a loop to open each SST_TROPATL, and squeeze bcs 3D variable
#take lon, lat, time
#should maybe either conduct a spatial mean between -8 and 8 
# or more restrained 


mean_chl_tropatl_filtered=np.zeros((TIME_SST.size,LON_CHL.size))
mean_chl_tropatl=np.zeros((TIME_SST.size,LON_CHL.size))
#transform matlab to python datetime
date_chl=[]


for k in range(0,TIME_SST.size):
    #I filter the data based on the lat 
    #option 1 conduct the mean between -5 and 5 
    # chl_tropatl=np.squeeze(CHL_TROPATL[k,:,:])
    # chl_tropatl_filtered=np.squeeze(CHL_TROPATL_FILTERED[k,:,:])
    
    chl_tropatl=np.squeeze(CHL_TROPATL[k,102:139,:])     #for 1N and 1S
    chl_tropatl_filtered=np.squeeze(CHL_TROPATL_FILTERED[k,102:139,:])
    
    mean_chl_tropatl[k]=np.nanmean(chl_tropatl,axis=0) #for mean based on row= axis=1
    mean_chl_tropatl_filtered[k]=np.nanmean(chl_tropatl_filtered,axis=0) #for mean based on row= axis=1
    python_date = datetime.fromordinal(int(TIME_SST[k])) + timedelta(days=TIME_SST[k]%1) - timedelta(days = 366)
    date_chl.append(python_date)


###############################
###### #######
###############################
    
jet = plt.get_cmap('Spectral')

   
    
#add contouring level 
#add contouring level 
levels = 20
min_contour_level = -0.38
max_contour_level = 0.38
contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)


width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = datetime(2021, 7, 1, 0, 0),datetime(2022, 4, 1, 0, 0)
fig = plt.figure(1, figsize=(3,5))
ax=fig.add_subplot(412)
p1 = plt.contourf(date_sst, LON_SST, mean_sst_tropatl_filtered.transpose(),contour_levels, cmap=cm, alpha=1,vmin=-0.38, vmax=0.38, extend='both')
p4=plt.plot(date_float22,df_float22.Lon, color='black')

ax.set_ylim([-40, 10])
ax.set_xlim([set_ylim_lower, set_ylim_upper])

loc = mdates.MonthLocator(interval=1)
ax.xaxis.set_major_locator(loc)
fmt = mdates.DateFormatter('%b\n%y')
ax.xaxis.set_major_formatter(fmt)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(b)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')
ax.set_ylabel('Longitude (°)', fontsize=7)
ax.get_xaxis().set_ticklabels([])

#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('SST anomaly(°C)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)

for p1 in p1.collections:
    p1.set_rasterized(True)
    


jet = plt.get_cmap('Spectral_r')


    
#add contouring level 
#add contouring level 
levels = 100
min_contour_level =22
max_contour_level = 30.5
contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)


width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = datetime(2021, 7, 1, 0, 0),datetime(2022, 4, 1, 0, 0)
fig = plt.figure(1, figsize=(3,5))
ax=fig.add_subplot(411)
p1 = plt.contourf(date_sst, LON_SST, mean_sst_tropatl.transpose(),contour_levels, cmap=jet, alpha=1,vmin=22, vmax=30.5, extend='both')
p4=plt.plot(date_float22,df_float22.Lon, color='black')

ax.set_ylim([-40, 10])
ax.set_xlim([set_ylim_lower, set_ylim_upper])

loc = mdates.MonthLocator(interval=1)
ax.xaxis.set_major_locator(loc)
fmt = mdates.DateFormatter('%b\n%y')
ax.xaxis.set_major_formatter(fmt)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(a)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')
ax.set_ylabel('Longitude (°)', fontsize=7)
ax.get_xaxis().set_ticklabels([])

#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('SST (°C)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)

for p1 in p1.collections:
    p1.set_rasterized(True)
    


#add contouring level 
levels = 20
min_contour_level = 0
max_contour_level = 1
contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)


jet=cmocean.cm.algae

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = datetime(2021, 7, 1, 0, 0),datetime(2022, 4, 1, 0, 0)
ax=fig.add_subplot(413)
p1 = plt.contourf(date_chl, LON_CHL, mean_chl_tropatl.transpose(),contour_levels, cmap=jet, alpha=1,vmin=0, vmax=1,  extend='both')
p4=plt.plot(date_float22,df_float22.Lon , color='black')

ax.set_ylim([-40, 10])
ax.set_xlim([set_ylim_lower, set_ylim_upper])

loc = mdates.MonthLocator(interval=1)
ax.xaxis.set_major_locator(loc)
fmt = mdates.DateFormatter('%b\n%Y')
ax.xaxis.set_major_formatter(fmt)
ax.get_xaxis().set_ticklabels([])

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(c)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')
ax.set_ylabel('Longitude (°)', fontsize=7)

#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('Chl-a \n (mg m$^{-3}$)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)
cbar.ax.tick_params(labelsize=7) 
ax.get_xaxis().set_ticklabels([])


y_min, y_max = ax.get_ylim()

# Calculate ymin and ymax fractions
ymin_fraction = (-28 - y_min) / (y_max - y_min)
ymax_fraction = (0 - y_min) / (y_max - y_min)

for date22 in date_list1:
    ax.axvline(x=date22, color='r', linestyle='--', linewidth=1, ymin=ymin_fraction, ymax=ymax_fraction)
    
for date22 in date_list2:
    ax.axvline(x=date22, color='b', linestyle='--', linewidth=1, ymin=ymin_fraction, ymax=ymax_fraction)

for p1 in p1.collections:
    p1.set_rasterized(True)
    
    
#add contouring level 
levels = 20
min_contour_level = -0.1
max_contour_level = 0.1
contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)


jet = plt.get_cmap('BrBG')

width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = datetime(2021, 7, 1, 0, 0),datetime(2022, 4, 1, 0, 0)
ax=fig.add_subplot(414)
p1 = plt.contourf(date_chl, LON_CHL, mean_chl_tropatl_filtered.transpose(),contour_levels, cmap=jet, alpha=1,vmin=-0.1, vmax=0.1,  extend='both')
p4=plt.plot(date_float22,df_float22.Lon , color='black')

ax.set_ylim([-40, 10])
ax.set_xlim([set_ylim_lower, set_ylim_upper])

loc = mdates.MonthLocator(interval=1)
ax.xaxis.set_major_locator(loc)
fmt = mdates.DateFormatter('%b\n%y')
ax.xaxis.set_major_formatter(fmt)


plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.2, 1.15, '(d)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')
ax.set_ylabel('Longitude (°)', fontsize=7)

#colorbar
cbar = plt.colorbar(p1)
cbar.ax.set_ylabel('Chl-a anomaly \n (mg m$^{-3}$)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)
cbar.ax.tick_params(labelsize=7) 



plt.xticks(fontsize=7)
for p1 in p1.collections:
    p1.set_rasterized(True)
    

    

    

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float//new_version_052024/" )
fig_name_pdf = ("Fig01" + ".png")
plt.savefig(fig_name_pdf,dpi=300,bbox_inches="tight")