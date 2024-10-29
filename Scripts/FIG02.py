#I will use pierre data
import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from datetime import timedelta
import xarray as xr
import cmocean
from scipy.interpolate import griddata
import matplotlib.dates as mdates
from pathlib import Path
from seawater import dpth
import calendar
import time
from datetime import date,datetime
import matplotlib.dates as dates

import gsw
from oceanpy import mixed_layer_depth

from scipy.signal import savgol_filter


cmpa=cmocean.cm.algae
cmoxy=cmocean.cm.balance

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



cm = nonlinear_colormap()


def contour_levels_func(min_contour_level, max_contour_level, levels):
    """Function to define contour levels for contourf"""
    distance_levels = max_contour_level / levels
    contour_levels = np.arange(min_contour_level, max_contour_level, distance_levels)
    return contour_levels


def gridding_func(pos_min_max, depth_min_max, pos_array, depth, param):
    grid_method = "linear"  # choose the gridding method here
    # method to do the regridding, can be nearest (for nearest neighbour) or linear

    xi = np.linspace(min(pos_min_max), max(pos_min_max), 5000)
    yi = np.linspace(min(depth_min_max), max(depth_min_max), 200)
    zi = griddata((pos_array, depth), param, (xi[None, :], yi[:, None]), method=grid_method)
    return xi, yi, zi

Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()

#data ctd between 8 and 8 
#it should be filtered using the longitude criteria 
#1.CTD
df_ = pd.read_csv(Path_to_data / "Ecopart_diagnostics_data_605.tsv",sep='\t', low_memory=False)

print(df_.columns.tolist()) #determine the variables we want 

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
# select only ascent profiles :
df_ = df_[df_['Profile'].str.contains('a')]    


df_ =df_[['RAWfilename', 'Profile', 'Latitude', 'Longitude', 'Date_Time', 'Pressure [dbar]',
          'Depth [m]', 'Potential density [kg/m3]', 'Temperature [degrees Celsius]', 
          'Practical salinity [psu]', 'Doxy [micromol/kg]', 'Chlorophyll-a [mg/m3]', 
          'bbp POC [mgC/m3]', 'bbp POC Koestner [mgC/m3]']]


# I extract the data
lon=np.array(df_['Longitude']);
lat=np.array(df_['Latitude'])
Date_Time=np.array(df_['Date_Time'])


pressure=-np.array(df_['Pressure [dbar]'])
chl=np.array(df_['Chlorophyll-a [mg/m3]'])
sal=np.array(df_['Practical salinity [psu]'])
temp=np.array(df_['Temperature [degrees Celsius]'])
bbp=np.array(df_['bbp POC [mgC/m3]'])
oxy=np.array(df_['Doxy [micromol/kg]'])


# # I convert the dates to float values (in seconds from 1970 1 1)
Date_Num=np.r_[0:lon.size]
# for i in Date_Num:
#     date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%SZ')
#     Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
#     #datetime.utcfromtimestamp(Date_Num[i])


for i in Date_Num:
    date_time_obj = datetime.strptime(Date_Time[i], '%Y-%m-%dT%H:%M:%S')
    Date_Num[i]= dates.date2num(date_time_obj)

CHL_filtered=np.array([]);pressure_chl=np.array([]);Date_Num_chl=np.array([])
SAL_filtered=np.array([]);pressure_sal=np.array([]);Date_Num_sal=np.array([])
TEMP_filtered=np.array([]);pressure_temp=np.array([]);Date_Num_temp=np.array([])
BBP_filtered=np.array([]);pressure_bbp=np.array([]);Date_Num_bbp=np.array([])
OXY_filtered=np.array([]);pressure_oxy=np.array([]);Date_Num_oxy=np.array([])


# plt.scatter(Date_Time,pressure,chl)
# plt.xticks(rotation=90,fontsize=7)


list_dates=np.unique(Date_Num)


for i in range(0,list_dates.size):
    sel=Date_Num==list_dates[i];x=Date_Num[sel];y=pressure[sel]
    z=chl[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
        z=savgol_filter(z,5,1)
        CHL_filtered = np.concatenate((CHL_filtered, z));Date_Num_chl = np.concatenate((Date_Num_chl, x2));pressure_chl = np.concatenate((pressure_chl, y2))
    
    z=sal[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
         z=savgol_filter(z,5,1)
         SAL_filtered = np.concatenate((SAL_filtered, z));Date_Num_sal = np.concatenate((Date_Num_sal, x2));pressure_sal = np.concatenate((pressure_sal, y2))
    
    
    z=temp[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
         z=savgol_filter(z,5,1)
         TEMP_filtered = np.concatenate((TEMP_filtered, z));Date_Num_temp = np.concatenate((Date_Num_temp, x2));
         pressure_temp = np.concatenate((pressure_temp, y2))
    
    z=bbp[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
         z=savgol_filter(z,5,1)
         BBP_filtered = np.concatenate((BBP_filtered, z));Date_Num_bbp= np.concatenate((Date_Num_bbp, x2));pressure_bbp = np.concatenate((pressure_bbp, y2))

    z=oxy[sel];sel2=~np.isnan(z);z=z[sel2];x2=x[sel2];y2=y[sel2]
    if sum(sel2)>0:
         z=savgol_filter(z,5,1)
         OXY_filtered = np.concatenate((OXY_filtered, z));Date_Num_oxy= np.concatenate((Date_Num_oxy, x2));pressure_oxy = np.concatenate((pressure_oxy, y2))



# I define the x and y arrays for the contourf plot
x_chl = np.linspace(Date_Num_chl.min(),Date_Num_chl.max(),100)
y_chl = np.linspace(pressure_chl.min(),pressure_chl.max(),50)
x_chl_g,y_chl_g=np.meshgrid(x_chl,y_chl)


x_sal = np.linspace(Date_Num_sal.min(),Date_Num_sal.max(),100)
y_sal = np.linspace(pressure_sal.min(),pressure_sal.max(),50)
x_sal_g,y_sal_g=np.meshgrid(x_sal,y_sal)

x_temp = np.linspace(Date_Num_temp.min(),Date_Num_temp.max(),100)
y_temp = np.linspace(pressure_temp.min(),pressure_temp.max(),50)
x_temp_g,y_temp_g=np.meshgrid(x_temp,y_temp)


x_bbp = np.linspace(Date_Num_bbp.min(),Date_Num_bbp.max(),100)
y_bbp = np.linspace(pressure_bbp.min(),pressure_bbp.max(),50)
x_bbp_g,y_bbp_g=np.meshgrid(x_bbp,y_bbp)



x_oxy = np.linspace(Date_Num_oxy.min(),Date_Num_oxy.max(),100)
y_oxy = np.linspace(pressure_oxy.min(),pressure_oxy.max(),50)
x_oxy_g,y_oxy_g=np.meshgrid(x_oxy,y_oxy)


# I interpolate
CHL_interp = griddata((Date_Num_chl,pressure_chl), CHL_filtered, (x_chl_g, y_chl_g), method="nearest")
SAL_interp = griddata((Date_Num_sal,pressure_sal), SAL_filtered, (x_sal_g, y_sal_g), method="nearest")
TEMP_interp = griddata((Date_Num_temp,pressure_temp), TEMP_filtered, (x_temp_g, y_temp_g), method="nearest")
BBP_interp = griddata((Date_Num_bbp,pressure_bbp), BBP_filtered, (x_bbp_g, y_bbp_g), method="nearest")
OXY_interp = griddata((Date_Num_oxy,pressure_oxy), OXY_filtered, (x_oxy_g, y_oxy_g), method="nearest")


# Calculate the mean of 'chl' after removing the specified date range
mean_chl_after_removal = np.nanmean(CHL_interp[94:100], axis=0)


#5-The mixed layer depth 
#############################
#transform pressure to depth 
#compute absolute salinity and conservative temperature for the mixed layer based on CTD input 
#Variables


floa=np.array(df_['Profile'])
depth = np.array(dpth(np.array(pressure),lat)) #transform pressure into depth


#computing  potential density  using salinity and temperature
abs_psal_tmp = gsw.SA_from_SP(sal, pressure, lon, lat)  # compute absolute salinity
cons_tmp = gsw.CT_from_t(abs_psal_tmp, temp, pressure)    # compute conservative temperature
dens_tmp = gsw.density.sigma0(abs_psal_tmp, cons_tmp) #compute density from the calculated variables
dens_u=dens_tmp+1000 # unité kg


#create an index 
#if k-1 =k --> same profile other it's the following profile and we add 1
v=len(floa)
eco_index=[1]
prof_id=np.array(floa)
i=1
 
#create new indexes so I can have each profile alone based on lon,lat?  
for k in range(1,v):
     if ((lon[k-1])==(lon[k]) and (lat[k-1])==(lat[k])):
         i=i
         eco_index.append(i)
     else:
            i=i+1
            eco_index.append(i)
            
eco_INDEX=np.array(eco_index)

n=np.unique(eco_INDEX)


# #computing mld for each profile
# ####################
 
float_final=[] 
# mld_x=np.zeros((n.size))

# for d in range(1,max(eco_INDEX)+1):
#     float_name=Date_Num[eco_INDEX==d]
#     float_name=float_name[0]
    
#     depth_tmp=-depth[eco_INDEX==d]
#     temp_tmp=cons_tmp[eco_INDEX==d]
#     pres_tmp=pressure[eco_INDEX==d]
#     dens_tmp1=dens_tmp[eco_INDEX==d]
#     sal_tmp=abs_psal_tmp[eco_INDEX==d]
    
#     depth_tmp=depth_tmp[~np.isnan(temp_tmp)] #remove nan values
#     temp_tmp = temp_tmp[~np.isnan(temp_tmp)]
#     sal_tmp=sal_tmp[~np.isnan(sal_tmp)]
    
#     if len(temp_tmp) == 0:
#         mld_x[d-1]=np.nan
#         float_final.append(float_name)
#     else:
#         mld,_=  mixed_layer_depth(depth_tmp,temp_tmp,using_temperature='yes')
#         mld_x[d-1]=mld #I have a list of mixed layer for correlation 
#         float_final.append(float_name)
# stop       
        
# #I will fill the two nan cells with the value of the cells close to it 
# mld_x[65]=13.5
# mld_x[66]=13.7    

# # #save mld
# # dff=pd.DataFrame(mld_x)
# # dff.to_csv(r"/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/mld.csv", index=False)
mld_x = pd.read_csv(r"/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/mld.csv")

for d in range(1,max(eco_INDEX)+1):
    float_name=Date_Num[eco_INDEX==d]
    float_name=float_name[0]
    
    float_final.append(float_name)
float_final.sort()
    
### climatology thermocline 
thermocline_x=pd.read_csv(r"/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/traj_float2022_with_z20.csv")
## remove duplicates

# Remove duplicate rows based on 'lon' and 'lat' columns
thermocline_x = thermocline_x.drop_duplicates(subset=['Lon', 'Lat'])

thermocline_x = thermocline_x.reset_index(drop=True)

# Apply a rolling mean with a window size of 3
thermocline_x2 = thermocline_x.rolling(window=3).mean()


mld_time = pd.read_csv(r"/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/mld_time.csv")

# Convert Unix timestamp to datetime object
mld_time['datetime_column'] = pd.to_datetime(mld_time['0'], unit='s')
# Convert datetime object to numeric representation using dates.date2num()
mld_time['numeric_representation'] = dates.date2num(mld_time['datetime_column'])

float_final=mld_time['numeric_representation']
#stop

#date I need to remove 2022-01-05  -- 2022-01-11date_object.timestamp()

date_time_obj = datetime.strptime('2022-01-05', '%Y-%m-%d')
date_min= dates.date2num(date_time_obj)
date_time_obj = datetime.strptime('2022-01-11', '%Y-%m-%d')
date_max = dates.date2num(date_time_obj)



#nan between this period for the chl and the bbp 

CHL_interp[(x_chl_g>date_min) & (x_chl_g<date_max) ]=np.nan
BBP_interp[(x_chl_g>date_min) & (x_chl_g<date_max) ]=np.nan


# Convert numeric dates to datetime objects
reference_date = datetime(1970, 1, 1)  # Choose your reference date
real_dates = [reference_date + timedelta(days=date1) for date1 in x_chl]


# Calculate max chlorophyll concentration between 0-100m for each date
# Assuming y_chl represents the depths and CHL_interp is the chlorophyll data
depth_range_mask = (y_chl <= -20) & (y_chl >= -100)
chl_data_within_range = CHL_interp[depth_range_mask, :]

# Find the depth of the maximum chlorophyll concentration for each time point
max_chl_depth_indices = np.argmax(chl_data_within_range, axis=0)
max_chl_depths = y_chl[depth_range_mask][max_chl_depth_indices]



#### plots
levels = 15
min_contour_level = 0
max_contour_level = 1.5

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels) 



jet=cmocean.cm.algae
####### I plot
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -200, 0
# Set x-axis limit until the last day of data


fig = plt.figure(1, figsize=(3,5))
ax=fig.add_subplot(514)
ax.set_ylim([set_ylim_lower, set_ylim_upper])
ax_1 = plot2 = plt.contourf(real_dates, y_chl,CHL_interp,contour_levels,cmap=jet,vmin=0,vmax=1.5)#, cmap=cmhot)
# Add the max chlorophyll line
# Plot the depth of maximum chlorophyll concentration as a line
#plt.plot(real_dates, max_chl_depths, color='red', linestyle='-', linewidth=0.8, label='Max Chl Depth (0-100m)')

# Plot contour line for temperature at 20 degrees Celsius
#contour_line = plt.contour(real_dates, y_temp, TEMP_interp, levels=[20], colors='black', linewidths=0.5)
ax.set_xlim(real_dates[0], real_dates[-1])
# draw colorbar
cbar = fig.colorbar(plot2)
cbar.ax.set_ylabel('Chl \n (mg $m^{-3}$)', fontsize=7)
plt.ylabel('Depth (m)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
ax.text(-0.3, 1.15, '(d)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')
ax.get_xaxis().set_ticklabels([])
plt.xticks(rotation=90,fontsize=7)
xticklabels=[]
# I add the grid
# plt.grid(color='k', linestyle='dashed', linewidth=0.5)
# ax.set_xticklabels(xticklabels)
# plt.xticks(rotation=90,fontsize=7)
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter(' %b %y')) # %d

plt.xticks(rotation=90,fontsize=7)



lon158_2=np.squeeze(np.array(x_chl))
b=lon158_2[10]*30
line_1=lon158_2[10:45]
Y1=-30*line_1+b
plt.plot(line_1, Y1,'--r')


#equation 2
b=lon158_2[22]*30
line_2=lon158_2[22:66]
Y2=-30*line_2+b
plt.plot(line_2, Y2,'--r')


#first event 
lon158_2=np.squeeze(np.array(x_chl))
b=lon158_2[77]*30
line_4=lon158_2[77:100]
Y4=-30*line_4+b
plt.plot(line_4, Y4,'--b')


#equation 2
b=lon158_2[60]*30  #speed of 28 is better
line_3=lon158_2[60:87]
Y3=-30*line_3+b
plt.plot(line_3, Y3,'--b')


#### plots
levels = 100
min_contour_level = 80
max_contour_level = 210

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels) 

####### I plot
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -200, 0

ax=fig.add_subplot(513)
ax.set_ylim([set_ylim_lower, set_ylim_upper])
ax_1 = plot2 = plt.contourf(real_dates, y_oxy,OXY_interp,contour_levels,vmin=80,vmax=210)#, cmap=cmhot)
#contour_line = plt.contour(real_dates, y_temp, TEMP_interp, levels=[20], colors='black', linewidths=0.5)
ax.set_xlim(real_dates[0], real_dates[-1])
# draw colorbar
cbar = fig.colorbar(plot2)
cbar.ax.set_ylabel('Oxygen\n' r'($\mu$mol/kg)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.ylabel('Depth (m)', fontsize=7)
ax.text(-0.3, 1.15, '(c)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')

# #I set xticks
# nxticks=9
# xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
# xticklabels=[]
# for i in xticks:
#     xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B %y'))
# ax.set_xticks(xticks)
ax.get_xaxis().set_ticklabels([])
plt.xticks(rotation=90,fontsize=7)
# I add the grid
#plt.grid(color='k', linestyle='dashed', linewidth=0.5)







#### TEMP
levels = 100
min_contour_level =14
max_contour_level = 27

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels) 




####### I plot
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -200, 0

ax=fig.add_subplot(512)
ax.set_ylim([set_ylim_lower, set_ylim_upper])
ax_1 = plot2 = plt.contourf(real_dates, y_temp,TEMP_interp,contour_levels,vmin=14,vmax=27, extend='both')#, cmap=cmhot)
contour_line = plt.contour(real_dates, y_temp, TEMP_interp, levels=[20], colors='black', linewidths=0.5)
#contour_line = plt.contour(real_dates, y_temp, TEMP_interp, levels=[20], colors='black', linewidths=0.5)
ax.set_xlim(real_dates[0], real_dates[-1])
# draw colorbar
cbar = fig.colorbar(plot2)
cbar.ax.set_ylabel('Temperature \n (°C)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)
plt.plot(float_final,thermocline_x2.Z20,'b')
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.ylabel('Depth (m)', fontsize=7)
ax.text(-0.3, 1.15, '(b)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')

# #I set xticks
# nxticks=9
# xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
# xticklabels=[]
# for i in xticks:
#     xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B %y'))
# ax.set_xticks(xticks)
ax.get_xaxis().set_ticklabels([])
plt.xticks(rotation=90,fontsize=7)
# I add the grid
#plt.grid(color='k', linestyle='dashed', linewidth=0.5)


plt.xticks(fontsize=7)
plt.yticks(fontsize=7)



####SALINITY
levels = 150
min_contour_level =33
max_contour_level = 37

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels) 


####### I plot
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -200, 0

ax=fig.add_subplot(511)
ax.set_ylim([set_ylim_lower, set_ylim_upper])
ax.set_xlim(real_dates[0], real_dates[-1])
ax_1 = plot2 = plt.contourf(real_dates, y_sal,SAL_interp,contour_levels,vmin=33,vmax=37, extend='both')#, cmap=cmhot)
#contour_line = plt.contour(real_dates, y_temp, TEMP_interp, levels=[20], colors='black', linewidths=0.5)
# draw colorbar
cbar = fig.colorbar(plot2)
cbar.ax.set_ylabel('Salinity', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.ylabel('Depth (m)', fontsize=7)
plt.plot(float_final,- mld_x,'r')
ax.text(-0.3, 1.15, '(a)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')

# #I set xticks
# nxticks=9
# xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
# xticklabels=[]
# for i in xticks:
#     xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %b %y'))
# ax.set_xticks(xticks)
# ax.set_xticklabels(xticklabels)
ax.get_xaxis().set_ticklabels([])
plt.xticks(rotation=90,fontsize=7)
# I add the grid
#plt.grid(color='k', linestyle='dashed', linewidth=0.5)




# ####bbp
levels = 200
min_contour_level =0
max_contour_level = 30

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels) 


####### I plot
width, height = 0.8, 0.7
set_ylim_lower, set_ylim_upper = -200, 0

ax=fig.add_subplot(515)
ax.set_ylim([set_ylim_lower, set_ylim_upper])
ax_1 = plot2 = plt.contourf(real_dates,  y_bbp,(BBP_interp),contour_levels,vmin=0,vmax=30, extend='both')#, cmap=cmhot)


ax.set_xlim(real_dates[0], real_dates[-1])
# draw colorbar
cbar = fig.colorbar(plot2)
cbar.ax.set_ylabel('BBP POC \n (mgC m$^{-3}$)', fontsize=7)
plt.xlabel('Date', fontsize=7)
plt.ylabel('Depth (m)', fontsize=7)
cbar.ax.tick_params(labelsize=7) 
cbar.ax.locator_params(nbins=5)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.ylim(-200, 0)
ax.text(-0.3, 1.15, '(e)', transform=ax.transAxes,fontsize=7, fontweight='bold', va='top', ha='right')
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=7)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter(' %b %y')) # %d

plt.xticks(rotation=90,fontsize=7)




os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("Fig02" + ".png")
plt.savefig(fig_name_pdf,dpi=300,bbox_inches="tight")

