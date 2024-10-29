# trichodemsimum section
#interpolate the zooplankton 
import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'


from pathlib import Path
import pandas as pd
import numpy as np
import cmocean
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from seawater import dpth
import calendar
import time
from datetime import date,datetime
from scipy.signal import savgol_filter
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import itertools

from tabulate import tabulate

def contour_levels_func(min_contour_level, max_contour_level, levels):
    """Function to define contour levels for contourf"""
    distance_levels = max_contour_level / levels
    contour_levels = np.arange(min_contour_level, max_contour_level, distance_levels)
    return contour_levels


def gridding_func(pos_min_max, depth_min_max, pos_array, depth, param):
    grid_method = "linear"  # choose the gridding method here
    # method to do the regridding, can be nearest (for nearest neighbour) or linear

    xi = np.linspace(min(pos_min_max), max(pos_min_max), 50)
    yi = np.linspace(min(depth_min_max), max(depth_min_max), 200)
    zi = griddata((pos_array, depth), param, (xi[None, :], yi[:, None]), method=grid_method)
    return xi, yi, zi

#defining the path
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()

def rounditup(x, precision, method):
    import math
    if method == "floor":
        return math.floor(x / precision) * precision
    elif method == "ceiling":
        return math.ceil(x / precision) * precision
    else:
        return "give the parameter floor or ceiling"
#########################
#1- ecopart data
#########################
###open the mask csv
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/" )
mask=pd.read_csv("mask.csv")


M158_lon=pd.read_csv("Time_float50.csv", sep=',')
M158_depth=pd.read_csv("DEPTH_float50.csv", sep=',')
lon158_2=np.squeeze(np.array(M158_lon))
depth158=np.squeeze(np.array(M158_depth))


## create the matrix for time and depth in case I use them 
time_repeat=np.tile(lon158_2,(depth158.size,1))
depth_repeat=np.transpose([depth158]*lon158_2.size)


### ecopart
df_ecopart1 = pd.read_csv(Path_to_data / "full_equator_rough.csv", sep=',')

parameter_dic={"Cyanobacteria"}

df_ecopart1['depth_bin'] = df_ecopart1["depth"].apply(lambda x: rounditup(x, 30, "floor"))

# fig, axes = plt.subplots(1, 1, figsize=(9, 8), sharex=True)

df_filtered = df_ecopart1[df_ecopart1["taxon"] == "Cyanobacteria"]
 
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

 # and the sampled volume
volumes = pd.read_csv(Path_to_data/'volumes_float.csv')

volumes=volumes.rename(columns={"sample_id": "profile","depth_bin":"depth"})

 # add the volume and compute the concentrations
obs = pd.merge(df_filtered, volumes, how="left", on=['profile', 'depth'])


#ADD  lon lat and date  
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()


profiles = pd.read_csv(Path_to_data/"list_of_profiles.csv")

full_final=pd.merge(obs, profiles, how="left", on=['profile'])

# I extract the data
lon=np.array(full_final['lon']);
lat=np.array(full_final['lat'])
Date_Time=np.array(full_final['date'])

#Date_Time.replace('.','')


pressure=-np.array(full_final['depth'])
#abund=np.log(np.array(df_1['conc']))

abund=(np.array(full_final['conc']))

abund[abund==0] = 1
abund=np.log(abund)


Date_Num=np.r_[0:lon.size]
for i in Date_Num:
   date_time_obj = datetime.strptime(str(int(Date_Time[i])), '%Y%m%d')
   Date_Num[i] = calendar.timegm(date_time_obj.timetuple())




#interpolation method 1


# I define the x and y arrays for the plot
x_date = np.linspace(Date_Num.min(),Date_Num.max(),100)
y_pressure = np.linspace(pressure.min(),pressure.max(),200)
x_date_g,y_pressure_g=np.meshgrid(x_date,y_pressure)




# I interpolate
# I interpolate
zi1 = griddata((Date_Num,pressure), abund, (x_date_g, y_pressure_g), method="nearest",rescale=True)


#mean or median of the data to try detecting anomaly 
#trying map
mean_zi=np.mean(zi1,1) #200 levels of depth 



#value-mean to conduct anomaly 
ZI1= pd.DataFrame(zi1) 
New_value_PART=ZI1.sub(mean_zi,axis=0)



levels = 20
min_contour_level = 0
max_contour_level = 8

contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)


set_ylim_lower, set_ylim_upper =  -1000, pressure.max()
fig = plt.figure(1, figsize=(5,3))
ax=fig.add_subplot(1,1,1)

xlim=(Date_Num.min(), Date_Num.max())
plt.title("Cyanobacteria")

ax_1 = plot2 = plt.contourf(x_date, y_pressure, zi1,contour_levels, cmap='viridis',extend='both')
cbar = fig.colorbar(plot2)
cbar.ax.set_ylabel('Log', fontsize=7)
#ax_2 = plot2 = plt.contourf(x_date, y_pressure, mask,levels=5, cmap=cm,alpha=0.3)
# draw colorbar

ax.set_ylim([-400, 0])
#I set xticks
nxticks=10
xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
xticklabels=[]
for i in xticks:
        xticklabels.append(datetime.utcfromtimestamp(i).strftime('%d %B'))
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)
plt.xticks(rotation=90,fontsize=7)
ax.set_ylabel('Depth')



os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024//")
fig_name_pdf = ("Supp_trichodesmium" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")

plt.close()
