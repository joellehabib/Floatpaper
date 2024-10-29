"""
Script to plot the uvp histogram data as a contour plot. Trying out pycharm.
"""


import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'

from matplotlib.ticker import FormatStrFormatter
from pathlib import Path
import matplotlib.dates as dates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings # these two lines to remove the annoying warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)
import itertools
import matplotlib.cm as cm
from scipy.interpolate import griddata
from seawater import dpth
import calendar
import time
from datetime import date,datetime
from scipy.signal import savgol_filter
import seaborn as sns
import cmocean
import matplotlib.dates as mdates
from datetime import timedelta

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

Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

#open mask 
###open the mask csv
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )

M158_lon=pd.read_csv("Time_float100.csv", sep=',')
M158_depth=pd.read_csv("DEPTH_float100.csv", sep=',')
M158_Cflux=pd.read_csv("FluxC_float100.csv", sep=',')
M158_Mip=pd.read_csv("Mip_float100.csv", sep=',')
M158_Map=pd.read_csv("Map_float100.csv", sep=',')

lon158=np.squeeze(np.array(M158_lon))
depth158=np.squeeze(np.array(M158_depth))

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/" )

mask=pd.read_csv("mask100.csv")
mask[mask==0]=np.nan
## create the matrix for time and depth in case I use them 
time_repeat=np.tile(lon158,(depth158.size,1))
depth_repeat=np.transpose([depth158]*lon158.size)



Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/pierre").expanduser()


#open the files where I have the different particle 
# first I need to import the csv indiv_float that contains object id in order to compute the concentration of each cluster considered as species

indiv = pd.read_csv(Path_to_data/"clusters_concentrations.csv")

# Compute abundance and biovolume per taxa and bin for either the rough or medium taxonomic definiton
indiv_binned = indiv.groupby(['Cluster', 'depth', 'profile'], group_keys=False).agg(n = pd.NamedAgg(column = 'Cluster', aggfunc = 'count'),
                                                                                         vol_sph = pd.NamedAgg(column = "vol_sph", aggfunc = 'sum'),
                                                                                         perim = pd.NamedAgg(column = "perim", aggfunc = 'mean'),
                                                                                         circ = pd.NamedAgg(column = "circ", aggfunc = 'mean'),
                                                                                         mean = pd.NamedAgg(column = "mean", aggfunc = 'mean'),
                                                                                         kurt = pd.NamedAgg(column = "kurt", aggfunc = 'mean'),
                                                                                         esd = pd.NamedAgg(column = "esd", aggfunc = 'mean'),
                                                                                         fractal = pd.NamedAgg(column = "fractal", aggfunc = 'mean'),
                                                                                         conc=pd.NamedAgg(column = "conc", aggfunc = 'mean'),
                                                                                         vol=pd.NamedAgg(column = "watervolume", aggfunc = 'mean'))

indiv_binned.reset_index(inplace = True) # to keep a column with exports_groups and depth_bin values

Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

# and the sampled volume
volumes = pd.read_csv(Path_to_data/'volumes_float.csv')

volumes=volumes.rename(columns={"sample_id": "profile","depth_bin":"depth"})

# add the volume and compute the concentrations
obs = pd.merge(indiv_binned, volumes, how="left", on=['profile', 'depth'])

#ADD  lon lat and date  
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()


profiles = pd.read_csv(Path_to_data/"list_of_profiles.csv")

full_final=pd.merge(obs, profiles, how="left", on=['profile'])

# Add subplot labels
labels = ['(a)', '(b)', '(c)', '(d)', '(e)']

parameter_dic={"cluster 1":[0, 2],
               "cluster 2":[0, 0.5], 
               "cluster 3":[0,0.2],
               "cluster 4":[0, 0.3],
               "cluster 5":[0, 0.3]}

flux1=[];flux2=[];flux3=[];stdflux1=[];flux2=[];flux1_std=[];flux2_std=[];flux3_std=[];

# full_final = full_final.dropna()
#Interpolation based on time 
n=1
for parameter in parameter_dic:
    df_1=full_final[(full_final["Cluster"]==parameter)]
    df_1 = df_1.reset_index()
    df_1.drop(['index'], inplace=True, axis=1)
    df_1 = df_1[df_1['date'].notna()]

    
    # I extract the data
    lon=np.array(df_1['lon']);
    lat=np.array(df_1['lat'])
    Date_Time=np.array(df_1['date'])
    
    #Date_Time.replace('.','')


    pressure=-np.array(df_1['depth'])
    #abund=np.log(np.array(df_1['conc']))
   
    abund=(np.array(df_1['conc']))
    
    abund[abund==0] = 1
    abund=np.log(abund)
   

 
   
    Date_Num=np.r_[0:lon.size]
    Date_Num2=np.r_[0:lon.size]
    for i in Date_Num:
       date_time_obj = datetime.strptime(str(int(Date_Time[i])), '%Y%m%d')
       Date_Num[i]= dates.date2num(date_time_obj)
       #Date_Num[i] = calendar.timegm(date_time_obj.timetuple())
    
    
    
    
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
    
    
   
    x_date=np.squeeze(np.array(x_date))
    x_date=x_date.astype(float)
    
    # Convert numeric dates to datetime objects
    reference_date = datetime(1970, 1, 1)  # Choose your reference date
    real_dates = [reference_date + timedelta(days=date) for date in x_date]



    ### anomaly 1
    levels = 15
    min_contour_level = 0
    max_contour_level = 8

    contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)
    fig = plt.figure(1, figsize=(5,8))
    ax=fig.add_subplot(5,1,n)
    
    
    jet = plt.get_cmap('BrBG')
    p1 = plt.contourf(real_dates, y_pressure ,ZI1, contour_levels, cmap='viridis', alpha=1,extend='both')
    # ax_2 = plot2 = plt.contourf(x_date, y_pressure, mask,levels=5, cmap=cm,alpha=0.3)
    ax.set_ylim([-1000, 0])
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.ylabel('Depth (m)')
    cbar = fig.colorbar(p1)
    cbar.ax.set_ylabel('Log', fontsize=8)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  # Set format to one decimal place
    #plt.title(parameter, fontsize=8)
    

  # Add subplot labels outside the plots
    ax.annotate(labels[n-1], xy=(-0.1, 1), xycoords='axes fraction', fontsize=12, fontweight='bold', va='top', ha='right')
    
    
    nxticks=10
    xticks=np.linspace(Date_Num.min(),Date_Num.max(),nxticks)
    xticklabels=[]
    
    
    if n!=5:
        ax.get_xaxis().set_ticklabels([])
    else:
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %y')) 
        plt.xticks(rotation=90,fontsize=7)
   
    

    n+=1


os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("Fig6" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")
    
plt.close()


