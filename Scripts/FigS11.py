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
mask=pd.read_csv("mask_final.csv")

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )
M158_lon=pd.read_csv("Time_float100.csv", sep=',')
M158_depth=pd.read_csv("DEPTH_float100.csv", sep=',')
lon158_2=np.squeeze(np.array(M158_lon))
depth158=np.squeeze(np.array(M158_depth))


## create the matrix for time and depth in case I use them 
time_repeat=np.tile(lon158_2,(depth158.size,1))
depth_repeat=np.transpose([depth158]*lon158_2.size)


### ecopart
df_ecopart1 = pd.read_csv(Path_to_data / "full_equator_rough.csv", sep=',')



parameter_dic={"Other_Rhizaria",
               "Copepoda",
               "Phaeodaria",
               "Collodaria"}



df_ecopart1['depth_bin'] = df_ecopart1["depth"].apply(lambda x: rounditup(x, 30, "floor"))

# Filter the dataframe based on the values in parameter_dic
df_filtered = df_ecopart1[df_ecopart1["taxon"].isin(parameter_dic)]




Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

# and the sampled volume
volumes = pd.read_csv(Path_to_data/'volumes_float.csv')

volumes=volumes.rename(columns={"sample_id": "profile","depth_bin":"depth","volume_L":"vol"})

# add the volume and compute the concentrations
obs = pd.merge(df_filtered, volumes, how="left", on=['profile', 'depth'])

#ADD  lon lat and date  
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()


profiles = pd.read_csv(Path_to_data/"list_of_profiles.csv")

full_final=pd.merge(obs, profiles, how="left", on=['profile'])



#open the mask I created for the different events 

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/" )

###open the mask csv
mask=pd.read_csv("mask_final.csv")



#I will create bins of 50m to stack the data 

# bin the data per 50m
def rounditup(x, precision, method):
    import math
    if method == "floor":
        return math.floor(x / precision) * precision
    elif method == "ceiling":
        return math.ceil(x / precision) * precision
    else:
        return "give the parameter floor or ceiling"


full_final.loc[:, 'depth_part'] = 99999
full_final['depth_part'] = full_final["depth"].apply(lambda x: rounditup(x, 50, "floor") + 0)

#calculate the number of paticles bcs my table has concentrations
full_final['n_part']=full_final['conc']*full_final['vol']


# Now compute number of particles per profile, size class, taxo group and depth bin
# I modified this, I summed the 
df = full_final.groupby(['profile', 'depth_part','taxon']) \
    .agg({'esd': 'mean', 'vol_sph': 'sum', 'n_part': 'sum', 'vol':'sum'}) \
    .rename(columns={'esd': 'esd_mean', 'vol_sph': 'vol_sph_tot', 'n_part': 'n','vol':'watervolume'}) \
    .reset_index()


df["conc"] = df["n"] / df["watervolume"]  # concentration in n/m3
df["vol_sph"] = df["vol_sph_tot"] / df["watervolume"]  # biovolume concentration in mm3/m3


# Add the date
list_prof_500 = pd.read_csv(Path_to_data/"list_of_profiles.csv")

df = df.merge(list_prof_500, on=['profile'])

ref_date = datetime(2021, 7, 1) #beginning of my timeseries
df.head()


width, height = 0.8, 0.7
fig = plt.figure(1, figsize=(10, 5))

flux_0=np.zeros((41,4))
flux_1=np.zeros((41,4))
flux_2=np.zeros((41,4))

std_0=np.zeros((41,4))
std_1=np.zeros((41,4))
std_2=np.zeros((41,4))
  
for mas in range(0,3):
    parameter_dic={"Other_Rhizaria":[0, 2],
                   "Copepoda":[0, 0.5], 
                   "Phaeodaria":[0,0.2]}
    
    title_dic=("outside-between mask",
              "event 1",
              "event 2")

    flux1=[];flux2=[];flux3=[];stdflux1=[];flux2=[];flux1_std=[];flux2_std=[];flux3_std=[];

    # full_final = full_final.dropna()
    #Interpolation based on time 
    n=0 
    
    ax=fig.add_subplot(1,3,mas+1) 
    
    
    for parameter in parameter_dic:
        df_1=df[(df["taxon"]==parameter)]
        df_1 = df_1.reset_index()
        df_1.drop(['index'], inplace=True, axis=1)
        df_1 = df_1[df_1['date'].notna()]
        
        
        # I extract the data
        lon=np.array(df_1['lon']);
        lat=np.array(df_1['lat'])
        Date_Time=np.array(df_1['date'])
        
        #Date_Time.replace('.','')


        pressure=-np.array(df_1['depth_part'])
        #abund=np.log(np.array(df_1['conc']))
       
        abund=(np.array(df_1['conc']))
        
        abund[abund==0] = 1
        abund=np.log(abund)
        #abund[np.isnan(abund)] = 0
       
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
        zi2 = griddata((Date_Num,pressure), abund, (x_date_g, y_pressure_g), method="nearest",rescale=True)
        
         #round the depth every 50 m is one depth 
        y_pressure01 = pd.DataFrame(y_pressure, columns =["depth"])
        y_pressure01['depth_equiv']= y_pressure01["depth"].apply(lambda x: rounditup(x, 50, "floor") + 0)
        
        #I will take the different masks of 50,200
        zi1=zi2
        
        
        zi1[mask!=mas]=np.nan
        mean_flux1=np.nanmean(zi1,1)
        std_flux1=np.nanstd(zi1,1)
        
        
        # plt.scatter(mean_flux1,y_pressure)
        # plt.ylim(-600,0)
        
       
        #combine flux and pressure 
        y_pressure01['flux1']=mean_flux1
        y_pressure01['depth_equiv']= y_pressure01['depth_equiv']
        
        def std(x): return np.std(x)
        
        df01 = y_pressure01.groupby(['depth_equiv']) \
            .agg({'flux1': ['mean',std]}) \
            .rename(columns={'mean': 'mean_flux', 'std': 'std_flux'}) \
           .set_axis(['mean_flux', 'std_flux'], axis=1)\
           .reset_index()
        
        if mas==0:
            flux_0[:,n]=df01.mean_flux
            std_0[:,n]=df01.std_flux
            n+=1
        if mas==1:
            flux_1[:,n]=df01.mean_flux
            std_1[:,n]=df01.std_flux
            n+=1
        if mas==2:
            flux_2[:,n]=df01.mean_flux
            std_2[:,n]=df01.std_flux
            n+=1
        
        
       
        plt.plot(df01.mean_flux,df01.depth_equiv, label=parameter)
        plt.fill_betweenx(df01.depth_equiv,df01.mean_flux+df01.std_flux,df01.mean_flux-df01.std_flux, alpha=0.1)
        plt.ylim(-1000,0)
        plt.axhline(y=-351, color='r', linestyle='--')
        plt.xlim(0, 8)
        #plt.legend()
        plt.title(title_dic[mas])
        ax.set_xlabel('concentration (ind.m$^{-3}$)')
        
        
        if mas==2:
            plt.legend()
        else:
            print('no legend')
            

        
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("supp_S11" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")

plt.close()
            