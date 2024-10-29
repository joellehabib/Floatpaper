"""
Script to plot the uvp histogram data as a contour plot. Trying out pycharm.
"""


import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'


from pathlib import Path

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
import itertools
from datetime import date,datetime
from scipy.signal import savgol_filter
import seaborn as sns
import cmocean


cm=cmocean.cm.balance


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
cm2=cmocean.cm.balance

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

volumes.loc[:, 'depth_part'] = 99999
volumes['depth_part'] = volumes["depth"].apply(lambda x: rounditup(x, 50, "floor") + 0)



full_final.loc[:, 'depth_part'] = 99999
full_final['depth_part'] = full_final["depth"].apply(lambda x: rounditup(x, 50, "floor") + 0)

#calculate the number of paticles bcs my table has concentrations
full_final['n_part']=full_final['conc']*full_final['vol']


# Now compute number of particles per profile, size class, taxo group and depth bin
# I modified this, I summed the 
df = full_final.groupby(['profile', 'depth_part','Cluster']) \
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
fig = plt.figure(1, figsize=(11, 5))

flux_0=np.zeros((41,5))
flux_1=np.zeros((41,5))
flux_2=np.zeros((41,5))

std_0=np.zeros((41,5))
std_1=np.zeros((41,5))
std_2=np.zeros((41,5))
  
for mas in range(0,3):
    parameter_dic={"cluster 1":[0, 2],
                   "cluster 2":[0, 0.5], 
                   "cluster 3":[0,0.2],
                   "cluster 4":[0, 0.3],
                   "cluster 5":[0, 0.3]}
    
    parameter_leg=["BDP",
                   "FP", 
                   "BPP",
                   "SDP",
                   "SPP"]
    
    title_dic=("outside-between mask",
              "event 1",
              "event 2")

    flux1=[];flux2=[];flux3=[];stdflux1=[];flux2=[];flux1_std=[];flux2_std=[];flux3_std=[];

    # full_final = full_final.dropna()
    #Interpolation based on time 
    n=0 
    
    ax=fig.add_subplot(1,3,mas+1) 
    
    
    for idx, parameter in enumerate(parameter_dic):
        df_1=df[(df["Cluster"]==parameter)]
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
        
        # abund[abund==0] = 1
        # abund=np.log(abund)
        abund[np.isnan(abund)] = 0
       
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
        
        
       
        plt.plot(df01.mean_flux,df01.depth_equiv, label=parameter_leg[idx])
        plt.fill_betweenx(df01.depth_equiv,df01.mean_flux+df01.std_flux,df01.mean_flux-df01.std_flux, alpha=0.1)
        plt.ylim(-1000,0)
        if mas==0:
            plt.xlim(0, 900)
        else:
            plt.xlim(0, 900)
        #plt.legend()
        plt.title(title_dic[mas])
        ax.set_xlabel('concentration (#.m$^{-3}$)')
        if mas==2:
            plt.legend()
        else:
            print('no legend')


stop

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("Fig_profile_clusters_100mbin_1000m" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")

plt.close()

stop
#calculate the transfer efficiency or the attenuation flux b 
# Initialize an empty list to store the results
result_list = []
z0=100

from sklearn.linear_model import LinearRegression
#38  and  20

# Filter values between 100 and 1000
filtered_flux = flux_0[(df01.depth_equiv <= -100) & (df01.depth_equiv >= -1000)]
z=df01.depth_equiv[(df01.depth_equiv <= -100) & (df01.depth_equiv >= -1000)]

# Loop through each column
for i in range(filtered_flux.shape[1]):
    # Create a linear regression model
    model = LinearRegression()

    # Fit the model
    x=np.log(z.values/-z0).reshape(-1, 1)
    y=np.log(filtered_flux[:,i])
    model.fit(x,y)

    #add regression line 

    # Print the slope and intercept
    print(f'Slope (m): {model.coef_[0]}')

    r_sq = model.score(x, y)
    print(f"coefficient of determination: {r_sq}")


    y_pred = np.exp((model.predict(x)))
    print(f"predicted response:\n{y_pred}")

    plt.plot(y_pred,z.values,color='k')
    plt.plot(filtered_flux[:,i], z.values)
    # Append the results to the list
    result_list.append(model.coef_[0])

   


#calculate the %
#percentage by depth 
sum_0=np.sum(flux_0,1)


perc_flux0=np.zeros((41,5))
perc_flux1=np.zeros((41,5))
perc_flux2=np.zeros((41,5))



for n in range(0,5):
    perc_flux0[:,n]=flux_0[:,n]*100/sum_0
    perc_flux1[:,n]=flux_1[:,n]*100/sum_0
    perc_flux2[:,n]=flux_2[:,n]*100/sum_0
    

# Create a line plot for each cluster
# Create subplots for line plots
fig_line, axs_line = plt.subplots(3, 1, figsize=(12, 18), sharex=True)

# Plot line for each cluster in perc_flux0
for n in range(5):
    axs_line[0].plot(perc_flux0[:, n],df01.depth_equiv, label=f'Cluster {n + 1}', alpha=0.7)

axs_line[0].set_ylabel('Percentage')
axs_line[0].set_title('Line Plot of Percentage abundance for Clusters (Event 0)')
axs_line[0].legend(loc='upper left', bbox_to_anchor=(1, 1))
axs_line[0].set_ylim(0, -1000)  # Set x-axis limits
# Plot line for each cluster in perc_flux1
for n in range(5):
    axs_line[1].plot(perc_flux1[:, n],df01.depth_equiv, label=f'Cluster {n + 1}', alpha=0.7)

axs_line[1].set_ylabel('Percentage')
axs_line[1].set_ylim(0, -1000)  # Set x-axis limits

# Plot line for each cluster in perc_flux2
for n in range(5):
    axs_line[2].plot(perc_flux2[:, n],df01.depth_equiv, label=f'Cluster {n + 1}', alpha=0.7)

axs_line[2].set_xlabel('Depth Equivalent (m)')
axs_line[2].set_ylabel('Percentage')
axs_line[2].set_title('Line Plot of Percentage abundance for Clusters (Event 2)')
axs_line[2].set_ylim(0, -1000)  # Set x-axis limits



os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("percentage_abundance" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")

plt.close()






sum_0=np.sum(flux_0,1)


perc_flux0=np.zeros((41,5))
perc_flux1=np.zeros((41,5))
perc_flux2=np.zeros((41,5))



for n in range(0,5):
    perc_flux0[:,n]=flux_0[:,n]*100/sum_0
    perc_flux1[:,n]=flux_1[:,n]*100/sum_0
    perc_flux2[:,n]=flux_2[:,n]*100/sum_0
    


fig = plt.figure(1, figsize=(11, 5))
ax=fig.add_subplot(131) 
for n in range(5):
    ax.plot(perc_flux0[:, n],df01.depth_equiv, alpha=0.7)


ax.set_title('outside-between mask')
ax.set_ylim(-1000,0)  # Set x-axis limits
ax.set_xlabel('Percentage')  # Set x-axis limits
ax.set_xlim(0, 70) 

ax=fig.add_subplot(132) 
for n in range(5):
    ax.plot(perc_flux1[:, n],df01.depth_equiv, alpha=0.7)


ax.set_title('event 1')
ax.set_ylim(-1000,0)  # Set x-axis limits
ax.set_xlim(0, 70) 
ax.set_xlabel('Percentage')  # Set x-axis limits

ax=fig.add_subplot(133) 
for n in range(5):
    ax.plot(perc_flux2[:, n],df01.depth_equiv, alpha=0.7)


ax.set_title('event 2')
ax.set_ylim(-1000,0)  # Set x-axis limits
ax.set_xlim(0, 70) 
ax.set_xlabel('Percentage')  # Set x-axis limits



os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("percentage_abundance" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")

plt.close()
    
    