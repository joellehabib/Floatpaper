"""
Figure 5: Contour plots of log-concentration per particle cluster over depth and time.

"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import griddata
# datetime not required directly; using pandas + mdates for conversions

def contour_levels_func(min_contour_level: float, max_contour_level: float, levels: int) -> np.ndarray:
    """Return evenly spaced contour levels between min and max."""
    step = (max_contour_level - min_contour_level) / levels
    return np.arange(min_contour_level, max_contour_level, step)


def gridding_func(pos_min_max, depth_min_max, pos_array, depth, param):
    """Grid scattered points (pos_array, depth) onto a regular grid."""
    xi = np.linspace(min(pos_min_max), max(pos_min_max), 50)
    yi = np.linspace(min(depth_min_max), max(depth_min_max), 200)
    zi = griddata((pos_array, depth), param, (xi[None, :], yi[:, None]), method="nearest")
    return xi, yi, zi

#defining the path

Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

## Load base float/time/depth data (kept paths for backwards compatibility)
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot")
M158_lon = pd.read_csv("Time_float100.csv", sep=",")
M158_depth = pd.read_csv("DEPTH_float100.csv", sep=",")

lon158 = np.squeeze(np.array(M158_lon))
depth158 = np.squeeze(np.array(M158_depth))



Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/pierre").expanduser()


"""Open the file that contains object IDs to compute concentration per cluster"""
indiv = pd.read_csv(Path_to_data / "clusters_concentrations.csv")

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
# Sampled volume per profile/depth
volumes = pd.read_csv(Path_to_data / "volumes_float.csv")

volumes=volumes.rename(columns={"sample_id": "profile","depth_bin":"depth"})

# add the volume and compute the concentrations
obs = pd.merge(indiv_binned, volumes, how="left", on=['profile', 'depth'])

#ADD  lon lat and date  
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()
profiles = pd.read_csv(Path_to_data / "list_of_profiles.csv")

full_final=pd.merge(obs, profiles, how="left", on=['profile'])

# Add subplot labels
labels = ['(a)', '(b)', '(c)', '(d)', '(e)']

parameter_dic={"cluster 1":[0, 2],
               "cluster 2":[0, 0.5], 
               "cluster 3":[0,0.2],
               "cluster 4":[0, 0.3],
               "cluster 5":[0, 0.3]}

## Interpolation based on time
n = 1
for parameter in parameter_dic:
    df_1 = full_final[full_final["Cluster"] == parameter].reset_index(drop=True)
    df_1 = df_1[df_1["date"].notna()]

    
    # I extract the data
    lon = df_1["lon"].to_numpy()
    lat = df_1["lat"].to_numpy()
    date_vals = df_1["date"].astype(int).astype(str)

    # Depth is provided positive downward; use negative for plotting with 0 at surface
    pressure = -df_1["depth"].to_numpy()

    # Log-transform concentrations (avoid log(0))
    abund = df_1["conc"].to_numpy()
    abund[abund == 0] = 1
    abund = np.log(abund)
   

 
   
    # Vectorized datetime parsing and conversion to Matplotlib date numbers
    dt = pd.to_datetime(date_vals, format="%Y%m%d", errors="coerce")
    # Convert pandas datetime Series to Matplotlib date numbers without to_pydatetime()
    Date_Num = dt.map(mdates.date2num).to_numpy()
    
    
    
    
    #interpolation method 1
    
    
    # I define the x and y arrays for the plot
    x_date = np.linspace(np.nanmin(Date_Num), np.nanmax(Date_Num), 100)
    y_pressure = np.linspace(np.nanmin(pressure), np.nanmax(pressure), 200)
    x_date_g, y_pressure_g = np.meshgrid(x_date, y_pressure)




    # I interpolate
    # I interpolate
    zi1 = griddata((Date_Num, pressure), abund, (x_date_g, y_pressure_g), method="nearest", rescale=True)
    
    #mean or median of the data to try detecting anomaly 
    #trying map
    # Keep original field (non-anomaly) for plotting
    ZI1 = pd.DataFrame(zi1)
    
    # Convert Matplotlib date numbers back to datetime for plotting
    real_dates = mdates.num2date(x_date)



    ### anomaly 1
    levels = 15
    min_contour_level = 0
    max_contour_level = 8

    contour_levels = contour_levels_func(min_contour_level, max_contour_level, levels)
    fig = plt.figure(1, figsize=(5,8))
    ax=fig.add_subplot(5,1,n)
    
    
    p1 = ax.contourf(real_dates, y_pressure, ZI1, contour_levels, cmap="viridis", alpha=1, extend="both")
    ax.set_ylim([-1000, 0])
    ax.tick_params(axis='both', labelsize=7)
    ax.set_ylabel('Depth (m)')
    cbar = fig.colorbar(p1, ax=ax)
    cbar.ax.set_ylabel('Log', fontsize=8)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  # Set format to one decimal place
    #plt.title(parameter, fontsize=8)
    

  # Add subplot labels outside the plots
    ax.annotate(labels[n-1], xy=(-0.1, 1), xycoords='axes fraction', fontsize=12, fontweight='bold', va='top', ha='right')
    
    
    if n!=5:
        ax.get_xaxis().set_ticklabels([])
    else:
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %y')) 
        plt.xticks(rotation=90,fontsize=7)
   
    

    n+=1


os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/")
fig_name_pdf = ("Fig5" + ".png")
plt.tight_layout(pad=1.0)
plt.savefig(fig_name_pdf, dpi=300, bbox_inches="tight")
print(f"Figure saved to: {os.path.abspath(fig_name_pdf)}")

# Display the figure on screen
plt.show()


