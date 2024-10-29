import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.dates as mdates

### open the float interpolation data
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )

M158_lon=pd.read_csv("Time_float100.csv", sep=',')
M158_depth=pd.read_csv("DEPTH_float100.csv", sep=',')
M158_Cflux=pd.read_csv("FluxC_float100.csv", sep=',')
M158_MiP=pd.read_csv("Mip_float100.csv", sep=',')
M158_MaP=pd.read_csv("Map_float100.csv", sep=',')

lon158_2=np.squeeze(np.array(M158_lon))
depth158=np.squeeze(np.array(M158_depth))

lon158=mdates.num2date(lon158_2)

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/" )

###open the mask csv
mask=pd.read_csv("mask_final.csv")

#open file of mld 
mld=pd.read_csv("mld.csv", sep=',')
mld_time=pd.read_csv("mld_time.csv", sep=',')

## create the matrix for time and depth in case I use them 
time_repeat=np.tile(lon158_2,(depth158.size,1))
depth_repeat=np.transpose([depth158]*lon158_2.size)

# Define the bin size (20m)
bin_size = 30

# Round the depth values to the nearest multiple of bin_size
rounded_depth_repeat = np.round(depth158 / bin_size) * bin_size

def calculate_flux_and_plot(dataframe, label, color, depth_unique):
    # Event 1
    FLUX_event1 = dataframe.copy(deep=True)
    FLUX_event1[mask != 1] = np.nan
    FLUX_event1['Rounded_Depth'] = rounded_depth_repeat
    mean_flux_by_depth1 = FLUX_event1.groupby('Rounded_Depth').mean()
    mean_flux1 = np.nanmean(mean_flux_by_depth1, axis=1)
    std_flux1 = np.nanstd(mean_flux_by_depth1, axis=1)
    yerr1_0 = mean_flux1 - std_flux1
    yerr1_1 = mean_flux1 + std_flux1

    # Event 2
    FLUX_event2 = dataframe.copy(deep=True)
    FLUX_event2[mask != 2] = np.nan
    FLUX_event2['Rounded_Depth'] = rounded_depth_repeat
    condition = (FLUX_event2['Rounded_Depth'] <= -750) & (FLUX_event2['Rounded_Depth'] >= -800)
    for column in range(100):  # Assuming columns are named from '0' to '99'
        FLUX_event2.loc[condition & (FLUX_event2[str(column)] > 50), str(column)] = np.nan
    mean_flux_by_depth2 = FLUX_event2.groupby('Rounded_Depth').mean()
    mean_flux2 = np.nanmean(mean_flux_by_depth2, axis=1)
    std_flux2 = np.nanstd(mean_flux_by_depth2, axis=1)
    yerr2_0 = mean_flux2 - std_flux2
    yerr2_1 = mean_flux2 + std_flux2

    # Between events
    FLUX_between_event = dataframe.copy(deep=True)
    FLUX_between_event[mask != 0] = np.nan
    FLUX_between_event['Rounded_Depth'] = rounded_depth_repeat
    mean_flux_by_depth3 = FLUX_between_event.groupby('Rounded_Depth').mean()
    mean_flux_between = np.nanmean(mean_flux_by_depth3, axis=1)
    std_flux_between = np.nanstd(mean_flux_by_depth3, axis=1)
    yerr3_0 = mean_flux_between - std_flux_between
    yerr3_1 = mean_flux_between + std_flux_between

    # Plotting
    plt.plot(mean_flux1, depth_unique, 'b', label=f'{label} ev1')
    plt.fill_betweenx(depth_unique, yerr1_1, yerr1_0, color='b', alpha=0.1)
    plt.plot(mean_flux2, depth_unique, 'r', label=f'{label} ev2')
    plt.fill_betweenx(depth_unique, yerr2_1, yerr2_0, color='r', alpha=0.1)
    plt.plot(mean_flux_between, depth_unique, 'g', label=f'{label} between')
    plt.fill_betweenx(depth_unique, yerr3_1, yerr3_0, color='g', alpha=0.1)

depth_unique = np.unique(rounded_depth_repeat)

# Plotting
fig = plt.figure(1, figsize=(6, 4))
ax1 = fig.add_subplot(131)



# Plot for Cflux
ax1.set_title('Cflux')
calculate_flux_and_plot(M158_Cflux, 'Cflux', 'b', depth_unique)
ax1.set_ylim(-2000, 0)
ax1.set_ylabel('Depth (m)', fontsize=10)
ax1.set_xlim(0, 50)


ax2 = fig.add_subplot(132)
# Plot for MiP
ax2.set_title('MiP')
calculate_flux_and_plot(M158_MiP, 'MiP', 'r', depth_unique)
ax2.set_ylim(-2000, 0)
ax2.get_yaxis().set_ticklabels([])
ax2.set_xlim(0, 90)

ax3 = fig.add_subplot(133)
# Plot for MaP
ax3.set_title('MaP')
calculate_flux_and_plot(M158_MaP, 'MaP', 'g', depth_unique)
ax3.set_ylim(-2000, 0)
ax3.get_yaxis().set_ticklabels([])
ax3.set_xlim(0, 1)




plt.tight_layout()
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("Fig02_supplementary" + ".png")
plt.savefig(fig_name_pdf, dpi=300, bbox_inches="tight")
plt.show()
