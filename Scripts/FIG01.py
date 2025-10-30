import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import xarray as xr
import cmocean
import matplotlib.dates as mdates

cm = cmocean.cm.balance

def contour_levels_func(min_contour_level, max_contour_level, levels):
    """Function to define contour levels for contourf"""
    return np.linspace(min_contour_level, max_contour_level, levels)


def matlab_to_datetime(matlab_times):
    """Vectorized conversion of MATLAB datenum to Python datetime"""
    # MATLAB datenum: days since January 0, 0000
    # Python datetime: days since January 1, 1970 (Unix epoch)
    # The offset is 719529 days (difference between MATLAB and Unix epoch)
    days_since_epoch = matlab_times - 719529
    return pd.to_datetime(days_since_epoch, unit='D', origin='unix')
           
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/")

# Load trajectory files
df_158 = pd.read_csv("traj_M158.csv", sep=',')
df_181 = pd.read_csv("traj_M181.csv", sep=',')
df_float22 = pd.read_csv("traj_float2022.csv", sep=',')

# Vectorized datetime conversion (much faster than loops)
date_158 = matlab_to_datetime(df_158.Time.values)
date_181 = matlab_to_datetime(df_181.Time.values)
date_float22 = matlab_to_datetime(df_float22.Time.values)

# Load SST data
dataset = xr.open_dataset('OI_SST_TROPATL.nc')

LAT_SST = dataset.variables['LAT'].values
LON_SST = dataset.variables['LON'].values
SST_TROPATL = dataset.variables['SST_TROPATL'].values
SST_TROPATL_FILTERED = dataset.variables['SST_TROPATL_FILTERED'].values
TIME_SST = dataset.variables['TIME'].values

# Vectorized computation of spatial means (much faster than loops)
# Filter latitude range 28:35 and compute mean along latitude axis
# Check array dimensions and adjust indexing accordingly
if SST_TROPATL.ndim == 3:
    mean_sst_tropatl = np.nanmean(SST_TROPATL[:, 28:35, :], axis=1)
    mean_sst_tropatl_filtered = np.nanmean(SST_TROPATL_FILTERED[:, 28:35, :], axis=1)
elif SST_TROPATL.ndim == 2:
    mean_sst_tropatl = SST_TROPATL
    mean_sst_tropatl_filtered = SST_TROPATL_FILTERED
else:
    raise ValueError(f"Unexpected SST array dimensions: {SST_TROPATL.shape}")

# Vectorized datetime conversion
date_sst = matlab_to_datetime(TIME_SST)

# Define export event dates
date_list1 = [datetime(2021, 8, 8), datetime(2021, 9, 8)]
date_list2 = [datetime(2021, 12, 13), datetime(2022, 1, 26)]

# Load CHLA data
dataset = xr.open_dataset('CHLA_TROPATL.nc')

LAT_CHL = dataset.variables['LAT'].values
LON_CHL = dataset.variables['LON'].values
CHL_TROPATL = dataset.variables['CHLA_TROPATL'].values
CHL_TROPATL_FILTERED = dataset.variables['CHLA_TROPATL_FILTERED'].values
TIME_CHL = dataset.variables['TIME'].values

# Vectorized computation of spatial means (latitude range 102:139 for 1N to 1S)
# Check array dimensions and adjust indexing accordingly
if CHL_TROPATL.ndim == 3:
    mean_chl_tropatl = np.nanmean(CHL_TROPATL[:, 102:139, :], axis=1)
    mean_chl_tropatl_filtered = np.nanmean(CHL_TROPATL_FILTERED[:, 102:139, :], axis=1)
elif CHL_TROPATL.ndim == 2:
    mean_chl_tropatl = CHL_TROPATL
    mean_chl_tropatl_filtered = CHL_TROPATL_FILTERED
else:
    raise ValueError(f"Unexpected CHL array dimensions: {CHL_TROPATL.shape}")

# Vectorized datetime conversion
date_chl = matlab_to_datetime(TIME_CHL)



# Helper function to reduce plotting code repetition
def setup_subplot(ax, hide_xticklabels=True):
    """Configure common subplot settings"""
    ax.set_ylim([-40, 10])
    ax.set_xlim([datetime(2021, 7, 1), datetime(2022, 4, 1)])
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%y'))
    ax.tick_params(axis='both', labelsize=7)
    ax.set_ylabel('Longitude (°)', fontsize=7)
    if hide_xticklabels:
        ax.set_xticklabels([])
    return ax


# Create figure with all subplots
fig, axes = plt.subplots(4, 1, figsize=(3, 5), sharex=True)

# (a) SST absolute values
contour_levels = contour_levels_func(22, 30.5, 100)
p1 = axes[0].contourf(date_sst, LON_SST, mean_sst_tropatl.T, 
                      contour_levels, cmap='Spectral_r', vmin=22, vmax=30.5, extend='both')
axes[0].plot(date_float22, df_float22.Lon, color='black', linewidth=0.8)
setup_subplot(axes[0])
axes[0].text(-0.2, 1.15, '(a)', transform=axes[0].transAxes, 
            fontsize=7, fontweight='bold', va='top', ha='right')
cbar = plt.colorbar(p1, ax=axes[0])
cbar.ax.set_ylabel('SST (°C)', fontsize=7)
cbar.ax.tick_params(labelsize=7)
cbar.ax.locator_params(nbins=5)

# (b) SST anomaly
contour_levels = contour_levels_func(-0.38, 0.38, 20)
p2 = axes[1].contourf(date_sst, LON_SST, mean_sst_tropatl_filtered.T,
                      contour_levels, cmap=cm, vmin=-0.38, vmax=0.38, extend='both')
axes[1].plot(date_float22, df_float22.Lon, color='black', linewidth=0.8)
setup_subplot(axes[1])
axes[1].text(-0.2, 1.15, '(b)', transform=axes[1].transAxes,
            fontsize=7, fontweight='bold', va='top', ha='right')
cbar = plt.colorbar(p2, ax=axes[1])
cbar.ax.set_ylabel('SST anomaly (°C)', fontsize=7)
cbar.ax.tick_params(labelsize=7)
cbar.ax.locator_params(nbins=5)

# (c) Chl-a absolute values
contour_levels = contour_levels_func(0, 1, 20)
p3 = axes[2].contourf(date_chl, LON_CHL, mean_chl_tropatl.T,
                      contour_levels, cmap=cmocean.cm.algae, vmin=0, vmax=1, extend='both')
axes[2].plot(date_float22, df_float22.Lon, color='black', linewidth=0.8)
setup_subplot(axes[2])
axes[2].text(-0.2, 1.15, '(c)', transform=axes[2].transAxes,
            fontsize=7, fontweight='bold', va='top', ha='right')
cbar = plt.colorbar(p3, ax=axes[2])
cbar.ax.set_ylabel('Chl-a \n (mg m$^{-3}$)', fontsize=7)
cbar.ax.tick_params(labelsize=7)
cbar.ax.locator_params(nbins=5)

# Add vertical lines for export events
y_min, y_max = axes[2].get_ylim()
ymin_fraction = (-28 - y_min) / (y_max - y_min)
ymax_fraction = (0 - y_min) / (y_max - y_min)
for date in date_list1:
    axes[2].axvline(x=date, color='r', linestyle='--', linewidth=1, 
                   ymin=ymin_fraction, ymax=ymax_fraction)
for date in date_list2:
    axes[2].axvline(x=date, color='b', linestyle='--', linewidth=1,
                   ymin=ymin_fraction, ymax=ymax_fraction)

# (d) Chl-a anomaly
contour_levels = contour_levels_func(-0.1, 0.1, 20)
p4 = axes[3].contourf(date_chl, LON_CHL, mean_chl_tropatl_filtered.T,
                      contour_levels, cmap='BrBG', vmin=-0.1, vmax=0.1, extend='both')
axes[3].plot(date_float22, df_float22.Lon, color='black', linewidth=0.8)
setup_subplot(axes[3], hide_xticklabels=False)
axes[3].text(-0.2, 1.15, '(d)', transform=axes[3].transAxes,
            fontsize=7, fontweight='bold', va='top', ha='right')
cbar = plt.colorbar(p4, ax=axes[3])
cbar.ax.set_ylabel('Chl-a anomaly \n (mg m$^{-3}$)', fontsize=7)
cbar.ax.tick_params(labelsize=7)
cbar.ax.locator_params(nbins=5)

plt.tight_layout()

# Save figure
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float//new_version_052024/")
fig_name_pdf = "Fig01.png"
# plt.savefig(fig_name_pdf, dpi=300, bbox_inches="tight")