import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from datetime import timedelta
import xarray as xr
from matplotlib.dates import DateFormatter
import matplotlib.dates as dates
from scipy.interpolate import griddata
import matplotlib.dates as mdates
from pathlib import Path
from seawater import dpth
from datetime import datetime,date,timezone
import calendar
import time
from datetime import date,datetime
from pyhdf.SD import SD
import gsw
from oceanpy import mixed_layer_depth

from scipy.signal import savgol_filter

#function to calculate  the  day of the year
def day_of_year(date):
    # Calculate the day of the year
    start_of_year = datetime(date.year, 1, 1, tzinfo=timezone.utc)
    day_number = (date - start_of_year).days + 1
    return day_number


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

### open the float interpolation data
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )

lon=pd.read_csv("Time_float100.csv", sep=',')
depth=pd.read_csv("DEPTH_float100.csv", sep=',')
Cflux=pd.read_csv("FluxC_float100.csv", sep=',')
Mip=pd.read_csv("Mip_float100.csv", sep=',')
Map=pd.read_csv("Map_float100.csv", sep=',')

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/Article/" )
###open the mask csv
mask=pd.read_csv("mask_final.csv")

# Get the last row of the DataFrame
last_row = mask.iloc[-1]

# Replace each column's values with the value from the last row
for column in mask.columns:
    mask[column] = last_row[column]


os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )
longitude=pd.read_csv("Lon_float100.csv", sep=',')
latitude=0

Cflux=np.squeeze(np.array(Cflux))
lon158_01=np.squeeze(np.array(lon))
depth158=np.squeeze(np.array(depth))

lon158=mdates.num2date(lon158_01)
#storing array
day_number01=np.zeros((lon.size,1))

# Convert all datetime objects to UTC by adding timezone info
datetime_list_utc = [dt.replace(tzinfo=timezone.utc) for dt in lon158]


# Check if the list contains only datetime objects
if all(isinstance(dt, datetime) for dt in datetime_list_utc):
    # Convert all datetime objects to UTC by adding timezone info
    datetime_list_utc01 = [dt.replace(tzinfo=timezone.utc) for dt in datetime_list_utc]
else:
    print("The list contains non-datetime objects.")



# Function to find the nearest day of the year
def find_nearest_day_of_year(target_day, day_list):
    return min(day_list, key=lambda x: abs(x - target_day))


###the doy list  of NPP

monthly_or_8d = ""



plot_map = 0

if plot_map != 1:
    plot_timeseries = 1
else:
    plot_timeseries = 0
    


doy_list = ["001", "009", "017", "025", "033", "041", "049", "057",
            "065", "073", "081", "089", "097", "105", "113", "121", "129", "137",
            "145", "153", "161", "169", "177", "185", "193", "201", "209", "217",
            "225", "233", "241", "249", "257", "265", "273",
             "281", "289", "297", "305", "313", "321", "329", "337", "345", "353", "361"]



corresponding_value21=np.zeros((lon.size,1))
near_doy=np.zeros((lon.size,1))

n1=0
n=0

# Iterate through each datetime in the list
for dt in datetime_list_utc01:
    # Calculate the day of the year for the current datetime
    target_day = day_of_year(dt)
    
    
    # Find the nearest day of the year from the predefined list
    nearest_doy = find_nearest_day_of_year(target_day, [int(doy) for doy in doy_list])
    near_doy[n1,0]=nearest_doy
    n1+=1
    if target_day>180: # for  2021
        
        os.chdir("/Users/joellehabib/GIT/joelle/NPP/2021/" )
        # Construct the filename for the corresponding day # I added next iter so the output would be
        #different than list
        filename=next(iter(sorted(glob.glob("eppley.*" + str(nearest_doy) + ".hdf"))))
    
        
        # Check if the file exists
        if os.path.isfile(filename):
            # Open the file and do further processing
            data_from_file = SD(filename)
            d = data_from_file.select("npp")
            d = d.get()
            
            # figuring out dimensions and spacing of the dataset

            dimensions = data_from_file.datasets()['npp'][1]

            lat_dim = dimensions[0]
            lon_dim = dimensions[1]
            
              
            spacing = lon_dim/360.0
            
            lat_range_timeseries=[-24, -7]
            lon_range_timeseries=[-1, 1]

            lat_range = lat_range_timeseries
            lon_range = lon_range_timeseries

            lat_point_range = range((90 - lat_range[1])*int(spacing), (90 - lat_range[0])*int(spacing))
            lon_point_range = range((lon_range[0]+180)*int(spacing), (lon_range[1]+180)*int(spacing))
            
            lats = [(90 - i/spacing) - 1/(spacing*2) for i in lat_point_range]
            lons = [(j/spacing - 180) - 1/(spacing*2) for j in lon_point_range]
            data_ = d[(90 - lat_range[1])*int(spacing):(90 - lat_range[0])*int(spacing), (lon_range[0]+180)*int(spacing):(lon_range[1]+180)*int(spacing)].T
            data_ =  np.where(data_==-9999, np.nan, data_)
            
            # I transform longitude into an array, I should put an  index toc hange for each file
            target_lon = np.squeeze(np.array(longitude.iloc[n]))

            # Find the index in lons that is closest to the target_lon
            index_lon = np.argmin(np.abs(np.array(lons) - target_lon))

            # Now, you can access the corresponding data value in data_
            corresponding_value21[n,0] = data_[index_lon,5] #I chose the closest 0.08        
            
            #store values of 2021 in an  array
            n+=1
            
        
        else:
                print(f"File not found for day of year {nearest_doy:03d}")
    
    else:
        
        os.chdir("/Users/joellehabib/GIT/joelle/NPP/2022/" )
        print("Doing 2022 now")
        # Construct the filename for the corresponding day
        filename=next(iter(sorted(glob.glob("eppley.*" + str(nearest_doy) + ".hdf"))))
   
        
        # Check if the file exists
        if os.path.isfile(filename):
            # Open the file and do further processing
            data_from_file = SD(filename)
            d = data_from_file.select("npp")
            d = d.get()
            
            # figuring out dimensions and spacing of the dataset

            dimensions = data_from_file.datasets()['npp'][1]

            lat_dim = dimensions[0]
            lon_dim = dimensions[1]
            
              
            spacing = lon_dim/360.0
            
            lat_range_timeseries=[-24, -7]
            lon_range_timeseries=[-1, 1]

            lat_range = lat_range_timeseries
            lon_range = lon_range_timeseries

            lat_point_range = range((90 - lat_range[1])*int(spacing), (90 - lat_range[0])*int(spacing))
            lon_point_range = range((lon_range[0]+180)*int(spacing), (lon_range[1]+180)*int(spacing))
            
            lats = [(90 - i/spacing) - 1/(spacing*2) for i in lat_point_range]
            lons = [(j/spacing - 180) - 1/(spacing*2) for j in lon_point_range]
            data_ = d[(90 - lat_range[1])*int(spacing):(90 - lat_range[0])*int(spacing), (lon_range[0]+180)*int(spacing):(lon_range[1]+180)*int(spacing)].T
            data_ =  np.where(data_==-9999, np.nan, data_)
            
            # I transform longitude into an array, I should put an  index toc hange for each file
            target_lon = np.squeeze(np.array(longitude.iloc[n]))

            # Find the index in lons that is closest to the target_lon
            index_lon = np.argmin(np.abs(np.array(lons) - target_lon))

            # Now, you can access the corresponding data value in data_
            corresponding_value21[n,0] = data_[index_lon,5] #I chose the closest 0.08   
            
            n+=1
        else:
                print(f"File not found for day of year {nearest_doy:03d}")





#### lagrangian data
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

df_float1 = pd.read_csv(Path_to_data/"final_data_an37eq.csv", delimiter=',' ) #1411

 
#remove these columns reduce from 1407 to 239 variables
#consider looking at them later
df_float1.drop(list(df_float1.filter(regex = 'percentile|std|signif|slope|pvalue|_R|Wind')), axis = 1, inplace = True)  #remove columns with the string percentile 

df_float1.rename(columns={'Profile':'profile'}, inplace=True)


df_globekman =df_float1.filter(regex=('Lagr_Chl_GlobEkmanDt_Oceancolour_015|Date_Num|Longitude|Latitude'))

# Remove duplicate rows based on 'lon' and 'lat' columns
df_globekman = df_globekman.drop_duplicates(subset=['Longitude','Latitude'])
df_globekman = df_globekman.reset_index(drop=True)

TIME_lag=df_globekman.Date_Num
CHL_lag=df_globekman.Lagr_Chl_GlobEkmanDt_Oceancolour_015daysBackward_mean


date_chl2=[]
for k in range(TIME_lag.size):
    python_date = datetime.fromordinal(int(TIME_lag[k])) + timedelta(days=TIME_lag[1]%1) - timedelta(days = 366)
    date_chl2.append(python_date)


##### extract satellite chlorophyll data 
##################
###CHLA
###################



os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/" )

    
dataset = xr.open_dataset('CHLA_TROPATL.nc')

#I extract lon/lat
LAT_SST = dataset.variables['LAT'].values
LON_CHL = dataset.variables['LON'].values
CHL_TROPATL = np.squeeze(dataset.variables['CHLA_TROPATL'].values)
TIME_SST=dataset.variables['TIME'].values


mean_sst_tropatl_filtered=np.zeros((TIME_SST.size,LON_CHL.size))
mean_sst_tropatl=np.zeros((TIME_SST.size,LON_CHL.size))


#lat (241) lon(1202)  time(4018)

#find the dates and then find the lon and the lat 

#open the file of lon lat for float 
### ecopart
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()
df_profile = pd.read_csv(Path_to_data / "list_of_profiles.csv", sep=',')
float_lon=df_profile.lon; float_lat=df_profile.lat; float_date=df_profile.date


# Convert the Series to a datetime format
datetime_series = float_date.apply(lambda x: pd.to_datetime(str(int(x)), format='%Y%m%d'))

# Convert the entire series to pandas datetime
date_series = pd.to_datetime(datetime_series)

# Convert to Unix timestamp (seconds since 1970-01-01)
timestamp_seconds = (date_series - pd.to_datetime('1970-01-01')).dt.total_seconds()

# Convert to MATLAB time
matlab_time = timestamp_seconds / (24 * 60 * 60) + 719529


date_chl=[]
for k in range(0,TIME_SST.size):
    python_date = datetime.fromordinal(int(TIME_SST[k])) + timedelta(days=TIME_SST[1]%1) - timedelta(days = 366)
    date_chl.append(python_date)


chl_sat=[]

for k in range(0,float_date.size):
   float_lon1=float_lon[k]; float_lat1=float_lat[k]; float_date1=matlab_time[k];
   
   
   # Find the index of the closest latitude in lat_sat snd long and time
   index_lat = np.argmin(np.abs(LAT_SST - float_lat1))
   index_lon = np.argmin(np.abs(LON_CHL - float_lon1))
   
   # Find the index where TIME_SST is equal to float_date1
   index_time = np.where( TIME_SST== float_date1)[0]

    
    #extract chl_troptl value 

   chl_tropatl=np.squeeze(CHL_TROPATL[index_time,index_lat,index_lon])
   chl_sat.append(chl_tropatl)




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
y_chl = np.linspace(pressure_chl.min(),pressure_chl.max(),100)
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


#date I need to remove 2022-01-05  -- 2022-01-11date_object.timestamp()

date_time_obj = datetime.strptime('2022-01-05', '%Y-%m-%d')
date_min= dates.date2num(date_time_obj)
date_time_obj = datetime.strptime('2022-01-11', '%Y-%m-%d')
date_max = dates.date2num(date_time_obj)



#nan between this period for the chl and the bbp 

CHL_interp[(x_chl_g>date_min) & (x_chl_g<date_max) ]=np.nan
# BBP_interp[(x_chl_g>date_min) & (x_chl_g<date_max) ]=np.nan


# Convert numeric dates to datetime objects
reference_date = datetime(1970, 1, 1)  # Choose your reference date
real_dates = [reference_date + timedelta(days=date1) for date1 in x_chl]


# Calculate the mean of 'chl' after removing the specified date range
mean_chl_after_removal = np.nanmean(CHL_interp[94:100], axis=0)
mean_sal_after_removal = np.nanmean(SAL_interp[46:], axis=0)
mean_temp_after_removal = np.nanmean(TEMP_interp[46:], axis=0)


fig = plt.figure(1, figsize=(6, 10))  # Adjust figure size as needed
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)  # Adjust spacing between subplots

ax=fig.add_subplot(513)
p1=plt.plot(date_chl2,CHL_lag,'r',label='Chl lag 15 days')
p2=plt.plot(real_dates,mean_chl_after_removal,label='chl in situ')
ax.legend(fontsize=7)
plt.title("Surface chl")
plt.ylabel('mg $m^{-3}$')
ax.text(-0.1, 1.15, '(c)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')


ax.get_xaxis().set_ticklabels([])
ax=fig.add_subplot(511)
plt.plot(real_dates,mean_sal_after_removal)
plt.title("Surface salinity")
ax.get_xaxis().set_ticklabels([])
ax.text(-0.1, 1.15, '(a)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')


ax=fig.add_subplot(512)
plt.plot(real_dates,mean_temp_after_removal)
plt.title("Surface temperature")

ax.get_xaxis().set_ticklabels([])
plt.ylabel('°C')
ax.text(-0.1, 1.15, '(b)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')

##ADD the carbon flux in the surface as well as MIP and MAP
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )

M158_lon=pd.read_csv("Time_float100.csv", sep=',')
M158_depth=pd.read_csv("DEPTH_float100.csv", sep=',')
M158_Cflux=pd.read_csv("FluxC_float100.csv", sep=',')
M158_Mip=pd.read_csv("Mip_float100.csv", sep=',')
M158_Map=pd.read_csv("Map_float100.csv", sep=',')


M158_Cflux=np.array(M158_Cflux)
M158_Mip=np.array(M158_Mip)
M158_Map=np.array(M158_Map)

lon158=np.squeeze(np.array(M158_lon))
depth158=np.squeeze(np.array(M158_depth))

lon158=mdates.num2date(lon158)


#calculate mean 
mean_cflux_after_removal = np.nanmean(M158_Cflux[188:190], axis=0)
mean_Mip_after_removal = np.nanmean(M158_Mip[189:], axis=0)
mean_Map_after_removal = np.nanmean(M158_Map[189:], axis=0)


# ax=fig.add_subplot(717)
# plt.plot(lon158,mean_cflux_after_removal)
# plt.title("Carbon flux at 100 m")
# plt.xticks(rotation=90,fontsize=7)
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter(' %b %y')) # %d

# plt.ylabel('mgC $m^{-2}$ $day^{-1}$')
# ax.text(-0.1, 1.3, '(g)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
# # r = np.corrcoef(mean_cflux_after_removal,mean_chl_after_removal)

# from scipy.stats import pearsonr
# pearsonr(mean_cflux_after_removal,mean_chl_after_removal)


ax=fig.add_subplot(514)
plt.plot(lon158,mean_Mip_after_removal)
plt.title("Mip abundance")
plt.xticks(rotation=90,fontsize=7)
ax.get_xaxis().set_ticklabels([])
plt.ylabel('# $L^{-1}$')
ax.text(-0.1, 1.15, '(d)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')

ax=fig.add_subplot(515)
plt.plot(lon158,mean_Map_after_removal)
plt.title("Map abundance")
plt.xticks(rotation=90,fontsize=10)
plt.ylabel('# $L^{-1}$')
ax.text(-0.1, 1.15, '(e)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
plt.xticks(rotation=90,fontsize=7)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter(' %b %y')) # %d

# ax=fig.add_subplot(714)
# p1=plt.plot(lon158,corresponding_value21)

# ax.legend(fontsize=7)
# plt.title("NPP")
# plt.ylabel('mg $m^{-2}$ $day^{-1}$')
# ax.get_xaxis().set_ticklabels([])
# ax.text(-0.1, 1.15, '(d)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
# ax.get_xaxis().set_ticklabels([])





stop
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("S03_v2" + ".png")
plt.savefig(fig_name_pdf,dpi=300,bbox_inches="tight")

plt.close()