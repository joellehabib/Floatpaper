#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 19:09:24 2024

@author: joellehabib
"""

import os
import glob
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'
from pyhdf.SD import SD
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
import calendar
import time
from datetime import date,datetime
from datetime import timezone

import gsw
from oceanpy import mixed_layer_depth

from scipy.signal import savgol_filter

#function to calculate  the  day of the year
def day_of_year(date):
    # Calculate the day of the year
    start_of_year = datetime(date.year, 1, 1, tzinfo=timezone.utc)
    day_number = (date - start_of_year).days + 1
    return day_number

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


df_globekman =df_float1.filter(regex=('Eul_Chl_Oceancolour_000|Date_Num|Longitude|Latitude'))

# Remove duplicate rows based on 'lon' and 'lat' columns
df_globekman = df_globekman.drop_duplicates(subset=['Longitude','Latitude'])
df_globekman = df_globekman.reset_index(drop=True)

TIME_lag=df_globekman.Date_Num
CHL_lag=df_globekman.Eul_Chl_Oceancolour_000daysBackward_mean


date_chl2=[]
for k in range(TIME_lag.size):
    python_date = datetime.fromordinal(int(TIME_lag[k])) + timedelta(days=TIME_lag[1]%1) - timedelta(days = 366)
    date_chl2.append(python_date)
    
    
#extract temperature 

os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/DATA_PETER/" )

    
dataset = xr.open_dataset('OI_SST_TROPATL.nc')

#I extract lon/lat
LAT_SST = dataset.variables['LAT'].values
LON_CHL = dataset.variables['LON'].values
SST_TROPATL = np.squeeze(dataset.variables['SST_TROPATL'].values)
SST_TROPATL_FILTERED=np.squeeze(dataset.variables['SST_TROPATL_FILTERED'].values)

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


# Function to convert MATLAB datenum to Python datetime
def matlab_to_datetime(matlab_datenum):
    days = matlab_datenum - 366
    return datetime.fromordinal(int(days)) + timedelta(days=days%1)

# Convert TIME_SST from MATLAB datenum to datetime
TIME_SST_datetime = [matlab_to_datetime(t) for t in TIME_SST]


chl_sat=[]
chl_anomaly=[]
sst_climato=[]
sst_climato_anomaly=[]
anomaly_chl_time=[]
anomaly_chl_points=[]


for k in range(0,float_date.size):
   float_lon1=float_lon[k]; float_lat1=float_lat[k]; float_date2=matlab_time[k];
   float_date1 = matlab_to_datetime(matlab_time[k])  # Convert MATLAB time to datetime
   
   # Find the index of the closest latitude in lat_sat snd long and time
   index_lat = np.argmin(np.abs(LAT_SST - float_lat1))
   index_lon = np.argmin(np.abs(LON_CHL - float_lon1))
   
   # Find the index where TIME_SST is equal to float_date1
   index_time = np.where( TIME_SST== float_date2)[0]
   
   # Find the indices where TIME_SST has the same day and month as float_date1
   float_day = float_date1.day
   float_month = float_date1.month
   indices_time_climatology = [i for i, t in enumerate(TIME_SST_datetime) if t.month == float_month and t.day == float_day and t.year < float_date1.year]

   
   
    
    #extract chl_troptl value 

   chl_tropatl=np.squeeze(SST_TROPATL[index_time,index_lat,index_lon])
   chl_tropatl_anomaly=np.squeeze(SST_TROPATL_FILTERED[index_time,index_lat,index_lon])
   chl_sat.append(chl_tropatl)
   
   chl_anomaly.append(chl_tropatl_anomaly)
   
   anomaly_chl_time1=np.squeeze(TIME_SST[indices_time_climatology])
   anomaly_chl_time.append(anomaly_chl_time1)
   
   SST_climato=np.squeeze(SST_TROPATL[indices_time_climatology,index_lat,index_lon])
   
   SST_climato_mean=np.mean(SST_climato)
   sst_climato.append(SST_climato_mean)
   
   SST_climato_ano=np.squeeze(SST_TROPATL_FILTERED[indices_time_climatology,index_lat,index_lon])
   anomaly_chl_points.append(SST_climato_ano)
   SST_climato_mean_ano=np.mean(SST_climato_ano)
   sst_climato_anomaly.append(SST_climato_mean_ano)
   
   

# Convert to a single NumPy array and then flatten it
flattened_anomaly_chl_time = np.concatenate(anomaly_chl_time)
flattened_anomaly_chl_points = np.concatenate(anomaly_chl_points)


time_final=[]
for k in range(0,flattened_anomaly_chl_time.size):
     time1 = matlab_to_datetime(flattened_anomaly_chl_time[k])  
     time_final.append(time1)

# Function to generate date ranges from 2011 to 2022
def generate_date_ranges(start_year, end_year):
    ranges = []
    for year in range(start_year, end_year + 1):
        if year % 2 == 1:  # Odd years
            ranges.append((datetime(year, 7, 1), datetime(year + 1, 3, 31)))
    return ranges

# Generate date ranges from 2011 to 2022
date_ranges = generate_date_ranges(2011, 2020)

# Create a dictionary to hold the data for each year range
data = {f'Year{i+1}': [] for i in range(len(date_ranges))}

# Function to check if a date falls within a given range
def falls_in_range(date, start, end):
    return start <= date <= end

# Iterate over the list of datetimes and categorize them
for dt in time_final:
    for i, (start, end) in enumerate(date_ranges):
        if falls_in_range(dt, start, end):
            data[f'Year{i+1}'].append(dt)
            break  # Once the date is added to the correct range, move to the next date

# Create a DataFrame from the dictionary, filling missing values with NaN
df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))

# Display the DataFrame
print(df)




# Function to generate date ranges from 2011 to 2020
def generate_date_ranges(start_year, end_year):
    ranges = []
    for year in range(start_year, end_year + 1):
        if year % 2 == 1:  # Odd years
            ranges.append((datetime(year, 7, 1), datetime(year + 1, 3, 31)))
    return ranges

# Generate date ranges from 2011 to 2020
date_ranges = generate_date_ranges(2011, 2020)

# Create a dictionary to hold the data for each year range
data = {f'Year{i+1}': [] for i in range(len(date_ranges))}
indices = {f'Year{i+1}': [] for i in range(len(date_ranges))}

# Function to check if a date falls within a given range
def falls_in_range(date, start, end):
    return start <= date <= end

# Iterate over the list of datetimes and categorize them
for idx, dt in enumerate(time_final):
    for i, (start, end) in enumerate(date_ranges):
        if falls_in_range(dt, start, end):
            data[f'Year{i+1}'].append(dt)
            indices[f'Year{i+1}'].append(idx)
            break  # Once the date is added to the correct range, move to the next date

# Create DataFrames from the dictionaries, filling missing values with NaN
date_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))
index_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in indices.items()]))


# Extract the values from flattened_anomaly_chl_points using the indices in index_df
chl_data = {f'Year{i+1}': [flattened_anomaly_chl_points[idx] for idx in index_df[f'Year{i+1}'].dropna().astype(int)] for i in range(len(date_ranges))}

# Create a DataFrame for the chlorophyll data
chl_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in chl_data.items()]))



# Replace the year in date_df based on the conditions
def replace_year(date):
    if pd.isna(date):
        return date
    if date.month >= 7:  # July to December
        return datetime(2021, date.month, date.day)
    else:  # January to March
        return datetime(2022, date.month, date.day)

# Apply the replace_year function to each column in date_df
for column in date_df.columns:
    date_df[column] = date_df[column].apply(replace_year)



   
dataset = xr.open_dataset('CHLA_TROPATL.nc')

#I extract lon/lat
LAT_SST = dataset.variables['LAT'].values
LON_CHL = dataset.variables['LON'].values
CHL_TROPATL = np.squeeze(dataset.variables['CHLA_TROPATL'].values)
CHL_TROPATL_FILTERED=np.squeeze(dataset.variables['CHLA_TROPATL_FILTERED'].values)
TIME_SST=dataset.variables['TIME'].values

chloro_anomaly=[]
chl_climato=[]
chl_climato_ano=[]

chl_bandpassfilter_fullyears=[]

anomaly_chl_time123=[]
anomaly_chl_points12=[]

for k in range(0,float_date.size):
   float_lon1=float_lon[k]; float_lat1=float_lat[k]; float_date2=matlab_time[k];
   float_date1 = matlab_to_datetime(matlab_time[k])  # Convert MATLAB time to datetime
   
   
   # Find the index of the closest latitude in lat_sat snd long and time
   index_lat = np.argmin(np.abs(LAT_SST - float_lat1))
   index_lon = np.argmin(np.abs(LON_CHL - float_lon1))
   
   # Find the index where TIME_SST is equal to float_date1
   index_time = np.where( TIME_SST== float_date2)[0]
   
   # Find the indices where TIME_SST has the same day and month as float_date1
   float_day = float_date1.day
   float_month = float_date1.month
   indices_time_climatology = [i for i, t in enumerate(TIME_SST_datetime) if t.month == float_month and t.day == float_day and t.year < float_date1.year]

    #extract chl_troptl value 
   chloro_tropatl_anomaly=np.squeeze(CHL_TROPATL_FILTERED[index_time,index_lat,index_lon])
   chloro_anomaly.append(chloro_tropatl_anomaly)
   
   
   anomaly_chl_time12=np.squeeze(TIME_SST[indices_time_climatology])
   anomaly_chl_time123.append(anomaly_chl_time12)
   
   CHL_climato=np.squeeze(CHL_TROPATL[indices_time_climatology,index_lat,index_lon])
   CHL_climato_mean=np.mean(CHL_climato)
   chl_climato.append(CHL_climato_mean)
   
   CHL_climato_ano=np.squeeze(CHL_TROPATL_FILTERED[indices_time_climatology,index_lat,index_lon])
   anomaly_chl_points12.append(CHL_climato_ano)
   CHL_climato_mean_ano=np.mean(CHL_climato_ano)
   chl_climato_ano.append(CHL_climato_mean_ano)

# Convert to a single NumPy array and then flatten it
flattened_anomaly_chl_time2 = np.concatenate(anomaly_chl_time123)
flattened_anomaly_chl_points2 = np.concatenate(anomaly_chl_points12)


time_final=[]
for k in range(0,flattened_anomaly_chl_time.size):
     time1 = matlab_to_datetime(flattened_anomaly_chl_time2[k])  
     time_final.append(time1)


# Iterate over the list of datetimes and categorize them
for dt in time_final:
    for i, (start, end) in enumerate(date_ranges):
        if falls_in_range(dt, start, end):
            data[f'Year{i+1}'].append(dt)
            break  # Once the date is added to the correct range, move to the next date

# Create a DataFrame from the dictionary, filling missing values with NaN
df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))


# Create a dictionary to hold the data for each year range
data = {f'Year{i+1}': [] for i in range(len(date_ranges))}
indices = {f'Year{i+1}': [] for i in range(len(date_ranges))}


# Iterate over the list of datetimes and categorize them
for idx, dt in enumerate(time_final):
    for i, (start, end) in enumerate(date_ranges):
        if falls_in_range(dt, start, end):
            data[f'Year{i+1}'].append(dt)
            indices[f'Year{i+1}'].append(idx)
            break  # Once the date is added to the correct range, move to the next date

# Create DataFrames from the dictionaries, filling missing values with NaN
date_df1 = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))
index_df1 = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in indices.items()]))


# Extract the values from flattened_anomaly_chl_points using the indices in index_df
chl_data = {f'Year{i+1}': [flattened_anomaly_chl_points2[idx] for idx in index_df[f'Year{i+1}'].dropna().astype(int)] for i in range(len(date_ranges))}

# Create a DataFrame for the chlorophyll data
chl_df1 = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in chl_data.items()]))




# Apply the replace_year function to each column in date_df
for column in date_df.columns:
    date_df1[column] = date_df[column].apply(replace_year)




## plots 

fig = plt.figure(1, figsize=(6, 10))  # Adjust figure size as needed
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)  # Adjust spacing between subplots

ax=fig.add_subplot(511)
p1=plt.plot(date_chl2,chl_sat,'b')
p2=plt.plot(date_chl2,sst_climato,'r',label='climato')
ax.legend(fontsize=7)
ax.legend(fontsize=7)
plt.title("Surface temperature")
plt.ylabel('°C')
ax.text(-0.1, 1.15, '(a)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter(' %b %y')) # %d
ax.get_xaxis().set_ticklabels([])

ax=fig.add_subplot(512)
p1=plt.plot(date_chl2,CHL_lag,'b')
p2=plt.plot(date_chl2,chl_climato,'r')
ax.legend(fontsize=7)
plt.title("Surface chl")
plt.ylabel('mg $m^{-3}$')
ax.text(-0.1, 1.15, '(b)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
ax.get_xaxis().set_ticklabels([])
plt.xticks(rotation=90,fontsize=8)

ax=fig.add_subplot(513)
p1=plt.plot(date_chl2,chl_anomaly,'b')
plt.plot(date_df,chl_df,color='red', alpha=0.2)
#p2=plt.plot(date_chl2,sst_climato_anomaly,'r')
ax.legend(fontsize=7)
plt.title("bandpass filtered temperature")

plt.ylabel('°C')
ax.text(-0.1, 1.15, '(c)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
ax.get_xaxis().set_ticklabels([])

ax=fig.add_subplot(514)
p1=plt.plot(date_chl2,chloro_anomaly,'b')
#plt.scatter(updated_time_final2,flattened_anomaly_chl_points2,color='red', alpha=0.1)
plt.plot(date_df1,chl_df1,color='red', alpha=0.2)
plt.title("bandpass filtered chl")
plt.ylabel('mg $m^{-3}$')
ax.text(-0.1, 1.15, '(d)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
ax.get_xaxis().set_ticklabels([])
ax.legend(fontsize=7)




ax=fig.add_subplot(515)
p1=plt.plot(lon158,corresponding_value21)

ax.legend(fontsize=7)
plt.title("NPP")
plt.ylabel('mg $m^{-2}$ $day^{-1}$')
ax.get_xaxis().set_ticklabels([])
ax.text(-0.1, 1.15, '(e)', transform=ax.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')
plt.xticks(rotation=90,fontsize=7)
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter(' %b %y')) # %d



 
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("S2_anomaly1" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")








