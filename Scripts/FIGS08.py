# Script to conduc correlation between float data and lagrangian data
#should also add the physical parameters of the float
# 1. Setup
# import needed libraries
import os
from pathlib import Path
import pandas as pd
import numpy as np
from seawater import dpth
from oceanpy import mixed_layer_depth
import gsw

# for plots
import seaborn as sns
import matplotlib.pyplot as plt


#function  determining the pvalue
def r_pvalues(df):
    from scipy.stats import spearmanr
    cols = pd.DataFrame(columns=df.columns)
    p = cols.transpose().join(cols, how='outer')
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            p[r][c] = round(spearmanr(tmp[r], tmp[c])[1], 4) #you can change between pearson or spearman 
    return p

############################################
#1- Lagrangian dataset
###########################################
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

df_float1 = pd.read_csv(Path_to_data/"final_data_an37eq.csv", delimiter=',' ) #1411

 
#remove these columns reduce from 1407 to 239 variables
#consider looking at them later
df_float1.drop(list(df_float1.filter(regex = 'percentile|std|signif|slope|pvalue|_R|Wind')), axis = 1, inplace = True)  #remove columns with the string percentile 

df_float1.rename(columns={'Profile':'profile'}, inplace=True)

############################################
#2- carbon dataset and physical parameters
###########################################
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()

#data ctd between 8 and 8 
#it should be filtered using the longitude criteria 
#1.CTD
df_ = pd.read_csv(Path_to_data / "Ecopart_diagnostics_data_605.tsv",sep='\t', low_memory=False)

print(df_.columns.tolist()) #determine the variables we want 

#I select only the profiles data, which contain 'ASC' in the filename, and I exclude the parkings
# select only ascent profiles :
df_ = df_[df_['Profile'].str.contains('a')]    


df_ =df_[['RAWfilename', 'Profile', 'Latitude', 'Longitude', 'Date_Time', 
          'Pressure [dbar]', 'Vol [L] (sampled for this depth bin)', 'Flux_mgC_m2',
          'MiP_abun', 'MaP_abun',
          'Depth [m]','Chlorophyll-a [mg/m3]']]




#take clusters
#defining the path
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/pierre").expanduser()


#open the files where I have the different particle 
# first I need to import the csv indiv_float that contains object id in order to compute the concentration of each cluster considered as species

indiv_binned = pd.read_csv(Path_to_data/"clusters_concentrations.csv")

#rename the depth 
#rename pressure to pressure bin
df_.rename(columns={'Pressure [dbar]':'depth', 'Profile':'profile'}, inplace=True)



Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

# and the sampled volume
volumes = pd.read_csv(Path_to_data/'volumes_float.csv')

volumes=volumes.rename(columns={"sample_id": "profile","depth_bin":"depth"})

# add the volume and compute the concentrations
obs = pd.merge(indiv_binned, volumes, how="left", on=['profile', 'depth'])


# Replace '0008a' with '0008a_WMO6904139_recovery' in the 'profile' column
obs['profile'] = obs['profile'].replace('0008a', '0008a_WMO6904139_recovery')


#ADD  lon lat and date  
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()


profiles = pd.read_csv(Path_to_data/"list_of_profiles.csv")

full_final=pd.merge(obs, profiles, how="left", on=['profile'])

# Drop rows where 'lat' column has NaN values
full_final = full_final.dropna(subset=['lat'])


#create column for cluster 1, 2, 3,4,5

# Pivot the DataFrame to create new columns for each cluster and their respective concentrations
pivot_df = full_final.pivot_table(index=full_final.index, columns='Cluster', values='conc', aggfunc='first')

# Rename the columns to include the prefix 'cluster_' and the cluster number
pivot_df.columns = ['cluster_' + col.split(' ')[1] for col in pivot_df.columns]

# Concatenate the pivoted DataFrame with the original DataFrame
result_df = pd.concat([full_final, pivot_df], axis=1)



#combine Cflux with cluster

#merge metadata and ecopart data based on the Rawfilename 
df_grouped = df_.merge(result_df, how="left", on=["profile","depth"])



# Remove 'lat', 'lon', and 'time' columns
columns_to_remove = ['lat', 'lon', 'date','Date_Time','watervolume','fractal','esd','kurt','mean','circ','perim','vol_sph']
df_grouped = df_grouped.drop(columns=columns_to_remove)


# I will do the mean between 0-100m 
deptha1=100
deptha2=0


df_mean=df_grouped[(df_grouped.depth>=deptha2-10) & (df_grouped.depth<=deptha1+10)]

# calculate mean based on taxon, lat, lon, profile 
# Calculate mean based on 'lat', 'lon', 'profile', and 'taxon'
mean_df = df_mean.groupby(['Latitude', 'Longitude', 'profile', 'Cluster','RAWfilename']).mean().reset_index()



#merge 
df_grouped_lagr= mean_df.merge(df_float1, how="left", on=["profile","Longitude","Latitude"])

#remove these columns
df_grouped_lagr.drop(list(df_grouped_lagr.filter(regex = 'percentile|std|signif|slope|pvalue|_R|Wind')), axis = 1, inplace = True)  #remove columns with the string percentile 


#remove the int columns  before calculating the p value
#drop string columns      
df_grouped_lagr.drop(['Date_Num', 'Cruise', 'Glob_dt_or_nrt', 
           'Glob_Ekman_dt_or_nrt','RAWfilename','profile',"Cluster"],
          inplace=True, axis=1)



df_globekman =df_grouped_lagr.filter(regex=('MiP|MaP|Flux|Eul_Chl_Oceancolour_|cluster_|Lagr_Chl_GlobEkmanDt_Oceancolour_|Chloro'))


#here I should put a loop on the variable so I can plot each variable alone

co=df_globekman.corr(method='spearman') #montonic function instead of only linear


#use custom function to calculate p-values
Pvalues=r_pvalues(df_globekman)

#Replace p values >0.05  with nan 
co1=co
co1[Pvalues>0.05]=np.nan

co2=co1[['Flux_mgC_m2','MiP_abun','MaP_abun','cluster_1','cluster_2','cluster_3','cluster_4']]





co2.drop(['Flux_mgC_m2','MiP_abun','MaP_abun','cluster_1','cluster_2','cluster_3','cluster_4','cluster_5'], inplace=True, axis=0) #drop lines I don't need

#draw subplots of correlation between Map and the chlorophyll 
list_chl=['Lag','Eul']


# # # Find the index where the word 'lag' is present
# lag_index = [idx for idx in co2.index if 'Lag' in str(idx)]

# # Extract the values of the 'Map' column corresponding to the indices in lag_index
# map_values_lag_index = co2.loc[lag_index, 'MaP_abun']


# Extracting the first part of each index name up to the first occurrence of "days"

import re

# Creating a figure and subplots
fig, axes = plt.subplots(nrows=1, ncols=len(list_chl), figsize=(10,4))

for i, chl_type in enumerate(list_chl):

    # Searching for chl_type in lag_index
    chl_indices = [idx for idx in co2.index if chl_type in idx]
    
    # Extracting the corresponding map values
    #map_values_chl_index = co2.loc[chl_indices, 'MaP_abun']
    map_values_chl_index = co2.loc[chl_indices, 'Flux_mgC_m2']
    
    
    # Extracting the corresponding days before for the x-axis
    days_before_list_chl = []
    for idx in chl_indices:
        match = re.search(r'(\d+)days', idx)
        if match:
            days_before = int(match.group(1))
            days_before_list_chl.append(days_before)
        else:
            days_before_list_chl.append(None)
    
    # Plotting the bar plot in the current subplot
    axes[i].bar(days_before_list_chl, map_values_chl_index)
    axes[i].set_xlabel('Day')
    axes[i].set_ylabel('Corr Value')
    axes[i].set_title(f'Corr value of MiP abund and {chl_type} chl')


plt.tight_layout()
stop
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/Lagrangian/" )
fig_name_png = "Mip_vs_chloro_correlation.png"
plt.savefig(fig_name_png,dpi=300,bbox_inches="tight")
plt.close()