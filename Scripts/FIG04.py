#Calculate flux for the 2 events and between the events

import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.dates as mdates

plt.rcParams.update(plt.rcParamsDefault)
### open the float interpolation data
os.chdir("/Users/joellehabib/GIT/TRATLEQ/Data/Float_data/matrice_plot" )

M158_lon=pd.read_csv("Time_float100.csv", sep=',')
M158_depth=pd.read_csv("DEPTH_float100.csv", sep=',')
M158_Cflux=pd.read_csv("FluxC_float100.csv", sep=',')
# M158_Cflux=pd.read_csv("Mip_float.csv", sep=',')
# M158_Cflux=pd.read_csv("Map_float.csv", sep=',')



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


#profile event 1
#estimate event number 1
FLUX_event1=M158_Cflux.copy(deep=True)
FLUX_event1[mask!=1]=np.nan

# Define the bin size (20m)
bin_size = 20

# Round the depth values to the nearest multiple of bin_size
rounded_depth_repeat = np.round(depth158 / bin_size) * bin_size
# Add the rounded depth values as a new column to FLUX_event1
FLUX_event1['Rounded_Depth'] = rounded_depth_repeat


# Calculate the mean based on the depth groups
mean_flux_by_depth1 = FLUX_event1.groupby('Rounded_Depth').mean()

# If you want to print or display the mean values
print(mean_flux_by_depth1)

mean_flux1=np.nanmean(mean_flux_by_depth1,1)
std_flux1=np.nanstd(mean_flux_by_depth1,1)

yerr1_0=mean_flux1-std_flux1
yerr1_1=mean_flux1+std_flux1


#estimate event number 2
FLUX_event2=M158_Cflux.copy(deep=True)
FLUX_event2[mask!=2]=np.nan


FLUX_event2['Rounded_Depth'] = rounded_depth_repeat



# Condition to filter rows where Rounded_Depth is -750
condition = (FLUX_event2['Rounded_Depth'] <= -700) & (FLUX_event2['Rounded_Depth'] >= -800)


# Iterate through columns from '0' to '99' and replace values with NaN if flux > 50
for column in range(100):  # Assuming columns are named from '0' to '99'
    FLUX_event2.loc[condition & (FLUX_event2[str(column)] > 50), str(column)] = None  # Replace with np.nan if not None

print(FLUX_event2)

# Calculate the mean based on the depth groups
mean_flux_by_depth2 = FLUX_event2.groupby('Rounded_Depth').mean()


mean_flux2=np.nanmean(mean_flux_by_depth2,1)
std_flux2=np.nanstd(mean_flux_by_depth2,1)

yerr2_0=mean_flux2-std_flux2
yerr2_1=mean_flux2+std_flux2


#estimate between two events
#put 0 instead of nan
FLUX_between_event=M158_Cflux.copy(deep=True)
FLUX_between_event[mask!=0]=np.nan

#FLUX_between_event[(time_repeat<18850)]=np.nan
FLUX_between_event['Rounded_Depth'] = rounded_depth_repeat

# Calculate the mean based on the depth groups
mean_flux_by_depth3 = FLUX_between_event.groupby('Rounded_Depth').mean()

mean_flux_between=np.nanmean(mean_flux_by_depth3,1)
std_flux_between=np.nanstd(mean_flux_by_depth3,1)

yerr3_0=mean_flux_between-std_flux_between
yerr3_1=mean_flux_between+std_flux_between


depth_unique=np.unique(rounded_depth_repeat)



P_values=np.zeros((mean_flux_by_depth3.shape[0],1))
# Lists to store significance markers
significance_markers = []
positions = []
tukey_results= []


from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
# Initialize variables
num_rows = mean_flux_by_depth3.shape[0]
alpha = 0.07


# Initialize a Pandas DataFrame to store the Tukey pairwise results
tukey_pairwise_columns = ['Row Index', 'group1', 'group2', 'reject']
tukey_pairwise_table = pd.DataFrame(columns=tukey_pairwise_columns)


for row_index in range(num_rows):
    # Extract data for each row
    a1 = mean_flux_by_depth3.iloc[row_index].dropna().values
    a2 = mean_flux_by_depth2.iloc[row_index].dropna().values
    a3 = mean_flux_by_depth1.iloc[row_index].dropna().values
    

    # One-way ANOVA test
    #_, p_anova = f_oneway(a1[4:a2.size+5], a2, a3)
    _, p_anova = f_oneway(a1, a2, a3)
    
    P_values[row_index]=p_anova
    

    # Print results
    # print(f"Row {row_index + 1} - One-way ANOVA p-value: {p_anova}")
    

    # Check if the p-value is less than the significance level
    if p_anova < alpha:
        significance_markers.append('+')
        #print("There is a significant difference among the groups.")
    else:
        #print("There is no significant difference among the groups.")
        significance_markers.append('-')

    positions.append(row_index + 0.5)
    
    # #If significant, perform Tukey post hoc test
    # if p_anova < alpha:
    tukey_result = pairwise_tukeyhsd(
    np.concatenate([a1, a2, a3]),
    np.concatenate([np.full(a1.size, 'out-betw'), np.full(a2.size, 'Event2'), np.full(a3.size, 'Event1')]))
    tukey_results.append(tukey_result)
    #print(tukey_result)
    
    
    # Extract relevant information from the Tukey test result
    tukey_df = pd.DataFrame(tukey_result._results_table.data[1:], columns=tukey_result._results_table.data[0])
    tukey_df['Row Index'] = row_index
        
        # Filter relevant columns and add the information to the tukey_pairwise_table
    tukey_pairwise_table = pd.concat([tukey_pairwise_table, tukey_df[['Row Index', 'group1', 'group2', 'reject']].copy()], ignore_index=True)

# Map True to '+' and False to '-'
tukey_pairwise_table['reject'] = tukey_pairwise_table['reject'].map({True: '+', False: '-'})

# Display the Tukey pairwise results table
print("\nTukey Pairwise Results Table:")
print(tukey_pairwise_table)

# Grouping by 'group1' and 'group2' and aggregating the 'reject' values
grouped_tukey_table = tukey_pairwise_table.groupby(['group1', 'group2'], as_index=False)['reject'].agg(lambda x: ''.join(x))

# Display the grouped Tukey pairwise results table
print("\nGrouped Tukey Pairwise Results Table:")
print(grouped_tukey_table)

#martin curve


from sklearn.linear_model import LinearRegression
z=100
#plor loglog transformed data to check fit
d=depth_unique[depth_unique<0]
f=mean_flux2[depth_unique<0]

# Create a linear regression model
model = LinearRegression()

# Fit the model
x=np.log(-d/z).reshape(-1, 1)
y=np.log(f)
model.fit(x,y)

y_pred = np.exp((model.predict(x)))


#now I want to plot the profiles, with 

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(7, 8), sharex=False, gridspec_kw={'width_ratios': [3, 1],'height_ratios': [2, 2]})

# Plot '+' and '-' based on Tukey pairwise results


# Plot for the first subplot [-70, 0]
ax1 = axs[0,0]
## plot the profiles
p1=ax1.plot(mean_flux1,depth_unique,'b', label='ev1') 
ax1.fill_betweenx(depth_unique,yerr1_1,yerr1_0, color='C0', alpha=0.1)
ax1.set_ylim(-100, 0)
ax1.set_xlim(0, 200)  # Adjust as needed
#flux 2
p2=ax1.plot(mean_flux2,depth_unique,'r', label='ev2')
ax1.fill_betweenx(depth_unique,yerr2_1,yerr2_0, color='r', alpha=0.1)

#flux 3
p3=ax1.plot(mean_flux_between,depth_unique,'g', label='outside-between')
ax1.fill_betweenx(depth_unique,yerr3_1,yerr3_0, color='g', alpha=0.1)
ax1.set_ylabel('Depth (m)', fontsize=10)
ax1.tick_params(axis='x', labelsize=10)  # Add this line for x-ticks
ax1.tick_params(axis='y', labelsize=10)  # Add this line for y-ticks
ax1.text(-0.2, 1, '(a)', transform=ax1.transAxes, fontsize=11, fontweight='bold', va='top', ha='right')
ax1.set_xlim(0, 200)  # Set independent x-limits



# Plot for the second subplot [-200, -70]
ax2 = axs[1,0]

## plot the profiles
p1=ax2.plot(mean_flux1,depth_unique,'b', label='ev1') 
ax2.fill_between(depth_unique,yerr1_1,yerr1_0, color='C0', alpha=0.1)
ax2.set_ylim(-2000, -100)
ax2.set_xlim(0, 50)  # Adjust as needed

#flux 2
p2=ax2.plot(mean_flux2,depth_unique,'r', label='ev2')
ax2.fill_betweenx(depth_unique,yerr2_1,yerr2_0, color='r', alpha=0.1)

#flux 3
p3=ax2.plot(mean_flux_between,depth_unique,'g', label='out-betw')
ax2.fill_betweenx(depth_unique,yerr3_1,yerr3_0, color='g', alpha=0.1)

p4=ax2.plot(y_pred,d,'k', label='Martin curve')


ax2.legend(fontsize=7)
ax2.set_ylabel('Depth (m)', fontsize=10)
ax2.set_xlabel('Carbon flux (mgC m$^{-3}$)', fontsize=10)
ax2.tick_params(axis='x', labelsize=10)  # Add this line for x-ticks
ax2.tick_params(axis='y', labelsize=10)  # Add this line for y-ticks
ax2.text(-0.2, 1, '(b)', transform=ax2.transAxes,fontsize=10, fontweight='bold', va='top', ha='right')






# Create a subplot for displaying Tukey pairwise results on the right
ax3 = axs[0,1]

# Assign colors to '+' and '-' values
colors = {'+': 'black', '-': 'white'}



# Plot horizontal bars for each group pair
for index, row in grouped_tukey_table.iterrows():
    reject_values = row['reject']
    
    for i, reject in enumerate(reject_values):
        if i > 95:  # Ensure valid range
            if reject == '-':
                # Use empty circle 'o' for '-'
                ax3.scatter(index, i, marker='o', s=50, facecolors='none', edgecolors='black', linewidths=1)
            else:
                # Use '+' for other values
                ax3.scatter(index, i, marker='+', s=50, color=colors[reject], edgecolors='black', linewidths=1)
        else:
            print('no')



# Customize the plot
ax3.set_xticks(range(len(grouped_tukey_table)))
ax3.set_xticklabels([f"{row['group1']} vs {row['group2']}" for _, row in grouped_tukey_table.iterrows()])
ax3.set_title('Tukey Pairwise Results')
ax3.invert_xaxis()  # Invert y-axis to have the first group pair at the top
# Rotate x-axis tick labels and center them
ax3.set_xticklabels([], ha='center')

# Choose 5 specific y-tick positions
y_ticks_positions = [96,97,98, 99, 100, 101]

# Set corresponding y-tick labels
y_ticks_labels = [-100,-80, -60, -40, -20,0]

# Set y-tick positions and labels
ax3.set_yticks(y_ticks_positions)
ax3.set_yticklabels(y_ticks_labels)

# Set xlim to start a bit earlier
ax3.set_xlim(-0.5, len(grouped_tukey_table) - 0.5)  # Adjust the values as needed

# Adjust y-axis limits to make markers more visible
ax3.set_ylim(95.8, 101.3)  # Adjust the values as needed



# Create a subplot for displaying Tukey pairwise results on the right
ax4 = axs[1,1]

# Assign colors to '+' and '-' values
colors = {'+': 'black', '-': 'white'}




# Plot horizontal bars for each group pair
for index, row in grouped_tukey_table.iterrows():
    reject_values = row['reject']
    
    for i, reject in enumerate(reject_values):
        if i<95:
            #ax4.scatter(index, i, marker='+', s=50, color=colors[reject], edgecolors='black', linewidths=1)
            #ax4.scatter(index, i, marker=',' if reject == '-' else '+', s=50, color=colors[reject], edgecolors='black', linewidths=1)
            ax4.scatter(index, i, marker='+', s=50, color=colors[reject], edgecolors='black', linewidths=1)
            
        else:
            print('no')


# Customize the plot
# Set background color to white





ax4.set_xticks(range(len(grouped_tukey_table)))
ax4.set_xticklabels([f"{row['group1']} vs {row['group2']}" for _, row in grouped_tukey_table.iterrows()])

ax4.invert_xaxis()  # Invert y-axis to have the first group pair at the top
# Rotate x-axis tick labels and center them
ax4.set_xticklabels(ax4.get_xticklabels(), rotation=90, ha='center')

# Choose 8 specific y-tick positions
y_ticks_positions = [ 0, 12, 24, 36, 50, 63, 75, 88, 94]

# Set corresponding y-tick labels
y_ticks_labels = [-2000,-1750, -1500, -1250, -1000,-750,-500,-250,-100]

# Set y-tick positions and labels
ax4.set_yticks(y_ticks_positions)
ax4.set_yticklabels([f'{label}' for label in y_ticks_labels[:9]])  # Select only the first 8 labels

# Set xlim to start a bit earlier
ax4.set_xlim(-0.5, len(grouped_tukey_table) - 0.5)  # Adjust the values as needed

# Adjust y-axis limits to make markers more visible
ax4.set_ylim(0, 95)  # Adjust the values as needed



# Adjust layout to make the size of the third and fourth plots the same as the first two plots
plt.tight_layout(pad=1.5)






os.chdir("/Users/joellehabib/GIT/TRATLEQ/Plots/Article_float/new_version_052024/" )
fig_name_pdf = ("Fig_04_ver03" + ".png")
plt.savefig(fig_name_pdf, dpi=300,bbox_inches="tight")


plt.close()

