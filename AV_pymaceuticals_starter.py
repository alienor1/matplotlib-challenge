#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Dependencies and Setup
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Hide warning messages in notebook
import warnings
warnings.filterwarnings('ignore')

# File to Load (Remember to Change These)
mouse_drug_data_to_load = "data/mouse_drug_data.csv"
clinical_trial_data_to_load = "data/clinicaltrial_data.csv"

# Read the Mouse and Drug Data and the Clinical Trial Data
got_mouse_data = pd.read_csv(mouse_drug_data_to_load)

# Preview Mouse and Drug Data
got_mouse_data.head()


# In[2]:


# Preview Clinical Trial Data
got_clinical_data = pd.read_csv(clinical_trial_data_to_load)
got_clinical_data.head()


# In[3]:


# Combine the data into a single dataset
combined_data = pd.merge(got_clinical_data, got_mouse_data,
                         on="Mouse ID", how = "left")

# Display the data table for preview
combined_data.head()


# ## Tumor Response to Treatment

# In[4]:


# Store the Mean Tumor Volume Data Grouped by Drug and Timepoint 
tumor_group = combined_data.groupby(["Drug","Timepoint"])
tumor_group_mean = tumor_group["Tumor Volume (mm3)"].mean()

# Convert to DataFrame
tumor_group_table = pd.DataFrame(tumor_group_mean).reset_index()

# Preview DataFrame
tumor_group_table


# In[5]:


# Store the Standard Error of Tumor Volumes Grouped by Drug and Timepoint
tumor_group_sem = tumor_group["Tumor Volume (mm3)"].sem()

# Convert to DataFrame
tumor_group_sem_table = pd.DataFrame(tumor_group_sem).reset_index()

# Preview DataFrame
tumor_group_sem_table.head()


# In[6]:


# Minor Data Munging to Re-Format the Data Frames
tumor_group_mean = tumor_group_mean.reset_index()
tumor_vols_pivot_mean = tumor_group_mean.pivot(index="Timepoint", 
                                               columns="Drug")["Tumor Volume (mm3)"]

# Preview that Reformatting worked
tumor_vols_pivot_mean.head()


# In[7]:


# Generate the Plot (with Error Bars)
Capomulin_error = tumor_group_sem_table.loc[tumor_group_sem_table["Drug"] 
                                            == "Capomulin", "Tumor Volume (mm3)"]
Infubinol_error = tumor_group_sem_table.loc[tumor_group_sem_table["Drug"] 
                                            == "Infubinol", "Tumor Volume (mm3)"]
Ketapril_error = tumor_group_sem_table.loc[tumor_group_sem_table["Drug"] 
                                           == "Ketapril", "Tumor Volume (mm3)"]
Placebo_error = tumor_group_sem_table.loc[tumor_group_sem_table["Drug"] 
                                          == "Placebo", "Tumor Volume (mm3)"]

# Set time
Time = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]

plt.errorbar(Time, tumor_vols_pivot_mean["Capomulin"] , 
             yerr= Capomulin_error, label= "Capomulin", marker= "o", color="red", linestyle="--", linewidth=0.5)
plt.errorbar(Time, tumor_vols_pivot_mean["Infubinol"] , 
             yerr= Infubinol_error, label= "Infubinol", marker= "^", color="blue", linestyle="--", linewidth=0.5)
plt.errorbar(Time, tumor_vols_pivot_mean["Ketapril"] , 
             yerr= Ketapril_error, label= "Ketapril", marker= "d", color="green", linestyle="--", linewidth=0.5)
plt.errorbar(Time, tumor_vols_pivot_mean["Placebo"] , 
             yerr= Placebo_error , label= "Placebo", marker= "s", color="black", linestyle="--", linewidth=0.5)

# Add labels, legend, title and grid
plt.xlabel("Time (Days)")
plt.ylabel("Tumor Volume (mm3)")
plt.legend(loc="upper left")
plt.title("Tumor Response to Treatment")
plt.grid(axis='y')

# Save the Figure
plt.savefig("charts_images/treatment_graph.png")

# Show the Figure
plt.show()


# ## Metastatic Response to Treatment

# In[8]:


# Store the Mean Met. Site Data Grouped by Drug and Timepoint 
meta_group_data = combined_data.groupby(["Drug", "Timepoint"])
meta_group_mean = meta_group_data["Metastatic Sites"].mean()

meta_group_table = pd.DataFrame(meta_group_mean)

# Preview DataFrame
meta_group_table.head()


# In[9]:


# Store the Standard Error associated with Met. Sites Grouped by Drug and Timepoint 
meta_group_sem = meta_group_data["Metastatic Sites"].sem()
# Convert to DataFrame
meta_group_sem_table = pd.DataFrame(meta_group_sem)

# Preview DataFrame
meta_group_sem_table.head()


# In[10]:


# Minor Data Munging to Re-Format the Data Frames
meta_group_mean = meta_group_mean.reset_index()
meta_vols_pivot_mean = meta_group_mean.pivot(index="Timepoint", 
                                               columns="Drug")["Metastatic Sites"]

# Preview that Reformatting worked
meta_vols_pivot_mean.head()


# In[11]:


# Generate the Plot (with Error Bars)
plt.errorbar(Time, meta_vols_pivot_mean["Capomulin"] , 
             yerr= Capomulin_error, label= "Capomulin", marker= "o", color="red", linestyle="--", linewidth=0.5)
plt.errorbar(Time, meta_vols_pivot_mean["Infubinol"] , 
             yerr= Infubinol_error, label= "Infubinol", marker= "^", color="blue", linestyle="--", linewidth=0.5)
plt.errorbar(Time, meta_vols_pivot_mean["Ketapril"] , 
             yerr= Ketapril_error, label= "Ketapril", marker= "d", color="green", linestyle="--", linewidth=0.5)
plt.errorbar(Time, meta_vols_pivot_mean["Placebo"] , 
             yerr= Placebo_error , label= "Placebo", marker= "s", color="black", linestyle="--", linewidth=0.5)

# Add labels, legend, title and grid
plt.xlabel("Treatment Duration (Days)")
plt.ylabel("Met. Sites")
plt.axis([-2, 47.5, -0.3, 3.75])
plt.legend(loc="upper left")
plt.title("Metastatic Spread During Treatment")
plt.grid(axis='y')

# Save the Figure
plt.savefig("charts_images/metastatic_spread_graph.png")

# Show the Figure
plt.show()


# ## Survival Rates

# In[12]:


# Store the Count of Mice Grouped by Drug and Timepoint (W can pass any metric)
survival_count = combined_data.groupby(["Drug", "Timepoint"]).count()["Tumor Volume (mm3)"]
survival_count 

# Convert to DataFrame
survival_count_table = pd.DataFrame(survival_count).reset_index()
survival_count_table=survival_count_table.rename(columns = {"Tumor Volume (mm3)":"Mouse Count"})

# Preview DataFrame
survival_count_table.head()


# In[13]:


# Minor Data Munging to Re-Format the Data Frames
survival_count_table = survival_count_table.reset_index()
survival_count_table_pivot = survival_count_table.pivot(index="Timepoint", 
                                               columns="Drug")["Mouse Count"]

# Preview the Data Frame
survival_count_table_pivot.head()


# In[ ]:





# In[14]:


# Generate the Plot
plt.errorbar(Time, survival_count_table_pivot["Capomulin"] , 
             yerr= Capomulin_error, label= "Capomulin", marker= "o", color="red", linestyle="--", linewidth=0.5)
plt.errorbar(Time, survival_count_table_pivot["Infubinol"] , 
             yerr= Infubinol_error, label= "Infubinol", marker= "^", color="blue", linestyle="--", linewidth=0.5)
plt.errorbar(Time, survival_count_table_pivot["Ketapril"] , 
             yerr= Ketapril_error, label= "Ketapril", marker= "d", color="green", linestyle="--", linewidth=0.5)
plt.errorbar(Time, survival_count_table_pivot["Placebo"] , 
             yerr= Placebo_error , label= "Placebo", marker= "s", color="black", linestyle="--", linewidth=0.5)

# Add labels, legend, title and grid
plt.xlabel("Time (Days)")
plt.ylabel("Survival Rate")
plt.legend(loc="lower left")
plt.title("Survival During Treatment")
plt.grid()

# Save the Figure
plt.savefig("charts_images/survival_graph.png")

# Show the Figure
plt.show()


# ## Summary Bar Graph

# In[15]:


# Calculate the percent changes for each drug
tumor_vols_pivot_mean_data = 100*(tumor_vols_pivot_mean.iloc[-1]/tumor_vols_pivot_mean.iloc[0]-1)


tumor_vols_pivot_mean_data = 100*(tumor_vols_pivot_mean.iloc[-1]/tumor_vols_pivot_mean.iloc[0]-1)
# Display the data to confirm
tumor_vols_pivot_mean_data


# In[16]:


# Store all Relevant Percent Changes into a Tuple
tuple_percent_changes = tuple(zip(tumor_vols_pivot_mean_data.index, tumor_vols_pivot_mean_data))
tuple_percent_changes_list = list(tuple_percent_changes)


# In[17]:


# Splice the data between passing and failing drugs
passing = tumor_vols_pivot_mean_data< 0

# Orient widths. Add labels, tick marks, etc. 
drug_list = ['Capomulin','Infubinol','Ketapril','Placebo']
change_list = [(tumor_vols_pivot_mean_data[drug])for drug in drug_list]
change_plt = plt.bar(drug_list,change_list,width=-1,align='edge',color=passing.map({True:'g',False:'r'}))

plt.ylim(-30,75)
plt.ylabel('% Tumor Volume Change')
plt.title('Tumor Change over 45 Day Treatment')
plt.grid()

# Use functions to label the percentages of changes
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        if height > 0:
            label_position = 2
        else:
            label_position = -8
        plt.text(rect.get_x() + rect.get_width()/2., label_position,
                '%d' % int(height)+'%',color='white',
                ha='center', va='bottom')
                
# Call functions to implement the function calls
autolabel(change_plt)

# Save the Figure
plt.savefig("charts_images/tumor_change.png")


# In[ ]:




