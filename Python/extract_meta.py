
# coding: utf-8

# ### Extract and Clean feeding variables
# 
# 

# In[7]:


# import modules
import pandas as pd
import numpy as np
import os
import re
# path of file 
directory = "data/excel_files/"
all_files = os.listdir(directory) 
# exlude non-xls files 
files = [file for file in all_files if file[-3:] == "xls"]


# In[8]:


# Extract  ---------------------------

# create df in long format but only if data is present
# thus, exclude if sheet2 no present
sheet_not_found = []
df_long = pd.DataFrame()
for file in files:
  id = file[0:3]
  temp_xl = pd.ExcelFile(f"{directory}/{file}")
  if "Sheet2" in temp_xl.sheet_names:
    temp_df = temp_xl.parse("Sheet2")
    temp_df.insert(loc=0, column='subject', value=id)
    df_long = df_long.append(temp_df, sort = True)
  else:
    sheet_not_found.append(id)

# how many files worked
print(df_long.subject.unique().size)
# which did not work? 
print(sheet_not_found)
# Note: I double checked those excel files and indeed the data is not there


# In[9]:


# drop some unused columns and change colnames
df_long = df_long.drop(columns = [
  "Unnamed: 10", 
  "merk en naam poedermelk",
  "opmerking",
  "opmerkingen",
  "welke ziekte",
  "ziek"
  ])
colnames = [
    "expressed_bf", 
    "bf", 
    "pet", 
    "formula", 
    "type_solid_food", 
    "subject", 
    "solid_food", 
    "week"
]
df_long.columns = colnames


# In[10]:


# define cleaner functions ---------------------------

# to learn about the trouble-values I catch them if they are not convertable to floats:
def catch_exceptions(values):
  problem_values = []
  for value in values:
    value_str = str(value)
    try:
       float(value_str)
    except ValueError:
      problem_values.append(value)
  return(problem_values)

# which values are not convertable (only one example is printed):
for col in ["bf", "expressed_bf", "solid_food", "formula"]:
    print(catch_exceptions(df_long.loc[:, col]))


# In[11]:


# disregard the rest non-information 
def clean_value(value):
    value = str(value)
    value = value.replace(",", ".")
    if value == "0-1":
        value = "0"
    # we can use the information sometimes for the x tot/a y pattern
    # for that I use the average. e.g. 5 tot 6 will become 5.5
    if " tot " in value or " a " in value:
        pattern = re.compile(r'(\d+)(\s\w+\s)(\d+)')
        match_groups = pattern.match(value).groups()
        first_n = float(match_groups[0])
        second_n = float(match_groups[2])
        # if the values are more than 2 counts apart, then I put NA
        if abs(first_n - second_n) <= 2:
          value = str(np.mean([first_n, second_n]))
        else:
          value = str(min(first_n, second_n))
    # filter some leftover non-information
    if value in [" ", "leeg", "kwijt", ".", "?", "nan", "`"] or len(value) > 4:
        value = np.nan
    return(value)

# to later compare if everything worked I clone the df 
df_long2 = df_long.copy()


# use cleaner ---------------------------

df_long.bf = df_long.bf.apply(clean_value)
df_long.expressed_bf = df_long.expressed_bf.apply(clean_value)
df_long.solid_food = df_long.solid_food.apply(clean_value)
df_long.formula = df_long.formula.apply(clean_value)
df_long.pet = df_long.pet.apply(clean_value)


# In[12]:


# check if we can now convert to numericals 
for col in ["bf", "expressed_bf", "solid_food", "formula"]:
    print(catch_exceptions(df_long.loc[:, col]))
# now everything is convertible to numerics where necessary


# In[13]:


# export back to csv 
df_long.to_csv("data/meta_variables/metadata_long.csv")

