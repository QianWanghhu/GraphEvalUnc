## Import Packages
import os #Navigating OW model and OW binary filepaths
import io #File inputs/outputs
import numpy as np 
import pandas as pd 
import random #Generating calibration seeds
import datetime

## Import Openwater features
import openwater.nodes as node_types #Nodes 
import openwater.discovery #Package that discovers the base algorithm binaries for OpenWater (written in Go). 
from openwater.discovery import discover 
from openwater.template import ModelFile # One of two classes for working with models
#ModelFile is used when working with existing models that are saved to disk
from openwater.template import run_simulation #Runs simulation by writing onto harddrive (SSD preferred)
from openwater import config #Deals with parameterising the OW model 
from openwater.results import OpenwaterResults #Imports results (reads .h5 files, breaksdown result by dimensions/variables, and summarises these results as tables)
from dsed import migrate #SedNet implementation of Source in OW

## Calibration Algorithm
import spotpy
from spotpy.parameter import Uniform #Method of specifying how to draw from parameter space. 
from spotpy.objectivefunctions import rmse #Root Mean Square Function, not used (we use NSE)
 
#Tell OpenWater where the base algorithm binaries are located on your PC
OpenWater_engine_binaries = r'D:\BaiduSyncdisk\Projects\OW Model Packaged\openwater-binaries-20240117\bin' 

#A list of available models should appear
openwater.discovery.OW_BIN = OpenWater_engine_binaries
discover()
model_fn_original = os.path.join(r'D:\BaiduSyncdisk\Projects\OW Model Packaged\OW Models Trimmed\MW', \
                                 'MW_RC_13_WQT_126001A' +'.h5')
model_original = ModelFile(model_fn_original) #Model object

## Part 1. Load the model
#Range of model
START='1986/07/01' #Here we specify the whole date range of the model (even though we are calibrating to only 2009-2018), OW doesn't like it if we only run for certain years
END='2023/06/30'
TIME_PERIOD = pd.date_range(START,END)
YEARS=len(TIME_PERIOD)/365.25

# MODEL= 'MW_RC_13_WQT_126001A_copy' #The name of the model that we use in python #Qian: don't think it is used in the code for now, so comment it out.
MODEL_NAME= 'MW_RC_13_WQT_126001A_copy' #The name of the actual .h5 file. Can be the same as the model name we call in python

model_directory = r'D:\BaiduSyncdisk\Projects\OW Model Packaged\OW Models Trimmed\MW' #Location of the model
time_string = datetime.datetime.now().strftime("on %B %d %Y") #Gets the current time and date and converts into str, we use this to record model runs [CURRENTLY NOT USED]

model_dir = os.path.join(model_directory) #The folder where the file is located
model_fn = os.path.join(model_dir,MODEL_NAME+'.h5') #The path for the model itself
new_output_fn = os.path.join(model_directory, MODEL_NAME + '_TEMP'+ '.h5') #Where we wan to store the output
##======================Run model with default parameter values===========================#
# Note that the default values might have been changed due to calibration and the oldest model verson used is 2024/06/25.
model = ModelFile(model_fn)
new_results = model.run(TIME_PERIOD, new_output_fn, overwrite = True)
new_results = OpenwaterResults(model_fn,new_output_fn,TIME_PERIOD)
SC = 'SC #104'
din =  new_results.time_series('InstreamFineSediment','upstreamMass', 'catchment')[SC] #Query the results
din = din*60*60*24/1000 #OW uses SI units for consistency, here we convert kg/s into kg/day # QIAN: need value to be tones
flow = new_results.time_series('StorageRouting', 'inflow', 'catchment')[SC] # Unit for flow is m^3/s.
# annual_sum = din.resample('AS-JUL').sum() #Resample based on WATER YEAR
annual_sum = din['2009-07-01':'2018-07-01 '] # years of interest
flow_filter = flow['2009-07-01':'2018-07-01 ']
flow_load = pd.concat([flow_filter, annual_sum], axis = 1)
flow_load.columns = ['Q (m3/s)', 'Load (t)']
flow_load.to_csv('../output/TN_TP_Sed/sediment.csv')