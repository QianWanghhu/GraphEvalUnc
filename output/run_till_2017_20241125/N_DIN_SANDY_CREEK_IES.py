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

## NSE Objective Function, SpotPy supports custom objective functions. 
def calculate_nse(observed, simulated):
    #Calculate Nash-Sutcliffe Efficiency
    observed = np.array(observed)
    simulated = np.array(simulated) 
    mean_observed = np.mean(observed)
    numerator = np.sum((observed - simulated) ** 2)
    denominator = np.sum((observed - mean_observed) ** 2)
    nse_value = (numerator / denominator) #NSE without the 1-
    return nse_value

class spot_setup(object):
     
    def __init__(self, master_input, obj_func=None):
        self.params = [] #Initialises a list that will hold parameters
        self.parameters_df = master_input #The input we specified, master_input, that contains the parameters
        # self.define_parameters() #Calls the function below # Not needed for Qian

    def define_parameters(self): #Normally, we need to define parameters line by line. We avoid this by iterating by rows over master_input
        for index, row in self.parameters_df.iterrows():
            self.params.append(spotpy.parameter.Uniform(row['identifier'], row['min'], row['max'], optguess = row['optguess'])) #Assumes UNIFORM distribution of parameters

    def parameters(self): #Generates parameter samples 
        return spotpy.parameter.generate(self.params)

    def simulation(self,parameters_vector, model_fn, model_dataframes): #Body of the class where we call and run OW

        ## Part 2. Parameterise the Model
        # Load the model
        model = ModelFile(model_fn)
        # The main worker here is a Parameteriser, which manages a list of steps
        p = config.Parameteriser()
        # And provide it to the model
        model._parameteriser = p
        #Call parameters from SpotPy from vector into list form     
        parameters_list = list(parameters_vector)  #List form
        ##Split and transform data
        master_input_dummy = pd.DataFrame({'Identifier': master_input['identifier'], \
                                           'value': parameters_list}) #Simply contains an identifier and value column 
        #We iterate over master_input_dummy, which is the parameter set SpotPy has given us
        for _, row in master_input_dummy.iterrows():
            string = row['Identifier']
            parameter_value = row['value']
            parts = string.split('-') #split identifier using delimiter -
            #Function flags
            if parts[0] == 'SednetDissolvedNutrientGeneration':  #The first part of the identifier string (part[0]) is the model
                specific_model_df = model_dataframes[parts[0]] #Find the appropreate model
                specific_model_df.loc[(specific_model_df['cgu'] == parts[3]) & (specific_model_df['constituent'] == parts[2]), parts[1]] = parameter_value
                #We locate the parameter, FU and constituent and write the parameter value into the dataframe
                #Recall from earlier: SedNet_EMC_parameters_df['Identifier'] = model_name + '-' + SedNet_parameters_df['catchment'] + '-' + SedNet_parameters_df['cgu'] + '-' + 'dissConst_EMC'
                #parts[3] is the parameter
                #parts[2] is the CGU/FU
                #parts[1] is the sub-catchment

                # Save back into dictionary
                model_dataframes[parts[0]] = specific_model_df
            elif parts[0] == 'ApplyScalingFactor':
                specific_model_df = model_dataframes[parts[0]]
                specific_model_df.loc[(specific_model_df['constituent'] == parts[2]), parts[1]] = parameter_value #Note that the parameters in ApplyScaling Factor have a different shape, hence .loc is different
                # Save back into dictionary
                model_dataframes[parts[0]] = specific_model_df
            else:
                print(f"DataFrame for '{parts[0]}' not found.")

        #Now we have taken the parameters from SpotPy and written them into a dataframe in the form that OW can accept, we WRITE these parameter values into the model file
        action = config.ParameterTableAssignment(model_dataframes['SednetDissolvedNutrientGeneration'], 'SednetDissolvedNutrientGeneration' ,dim_columns=['constituent','catchment', 'cgu'],complete=False)
        p._parameterisers.append(action)
        action = config.ParameterTableAssignment(model_dataframes['ApplyScalingFactor'], 'ApplyScalingFactor' ,parameter = 'scale',complete=False)
        p._parameterisers.append(action)
        #Run the model (overwrite=True)
        model.write()
        new_results = model.run(TIME_PERIOD, new_output_fn, overwrite = True)

        ## Part 3. Extracting Results:
        new_results = OpenwaterResults(model_fn,new_output_fn,TIME_PERIOD)
        din =  new_results.time_series('InstreamDissolvedNutrientDecay','incomingMassUpstream', 'catchment',constituent='N_DIN')['SC #104'] #Query the results
        din = din*60*60*24/1000 #OW uses SI units for consistency, here we convert kg/s into kg/day # QIAN: need value to be tones
        annual_sum = din.resample('AS-JUL').sum() #Resample based on WATER YEAR
        annual_sum = annual_sum['2009-07-01':'2018-07-01 '] #years of interest
        annual_sum = annual_sum.to_list()
        return annual_sum

    def evaluation(self):
        annual_eval #The pandas series (maybe also list?) that we are calibrating to. 
        return annual_eval

    def objectivefunction(self, simulation, evaluation, params=None): #Objective function. 
        like = calculate_nse(evaluation, simulation) #NORMALISED RMSE wrt. mean
        return like
    
#Tell OpenWater where the base algorithm binaries are located on your PC"
OpenWater_engine_binaries = r'C:\Users\u1066632\LOCAL_WORK\Projects\SandyCreek-IES\src\OW_Model_Packaged\openwater-binaries-20240117\bin' 

#A list of available models should appear
openwater.discovery.OW_BIN = OpenWater_engine_binaries
discover()
model_fn_original = os.path.join(r'C:\Users\u1066632\LOCAL_WORK\Projects\SandyCreek-IES\src\OW_Model_Packaged\OW Models Trimmed\MW', \
                                 'MW_RC_13_WQT_126001A' +'.h5')
model_original = ModelFile(model_fn_original) #Model object

#dissConst_EMC
#Specific logic for 'SednetDissolvedNutrientGeneration'
model_name = ['SednetDissolvedNutrientGeneration', 'SednetDissolvedNutrientGeneration', 'ApplyScalingFactor']
SedNet_parameters_df = model_original.parameters(model_name[0]) #df with parameters for all catchments and FU
#We filter out all rows where area = 0 (no such FU exists in a catchment)
#Unique identifier to index later
SedNet_EMC_parameters_df = SedNet_parameters_df
SedNet_EMC_parameters_df['Identifier'] = model_name[0] + '-' + SedNet_parameters_df['catchment'] \
    + '-' + SedNet_parameters_df['cgu'] + '-' + 'dissConst_EMC'

#dissConst_DWC
#Specific logic for 'SednetDissolvedNutrientGeneration' for model_name = 'SednetDissolvedNutrientGeneration'
SedNet_parameters_df = model_original.parameters(model_name[1]) #df with parameters for all catchments and FU
#We filter out all rows where area = 0 (no such FU exists in a catchment)
#Unique identifier to index later
SedNet_DWC_parameters_df = SedNet_parameters_df
SedNet_DWC_parameters_df['Identifier'] = model_name[1] + '_' + SedNet_parameters_df['catchment'] \
    + '-' + SedNet_parameters_df['cgu'] + '-' + 'dissConst_DWC'

## Delivery Ratios (Surface and Subsurface) for model_name = 'ApplyScalingFactor'
ApplyScalingFactor_parameters_df = model_original.parameters(model_name[2]) #df with parameters for all catchments and FU
#We filter for N_DIN and NLeached
ApplyScalingFactor_parameters_df = \
    ApplyScalingFactor_parameters_df[(ApplyScalingFactor_parameters_df['constituent'] == 'N_DIN') \
                                     | (ApplyScalingFactor_parameters_df['constituent'] == 'NLeached')]

#Unique identifier to index later
ApplyScalingFactor_parameters_df_input = ApplyScalingFactor_parameters_df
ApplyScalingFactor_parameters_df_input['Identifier'] = model_name[2] + '-' + \
    ApplyScalingFactor_parameters_df['catchment'] + '-' \
        + ApplyScalingFactor_parameters_df['constituent']+ '-'+'scale'

#Here we have two dataframes, each corresponding to a model within OW (recall each model inputs/outputs parameters in dataframes with different shapes and column names)
SedNet_parameters_df_SPOTPY = SedNet_parameters_df 
ApplyScalingFactor_parameters_df_SPOTPY = ApplyScalingFactor_parameters_df
model_dataframes = {'SednetDissolvedNutrientGeneration': SedNet_parameters_df_SPOTPY, \
                    'ApplyScalingFactor': ApplyScalingFactor_parameters_df_SPOTPY} #Dict with 2 dataframes as entries

# Here we import a dataframe that contains our parameter: default values, min and max bounds for the calibration algorithm. 
# SpotPy requires parameter inputs in this form.
load_input = pd.read_excel(r'C:\Users\u1066632\LOCAL_WORK\Projects\SandyCreek-IES\src\OW_Model_Packaged\default_starting_params.xlsx')
load_input['identifier'] = load_input['model'] + '-' + load_input['variable'] + '-' + load_input['constituent'] + '-' + load_input['cgu'] 
master_input = load_input[['identifier','optguess', 'min', 'max']]

## SpotPy Calibration Class
## Part 1. Load the model
#Range of model
START='1986/07/01' #Here we specify the whole date range of the model (even though we are calibrating to only 2009-2018), OW doesn't like it if we only run for certain years
END='2023/06/30'
TIME_PERIOD = pd.date_range(START,END)
YEARS=len(TIME_PERIOD)/365.25

# MODEL= 'MW_RC_13_WQT_126001A_copy' #The name of the model that we use in python #Qian: don't think it is used in the code for now, so comment it out.
MODEL_NAME= 'MW_RC_13_WQT_126001A_copy' #The name of the actual .h5 file. Can be the same as the model name we call in python

model_directory = r'C:\Users\u1066632\LOCAL_WORK\Projects\SandyCreek-IES\src\OW_Model_Packaged\OW Models Trimmed\MW' #Location of the model
time_string = datetime.datetime.now().strftime("on %B %d %Y") #Gets the current time and date and converts into str, we use this to record model runs [CURRENTLY NOT USED]

model_dir = os.path.join(model_directory) #The folder where the file is located
model_fn = os.path.join(model_dir,MODEL_NAME+'.h5') #The path for the model itself
new_output_fn = os.path.join(model_directory, MODEL_NAME + '_TEMP'+ '.h5') #Where we want to store the output

##======================Processing outputs===========================#
spot_obj=spot_setup(master_input) # Calibration class into callable object
print('Read Parameters')
parameters = pd.read_csv('parameters.csv')
parameters_vector = []
for i,j in parameters.iterrows():
    scaled_value = (j.upper - j.lower) * j.value/100 + j.lower 
    parameters_vector.append(scaled_value)
annual_loads = spot_obj.simulation(parameters_vector, model_fn, model_dataframes)
output_file = "output.txt"
with open(output_file, 'w') as f:
    f.write('---- CONSTITUENT LOADS ----  \n')
    for j in annual_loads:
        f.write(str(j) + '\n')     