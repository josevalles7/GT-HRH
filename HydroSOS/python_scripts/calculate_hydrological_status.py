"""
BASED ON STATUSCALCV2.R BY KATIE FACER-CHILDS
JOSE VALLES  26-08-2024 
"""
# Import libraries
import argparse
import os
import sys

# Import user-defined librarys 
# Define the path to the python HydroSOS Package
sys.path.append('.')

import HydroSOS_scripts.get_status_product as SOS
import HydroSOS_scripts.aggregate_variable as aggvar

# Define arguments 
parser = argparse.ArgumentParser(
                    prog='calculate_hydrological_status',
                    description='Calculates monthly hydrological status from observation point',
                    epilog='Jose Valles, DINAGUA, 26082024')

# 
parser.add_argument('input_directory', help='input directory, should ONLY contain .csv daily flow timeseries')        
parser.add_argument('output_csv_directory', help='directory files will be saved to as cat_{input_file}.csv') 
parser.add_argument('output_json_directory', help='directory files will be saved to as json files') 

args = parser.parse_args()

for f in os.listdir(args.input_directory):
    if f.endswith('.csv'):
        flowdata, stationid = aggvar.import_data(args.input_directory, filename = f)
        MONTHLY_DISCHARGE = aggvar.calculate_monthly(flowdata)
        MONTHLY_STATUS = SOS.monthly_status(MONTHLY_DISCHARGE)
        SOS.export_csv(MONTHLY_STATUS, 
                       output_directory = args.output_csv_directory,
                       filename = stationid)
        print("==================================================================")

SOS.csv_to_json(args.output_csv_directory, args.output_json_directory)