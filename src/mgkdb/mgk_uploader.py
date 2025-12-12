# -*- coding: utf-8 -*-
"""
Main script to handle uploading a set of runs with same metadata to the MGK database
                    
@author: Austin Blackmon, Dongyang Kuang, Venkitesh Ayyar
"""

import os
import argparse

from support.mgk_file_handling import get_suffixes, upload_to_mongo, isLinear, Global_vars, f_get_linked_oid, f_set_metadata, f_load_config, f_user_input_metadata
from support.mgk_login import f_login_dbase

def f_parse_args():
    #==========================================================
    # argument parser
    #==========================================================
    parser = argparse.ArgumentParser(description='Process input for uploading files')

    parser.add_argument('-T', '--target', help='Target run output folder')
    parser.add_argument('-SIM', '--sim_type', choices=['GENE','CGYRO','TGLF','GS2','GX'], type=str, help='Type of simulation', required=True)
    parser.add_argument('-A', '--authenticate', default = None, help='locally saved login info, a .pkl file')
    parser.add_argument('-C', '--config_file', default = None, help='Configuration file (.yaml) to avoid terminal prompts.')
    parser.add_argument('-D', '--default', default = False, action='store_true', help='Using default inputs for all.')
    parser.add_argument('-V', '--verbose', dest='verbose', default = False, action='store_true', help='output verbose')
    
    parser.add_argument('-Ex', '--extra', dest='extra', default = False, action='store_true', help='whether or not to include any extra files for each suffix')
    parser.add_argument('-L', '--large_files', dest='large_files', default = False, action='store_true', help='whether or not to include large files')
    
    return parser.parse_args()

def main_upload(target, default, sim_type, extra, authenticate, verbose, large_files, config_file):
    '''
    Upload a set of suffixes with common Metadata
    '''
    ### Initial setup 
    upload_folder = os.path.abspath(target)
    global_vars = Global_vars(sim_type)    

    ### Connect to database 
    login = f_login_dbase(authenticate)
    client, database = login.connect()
    with client:
        user = login.login['user']

        # manual_time_flag = not default
        manual_time_flag = False

        print(f'Scanning in {upload_folder} *******************\n')
        linear = isLinear(upload_folder, sim_type)                      
        all_suffixes = get_suffixes(upload_folder, sim_type)

        if not all_suffixes:
            print("Did not find any suffixes in the folder",upload_folder)
            return

        if not default:
            if config_file is not None:
                config_dict = f_load_config(config_file)

                ## to do : Add check for format of config file

                user_input = config_dict['user_input']    
                metadata_info = config_dict['metadata']
                shared_files = user_input['shared_files']

                if user_input['extra_files']:
                    ex_files = user_input['extra_files'].split(',')
                    global_vars.Docs_ex +=ex_files

                ## Add suffixes
                ip_suffixes = user_input['suffixes']
                if ip_suffixes: 
                    suffixes=ip_suffixes.split(',')
                    ## Check if input suffixes exist in the folder
                    incorrect_suffixes = [s1 for s1 in suffixes if s1 not in all_suffixes]
                    if incorrect_suffixes: 
                        print("The following suffixes provided don't exist in this folder")
                        print(f"Found in {upload_folder} these suffixes:\n {all_suffixes}")
                        print('Skipping this folder')
                else:    suffixes = None

                ## Adding shared files info
                run_shared = shared_files.split(',') if shared_files else None
                
                ## Add metadata info
                metadata = metadata_info
                metadata['CodeTag']['sim_type']=sim_type
                linked_id_strg = metadata['DBtag']['linkedObjectID']
                if linked_id_strg is not None:
                    metadata['DBtag']['linkedObjectID'] = f_get_linked_oid(database, linked_id_strg)

                no_prompts         = user_input['no_prompts']
                reupload_if_exists = user_input['reupload_if_exists']

            else: ## Get data through user input 
                if extra: # this will change the global variable
                    ex_files = input('Please type FULL file names to update, separated by comma.\n').split(',')
                    global_vars.Docs_ex +=ex_files

                print("Found in {} these suffixes:\n {}".format(upload_folder, all_suffixes))
                
                suffixes = input('Which run do you want to upload? Separate them by comma. \n Press q to skip. Press ENTER to upload ALL.\n')
                if suffixes == 'q':
                    print("Skipping the folder {}.".format(upload_folder))
                elif len(suffixes):
                    suffixes = suffixes.split(',')
                else:
                    suffixes = None                              
                
                run_shared = input('Do you want to upload any shared files for all suffixes? Please specify path relative to parent folder.\n Separate them by comma. Press Enter to skip.\n')
                if len(run_shared):
                    run_shared = run_shared.split(',')
                else:
                    run_shared = None

                ### Metadata inputs
                user_ip_dict = f_user_input_metadata(database)
                metadata = f_set_metadata(**user_ip_dict,user=user, sim_type=sim_type)

                no_prompts = False
                reupload_if_exists = False
        else:
            suffixes = None
            run_shared = None
            metadata = f_set_metadata(user=user, sim_type=sim_type)
            no_prompts = True
            reupload_if_exists = False
            
        upload_to_mongo(database, linear, metadata, upload_folder, suffixes, run_shared,
                        large_files, verbose, manual_time_flag, global_vars, no_prompts=no_prompts, reupload_if_exists=reupload_if_exists)

def main():

    ### Parse arguments 
    args = f_parse_args()
    input_args = vars(args)
    print(input_args)

    main_upload(**input_args)

## Runner 
if __name__=="__main__":
    main()

# Example command : 
## python mgk_uploader.py -A <fname.pkl> -T test_data/test_gene1_tracer_efit -SIM GENE -C template_user_input.yaml
