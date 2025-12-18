import os
import sys
from pathlib import Path

# Add src/mgkdb to Python path so mgk_uploader can find the 'support' module
sys.path.insert(0, str(Path(__file__).parent / "src" / "mgkdb"))

from src.mgkdb import mgk_uploader


# Example command : 
## python template_upload -A <fname.pkl> -T test_data/test_cgyro_multi_runs/ -SIM CGYRO 
## Note, this assumes a structure of the form:
## test_data/test_cgyro_multi_runs/runs/run1/
## test_data/test_cgyro_multi_runs/runs/run2_miller_KY_0.30_PYX0_0.00/
## test_data/test_cgyro_multi_runs/configs/run1.yaml (or run1_config.yaml)
## test_data/test_cgyro_multi_runs/configs/run2_miller_KY_0.30_PYX0_0.00.yaml
## Config files are matched to run directories by folder name, not position.

if __name__=="__main__":

    ## Loop over a set of runs to upload 
    ### Parse arguments 
    args = mgk_uploader.f_parse_args()
    drop_keys = ['target','config_file']
    input_args = {k:v for k,v in vars(args).items() if k not in drop_keys}

    ## Get list of run directories    
    target_abs = os.path.abspath(args.target)
    d1 = os.path.join(target_abs,'runs')
    sub_dirs = sorted([os.path.join(d1,entry.name) for entry in os.scandir(d1) if entry.is_dir()])

    ## Get list of config files
    d2 = os.path.join(target_abs,'configs')
    config_files_dict = {entry.name: os.path.abspath(os.path.join(d2,entry.name)) for entry in os.scandir(d2) if entry.is_file()}

    ## Match config files to run directories by folder name
    ## Upload one by one. Metadata will be different for each folder
    for run_dir in sub_dirs:
        # Extract folder name from run directory path
        run_folder_name = os.path.basename(run_dir)
        
        # Look for config file matching the folder name
        # Try {folder_name}.yaml first, then {folder_name}_config.yaml
        config_file = None
        if f"{run_folder_name}.yaml" in config_files_dict:
            config_file = config_files_dict[f"{run_folder_name}.yaml"]
        else:
            print(f'Warning: No config file found for {run_folder_name}, skipping...')
            continue
        
        print(f'Uploading dir {run_dir} with config file {config_file}')
        mgk_uploader.main_upload(**input_args, target=run_dir, config_file=config_file)
        # mgk_uploader.main_upload(**input_args, target=run_dir, config_file = None)
    print("Upload complete!")

