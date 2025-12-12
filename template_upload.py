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
## test_data/test_cgyro_multi_runs/runs/run2/
## test_data/test_cgyro_multi_runs/configs/config1.yaml
## test_data/test_cgyro_multi_runs/configs/config2.yaml

if __name__=="__main__":

    ## Loop over a set of runs to upload 
    ### Parse arguments 
    args = mgk_uploader.f_parse_args()
    drop_keys = ['target','config_file']
    input_args = {k:v for k,v in vars(args).items() if k not in drop_keys}

    ## Get list of run directories    
    d1 = os.path.join(args.target,'runs')
    sub_dirs = sorted([os.path.join(d1,entry.name) for entry in os.scandir(os.path.join(os.getcwd(),d1)) if entry.is_dir()])

    ## Get list of config files
    d2 = os.path.join(args.target,'configs')
    config_files = sorted([os.path.join(d2,entry.name) for entry in os.scandir(os.path.join(os.getcwd(),d2)) if entry.is_file()])

    ## Warning: Ensure that folder and config file lists are in appropriate order

    ## Upload one by one. Metadata will be different for each folder
    for run_dir,fle in zip(sub_dirs,config_files):
        print(f'Uploading dir {run_dir} with config file {fle}')
        mgk_uploader.main_upload(**input_args, target=run_dir, config_file=fle)
        # mgk_uploader.main_upload(**input_args, target=run_dir, config_file = None)
    print("Upload complete!")

