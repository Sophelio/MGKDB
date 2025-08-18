import json, bson
import numpy as np
from pyrokinetics import Pyro,template_dir
from pyrokinetics.databases.imas import pyro_to_imas_mapping
import idspy_toolkit as idspy
from idspy_dictionaries import ids_gyrokinetics_local as gkids
from pathlib import Path
import os 


def convert_to_json(obj,separate_real_imag = False):
    """
    This function to recursively goes through GyrokineticsLocal class, and writes to json compatible dictionary

    Parameters
    ----------
    obj : GyrokineticsLocal
        GyrokineticsLocal object with data loaded
    separate_real_imag : bool
        If true, move real and imaginary parts to separate dictionary keys. 
        If False, encode complex np.array

    Returns
    -------
    json compatible dictionary.

    """
    
    if '__dataclass_fields__' in dir(obj):
        tmpdict = {}
        for item in obj.__dataclass_fields__.keys():
            
            value =  eval(f"obj.{item}")
            if separate_real_imag and isinstance(value, np.ndarray) and 'complex' in str(value.dtype).lower():
                #print(value)
                tmpdict[item + "_real"] = convert_to_json(value.real)
                tmpdict[item + "_imaginary"] = convert_to_json(value.imag)
            else:
                tmpdict[item] = convert_to_json(value)
        return tmpdict
    elif isinstance(obj, np.ndarray):
        if 'complex' in str(obj.dtype).lower():
                return dict(
                    __ndarray_tolist_real__=obj.real.tolist(),
                    __ndarray_tolist_imag__=obj.imag.tolist(),
                    dtype=str(obj.dtype),
                    shape=obj.shape,
                )
        else:
            return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_to_json(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [convert_to_json(item) for item in obj]
    else:
        return obj

def f_check_linear(fname,gkcode):
    '''
    Use parameter file to check if a run is linear or non-linear
    '''
    assert os.path.exists(fname), "Cannot find file %s"%(fname)

    if gkcode =='GENE': 
        with open(fname,'r') as f: 
            for line in f: 
                val = line.split('=')
                
                if val[0].strip()=='nonlinear':
                    if (val[1].strip()=='T'):
                        return False ## run is non linear
                    elif (val[1].strip()=='F'):
                        return True 
                    else : 
                        print("Unknown entry in parameter file for field \"nonlinear\" ",line)
                        raise SystemError
        
        print('Could not find \"nonlinear\" field in parameters file')
        raise SystemError
        
    elif gkcode =='CGYRO':
        with open(fname,'r') as f:
            for line in f: 
                val = line.split('=')

                if val[0].strip()=='NONLINEAR_FLAG':
                    if (val[1].strip()=='1'):
                        return False ## run is non linear
                    elif (val[1].strip()=='0'):
                        return True 
                    else : 
                        print("Unknown entry in parameter file for field \"NONLINEAR_FLAG\" ",line)
                        raise SystemError

        print('Could not find \"NONLINEAR_FLAG\" field in parameters file')
        raise SystemError
    
    elif gkcode =='TGLF':
        with open(fname,'r') as f:
            for line in f: 
                val = line.split('=')

                if val[0].strip()=='USE_TRANSPORT_MODEL':
                    true_present  = [strg in val[1].strip() for strg in ['true','T','t']]
                    false_present = [strg in val[1].strip() for strg in ['false','F','f']]
                    if any(true_present): 
                        return False ## run is non linear
                    elif any(false_present):
                        return True 
                    else : 
                        print("Unknown entry in parameter file for field \"USE_TRANSPORT_MODEL\" ",line)
                        raise SystemError

        print('Could not find \"USE_TRANSPORT_MODEL\" field in parameters file')
        raise SystemError

    elif gkcode in ['GX']:
        with open(fname,'r') as f:
            for line in f:
                strg = line.strip().split('#')[0].split('=')
                if (len(strg) == 2 and strg[0].strip() == 'nonlinear_mode'):
                    val = strg[1].strip()
                    if val in ['true','True','t','T']:
                        linear = False
                    elif val in ['false','False','f','F']:
                        linear = True
                    else : 
                        print("Unknown entry in parameter file for field \"nonlinear_mode\" ",line)
                        raise SystemError
                    break
        return linear
    
    elif gkcode in ['GS2']:
        linear = True
        return linear

def update_key_values(data, mod_key, new_value):
    '''
    Module to iteratively scan for entries in a dictionary or sub dictionaries or sublists
    with a specific key and repalce it with a given value
    '''
    
    if isinstance(data, dict):
        for key, value in data.items():
            if key == mod_key:
                data[key] = new_value
            else:
                update_key_values(value, mod_key, new_value)
    elif isinstance(data, list):
        for item in data:
            update_key_values(item, mod_key, new_value)

def prune_imas_gk_dict(gk_dict, pyro, linear):
    ''' Remove 4d files for linear and non-linear runs  
    Also set value of max_repr_length using input values 
    '''

    if linear: # If linear, drop entries in ['linear']['wavevector'][0]['eigenmode'][i]['fields'][key] 
        if gk_dict['non_linear']['fields_4d'] is not None:   
            assert gk_dict['non_linear']['fields_4d']['phi_potential_perturbed_norm']==[], "phi_potential_perturbed_norm field in non_linear, fields_4d is not empty"

        keys_list = ['phi_potential_perturbed_norm','a_field_parallel_perturbed_norm','b_field_parallel_perturbed_norm']
        if gk_dict['linear']['wavevector'] !=[]: 
            for i in range(len(gk_dict['linear']['wavevector'][0]['eigenmode'])): ## For each particle species, delete fields
                if gk_dict['linear']['wavevector'][0]['eigenmode'][0]['fields']: 
                    for key in keys_list:
                        gk_dict['linear']['wavevector'][0]['eigenmode'][i]['fields'][key]=None

    else: # If non-linear, drop  ['non_linear']['fields_4d]
        assert (gk_dict['linear']['wavevector']==[]),"wavevector field in linear is not empty"

        gk_dict['non_linear']['fields_4d']=None


    ## Setup max_repr_length
    if pyro._gk_code=='GENE':
        prec = pyro.gk_input.data["info"]["PRECISION"] 
        prec_dict = {'SINGLE':32, 'DOUBLE':64}
        max_repr_length = prec_dict[prec]

        update_key_values(gk_dict, 'max_repr_length', max_repr_length)

    return gk_dict 


def create_gk_dict_with_pyro(fname,gkcode):
    '''
    Create gyrokinetics dictionary to be upload to database
    '''

    assert gkcode in ['GENE','CGYRO','TGLF','GS2', 'GX'], "invalid gkcode type %s"%(gkcode)

    pyro = Pyro(gk_file=fname, gk_code=gkcode)
    # linear = f_check_linear(fname,gkcode)
    linear = not pyro.numerics.nonlinear

    if gkcode=='TGLF':   
        quasi_linear = pyro.numerics.nonlinear
        linear = True
    else:      
        quasi_linear = False 
    
    if linear:
        pyro.load_gk_output(load_fields=True)
    else: # Loading fields for non-linear runs can take too long, so do not read them 
        pyro.load_gk_output(load_fields=False)

    gkdict = gkids.GyrokineticsLocal()
    # gkdict.max_repr_length = 32 

    idspy.fill_default_values_ids(gkdict)
    gkdict = pyro_to_imas_mapping(
            pyro,
            comment=f"Testing IMAS %s"%(gkcode),
            ids=gkdict
        )
    
    json_data = convert_to_json(gkdict)

    json_data = prune_imas_gk_dict(json_data, pyro, linear)

    ## Confirm IMAS size is less than 2MB
    bson_data = bson.BSON.encode(json_data)
    imas_size = len(bson_data)/1e6 # IMAS size in Megabytes
    if imas_size  > 2.0 : 
        print("Size of IMAS dict: ",imas_size,"MB")
        print("IMAS size is larger than 2MB. Need to check IMAS content")
        raise SystemError
    
    
    return json_data,quasi_linear

if __name__=="__main__":

    # data_dir = "test_data/test_gene1_tracer_efit/"
    # data_dir = "test_data/test_gene2_miller_general/"
    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_gene1_tracer_efit/'
    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_gene2_miller_general/'
    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_gene5_non_st_single_prec/'
    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_gene6_non_st_double_prec/'
    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_gene7_non_st_double_prec/'
    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/GENE_fingerprints_march_2025/pscans_hask/scanfiles0002/'
    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/GENE_fingerprints_march_2025/test/'
    # suffix='_0001'
    # fname = data_dir+'parameters{0}'.format(suffix)
    # gkcode="GENE"
    
    # data_dir='/Users/venkitesh_work/Downloads/downloaded_select_gene_NERSC/tracer_5f34a52cbafb0f9d07b05731/'
    # suffix = '_0002'
    # fname= data_dir+'parameters{0}'.format(suffix)

    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/upload_datasets/non_lin_J78697_x985/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/upload_datasets/test_nonlin_gene_tracer_efit/"
    # data_dir = "/Users/venkitesh_work/Downloads/test_data/scanfiles0000/"
    # suffix='_0001'
    # fname = data_dir+'parameters{0}'.format(suffix)
    # gkcode="GENE"

    # data_dir = "test_data/test_cgyro_multi_runs/run1/"
    # fname = data_dir+'input.cgyro'
    # gkcode="CGYRO"

    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_cgyro_miller/KY_0.30_PX0_0.00/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/CGYRO_linear_scan/ky_0.10/"
    # fname = data_dir+'input.cgyro'
    # gkcode="CGYRO"

    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/CGYRO_nonlinear1_template/run1/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/CGYRO_nonlinear2/run1/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/CGYRO_nonlinear3_no_apar_saved/run1/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/CGYRO_nonlinear6_runs_multi_ebelli/d_lti00_ge00/"
    data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/CGYRO_nonlinear7_oldversion_nov26_2024/run1/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/CGYRO_nonlinear8_aaron_r_0.35_it_x/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/smarts_workflow_test/cygro_it_15/rho0.35_it_x/"
    fname = data_dir+'input.cgyro'
    gkcode="CGYRO"

    # data_dir = "test_data/TGLF/TGLF_linear/"
    # data_dir = "test_data/TGLF/TGLF_transport/"
    # data_dir = "test_data/TGLF/TGLF_transport2/"
    # fname = data_dir+'input.tglf'
    # gkcode="TGLF"

    # data_dir = "test_data/GS2_linear/run1_template/"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/GS2_linear/run1_template"
    # data_dir = "/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/GS2_linear/run4_id17"
    # fname = data_dir+'gs2.in'
    # gkcode="GS2"


    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_gx/1_template_gx_linear/'
    # fname = data_dir+'gx.in'
    # gkcode="GX"

    # data_dir = '/Users/venkitesh_work/Documents/work/Sapient_AI/Data/mgkdb_data/pyro_tests_data/data/test_gx/2_nonlinear_gx/'
    # fname = data_dir + 'cyclone_base.in'
    # gkcode="GX"

    json_data, quasi_linear = create_gk_dict_with_pyro(fname,gkcode)

    with open(data_dir+"gyrokinetics.json", "w") as json_file:
        json.dump(json_data, json_file, indent=4)
        
    ## Testing read from code 
    # with open(data_dir+"gyrokinetics.json", 'r') as j:
    #      contents = json.loads(j.read())

    # for key in json_data.keys():
    #     a,b=json_data[key],contents[key]
    #     print(key,a==b)
    #     if (a!=b):
    #         print("Unequal", key,a,b)
