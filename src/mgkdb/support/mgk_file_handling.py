#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File handling script for formatting output files, getting file lists, and
reading and writing to database containing:
    get_file_list(out_dir,begin):       input GENE output directory and base filepath 
                                           (nrg, energy, etc), return full list of files 
                                           in directory
    get_suffixes(out_dir):            input GENE output directory, return list of run 
                                           suffixes in the directory
    gridfs_put(filepath):               input filepath, upload  file to database, and 
                                           return object_id of uploaded file
    gridfs_read(db_file):               input database filename, return contents of file
    upload_to_mongo   
    isLinear
@author: Austin Blackmon, Dongyang Kuang, Venkitesh Ayyar
"""

import sys
import numpy as np
from bson.objectid import ObjectId
import os
from pathlib import Path
import gridfs
import json
import yaml
from time import strftime
import pickle
from bson.binary import Binary

from .pyro_gk import create_gk_dict_with_pyro
from .ParIO import Parameters
from .diag_plot import diag_plot
from .mgk_post_processing import get_parsed_params, get_suffixes, get_diag_from_run

#=======================================================

class Global_vars():
    '''
    Object to store global variables
    '''
    def __init__(self, sim_type):

        self.set_vars(sim_type)
        #User specified files#
        self.Docs_ex = [] 
        self.update_docs_keys()
        self.troubled_runs = [] # a global list to collection runs where exception happens

    def set_vars(self, sim_type):
        if sim_type=="GENE":

            self.required_files = ['field', 'nrg', 'omega','parameters']
            self.Docs         = ['autopar', 'nrg', 'omega','parameters']
            ## Geometry files added later in processing

            #Large files#
            self.Docs_L = ['field', 'mom', 'vsp']

        elif sim_type=='CGYRO':

            self.required_files =  ['input.cgyro', "out.cgyro.time","out.cgyro.grids","out.cgyro.equilibrium", "bin.cgyro.geo"]
            
            self.Docs = ['input.cgyro', 'input.cgyro.gen', 'input.gacode', 'out.cgyro.info']    
            
            #Large files#
            self.Docs_L = []

        elif sim_type=='TGLF':
            self.required_files = ['input.tglf', 'out.tglf.run']    

            self.Docs = ['input.tglf', 'input.tglf.gen', 'out.tglf.run']    

            #Large files#
            self.Docs_L = []

        elif sim_type=='GX':

            self.required_files = ['gx.in','gx.out.nc']    

            self.Docs = ['gx.in']    

            #Large files#
            self.Docs_L = []

        elif sim_type=='GS2':
            self.required_files = ['gs2.in','gs2.out.nc']    

            self.Docs = ['gs2.in']    

            #Large files#
            self.Docs_L = []
        
        else : 
            print("Invalid simulation type",sim_type)
            raise SystemError 
        
        ### Keys used for filenames 
        self.Keys = [fname.replace('.','_') for fname in self.Docs]

    def update_docs_keys(self):

        self.all_file_docs = self.Docs + self.Docs_ex
        ## Drop any duplicates 
        self.all_file_docs = list(set(self.all_file_docs))

        self.all_file_keys =  [fname.replace('.','_') for fname in self.all_file_docs]
    
    def reset_docs_keys(self,sim_type):
        ## Reset values 
        self.set_vars(sim_type)
        self.update_docs_keys() ## Update extra files
        print("File names and their key names are reset to default!")

def f_load_config(config_file):
    '''
    Load config file instead of user prompts
    '''
    
    with open(config_file) as f:
        config_dict=yaml.load(f, Loader=yaml.SafeLoader)
    
    return config_dict 

def f_check_required_files(global_vars, fldr, suffix, sim_type):

    files_exist=True 

    for fname in global_vars.required_files:
        file = os.path.join(fldr,fname+suffix) if sim_type=='GENE' else os.path.join(fldr,suffix,fname)

        if not os.path.isfile(file) :
            print('Necessary file %s does not exist. Skipping this suffix %s'%(file,suffix))
            files_exist = False
            break
        elif (os.path.getsize(file) == 0):
            print('Necessary file %s is empty. Skipping this suffix %s'%(file,suffix))
            files_exist = False
            break
    
    return files_exist

def f_user_input_metadata(database):
    '''
    Create a dictonary of user inputs for metadata
    Used as keyword arguments to construct metadata dictionary
    '''

    user_ip = {} 
    print("Filling metadata.")
    skip_metadata = input("To skip entering Metadata, please enter 0. press any other key to continue\n")
    if skip_metadata=='0': 
        return user_ip

    print("Please provide input for metadata. Press Enter to skip that entry.\n")

    guidelines = """Confidence Guidelines:
    Confidence | Description
    -----------|----------------------------------------------------------------
    3          | Simulation explicitly checked for convergence on at least one variable, add comments with details on which variable
    2          | Simulation derived or restarted from more rigorously checked simulation in database
    1          | Simulation checked visually for sufficiency, e.g.:
               |   • Fluxes appear to saturate near end of simulation
               |   • Perturbed structures not piling up at edge of simulated domains
               |   • Flux contributions not dominating from wavenumbers at edge of simulated range
    0          | Simulation output was generated
    """

    # Prompt the user for confidence level, displaying the guidelines
    confidence = input(f"{guidelines}\nWhat is your confidence (0-3) for the run? Press ENTER to use default value 0\n") or "0"

    if int(confidence) not in range(4):
        print("Invalid input confidence value. Must be between 0 and 3. Using default confidence value of 0.")
        print(f"Confidence level set to: {confidence}")
        confidence = 0

    user_ip['confidence']= confidence 

    keywords = input('Any keywords to categorize this run. Press Enter to skip\n')
    user_ip['keywords'] = keywords
    
    comments = input('Any comments for data in this folder?Press Enter to skip.\n')
    user_ip['comments'] = comments

    linked_id_strg = input('Do you want to link this run to another existing run in MGKDB. Press Enter to skip\n')
    user_ip['linked_ID'] = f_get_linked_oid(database, linked_id_strg) if linked_id_strg else None

    archive = input('Is there a location where the data is archived? Press Enter to skip.\n')
    user_ip['archive_loc'] = archive

    restart = input('Is this run a restart starting from a different run? For yes -> Y . no -> any other key \n')

    user_ip['restart'] = (restart=='Y')

    if restart=='Y':
        user_ip['restart_timestep'] = int(input('What was the timestep of the previous run used to start this run?\n'))
        initial_run_info = input('Has the initial run been uploaded to this database. For yes -> Y .\n')
        run_oid = ObjectId(input('Please enter the ObjectID for that run\n')) if initial_run_info == 'Y' else None
        if f_check_id_exists(database, run_oid):
            user_ip['initial_run_oid'] =  run_oid 
        else:
            print(f"Entered object id {run_oid} doesn't exist\n")
    
    expt = input('Name of actual or hypothetical experiment? Eg: diiid, iter, sparc, etc. Press Enter to skip.\n')
    user_ip['expt'] = expt

    scenario_id = input('Scenario ID : shot ID or time or runID? Eg: 129913.1500ms . Press Enter to skip.\n')
    user_ip['scenario_runid'] = scenario_id

    git_hash = input('Do you have git-hash to store?Press Enter to skip.\n')
    user_ip['git_hash'] = git_hash

    platform = input('Platform on which this was run? Eg: perlmutter, summit, engaging, pc . Press Enter to skip.\n')
    user_ip['platform'] = platform

    exec_date = input('Execution date?Press Enter to skip.\n')
    user_ip['ex_date'] = exec_date

    workflow = input('Workflow type? Eg: portals, smarts, standalone, etc. Press Enter to skip.\n')
    user_ip['workflow_type'] = workflow

    print("Publication information should be uploaded with a separate script")

    return user_ip

def f_set_metadata(user=None,out_dir=None,suffix=None,keywords=None,confidence=0,comments='Uploaded with default settings.',time_upload=None,\
                   last_update=None, linked_ID=None, expt=None, scenario_runid=None, linear=None, quasiLinear=None, has_1dflux = None, sim_type=None,\
                   git_hash=None, platform=None, ex_date=None, workflow_type=None, archive_loc=None, restart=False, restart_timestep=0, initial_run_oid=None):

    metadata={
        'DBtag': { 
            'user': user,
            'run_collection_name': out_dir,
            'run_suffix': suffix,
            'keywords':keywords,
            'confidence': confidence,
            'comments': comments,
            'time_uploaded': time_upload,
            'last_updated': last_update,
            'linkedObjectID': linked_ID, 
            'archiveLocation': archive_loc,
            'IsRestart': {'restart': restart,'Timestep': restart_timestep,'InitObjectId': initial_run_oid}
        },
        'ScenarioTag': { 
                    'Name of actual of hypothetical experiment': expt,
                    'scenario_runid': scenario_runid,
            },
        'CodeTag': { 
                'sim_type': sim_type,
                'IsLinear': linear,
                'quasi_linear': quasiLinear,
                'Has1DFluxes': has_1dflux,
                'git_hash': git_hash,
                'platform': platform,
                'execution_date': ex_date,
                'workflow_type': workflow_type
            },
        'Publications': { 
                'doi': [] 
            }
    }
    
    return metadata

def get_omega(out_dir, suffix):
    try:
        with open(os.path.join(out_dir, 'omega'+suffix)) as f:
            val = f.read().split()
            
        val = [float(v) for v in val] 
        if len(val) < 3:
            val = val + [np.nan for _ in range(3-len(val))]
    
    except:
        print('Omega file not found. Fill with NaN')
        val = [np.nan for _ in range(3)]
        
    return val
        

def get_time_for_diag(run_suffix):
    option = input('Please enter the tspan information for {}\n 1: Type them in manually.\n 2: Use default settings.\n 3. Use default settings for rest.\n'.format(run_suffix))
    if option == '1':      
        tspan = input('Please type start time and end time, separated by comma.\n').split(',')
        tspan[0] = float(tspan[0])
        tspan[1] = float(tspan[1])
    elif option == '2':
        tspan = None
    else:
        tspan = -1
    
    return tspan

def get_diag_with_user_input(out_dir, suffix,  manual_time_flag):

    if manual_time_flag:
        tspan = get_time_for_diag(suffix)
        if tspan == -1:
            manual_time_flag = False
            Diag_dict = get_diag_from_run(out_dir, suffix, None)
        else:
            Diag_dict = get_diag_from_run(out_dir, suffix, tspan) 
    else:
        Diag_dict = get_diag_from_run(out_dir, suffix, None)
        
    return Diag_dict, manual_time_flag

def get_data(key, *args):
    '''
    Use to get data from default files with functions defined in func_dic
    '''
    return func_dic[key](*args)

def get_data_by_func(user_func, *args):
    '''
    user_func takes args and should return a dictionary having at least two keys: '_header_' and '_data_'
    an example is provided as below: get_data_from_energy()
    '''
    return user_func(*args)

def get_data_from_energy(db, filepath):
    '''
    read GENE energy output, parsed into header and datapart
    '''
    fs = gridfs.GridFS(db)
    if fs.exists({"filepath": filepath}):
        file = fs.find_one({"filepath": filepath}) # assuming only one
        contents = file.read().decode('utf8').split('\n')
        header = []
        data = []
        for line in contents:
            if '#' in line:
                header.append(line)
            else:
                d_str = line.split()
                if d_str:
                    data.append([float(num) for num in d_str])
        
#        data = np.array(data)
        return {'_header_': header[:-1], '_data_': np.array(data)}
    
    else:
        print("No entry in current database matches the specified filepath.")
        return None

def get_data_from_nrg(db, filepath):
    fs = gridfs.GridFS(db)
    if fs.exists({"filepath": filepath}):
        file = fs.find_one({"filepath": filepath}) # assuming only one
        contents = file.read().decode('utf8').split('\n')
        header = []
        data = []
        time = []
        count = 0
        for line in contents[:-1]: # last line is ''
            if count % 2 == 0:
#               print(count)
               time.append(float(line))
            else:
                d_str = line.split()
                if d_str:
                    data.append([float(num) for num in d_str])
            count += 1
        
#        data = np.array(data)
        return {'_header_': header, '_time_': np.array(time), '_data_': np.array(data)}
    
    else:
        print("No entry in current database matches the specified filepath.")
        return None

def isfloat(a):
    try:
        float(a)
        return True
    except ValueError:
        return False

def to_float(a):
    try:
        b = float(a)
    except ValueError:
        b = a
    return b

def get_data_from_parameters(db, filepath):
    fs = gridfs.GridFS(db)
    if fs.exists({"filepath": filepath}):
        file = fs.find_one({"filepath": filepath}) # assuming only one
        contents = file.read().decode('utf8').split('\n')
        summary_dict=dict()
        for line in contents:
            if '&' in line:
                category = line[1:]
                sub_dict = {}
            elif '=' in line:
                pars = line.split('=')
                sub_dict[pars[0].rstrip()] = to_float(pars[1]) 
            elif '/' in line:
                summary_dict[category] = sub_dict            
            else:
                continue
            
        return summary_dict
    
    else:
        print("No entry in current database matches the specified filepath.")
        return None
  

def get_data_from_tracer_efit(db, filepath):      
    fs = gridfs.GridFS(db)
    if fs.exists({"filepath": filepath}):
        file = fs.find_one({"filepath": filepath}) # assuming only one
        contents = file.read().decode('utf8').split('\n')
        header_dict = {}
        data = []
        for line in contents:
            if '=' in line:
                item = line.split('=')
#                if '_' in item[1] or ' \' ' in item[1]:
                if isfloat(item[1]):
                    header_dict[item[0]] = float(item[1])
                else:
                    header_dict[item[0]] = item[1]
                    
            elif '/' in line or '&' in line:
                continue

            else:
                d_str = line.split()
                if d_str:
                    data.append([float(num) for num in d_str])
        
#        data = np.array(data)
        return {'_header_': header_dict, '_data_': np.array(data)}
    
    else:
        print("No entry in current database matches the specified filepath.")
        return None
    
func_dic = {'energy': get_data_from_energy,
            'nrg': get_data_from_nrg,
            'parameters': get_data_from_parameters
            }        

def get_file_list(out_dir, begin):
    '''
    Get files from out_dir that begins with "begin"
    '''
    files_list = []
    
    #unwanted filetype suffixes for general list
    bad_ext = ('.ps','.png', '.jpg', '.dat~', '.h5')
    
#    print('Searching in {} with key {}'.format(out_dir, begin))
    #scan files in GENE output directory, ignoring files in '/in_par', and return list
    
#    files = next(os.walk(out_dir))[2]
    files = os.listdir(out_dir)
    for count, name in enumerate(files, start=0):
        if name.startswith(begin) and name.endswith(bad_ext) == False: #and not os.path.isdir('in_par'):
            file = os.path.join(out_dir, name)
            if file not in  files_list:
                files_list.append(file)
            
    # print('{} files found in {} beginning with {}.'.format(len(files_list), out_dir, begin) )
    return files_list     


# def get_file_list(out_dir,fname):


def gridfs_put(db, filepath, sim_type):

    fs = gridfs.GridFS(db)
    with open(filepath, 'rb') as file:
        dbfile = fs.put(file, encoding='UTF-8', 
                        filepath=filepath,
                        filename=os.path.basename(filepath),
                        simulation_type=sim_type,
                        metadata=None)

    return dbfile
    
def gridfs_read(db, query):
    fs = gridfs.GridFS(db)
    file = fs.find_one(query)
    contents = file.read()
    return(contents)

def Array2Dict_dim1(npArray, key_names=None):
    '''
    Convert a 1d numpy array to dictionary
    '''
    assert len(npArray.shape) == 1, "Dimension of input numpy array should be 1."
    
    arraydict = dict()
    
    if key_names:
        for i in range(len(npArray)):
            arraydict[key_names[i]] = npArray[i]
    
    else:
        for i in range(len(npArray)):
            arraydict[str(i)] = npArray[i]
    
    return arraydict

def Array2Dict_dim2(npArray, row_keys=None, col_keys=None):
    '''
    Convert a 2d numpy array to dictionary
    '''
    assert len(npArray.shape) == 2, "Dimension of input numpy array should be 2."
    
    arraydict = dict()
    
    nrows, ncols = np.shape(npArray)
    if row_keys and col_keys:
        for i in range(nrows):
            row_dict = {}
            for j in range(ncols):
                row_dict[col_keys[j]] = npArray[i,j]
            arraydict[row_keys[i]] = row_dict
    
    else:
        for i in range(nrows):
            row_dict = {}
            for j in range(ncols):
                row_dict[str(j)] = npArray[i,j]
            arraydict[str(i)] = row_dict
        
    return arraydict

def Rep_OID(dic):
    '''
    Check a dictionary tree and replace any 'ObjectId' string to ObjectId object
    '''
    for key, val in dic.items():
        if isinstance(val, str) and 'ObjectId' in val:
#            oid_str = val[8:-1]
            oid_str = val[val.find('(')+1: val.find(')')].strip()
            dic[key] = ObjectId(oid_str)

        elif isinstance(val, dict):
            dic[key] = Rep_OID(val)
    return dic

def Str2Query(s):
    '''
    Convert a string s to python dictionary for querying the database
    '''
    q_dict = json.loads(s)
    q_dict = Rep_OID(q_dict)
    
    return q_dict

def get_oid_from_query(db, collection, query):
    
    records_found = collection.find(query)
    
    oid_list = []
    
    for record in records_found:
        oid_list.append(record['_id'])
        
    return oid_list

def _npArray2Binary(npArray):
    """Utility method to turn an numpy array into a BSON Binary string.
    utilizes pickle protocol 2 (see http://www.python.org/dev/peps/pep-0307/
    for more details).

    Called by stashNPArrays.

    :param npArray: numpy array of arbitrary dimension
    :returns: BSON Binary object a pickled numpy array.
    """
    return Binary(pickle.dumps(npArray, protocol=2), subtype=128 )

def _binary2npArray(binary):
    """Utility method to turn a a pickled numpy array string back into
    a numpy array.

    Called by loadNPArrays, and thus by loadFullData and loadFullExperiment.

    :param binary: BSON Binary object a pickled numpy array.
    :returns: numpy array of arbitrary dimension
    """
    return pickle.loads(binary)

def gridfs_put_npArray(db, value, filepath, filename, sim_type):
    '''
    Write numpy array to file and then to DB
    '''
    MAX_FILE_SIZE = 10 * 1024 * 1024  # 10 MB limit

    fs = gridfs.GridFS(db)
    binary_data = _npArray2Binary(value)
    data_size = len(binary_data)

    if data_size > MAX_FILE_SIZE: ## Ensure file is not too big
        print(f"Binary data for '{filename}' has size {data_size / (1024 * 1024)} MB. This exceeds the size limit of {MAX_FILE_SIZE / (1024 * 1024)} MB")
        print(f"Ignoring upload of this diagnostic: {filename}")
        obj_id = None
    else: 
        obj_id=fs.put(binary_data,encoding='UTF-8',
                    filename = filename,
                    simulation_type = sim_type,
                    filepath = filepath)
    
    return obj_id  
    
    
def load(db, collection, query, projection={'Metadata':1, 'gyrokineticsIMAS':1, 'Diagnostics':1}, getarrays=True):
    """Preforms a search using the presented query. For examples, see:
    See http://api.mongodb.org/python/2.0/tutorial.html
    The basic idea is to send in a dictionaries which key-value pairs like
    mdb.load({'basename':'ag022012'}).

    :param query: dictionary of key-value pairs to use for querying the mongodb
    :returns: List of full documents from the collection
    """
    
    results = collection.find(query, projection)
    
    if getarrays:
        allResults = [_loadNPArrays(db, doc) for doc in results]
    else:
        allResults = [doc for doc in results]
    
    if allResults:
#        if len(allResults) > 1:
#            return allResults
#        elif len(allResults) == 1:
#            return allResults[0]
#        else:
#            return None
        return allResults
    else:
        return None
    
def _loadNPArrays(db, document):
    """Utility method to recurse through a document and gather all ObjectIds and
    replace them one by one with their corresponding data from the gridFS collection

    Skips any entries with a key of '_id'.

    Note that it modifies the document in place.

    :param document: dictionary like-document, storable in mongodb
    :returns: document: dictionary like-document, storable in mongodb
    """
    fs = gridfs.GridFS(db)
    for (key, value) in document.items():
        if isinstance(value, ObjectId) and key != '_id':
            document[key] = _binary2npArray(fs.get(value).read())
        elif isinstance(value, dict):
            document[key] = _loadNPArrays(db, value)
    return document

def query_plot(db, collection, query, projection = {'Metadata':1, 'Diagnostics':1}):
    data_list = load(db, collection, query, projection)
    print('{} records found.'.format(len(data_list)))
    
    data_to_plot = [diag_plot(da) for da in data_list]
    
    for i in range(len(data_to_plot)):
         data_to_plot.plot_all()    
    
    
def isLinear(folder_name, sim_type):
    linear = None

    #check file for 'nonlinear' value
    suffixes = get_suffixes(folder_name, sim_type)
    
    if len(suffixes):
        suffix = suffixes[0] #assuming all parameters files are of the same linear/nonlinear type
        print('Scanning parameters for deciding linear/Nonlinar.')
    else:
        suffix = ''

    if sim_type=='GENE': 
        fname = os.path.join(folder_name, 'parameters' + suffix)
        if os.path.isfile( fname ):
            par = Parameters()
            par.Read_Pars( fname )
            pars = par.pardict
            linear = not pars['nonlinear']
            return(linear)
            
        #check folder name for nonlin
        elif folder_name.find('nonlin') != -1:
            linear = False
            return(linear)
        
        #check folder name for linear
        elif folder_name.find('linear') != -1:
            linear = True 
            return(linear)

        else:
            assert linear is None, "Can not decide, please include linear/nonlin as the suffix of your data folder!"
        
    elif sim_type=='CGYRO':
        fname=os.path.join(folder_name, suffix, 'input.cgyro')
        assert os.path.isfile(fname),"File %s does not exist"%(fname)

        non_lin = None
        with open(fname,'r') as f:
            for line in f: 
                if line.split('=')[0].strip()=='NONLINEAR_FLAG':
                    non_lin = int(line.split('=')[1].strip()) # Remove blank space to just get 0 or 1
        
        assert non_lin is not None, "Didn't find NONLINEAR_FLAG in file(%s)"%(fname) 

        linear = False if non_lin else True
        return linear
    
    elif sim_type=='TGLF':

        fname = os.path.join(folder_name, suffix, 'input.tglf')
        assert os.path.isfile(fname),"File %s does not exist"%(fname)

        with open(fname,'r') as f:
            for line in f: 
                val = line.split('=')

                if val[0].strip()=='USE_TRANSPORT_MODEL':
                    true_present  = [strg in val[1].strip() for strg in ['true','T','t']]
                    false_present = [strg in val[1].strip() for strg in ['false','F','f']]
                    if any(true_present): 
                        linear = False ## run is non linear
                    elif any(false_present):
                        linear = True 
                    else : 
                        print("Unknown entry in parameter file for field \"USE_TRANSPORT_MODEL\" ",line)
                        raise SystemError
                    break
        return linear
    
    elif sim_type=='GX':

        in_files = [f for f in os.listdir(os.path.join(folder_name,suffix)) if f.endswith('.in')]
        
        if len(in_files) > 1:
            print("Expected exactly one .in file in the folder, found %d. Using first one."%(in_files))
        fname=os.path.join(folder_name,suffix,in_files[0])
        assert os.path.isfile(fname),"File %s does not exist"%(fname)

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
    
    elif sim_type=='GS2': # Only support linear GS2 for now
        return True

def isUploaded(out_dir,runs_coll):
    '''
    check if out_dir appears in the database collection.
    Assuming out_dir will appear no more than once in the database
    '''
    inDb = runs_coll.find({ "Metadata.DBtag.run_collection_name": out_dir })

    entries = list(inDb)
    uploaded = True if len(entries)>0 else False 
    
    # uploaded = False
    # for run in inDb:
    #     if run["Metadata"]["run_collection_name"] == out_dir: # seems redundent?
    #         uploaded = True
    #         break
    
    return uploaded

def not_uploaded_list(out_dir, runs_coll, write_to = None):
    '''
    Get all subfolders in out_dir that are not in the database yet
    '''
    not_uploaded = []
    for dirpath, dirnames, files in os.walk(out_dir):
        if str(dirpath).find('in_par') == -1 and str(files).find('parameters') != -1:
            if not isUploaded(dirpath, runs_coll):
                not_uploaded.append(dirpath)
    
    if write_to is not None and len(not_uploaded):
        with open(os.path.abspath(write_to),'w') as f:
            f.write('\n'.join(not_uploaded))
    
    return not_uploaded

def get_record(out_dir, runs_coll):
    '''
    Get a list of summary dictionary for 'out_dir' in the database
    '''
    inDb = runs_coll.find({ "Metadata.DBtag.run_collection_name": out_dir })
    record = []
    for run in inDb:
#        dic = dict()
#        for key, val in run.items():
#            dic[key] = val
#        record.append(dic)
        record.append(run)
    return record
   
def f_check_id_exists(db, _id):
    ''' Given an object ID, check if it exists in linear or nonlinear collections
    '''

    for collection in ['LinearRuns','NonlinRuns']:
        runs_coll =  getattr(db,collection)

        try: 
            record = runs_coll.find_one({"_id": _id})
        except Exception as e:
            # print(e)
            print("Invalid object ID",_id)
            return False
    
        if record is not None: return True

    print("Entry %s not found in database, please double check the id"%(_id))
    return False

def f_get_linked_oid(database, linked_id_string):
    '''
    Get linked ObjectID
    '''

    if linked_id_string is not None: 
        oid = ObjectId(linked_id_string)
        id_exists = f_check_id_exists(database, oid)

        if id_exists:
            print("Linked OID %s"%(oid))
            return oid
        else :
            return None

def download_file_by_path(db, filepath, destination, revision=-1, session=None):
    '''
    db: database name
    filepath: filepath stored in database, that is "db.fs.files['filepath']"
    destination: local path to put the file
    
    Attention: filename may correspond to multiple entries in the database
    '''
    fs = gridfs.GridFSBucket(db)
    records = db.fs.files.find({"filepath": filepath})
    count = 0
    for record in records:
        _id = record['_id']
        filename = record['filepath'].split('/')[-1]
        with open(os.path.join(destination, filename+'_mgk{}'.format(count) ),'wb+') as f:
            fs.download_to_stream(_id, f)
            count +=1
#            fs.download_to_stream_by_name(filename, f, revision, session)
        
    print("Download completed! Downloaded: {}".format(count))
    
def download_file_by_id(db, _id, destination, fname=None, session = None):
    '''
    db: database name
    _id: object_id
    destination: local path to put the file
    fname: name you want to call for the downloaded file
    '''

    fs = gridfs.GridFSBucket(db)
    if not fname:
        fname = db.fs.files.find_one(_id)['filename']
    if not os.path.exists(destination):
        Path(destination).mkdir(parents=True) 
    with open(os.path.join(destination, fname),'wb+') as f:   
        fs.download_to_stream(_id, f)
    print("Download completed!")
    
def download_dir_by_name(db, runs_coll, dir_name, destination):  
    '''
    db: database name
    dir_name: as appear in db.Metadata['run_collection_name']
    destination: destination to place files
    '''
    path = os.path.join(destination, dir_name.split('/')[-1])
    if not os.path.exists(path):    
        try:
            #os.mkdir(path)
            Path(path).mkdir(parents=True) 
        except OSError:
            print ("Creation of the directory %s failed" % path)
    #else:
    fs = gridfs.GridFSBucket(db)
    inDb = runs_coll.find({ "Metadata.DBtag.run_collection_name": dir_name })

    if 'generr' in inDb[0]['Files'].keys(): ## Fix for when 'generr' doesn't exist
        if inDb[0]['Files']['geneerr'] != 'None':    
            with open(os.path.join(path, 'geneerr.log'),'wb+') as f:
                fs.download_to_stream(inDb[0]['Files']['geneerr'], f, session=None)

    for record in inDb:
        '''
        Download 'files'
        '''
        for key, val in record['Files'].items():
            if val != 'None' and key not in ['geneerr']:
                filename = db.fs.files.find_one(val)['filename']
                with open(os.path.join(path, filename),'wb+') as f:
#                    fs.download_to_stream_by_name(filename, f, revision=-1, session=None)
                    fs.download_to_stream(val, f, session=None)
                record['Files'][key] = str(val)
        if 'generr' in record['Files'].keys():  ## Fix for when 'generr' doesn't exist 
            record['Files']['geneerr'] = str(record['Files']['geneerr'])
        
        '''
        Deal with diagnostic data
        '''
        diag_dict={}
        fsf=gridfs.GridFS(db)
        for key, val in record['Diagnostics'].items():
            if isinstance(val, ObjectId):
#                data = _loadNPArrays(db, val)
#                data = _binary2npArray(fsf.get(val).read()) # no need to store data
                record['Diagnostics'][key] = str(val)
#                data = _binary2npArray(fsf.get(val).read()) 
#                np.save( os.path.join(path,str(record['_id'])+'-'+key), data)
                diag_dict[key] = _binary2npArray(fsf.get(val).read())
            
        with open(os.path.join(path,str(record['_id'])+'-'+'diagnostics.pkl'), 'wb') as handle:
            pickle.dump(diag_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
                
        record['_id'] = str(record['_id'])
        with open(os.path.join(path, 'mgkdb_summary_for_run'+record['Metadata']['run_suffix']+'.json'), 'w') as f:
            json.dump(record, f)
           
    print ("Successfully downloaded to the directory %s " % path)


def download_runs_by_id(db, runs_coll, _id, destination):
    '''
    Download all files in collections by the id of the summary dictionary.
    '''
    
    fs = gridfs.GridFSBucket(db)
    record = runs_coll.find_one({ "_id": _id })
    try:
        dir_name = record['Metadata']['DBtag']['run_collection_name']
    except TypeError:
        print("Entry not found in database, please double check the id")
        raise SystemExit
        
    path = os.path.join(destination, dir_name.split('/')[-1])

    if not os.path.exists(path):
        try:
#            path = os.path.join(destination, dir_name.split('/')[-1])
            #os.mkdir(path)
            Path(path).mkdir(parents=True)
        except OSError:
            print ("Creation of the directory %s failed" % path)
    #else:
    '''
    Download 'files'
    '''
    for key, val in record['Files'].items():
        if val != 'None':
            filename = db.fs.files.find_one(val)['filename']
            #print(db.fs.files.find_one(val)).keys()
            with open(os.path.join(path, filename),'wb+') as f:
#                    fs.download_to_stream_by_name(filename, f, revision=-1, session=None)
                fs.download_to_stream(val, f, session=None)
            record['Files'][key] = str(val)
            
    '''
    Deal with diagnostic data
    '''
    fsf=gridfs.GridFS(db)
    diag_dict = {}
    for key, val in record['Diagnostics'].items():
        if isinstance(val, ObjectId):
#                data = _binary2npArray(fsf.get(val).read()) # no need to store data
            record['Diagnostics'][key] = str(val)
#            data = _binary2npArray(fsf.get(val).read()) 
#            np.save( os.path.join(path,str(record['_id'])+'-'+key), data)
            diag_dict[key] = _binary2npArray(fsf.get(val).read())
            
    with open(os.path.join(path,str(record['_id'])+'-'+'diagnostics.pkl'), 'wb') as handle:
        pickle.dump(diag_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)        

    #print(record)
    record['_id'] = str(_id)

    with open(os.path.join(path, 'mgkdb_summary_for_run'+record['Metadata']['DBtag']['run_suffix']+'.json'), 'w') as f:
        json.dump(record, f)
    print("Successfully downloaded files in the collection {} to directory {}".format( record['_id'],path) )   
    
def update_mongo(db, metadata, out_dir, runs_coll, linear, suffixes=None):

    '''
    only update file related entries, no comparison made before update
    '''
    
    sim_type = metadata['CodeTag']['sim_type']

    fs = gridfs.GridFS(db)
    if suffixes is None:
        suffixes = get_suffixes(out_dir, sim_type)  
        
    update_option = input('Enter options for update:\n 0: Files shared by all runs, usually do not have a suffix. \n 1: Unique files used per run. Specify the keywords and suffixes. \n ')
    if update_option == '0':
        files_to_update = input('Please type FULL file names to update, separated by comma.\n').split(',')
        keys_to_update = input('Please type key names for each file you typed, separated by comma.\n').split(',')

        updated = []
        print('Uploading files .......')
        # update the storage chunk
        for doc, field in zip(files_to_update, keys_to_update):
            
            file = os.path.join(out_dir, doc)
            assert os.path.exists(file), "File %s not found"%(file)
            
            # delete ALL history
            grid_out = fs.find({'filepath': file})
            for grid in grid_out:
                print('File with path tag:\n{}\n'.format(grid.filepath) )
                fs.delete(grid._id)
                print('deleted!')

            with open(file, 'rb') as f:
                _id = fs.put(f, encoding='UTF-8', filepath=file, filename=file.split('/')[-1])
            
            updated.append([field, _id])
        
        # update the summary dictionary  
        print('Updating Metadata')              
        for entry in updated:
            for suffix in suffixes:                    
                runs_coll.update_one({ "Metadata.DBtag.run_collection_name": out_dir, "Metadata.DBtag.run_suffix": suffix}, 
                                 {"$set":{'Files.'+entry[0]: entry[1], 
                                          "Metadata.DBtag.last_updated": strftime("%y%m%d-%H%M%S")}}
                                 )
        print("Update complete")
                
    elif update_option == '1':
        files_to_update = input('Please type filenames (without suffixes) for files to update, separated by comma.\n').split(',')
        print("suffixes availables are:{}".format(suffixes))
        runs_to_update = input('Please type which suffixes to update, separated by comma. If you need to update all runs, just hit ENTER. \n').split(',')      
        # affect_QoI = input('Will the file change QoIs/Diagnostics? (Y/N)')
        affect_QoI = True

#        updated = []
        # update the storage chunk
        print('Uploading files .......')
        if len(runs_to_update[0]) != 0:
            run_suffixes = runs_to_update
        else:
            run_suffixes = suffixes
        
        for doc in files_to_update:
            manual_time_flag = True
            for suffix in run_suffixes:
                if affect_QoI:
                    input_fname = f_get_input_fname(out_dir, suffix, sim_type)
                    GK_dict, quasi_linear = create_gk_dict_with_pyro(input_fname, sim_type)   

                    if sim_type in ['CGYRO','TGLF','GS2','GX']:
                        Diag_dict = {}
                    elif sim_type=='GENE': 
                        Diag_dict, manual_time_flag = get_diag_with_user_input(out_dir, suffix, manual_time_flag)

                    run = runs_coll.find_one({ "Metadata.DBtag.run_collection_name": out_dir, "Metadata.DBtag.run_suffix": suffix})
                    for key, val in run['Diagnostics'].items():
                        if val != 'None':
                            # print((key, val))
                            fs.delete(val)
                            # print('deleted!')

                    for key, val in Diag_dict.items():
                        Diag_dict[key] = gridfs_put_npArray(db, Diag_dict[key], out_dir, key, sim_type)

                    runs_coll.update_one({ "Metadata.DBtag.run_collection_name": out_dir, "Metadata.DBtag.run_suffix": suffix },
                            { "$set": {'gyrokineticsIMAS': GK_dict, 'Diagnostics':Diag_dict}}
                                 )

                file = os.path.join(out_dir, doc  + suffix)
                grid_out = fs.find({'filepath': file})
                for grid in grid_out:
                    print('File with path tag:\n{}\n'.format(grid.filepath) )
                    fs.delete(grid._id)
                    print('deleted!')
                
                with open(file, 'rb') as f:
                    _id = fs.put(f, encoding='UTF-8', filepath=file, filename=file.split('/')[-1])

                runs_coll.update_one({ "Metadata.DBtag.run_collection_name": out_dir, "Metadata.DBtag.run_suffix": suffix }, 
                                 { "$set": {'Files.'+ doc: _id, "Metadata.DBtag.last_updated": strftime("%y%m%d-%H%M%S")} }
                                 )
        print("Update complete")
    
    else:
        print('Invalid input. Update aborted.')
        pass
    
def remove_from_mongo(out_dir, db, runs_coll):
    #find all documents containing collection name
        
    inDb = runs_coll.find({ "Metadata.DBtag.run_collection_name": out_dir })        
    fs = gridfs.GridFS(db)
    for run in inDb:
        # delete the gridfs storage:
        for key, val in run['Files'].items():
#            print(val)
#            if (key in all_file_keys) and val != 'None':
##                print((key, val))
#                target_id = ObjectId(val)
#                print((key, target_id))
#                fs.delete(target_id)
#                print('deleted!')
            if val != 'None':
                print((key, val))
                fs.delete(val)
                print('deleted!')
#                if fs.exists(target_id):
#                    print("Deleting storage for entry \'{}\' deleted with id: {}").format(key, val)
#                    fs.delete(target_id)
#                    print("Deleted!")
                
        for key, val in run['Diagnostics'].items():
            if val != 'None':
                print((key, val))
                fs.delete(val)
                print('deleted!')
                
#        delete the header file
        runs_coll.delete_one(run)
        
def f_update_global_var(global_vars, out_dir, suffix, sim_type, is_linear, large_files, count):
    '''
    Update global variables for each suffix
    '''
    ## For GENE, add geometry file to required and main files 
    ## Also modify large files when encountering many species 

    if sim_type=='GENE': ## Add geometry file for GENE 
        par = Parameters()
        par_file = os.path.join(out_dir,'parameters'+suffix)
        par.Read_Pars(par_file)
        pars = par.pardict
        ## Get geometry from parameters file and add that to list of files to save
        if 'magn_geometry' in pars:
            geom_file = pars['magn_geometry'][1:-1]
            global_vars.required_files.append(geom_file)
            global_vars.Docs.append(geom_file)

        n_spec = pars['n_spec']
        
        if large_files:
            if 'name1' in pars and 'mom' in global_vars.Docs_L:
                # global_vars.Docs_L.pop(global_vars.Docs_L.index('mom'))
                global_vars.Docs_L.remove('mom')
                for i in range(n_spec): # adding all particle species
                    global_vars.Docs_L.append('mom_'+pars['name{}'.format(i+1)][1:-1])

        if not is_linear: ## No omega file needed for non-linear runs
            global_vars.required_files.remove('omega')
            # global_vars.Docs.remove('omega')

    if large_files: 
        global_vars.Docs +=global_vars.Docs_L
    global_vars.update_docs_keys()

def f_get_full_fname(sim_type, fldr, suffix, fname):
    '''
    Give full file name for sim_type using suffix and filename
    '''
    if sim_type=='GENE':
        full_fname = os.path.join(fldr,fname+suffix)
    else: 
        full_fname = os.path.join(fldr,suffix,fname)

    return full_fname

def upload_file_chunks(db, out_dir, sim_type, suffix=None, run_shared=None, global_vars=None):
    '''
    This function does the actual uploading of gridfs chunks and
    returns object_ids for the chunk. If a ValueError is raised due to file size,
    it returns the current state of the dictionaries along with an error flag.
    '''
    
    MAX_FILE_SIZE = 10 * 1024 * 1024  # 10 MB limit 

    _docs, _keys = global_vars.all_file_docs, global_vars.all_file_keys

    file_upload_dict = {}
    s_dict = {}
    error_occurred = False
    error_message = None

    try:
        for doc, key in zip(_docs, _keys):
            file = f_get_full_fname(sim_type, out_dir, suffix, doc)
            file_upload_dict[key] = {'full_fname': file, 'oid': None}

            if os.path.isfile(file):
                ## Ensure file is not too big
                file_size = os.path.getsize(file)

                if file_size > MAX_FILE_SIZE:
                    raise ValueError(
                        "Size of the file %s is %s MB and it exceeds size limit of %s MB" %
                        (file, file_size / (1024 * 1024), MAX_FILE_SIZE / (1024 * 1024))
                    )

                _id = gridfs_put(db, file, sim_type)
                file_upload_dict[key]['oid'] = _id
            else:
                print(f'{file} not found in {out_dir}')

        
        ### Add shared files 
        if isinstance(run_shared, list):
            for sh in run_shared:
                key = sh.replace('.', '_')
                file = os.path.join(out_dir, sh)
                s_dict[key] = {'full_fname': file, 'oid': None}

                if os.path.isfile(file):
                    ## Ensure file is not too big
                    file_size = os.path.getsize(file)

                    if file_size > MAX_FILE_SIZE:
                        raise ValueError(
                            "Size of the file %s is %s MB and it exceeds size limit of %s MB" %
                            (file, file_size / (1024 * 1024), MAX_FILE_SIZE / (1024 * 1024))
                        )

                    _id = gridfs_put(db, file, sim_type)
                    s_dict[key]['oid'] = _id
                else:
                    print(f'{file} not found in {out_dir}')

        return file_upload_dict, s_dict, error_occurred, error_message

    except ValueError as e:
        error_occurred = True
        error_message = str(e)
        return file_upload_dict, s_dict, error_occurred, error_message

def f_get_input_fname(out_dir, suffix, sim_type):
    ''''
    Get the name of the input file with suffix for the simluation type
    '''

    fname_dict = {'CGYRO':os.path.join(out_dir,suffix,'input.cgyro'),
                    'TGLF':os.path.join(out_dir,suffix,'input.tglf'),
                    'GENE':os.path.join(out_dir,'parameters{0}'.format(suffix)),
                    'GS2': os.path.join(out_dir,suffix,'gs2.in'),
                    'GX': os.path.join(out_dir,suffix,'gx.in')
                }

    return fname_dict[sim_type]

def upload_runs(db, metadata, out_dir, is_linear=True, runs_coll = None, suffixes=None, run_shared=None,
                large_files=False, verbose=True, manual_time_flag=True, global_vars=None):
    """
    Uploads simulation run data to the database, handling both linear and nonlinear runs.

    Parameters:
    - db: Database connection object.
    - metadata: Dictionary containing metadata for the run.
    - out_dir: Output directory containing simulation files.
    - is_linear: Boolean indicating if the run is linear (True) or nonlinear (False). Default: True.
    - runs_coll: Name of collection in Database (Either db.LinearRuns or db.NonlinRuns)
    - suffixes: List of suffixes for files to upload. If None, determined automatically.
    - run_shared: List of shared files to upload (optional).
    - large_files: Boolean to handle large file uploads. Default: False.
    - verbose: Boolean to print detailed output. Default: True.
    - manual_time_flag: Boolean to handle user-specified time spans for diagnostics. Default: True.
    - global_vars: Object containing global variables for the upload process.

    Returns:
    None
    """
    sim_type = metadata['CodeTag']['sim_type']

    # Update files dictionary
    if suffixes is None:
        suffixes = get_suffixes(out_dir, sim_type)

    if isinstance(run_shared, list):
        shared_not_uploaded = [True for _ in run_shared]
    else:
        shared_not_uploaded = [False]
    shared_file_dict = {}

    for count, suffix in enumerate(suffixes):
        try:
            print('='*40)
            print(f'Working on files with suffix: {suffix} in folder {out_dir}.......')

            f_update_global_var(global_vars, out_dir, suffix, sim_type, is_linear, large_files, count)

            uploaded_ids = {}

            files_exist = f_check_required_files(global_vars, out_dir, suffix, sim_type)
            assert files_exist, "Required files don't exist. Skipping folder"

            # Compute gyrokinetics IMAS using pyrokinetics package
            print("Computing gyrokinetics IMAS using pyrokinetics")
            input_fname = f_get_input_fname(out_dir, suffix, sim_type)
            GK_dict, quasi_linear = create_gk_dict_with_pyro(input_fname, sim_type)

            # Upload files to DB
            print('Uploading files ....')

            if count == 0:
                f_dict, s_dict,err_occured, err_msg = upload_file_chunks(db, out_dir, sim_type, suffix, run_shared, global_vars)
                shared_file_dict = {k: v['oid'] for k, v in s_dict.items()}
            else:
                f_dict, s_dict, err_occured, err_msg = upload_file_chunks(db, out_dir, sim_type, suffix, None, global_vars)

            files_dict = {k: v['oid'] for k, v in f_dict.items()}
            files_dict = {**files_dict, **shared_file_dict}
            uploaded_ids = {k:v for k,v in files_dict.items() if v is not None}
            if err_occured: 
                print('Error occured during input file upload')
                raise ValueError(err_msg)

            print('='*60)
            # Metadata dictionary
            time_upload = strftime("%y%m%d-%H%M%S")

            metadata['DBtag']['run_collection_name'] = out_dir
            metadata['DBtag']['run_suffix'] = '' + suffix
            metadata['DBtag']['time_uploaded'] = time_upload
            metadata['DBtag']['last_updated'] = time_upload
            metadata['CodeTag']['IsLinear'] = is_linear
            metadata['CodeTag']['quasi_linear'] = quasi_linear
            metadata['CodeTag']['Has1DFluxes'] = GK_dict['non_linear']['fluxes_1d']['particles_phi_potential'] != 0

            meta_dict = metadata

            # Handle diagnostics based on sim_type
            if sim_type in ['CGYRO', 'TGLF', 'GS2', 'GX']:
                Diag_dict = {}
            elif sim_type == 'GENE':
                print('='*60)
                # print('\n Working on diagnostics with user specified tspan .....\n')
                Diag_dict, manual_time_flag = get_diag_with_user_input(out_dir, suffix, manual_time_flag)
                print('='*60)

                if is_linear:
                    # Add omega info to Diag_dict for linear runs
                    omega_val = get_omega(out_dir, suffix)
                    Diag_dict['omega'] = {
                        'ky': omega_val[0],
                        'gamma': omega_val[1],
                        'omega': omega_val[2]
                    }

                for key, val in Diag_dict.items():
                    oid = gridfs_put_npArray(db, val, out_dir, key, sim_type)
                    Diag_dict[key] = oid ## Rewrite array with oid of stored file
                    if oid is not None: 
                        uploaded_ids[key] = oid

            # Combine dictionaries and upload
            run_data = {
                'Metadata': meta_dict,
                'Files': files_dict,
                'gyrokineticsIMAS': GK_dict,
                'Diagnostics': Diag_dict
            }
            
            main_record_oid = runs_coll.insert_one(run_data).inserted_id

            print(f'Files with suffix: {suffix} in folder {out_dir} uploaded successfully.')
            print('='*40)
            if verbose:
                print('A summary is generated as below:\n')
                print(run_data)

        except Exception as e1:
            print(e1)
            print(f"Skip suffix {suffix} in \n {out_dir} \n")
            global_vars.troubled_runs.append(out_dir + '##' + suffix)
            print('cleaning ......')
            fs = gridfs.GridFS(db)
            try:
                for key, _id in uploaded_ids.items():
                    fs.delete(_id)
                    print(f'{key}: {_id} deleted.')
            except Exception as e3:
                print(f"Error deleting files from gridfs with exception:\t {e3}")
                pass

        global_vars.reset_docs_keys(sim_type)


def upload_to_mongo(db, linear, metadata, out_dir, suffixes=None, run_shared=None,
                    large_files=False, verbose=False, manual_time_flag=False, global_vars=None, no_prompts=False, reupload_if_exists=False):
    """
    Wrapper function to upload simulation runs to MongoDB, handling both linear and nonlinear runs.

    Parameters:
    - db: Database connection object.
    - linear: Boolean indicating if the run is linear (True) or nonlinear (False).
    - metadata: Dictionary containing metadata for the run.
    - out_dir: Output directory containing simulation files.
    - suffixes: List of suffixes for files to upload. If None, determined automatically.
    - run_shared: List of shared files to upload (optional).
    - large_files: Boolean to handle large file uploads. Default: False.
    - verbose: Boolean to print detailed output. Default: True.
    - manual_time_flag: Boolean to handle user-specified time spans for diagnostics. Default: False.
    - global_vars: Object containing global variables for the upload process.
    - no_prompts: Autoupload with no prompts. Default = False
    - reupload_if_exists: Delete and reupload if existing folder name is present in DB. Default: False
    Returns:
    None
    """

    # Connect to the appropriate collection
    if metadata['CodeTag']['sim_type']=='TGLF': 
        runs_coll = db.LinearRuns
    else:
        runs_coll = db.LinearRuns if linear else db.NonlinRuns

    # Determine run type for printing
    run_type = 'linear' if linear else 'nonlinear'
    print(f'Upload {run_type} runs ******')

    # Check if folder is already uploaded
    if isUploaded(out_dir, runs_coll):
        print(f'Folder tag:\n {out_dir} \n exists in database')
        
        if no_prompts: 
            update='0' if reupload_if_exists else '1'
        else: 
            update = input(f'You can:\n 0: Delete and reupload folder? \n 1: Run an update (if you have updated files to add) \n Press any other keys to skip this folder.\n')
        
        if update == '0':
            # Delete and reupload
            print("Deleting {out_dir} and reuploading")
            remove_from_mongo(out_dir, db, runs_coll)
            upload_runs(db, metadata, out_dir, is_linear=linear, runs_coll=runs_coll, suffixes=suffixes, run_shared=run_shared,
                        large_files=large_files, verbose=verbose, manual_time_flag=manual_time_flag, global_vars=global_vars)
        elif update == '1':
            update_mongo(db, metadata, out_dir, runs_coll, linear)
        else:
            print(f'Run collection \'{out_dir}\' skipped.')
    else:
        print(f'Folder tag:\n{out_dir}\n not detected, creating new.\n')
        upload_runs(db, metadata, out_dir, is_linear=linear, runs_coll=runs_coll, suffixes=suffixes, run_shared=run_shared,
                    large_files=large_files, verbose=verbose, manual_time_flag=manual_time_flag, global_vars=global_vars)
