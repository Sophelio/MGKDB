#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: Austin Blackmon, Dongyang Kuang
"""


import sys

import numpy as np
import optparse as op
import matplotlib.pyplot as plt

from sys import path
from sys import exit
import os
import glob
import base64

from .fieldlib import fieldfile
from .ParIO import Parameters
# from .finite_differences import *
# import mgkdb.support.pydiag.utils.comm as comm
from .pydiag.utils import comm 
from .pydiag.data import datafiles

from .pydiag.utils.gkdb import GKDB_linear, GKDB_nonlin
from .putils.loader import Loader
from .data.data import Data
#from .putils.geom import Geometry
from .putils.run import Run
from .putils.simulation import Simulation
#from .putils.spatial_grid import SpatialGrid
#from .putils.vsp_grid import VspGrid
#from .data.base_file import createGENEfile

from .diagnostics.diag_flux_spectra import DiagFluxSpectra
from .diagnostics.diag_amplitude_spectra import DiagAmplitudeSpectra
from .diagnostics.diag_field_mom_snapshots import DiagFieldMomSnapshots
#=======================================================

def get_nspec(out_dir,suffix):
    #grab parameters dictionary from ParIO.py - Parameters()
    par = Parameters()
    par.Read_Pars(os.path.join(out_dir, 'parameters' + suffix))
    pars = par.pardict 
    
    #find 'n_spec' value in parameters dictionary
    nspec = pars['n_spec']
    
    return(nspec)
    
def get_nrg(out_dir, suffix):
    #modified from IFSedge/get_nrg.py
    par = Parameters()
    par.Read_Pars(os.path.join(out_dir, 'parameters' + suffix))
    pars = par.pardict 
    
    #initializations
    ncols=pars['nrgcols']
    time=np.empty(0,dtype='float')  # confused about the use of 0
    nrg0=np.empty((1,ncols))
    nrg1=np.empty((0,ncols),dtype='float')
    
    #grab 'n_spec' from 'parameters'
#    nspec = get_nspec(out_dir,suffix)
    nspec = pars['n_spec']
    
    #separate initializations for different 'n_spec' values
    if nspec<=2:
        nrg2=np.empty((0,ncols),dtype='float')
    if nspec<=3:
        nrg2=np.empty((0,ncols),dtype='float')
        nrg3=np.empty((0,ncols),dtype='float')
    
    
    #open 'nrg' file
    f=open(os.path.join(out_dir , 'nrg' + suffix),'r')
    nrg_in=f.read()

    #format 'nrg' file for reading
    nrg_in_lines=nrg_in.split('\n')
    for j in range(len(nrg_in_lines)):
        if nrg_in_lines[j] and j % (nspec+1) == 0:
            time=np.append(time,nrg_in_lines[j])
        elif nrg_in_lines[j] and j % (nspec+1) == 1:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrg1=np.append(nrg1,nrg0,axis=0)
        elif nspec>=2 and nrg_in_lines[j] and j % (nspec+1) ==2:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrg2=np.append(nrg2,nrg0,axis=0)
        elif nspec==3 and nrg_in_lines[j] and j % (nspec+1) ==3:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrg3=np.append(nrg3,nrg0,axis=0)

    #return 'time' and 'nrgx' arrays
    if nspec==1:
        return time,nrg1
    elif nspec==2:
        return time,nrg1,nrg2
    else:
        return time,nrg1,nrg2,nrg3

    
def get_parsed_params(filepath):
    par = Parameters()
    par.Read_Pars(filepath)
    pars = par.pardict 
    
    return pars

def get_suffixes(out_dir, sim_type):
    suffixes = []

    if sim_type=='GENE': ## scan all file with parameters 
        suffixes = [ os.path.basename(file).split('parameters')[-1] for file in glob.glob(os.path.join(out_dir,'parameters*')) if os.path.isfile(file)]
        if len(suffixes) > 1:
            suffixes.remove('')
            print(f"Found multiple suffixes: {suffixes}. Removing empty suffix.")
    elif sim_type in ['CGYRO','TGLF','GS2','GX']:  ## scan folders return as list 
        # suffixes = [os.path.basename(fldr) for fldr in os.listdir(out_dir) os.path.isdir(os.path.join(out_dir,fldr))]
        suffixes = [fldr for fldr in os.listdir(out_dir) if os.path.isdir(os.path.join(out_dir,fldr))]
    suffixes.sort()
    return suffixes

def get_gyrokinetics_from_dir(out_dir,user, linear):
    suffixes = get_suffixes(out_dir)
    numscan = len(suffixes)
    assert numscan>0, "files must have a suffix!"
    suffixes.sort()
    GK_list = []
    for suffix in suffixes:
        GK_list.append(get_gyrokinetics_from_run(out_dir, suffix, user, linear))
    
    return GK_list

def get_gyrokinetics_from_run(out_dir, suffix, user, linear, tspan = None):
    if isinstance(tspan, list):
        starttime = tspan[0]
        endtime = tspan[-1]
    else:
        starttime = -1
        endtime=-2
    
    current_pwd = os.getcwd()
    os.chdir(out_dir) # change dir or make it as an argument
    common = comm.CommonData(suffix, starttime, endtime)
    rundatafiles = datafiles.RunDataFiles(common)
    if linear:
        tst = GKDB_linear(common, rundatafiles, user)
    else:
        tst = GKDB_nonlin(common, rundatafiles, user)
    
    os.chdir(current_pwd) # get back to original dir
    
    return tst.gkdict

def get_diag_from_run(out_dir, suffix, t_span = None):
#    t_start = 0.0 # use start/end time in nrg files?
#    t_end = 100.0

    par0 = Parameters()
    par0.Read_Pars(os.path.join(out_dir , 'parameters' + suffix))
    pars0 = par0.pardict 

    if pars0['nonlinear'] == 'T' or pars0['nonlinear'] == True:
        nonlinear = True
    else:
        nonlinear = False

    if t_span is not None and len(t_span)==2:
        t_start = t_span[0]
        t_end = t_span[-1]
        
    else:
        print('Time span not speficied, searching it in NRG file\n')
        
        nspec = pars0['n_spec']
        if os.path.isfile(os.path.join(out_dir , 'nrg' + suffix) ):
            if nspec==1:
                tn,_ =get_nrg(out_dir, suffix)
                t_start = float(tn[0])
                t_end = float(tn[-1])
            elif nspec==2:
                tn,_,_ =get_nrg(out_dir, suffix)
                t_start = float(tn[0])
                t_end = float(tn[-1])
            elif nspec==3:
                tn,_,_,_ =get_nrg(out_dir, suffix)
                t_start = float(tn[0])
                t_end = float(tn[-1])
            else:
                exit("Not ready for nspec>3")
        else:
            t_start = 0.0
            t_end = 100.0
            print('NRG files not found for suffix {}. Using t_start = 0.0 and t_end = 100.0'.format(suffix))
        

    print('****** Diagnostics start at {} and end at {}. ******'.format(t_start, t_end))
    
    Diag_dict = {}
       
    #all is also included in the run object
    simulation=Simulation(out_dir, None, [suffix])
    run = simulation.runs[0]

    data = Data(run.parameters,suffix)
    
    #selected_diags= {'Flux Spectra':DiagFluxSpectra(avail_vars=data.av_vars, specnames=run.specnames),
    #                 'Amplitude Spectra':DiagAmplitudeSpectra(avail_vars=data.av_vars, specnames=run.specnames,parameters = run.parameters[0].pardict,spatialgrid = run.spatialgrid[0]),
    
    '''
    Creat a folder for saving plots if it does not exist
    '''
    if nonlinear:   #Only calculate spectra for nonlinear
        selected_diags = []
        diag_keys = []
        #selected_diags.append(DiagFieldMomSnapshots())
        #diag_keys.append('Field Mom Snapshots')
        selected_diags.append(DiagFluxSpectra())
        diag_keys.append('Flux Spectra')
        selected_diags.append(DiagAmplitudeSpectra())
        diag_keys.append('Amplitude Spectra')

        loader=Loader()
        step = 1
        loader.set_interval(selected_diags, simulation.data,run, t_start,t_end,step)
        for it, time in enumerate(loader.times):
            print(" time {}".format(time))
            for i_d, diag in enumerate(selected_diags):
                diag.execute(simulation.data, simulation.runs[0], loader.steps[it], loader.files[it],time)
        for i in range(len(diag_keys)): 
            print(diag_keys[i])
            Diag_dict[diag_keys[i]] = selected_diags[i].dict_to_mgkdb()

    '''
    Get appropriate time information for snapshots.  We want maximum of 30, minumum of 15
    Note: the sparse z-grid snapshots are only done for nonlinear
    '''
    ntime = len(simulation.data.av_times['mom'].times)
    if nonlinear:
        step = int(np.floor(ntime/30)) + 1
    else:
        step = 1
    '''
    Get snapshots and final fields and moments
    '''
    if nonlinear:
        print("Getting snapshots and final for fields and moments")
    else:
        print("Getting final for fields and moments")
    #print("Number of time points, step",ntime,step)
    selected_diags = []
    diag_keys = []
    selected_diags.append(DiagFieldMomSnapshots())
    diag_keys.append('Field Mom Snapshots')
 
    loader2=Loader()
    loader2.set_interval(selected_diags, simulation.data,run, t_start,t_end,step)
    for it, time in enumerate(loader2.times):
        if nonlinear or time == loader2.times[-1]:  #Only final time for linear
            print(" time {}".format(time))
            for i_d, diag in enumerate(selected_diags):
                diag.execute(simulation.data, simulation.runs[0], loader2.steps[it], loader2.files[it],time)
    for i in range(len(diag_keys)): 
        print(diag_keys[i])
        Diag_dict[diag_keys[i]] = selected_diags[i].dict_to_mgkdb()
    
    '''
    Grid
    '''
    Diag_dict['Grid'] = vars(run.spatialgrid)

    return Diag_dict
    
        
