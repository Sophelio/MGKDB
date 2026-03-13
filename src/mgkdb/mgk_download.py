# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 14:39:26 2020

@author: dykua, venkitesh ayyar

For downloading files from mgk_fusion in shell
"""

import sys
import os
import argparse
from pymongo import MongoClient
from bson.objectid import ObjectId
import gridfs

from mgkdb.support.mgk_file_handling import get_oid_from_query, Str2Query, download_dir_by_name, download_file_by_path, download_file_by_id, download_runs_by_id
from mgkdb.support.mgk_login import mgk_login,f_login_dbase

def f_parse_args():
    #==========================================================
    # argument parser
    #==========================================================
    parser = argparse.ArgumentParser(description='Process input for downloading files')

    parser.add_argument('-Q', '--query', default= None,
        help='MongoDB query as JSON. E.g. {"Metadata.DBtag.run_collection_name": "path"} or {"Metadata.ScenarioTag.scenario_runid": {"$regex": "193870.3400"}}')
    parser.add_argument('-T', '--target', default= None,help='run collection_name, i.e. gene output folder path')
    parser.add_argument('-F', '--file', default = None, help='filename to be downloaded if any')
    parser.add_argument('-C', '--collection', choices=['linear','nonlinear','files'], default='linear', type=str, help='collection name in the database')
    parser.add_argument('-OID', '--objectID', default = None, help = 'Object ID in the database')
    parser.add_argument('-A', '--authenticate', default = None, help='locally saved login info, a .pkl file')
    parser.add_argument('-D', '--destination', default = '', help = 'directory where files are downloaded to.')
    parser.add_argument('-S', '--saveas', default = 'specific_file', help = 'Name to save the file as')

    return parser.parse_args()

### Main 
def main_download(target, file, objectID, destination, saveas, query, authenticate, collection):

    OID = objectID
    op_fname = saveas


    ### Connect to database 
    login = f_login_dbase(authenticate)
    client, database = login.connect()
    with client:
        ## Dict to convert from argument to collection name in database
        collection_dict={'linear':'LinearRuns','nonlinear':'NonlinRuns','files':'fs.files'}
        collection_name =  getattr(database,collection_dict[collection])

        if query:
            print("working on query: {} ......".format(query))
            found = get_oid_from_query(database, collection_name, Str2Query(query))
            for oid in found:
                download_runs_by_id(database, collection_name, oid, destination)

        elif file:
            download_file_by_path(database, file, destination, revision=-1, session=None)   
            
        elif OID:
            if collection=='files': 
                download_file_by_id(database, ObjectId(OID), destination, op_fname, session = None)
            elif collection in ['linear','nonlinear']:
                download_runs_by_id(database, collection_name, ObjectId(OID), destination)
            else : 
                print("Invalid option for collection for OID",collection)
                raise SystemError
        elif target:
            download_dir_by_name(database, collection_name, target, destination)


def main():

    ### Parse arguments 
    args = f_parse_args()
    # print(args)

    main_download(**vars(args))


if __name__=="__main__":
    main()

'''
python src/mgkdb/mgk_download.py -Q '{"Metadata.ScenarioTag.scenario_runid": {"$regex": "193870.3400"}}' -A ../db_ro_cred.pkl -C nonlinear -D ./out

fileid = 6933167a9cadc1844e312496
python src/mgkdb/mgk_download.py -OID 6933167a9cadc1844e312496 -A db_credentials.pkl -C files -D ./out

'''

