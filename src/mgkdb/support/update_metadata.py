# Code to download all non-linear runs in database 

import sys
import os
import argparse
from sys import exit

import gridfs
import json 
from bson.objectid import ObjectId

from mgkdb.support.mgk_login import mgk_login,f_login_dbase

# Run this as : 
# python update_meta_code.py -A <fname.pkl> -C linear -m 3

def f_parse_args():
    parser = argparse.ArgumentParser(description='Update Metadata entries. 3 modes.  1: Append to publication list.\n2: Append to comments.\n 3: Update any specific entry. \n Modes 1 and 2 add entered values to existing entry.\nUse mode=3 with caution as you are rewriting previous entry.')

    parser.add_argument('-C', '--collection', choices=['linear','nonlinear'], default='linear', type=str, help='collection name in the database')
    parser.add_argument('-OID', '--objectID', default = '678a700bd978ec23d2ebfd18', help = 'Object ID in the database')
    parser.add_argument('-m', '--mode', type=int, choices=[1,2,3], default = 3, help = 'Choose mode of operation for updating Metadata. 1: Append to publication list.\n2: Append to comments.\n 3: Update any specific entry.')
    parser.add_argument('-A', '--authenticate', default = None, help='locally saved login info, a .pkl file')

    return parser.parse_args()

def f_update_metadata(data, key_name, template_d):
    '''
    Update the existing metadata key entry 
    key_name : key to update 
    template_d : dictionary containing template
    '''
    if isinstance(data, dict):
        new_dict ={}
        for key,value in data.items():
            new_val = f_update_metadata(value, key, template_d)
            new_dict[key] = new_val

        return new_dict 
    
    elif isinstance(data, list):
        ## Create dictionary using last value 
        new_val = f_update_metadata(data[-1],key_name, template_d)
        
        return new_val

    else:
        print("The existing entry for this key is:\t%s"%(data))
        ## if entry is string, just replace
        new_value = input('Please enter the new entry you want for %s. If you don\'t want to proceed, enter : none \n'%(key_name))

        return new_value

if __name__=="__main__":
    
    ### Parse arguments 
    args = f_parse_args()
    Args = vars(args)
    
    authenticate = Args['authenticate']
    collection_ip = Args['collection']
    
    ### Connect to database 
    login = f_login_dbase(authenticate)
    client, database = login.connect()

    #### Dict to convert from argument to collection name in database
    collection_dict={'linear':'LinearRuns','nonlinear':'NonlinRuns','files':'fs.files'}
    collection =  getattr(database,collection_dict[collection_ip])
    
    # Get Object IDs 
    all_ids = [r['_id'] for r in collection.find({},{'id':1})]
    
    oid = ObjectId(args.objectID)
    assert oid in all_ids,"Given oid entry %s doesn't exist in database"%(oid)
    
    ### Extract document for this oid
    document = collection.find_one({"_id":oid},{'Metadata':1,'_id':0})

    ### Check user credential
    uname = login.login['user']
    input_usr = document['Metadata']['DBtag']['user']

    if uname!=input_usr:
        cont = input(f"Data was input by another use {input_usr}. Do you wish to continue ? y or n \n")
        if not cont=='y': 
            raise SystemExit
        
    keys = document.get('Metadata').keys()

    ### Options : Append to publications list, Append to comments, Update any specific entry
    if args.mode==1: 
        user_ip = input('Enter the doi strings to append, separated by commas\n')
        entry_add = user_ip.split(',')

        ## Existing entry
        entry_val = document['Metadata']['Publications']['doi']
        if not entry_val: entry_val = [] 
        
        new_pub_list = entry_val + entry_add ## Join the two lists 

        assert isinstance(new_pub_list,list),f"Entry to add must be a list"

        fltr = {"_id": oid}
        update = {"$set": {"Metadata.Publications.doi": new_pub_list}}
        result = collection.update_one(fltr, update)
        print("Appended publication record")

    elif args.mode==2: # Option to append to comment string
        old_comment = document['Metadata']['DBtag']['comments']
        # assert isinstance(old_comment,str),f"Existing entry {old_comment} is not a string" 
        user_ip = input(f'Enter the string to append to current entry. Current entry is: \n {old_comment} \n')
        
        fltr = {"_id": oid}
        # Using $set to add appended string to comments
        update = {"$set": {"Metadata.DBtag.comments": old_comment+'.\n'+user_ip}}
        result = collection.update_one(fltr, update)
        print("Appended comment")
        
    elif args.mode==3: 
        top_key = input(f'Please enter the key in Metadata that you want to update. Allowed options are: \n {list(keys)}\n')
        assert top_key in keys,f"Invalid input key {top_key}"

        sub_keys = document.get('Metadata').get(top_key).keys()
        key_name = input(f'The key {top_key} has the following subkeys {list(sub_keys)}. \nPlease select which subkey you want to modify\n')
        assert key_name in sub_keys,f"Invalid input key {key_name}"

        ans = document.get('Metadata').get(top_key).get(key_name)
        print("The existing entry for this key is:\t%s"%(ans))
        print("You will be resetting the entire value for this key. Please use caution")          
        new_value = input('Please enter the new entry you want for %s. If you don\'t want to proceed, enter : none \n'%(key_name))
        confirm = input(f'Confirm changing entry to: \n{new_value}? \n Enter Y or N\n').strip().upper()

        ## Some entries need to be in a list 
        if top_key in ['Publications']:   new_value = [new_value]

        if (confirm!='Y' or (new_value == 'none')):
            print('Didn not receive confirmation. Aborting update')
            raise SystemExit
            
        fltr = {"_id": oid}
        update = {"$set": {f"Metadata.{top_key}.{key_name}": new_value}}
        result = collection.update_one(fltr, update)

        document = collection.find_one({"_id":oid},{'Metadata':1,'_id':0})
        print('Updated entry: %s '%(document.get('Metadata').get(top_key).get(key_name)))    


