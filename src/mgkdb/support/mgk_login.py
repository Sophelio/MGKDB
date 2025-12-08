# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 12:49:43 2019

@author: dykua

Functions for handling credentials
"""

import pymongo
import pickle
import os
import urllib.parse

class mgk_login(object):
    def __init__(self, server='mongodb03.nersc.gov', port='27017', dbname='mgk_fusion', 
                 user='user', pwd = 'pwd'):

        self.login = {'server': server.strip(),
                      'port': port.strip(),
                      'dbname': dbname.strip(),
                      'user': user.strip(),
                      'pwd': pwd.strip()
                      }
        
    def from_saved(self, file_path):
        with open(file_path.strip(), 'rb') as pkl:
            info = pickle.load(pkl)
        
        self.login = info
        
        
    def save(self, file_path):
        
        with open(file_path, 'wb') as pkl:
            pickle.dump(self.login, pkl, protocol = pickle.HIGHEST_PROTOCOL)
        print('Login info saved to {}'.format(file_path) )    
            
    def update(self, dict_to_update):
        for key, val in dict_to_update.items():
            if key in self.login:
                self.login[key] = val 
    #   self.login.update(dict_to_update)
                
    def connect(self):
        ## Old way for pymongo 3.*
        # database = MongoClient(self.login['server'].strip(), int(self.login['port']) )[self.login['dbname'].strip()]
        # database.authenticate(self.login['user'].strip(), self.login['pwd'].strip())

        # URL encode the password to handle special characters safely
        encoded_pwd = urllib.parse.quote_plus(self.login['pwd'])
        
        # Create connection string with encoded password
        mongodb_url = 'mongodb://{user}:{pwd}@{server}:{port}/{dbname}'.format(
            user=self.login['user'],
            pwd=encoded_pwd, 
            server=self.login['server'],
            port=self.login['port'],
            dbname=self.login['dbname']
        )
        
        client = pymongo.MongoClient(mongodb_url, directConnection=True)
        database = client[self.login['dbname']]
        return client, database
        
def f_login_dbase(authenticate):

    info = authenticate

    if info is None:
        O1 = input("You did not enter credentials for accessing the database.\n You can \n 0: Enter it manually. \n 1: Enter the full path of the saved .pkl file\n")
        if O1 == '0':
            O2 = input("Please enter the server location, port, database name, username, password in order and separated by comma.\n").split(',')
            login = mgk_login(server= O2[0], port= O2[1], dbname=O2[2], user=O2[3], pwd = O2[4])
            O2_1 = input("You can save it by entering a target path, press ENTER if you choose not to save it\n")
            if len(O2_1)>1:
                login.save(os.path.abspath(O2_1) )
            else:
                print('Info not saved!')
                pass
        elif O1 == '1':
            O2= input("Please enter the target path\n")
            login = mgk_login()
            login.from_saved(os.path.abspath(O2))
        
        else:
            exit("Invalid input. Abort")
                
    else:
        login = mgk_login()
        try:
            login.from_saved(os.path.abspath(info))
        except OSError:
            exit("The specified credential file is not found!")

    return login
 