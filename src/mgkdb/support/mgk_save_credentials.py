# -*- coding: utf-8 -*-
"""

@author: dykua, venkitesh ayyar

For downloading files from mgk_fusion in shell
"""

import sys
import os
from pymongo import MongoClient
from mgkdb.support.mgk_login import mgk_login


### Main 
def main():
    O2 = input("Please enter the server location, port, database name, username, password in order and separated by comma.\n").split(',')
    login = mgk_login(server= O2[0], port= O2[1], dbname=O2[2], user=O2[3], pwd = O2[4])
    O2_1 = input("Plese enter a target path to save the file in .pkl format. For example : ../db_credentials.pkl  \n")
    login.save(os.path.abspath(O2_1))

    client, database = login.connect()

    print("Successfully connected to database")
    print("Following collections exist in database",database.list_collection_names())
    print("For upload and download, you can use the credential file with the -A argument as: -A <path_to_credentials.pkl>")


if __name__=="__main__":
    main()

