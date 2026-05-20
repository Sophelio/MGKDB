# **MGKDB**

## Installation
The easiest way to install this package is after installing MongoDB Compass and mongosh, by running pip install MGKDB-fusion in the command line

After running pip you will need to input database credentials using mgk_save_credentials. For testing, input localhost, 27017, mgk_fusion, <username>, <password>

## Setting up the environment
A convenient method to access the database tools is using conda environments.
You can find the instructions to build the conda environment in the wiki [here](https://github.com/Sapientai/MGKDB/wiki/Setting-up-the-environment).

## Test DB  locally
Please refer to this [Wiki](https://github.com/Sapientai/MGKDB/wiki/Local-MGKDB) for step-by-step instructions on building MGKDB locally 

## Test DB on NERSC

To access MGKDB at NERSC, you will need two sets of access credentials: 
1. NERSC : Access to the NERSC computing infrastructure. NERSC account request can be placed using the instructions given [here](https://docs.nersc.gov/accounts/).
2. MGKDB at NERSC: Access to the MGKDB database at NERSC. Please email michoski@oden.utexas.edu to request access.

Users can interact with the database in 3 ways: 
### 1. Using Command Line tools
This repository provides command line tools can be used to download data from the database and upload data to the database.

### 2. Using MongoDB Compass GUI

### 3. Directly using the terminal
Another method to interact with the database is directly from the terminal
1. From a terminal, `ssh` to Perlmutter using \
   ```ssh -l <nersc_username> perlmutter.nersc.gov```
3. Connect to the database using : 
```
module load mongodb/4.0.28
mongo -u <db_username> -p <db_password> mongodb03.nersc.gov/mgk_fusion
```
3. Now, one can use `db` commands.

For the former two approaches, please refer to this [Wiki](https://github.com/Sapientai/MGKDB/wiki/MGKDB-at-NERSC) for step-by-step instructions.

### Caution:
* If you have Read/Write access to the database, be careful while using mongodb compass, you may accidently modify the database.  
