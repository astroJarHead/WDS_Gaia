### build_and_write_wds_table.py

__author__ = "Daphne Zakarian, Robert T. Zavala, Jr."
__version__ = "1.0"
__license__ = " This is a work of the U.S. Government and is not subject to copyright \
 protection in the United States. The U.S. Government retains all exclusive \
 rights to use, duplicate, distribute, disclose, or release this software \
 and all associated documentation and data.Use of appropriate byline or \
 credit is requested. "

### Copyright and License 
#
# This is a work of the U.S. Government and is not subject to copyright 
# protection in the United States. The U.S. Government retains all exclusive 
# rights to use, duplicate, distribute, disclose, or release this software 
# and all associated documentation and data.Use of appropriate byline or 
# credit is requested. 
#
###
#
# STEP 1 OF WDS TO GAIA DR3 MATCHING 
#
# Take downloaded ASCII copy of the WDS catalog and:
#
#    1) Convert the WDS -> astropy Table
#    2) Using Primary RA, DEC and component rho-2, theta-2 calculate 
#       secondary J2000.0 RA, DEC
#    3) Apply proper motions and get positions at Gaia DR3 epoch
#    4) Write out the Table as CSV, ECSV and VOtables for next steps.
# 
# To use:
#
# At python prompt with conda env gaiaWds activated:
# 
# >>> from build_and_write_wds_table import *
# >>> gaiaTable()
# 
# On a VM with 156,618 lines in a WDS copy the code took
# (output copied and pasted below on next four lines):
# 
#     Tables built and written as CSV, ECSV and VOTables.
#     Elapsed time for building and writing tables:  72.32840490341187  seconds.
#     Elapsed CPU for building and writing tables:  67.851002069  seconds.
#     Finished, thank you for using gaiaTable!


import pandas as pd
from astropy.table import Table, Column, MaskedColumn
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord, FK5
from astropy.io import ascii
import math
import numpy as np
import os

# to time task completion only
import time

def gaiaTable():
    
    """
    Original Created on Thu 09 Jun 11:56:21 MDT 2022
    @author: Daphne Zakarian, Bob Zavala

    Read in the WDS using pandas, turn into an astropy table, find secondary 
    RA and DEC coordinates using offset from primary. Next, apply proper  
    proper motions to primary and secondary to Gaia DR3 2016.0 epoch. 
    At the end, about tables in csv, ecsv and VOtable formats for use in 
    matching against Gaia DR3. 

    Final table will be used to query Gaia -- most importantly, we have primary 
    and secondary stars' coordinates in degrees at the DR3 epoch. 

    Changes began on Thu 30 May 09:08:00 MDT 2024 
    by Bob Zavala to include Discoverer and Components in order to 
    pass unique identifier information along. Coordinates and WDS ID alone
    do not uniquely identify the WDS entries. 

    Accounts for small number of WDS entries with large proper motions 
    expressed in time units of kilo-years. Retains number of measures in 
    the WDS and the WDS spectral type. 

    Name of file changed to explicitly state the code builds and writes a WDS 
    table. 

    Parameters
    ----------
    None:

    Requires
    --------
    1) conda environment gaiaWds activated
    2) ASCII file wdsweb_summ2.txt in current directory
       This is a copy of the WDS catalog

    Output
    ------
    1) Output files with prefix wdstab_new in csv, ecsv and VOtable 
       formats written in current directory

    To run
    ------
    1) Start python
    2) import build_and_write_wds_table
       build_and_write_wds_table.gaiaTable()
    OR:
    1) Start python
    2) from  build_and_write_wds_table
    3) gaiaTable()

    """

    this_path = os.getcwd()

    # Time the code
    startTime    = time.time()
    startCpuTime = time.process_time()

    print('\n')
    print('Welcome to the gaiaTable module! This should take a few (<5) minutes.')
    print('Starting to build the WDS tables for the Gaia queries.\n')
    print('**********')
    print('  Any ErfaWarnings from function pmsafe are normal and appear because the WDS lacks distances.')
    print('**********')

    # Read in Table: fixed width file to pandas to astropy 

    # manually list WDS txt file column widths to read in 
    # the formatted data. A change made by Daphane is to 
    # separate the DM column into two parts and the WDS 
    # 2000 arcsecond coordinates into RA and DEC portions. 
    ## IMPORTANT: colspecs are HALF-OPEN intervals: [,)
    ## SO THE RIGHT SIDE IS NOT INCLUSIVE!!!
    columns = ((0,10),(10,17),(17,23),(23,28),(28,32),(33,37),
            (37,42),(42,45),(45,51),(51,57),(57,64),(64,68),
            (69,80),(80,84),(84,88),(88,93),(93,98),(98,101),
            (101,106),(107,111),(112,121),(121,130))

    oldnames= ['WDS_Identifier', 'Discovr', 'Comp', 'Epoch-1', 'Epoch-2', '#', 
            'Theta-1', 'Theta-2', 'Rho-1', 'Rho-2', 'Mag-pri', 'Mag-sec',
            'SpecType','PMpri-RA', 'PMpri-DEC', 'PMsec-RA', 'PMsec-DEC', 'DM', 
            'Desig', 'Note', 'WDS_RA', 'WDS_DEC']

    # read fixed width file into pandas 
        # easier to work with this fixed-width-formatted (fwf) file 
        # in pandas than astropy
    wdspd = pd.read_fwf(this_path+"/wdsweb_summ2.txt",
                        colspecs=columns,header=None,names=oldnames,skiprows=3)

    # pandas table -> astropy table
    wdstab = Table.from_pandas(wdspd)

    print('\n')
    print('An astropy table has been created from the WDS data.')
    print('Converting precise WDS coordinates, and some other tasks. This takes a bit.')

    # Establish format for RA and DEC

    # make new column to store updated RA and DEC -- is there a better way to
        # make column dtype longer so I don't have to manually type 15 spaces???
    wdstab['RAprideg']=0.0
    wdstab['DECprideg']=0.0
    wdstab['RAhms']='               '
    wdstab['DECdms']='               '
    

    # assign units to the new columns for coordinates in degrees
    wdstab['RAprideg'].unit=u.deg
    wdstab['DECprideg'].unit=u.deg


    # loop through the ra and dec of each row and make new columns for coords in 2 formats:
        # hms/dms format and degrees
    # setup counter for no coordinates
    noPrecise = 0
    for row in wdstab:
        try:
            ra1 = row['WDS_RA']
            dec1 = row['WDS_DEC']
            rastr = ra1[:2] + "h" + ra1[2:4]+"m"+ra1[4:]+"s"
            decstr = dec1[:3] + "d" + dec1[3:5]+"m"+dec1[5:]+'s'
            
            
            # put the strings of ra and dec into the table
            row['RAhms']= rastr
            row['DECdms']= decstr
            
            # put the coordinates into degree form - should work with GAIA
            coo = ICRS(rastr,decstr)
            # put the new coordinates into a column with a float dtype
            row['RAprideg']= coo.ra.deg
            row['DECprideg']=coo.dec.deg
        # skip rows that don't have coordinates
        except ValueError:
            noPrecise+=1
            pass

    print('\n')
    print('I found ',noPrecise,' WDS entries without precise (arcsec) coordinates.\n')
        

    """ FIND J2000.0 COORDINATE OF SECONDARY  """
    # WDS contains RA and DEC of primary and relative offset to companions
    
    # make new columns for secondary coords
    wdstab['RAsecdeg']=0.0
    wdstab['DECsecdeg']=0.0
    wdstab['RAsecdeg'].unit=u.deg
    wdstab['DECsecdeg'].unit=u.deg

    # The loop that the vector algorithm below replaced took 6.5 minutes 
    # on medea to run. 

    #startTimeVec = time.time()
    #cpuTimeVec_st = time.process_time()
    pricoord=SkyCoord(wdstab[0:]['RAprideg'],wdstab[0:]['DECprideg'], frame='icrs')
    angle = wdstab[0:]['Theta-2']*u.deg
    sep = wdstab[0:]['Rho-2']*u.arcsec
    seccoord = pricoord.directional_offset_by(angle, sep)

    wdstab['RAsecdeg'] = seccoord.ra.deg
    wdstab['DECsecdeg'] = seccoord.dec.deg

    print('\n')
    print('J2000.0 Coordinates (RA/DEC) of secondaries found via offsets from primaries.\n')
    #cpuTimeVec_end = time.process_time()
    #endTimeVec = time.time()
    #elapVecTime = endTimeVec - startTimeVec
    #elapCpuVec = cpuTimeVec_end - cpuTimeVec_st
    #print('Elapsed time for seccoord using vector method: ',elapVecTime,' seconds.')
    #print('Elapsed CPU for seccoord using vector method: ',elapCpuVec,' seconds.')

    # Insert some missing units into the table
    # The for loop seems necessary and takes little time
    col_need_units = ['Epoch-1','Epoch-2','Theta-1','Theta-2','Rho-1','Rho-2','Mag-pri','Mag-sec',
                    'PMpri-RA','PMpri-DEC','PMsec-RA','PMsec-DEC']
    units_for_cols = [u.yr,u.yr,u.deg,u.deg,u.arcsec,u.arcsec,u.mag,u.mag,u.mas/u.yr,u.mas/u.yr,
                    u.mas/u.yr,u.mas/u.yr]
    counter = 3
    for a_col in col_need_units:
        # print(a_col)
        wdstab[a_col].unit = units_for_cols[counter-3]
        counter +=1

    # Change the '--' to 0.0 in the proper motion columns 
    # in order to use zero proper motions in PM calculations
    wdstab['PMpri-RA'].fill_value = 0.0
    wdstab['PMpri-RA'] = wdstab['PMpri-RA'].filled(0.0)
    wdstab['PMpri-DEC'].fill_value = 0.0
    wdstab['PMpri-DEC'] = wdstab['PMpri-DEC'].filled(0.0)
    wdstab['PMsec-RA'].fill_value=0.0
    wdstab['PMsec-RA'] = wdstab['PMsec-RA'].filled(0.0)
    wdstab['PMsec-DEC'].fill_value=0.0
    wdstab['PMsec-DEC'] = wdstab['PMsec-DEC'].filled(0.0)

    """ACCOUNT FOR ENTRIES WITH NOTE P: PM IN KILO-YEARS"""
    # Set the mask for replacing the '--' entries in the Note column, and replace the '--' 
    # with ''
    wdstab['Note'].mask = (wdstab['Note'] == '--') 
    note_col_filled = wdstab['Note'].filled('')

    # Find the indices (mask) where the Note = P as this indicates proper motions have time units of 
    # per 1000 years. Use a temporary column pm_kyr_mask (proper motion kilo-year mask)
    pm_kyr_mask = np.char.find(note_col_filled, 'P') >= 0

    # Apply the resulting Boolean mask to the original table
    p_flagged_rows = wdstab[pm_kyr_mask]

    # Convert proper motion in wdstab with a 'P' note by multiplying by 10.0. 
    # These are fast movers! 
    # Define the columns that require the unit conversion.
    columns_to_convert = ['PMpri-RA', 'PMpri-DEC', 'PMsec-RA', 'PMsec-DEC']

    # Loop through the column names and apply the vectorized division.
    # The boolean mask ensures we only update the rows where 'Note' contains 'P'
    for col_name in columns_to_convert:
        wdstab[col_name][pm_kyr_mask] *= 10.0
    
    """ ACCOUNT FOR PROPER MOTION """
    # Need to add units for PMpri-RA,PMpri-DEC,PMsec-RA,PMsec-DEC as 
    # these did not make it into the table.

    # Change units of pm to deg/year
    wdstab['PMpri-RAdeg'] = wdstab['PMpri-RA'].to(u.deg/u.year)
    wdstab['PMpri-DECdeg'] = wdstab['PMpri-DEC'].to(u.deg/u.year)
    wdstab['PMsec-RAdeg'] = wdstab['PMsec-RA'].to(u.deg/u.year)
    wdstab['PMsec-DECdeg'] = wdstab['PMsec-DEC'].to(u.deg/u.year)

    # Apply proper motions to primaries -> 2016.0
    priRA    = wdstab['RAprideg']
    priDEC   = wdstab['DECprideg']
    pmpriRA  = wdstab['PMpri-RAdeg']
    pmpriDEC = wdstab['PMpri-DECdeg']

    coord = SkyCoord(priRA, priDEC, unit='deg', 
                    pm_ra_cosdec = pmpriRA, pm_dec = pmpriDEC)
    newcoord = coord.apply_space_motion(dt=16*u.yr)

    # While newcoord above IS a SkyCoord (astropy) object the conversion
    # to degrees from the hms default transforms newcoord -> numpy array.
    # I need to convert this numpy array back to an astropy column 
    # AND THEN insert the proper-motion adjusted coordinates into the 
    # table. Arrgh. 
    # The other thing I could not do was create the RApri-prepped 
    # and DECpri-prepped columns, prefill with zeroes and insert the new 
    # newcoord columns. 
    new_prira_col = Column(newcoord[0:].ra.deg)
    new_pridec_col = Column(newcoord[0:].dec.deg)
    # Here is how I finally entered the PM applied coordinates. 
    # Except I still need to fix the 'nan' proper motion prepped coordinates
    wdstab['RApri-prepped'] = new_prira_col
    wdstab['DECpri-prepped'] = new_pridec_col

    # Apply proper motions to secondaries -> 2016.0
    secRA    = wdstab[0:]['RAsecdeg']
    secDEC   = wdstab[0:]['DECsecdeg']
    pmsecRA  = wdstab[0:]['PMsec-RAdeg']
    pmsecDEC = wdstab[0:]['PMsec-DECdeg']

    coord = SkyCoord(secRA, secDEC, unit='deg', 
                    pm_ra_cosdec = pmsecRA, pm_dec = pmsecDEC)
    newcoord = coord.apply_space_motion(dt=16*u.yr)

    new_secra_col = Column(newcoord[0:].ra.deg)
    new_secdec_col = Column(newcoord[0:].dec.deg)

    wdstab['RAsec-prepped']  = new_secra_col
    wdstab['DECsec-prepped'] = new_secdec_col

    # Insert more missing units into the table
    col_need_units = ['RApri-prepped','DECpri-prepped','RAsec-prepped','DECsec-prepped','RAsecdeg','DECsecdeg']

    the_unit = u.deg

    for a_col in col_need_units:
        #print(a_col)
        #print(wdstab[a_col].unit)
        wdstab[a_col].unit = the_unit
        #print(wdstab[a_col].unit)


    print('\n')
    print('Proper motions applied to 2016.0 for querying against Gaia DR3.')

    # Save this for the END. The DM identifications are not available 
    # for all entries and are not a helpful alias.  
    wdstab.remove_columns(['DM', 'Desig'])

    # Avoid astropy generating a warning on a column with the name '#'
    # and use a more explicit name of 'Num_Meas'
    wdstab.rename_column('#','Num_meas')

    # SAVE TABLE AS ECSV and CSV and as a VOTable
    # Two different ways of coding the ascii-write are 
    # shown as an example for how to write these and one is a 
    # little more concise.
    ascii.write(wdstab, this_path+'/wdstab_new.ecsv', format='ecsv', overwrite=True)
    wdstab.write(this_path+'/wdstab_new.csv', overwrite=True)
    wdstab.write(this_path+'/wdstab_new.vot', format = 'votable', overwrite=True)

    print('\n')
    print('Tables built and written as CSV, ECSV and VOTables.')

    endTime    = time.time()
    endCpuTime = time.process_time()

    totalElapTime = endTime - startTime
    totalCpuTime  = endCpuTime - startCpuTime

    print('Elapsed time for building and writing tables: ',totalElapTime,' seconds.')
    print('Elapsed CPU for building and writing tables: ',totalCpuTime,' seconds.')
    print('Finished, thank you for using gaiaTable!')

    print('\n')









