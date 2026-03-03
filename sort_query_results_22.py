from astropy.io import ascii
from astropy.table import vstack, Table, unique, QTable, Column
from astropy.coordinates import SkyCoord 
import astropy.units as u
from astropy import table, log
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Distance, Angle
from astropy.time import Time
import pickle

import os
from IPython.display import display
from multiprocessing import set_start_method

from math import sqrt, cos, radians

# numpy needed as arrays required to add Quantity
# objects with add_column method for astropy Class Table
import numpy as np

# to time task completion only
import time

# FOR DEBUG/TESTING
import pdb

def extendFlagRows(nonassoc_TriplePlus,group_number,stars,flag):
    
    """
    Extends the length of the rows for the nonassociated triple systems
    that records the group number and 'stars' id or component of the 
    triple system the reason for the flagging. This information will 
    eventually be used when writing the output files for the 
    nonassociated triple systems.

    INPUT:

        nonassoc_TriplePlus: existing rows of nonassociated 3+ systems

        group_number: the group number to add for a newly nonassociated
        triple+ system

        stars: the 'stars' component id to add for a newly nonassociated
        triple+ system

        flag: numerical integer code of 1,2 or 3 encoding the reason 
        the star is flagged as nonassociated

    OUTPUT:

        extendedTripleFlagRows: the input with the row added using the 
        numpy vstack method
    
    """
    
    # Add a new row to the nonsssoc_TriplePlus entries. 
    
    # Compose the new row as a single, one row array 
    new_row = np.array([group_number, stars, flag])

    # Using the vstack method of np.array extend the input array 
    # by new_row
    extendedTripleFlagRows = np.vstack((nonassoc_TriplePlus, new_row))

    return extendedTripleFlagRows
    
def extendAssocRows(assoc_TriplePlus,group_number,stars):
    
    """
    Extends the length of the rows for the associated triple systems
    that records the group number and 'stars' id or component of the 
    triple system that is associated. This information will eventually 
    be used when writing the output files for the associated triple 
    systems.

    INPUT:

        assoc_TriplePlus: existing rows of associated 3+ systems

        group_number: the group number to add for a newly associated
        triple+ system

        stars: the 'stars' component id to add for a newly associated
        triple+ system

    OUTPUT:

        extendedTripleAssocRows: the input with the row added using the 
        numpy vstack method
        
    """

    # Add a new row to the nonsssoc_TriplePlus entries. 

    # Extend the rows of the associated 3+ systems for eventual
    # writing out as a file. 

    # Compose the new row as a single, one row array 
    new_row = np.array([group_number, stars])

    # Using the vstack method of np.array extend the input array 
    # by new_row
    extendedTripleAssocRows = np.vstack((assoc_TriplePlus, new_row))

    return extendedTripleAssocRows

def writeTriples(assoc_TriplePlus,nonassoc_TriplePlus,WDS_grouped,save_path):

    """
    Writes the output files for the associated and nonassociated triple+
    WDS systems. Loops are used to prepare the data for writing as vector 
    methods could not after much effort be made to work. This seems to be 
    due to the two indices of groups and stars in the grouped data. 

    Six about files are created: three each for the associated and 
    nonassociated systems using csv, ecsv and votable formats

    INPUT:

        assoc_TriplePlus: indices identifying group number and 
        components for the associated systems

        nonassoc_TriplesPlus: indices identifying group number and 
        components for the nonassociated systems and an integer for 
        the flag

        WDS_grouped: the WDS catalog in the grouped format

        save_path: the path where the data are to written

    OUTPUT:

        There are six output files written to save_path.
        For the associated triple+ systems:

        WDS_assoc_triples.csv
        WDS_assoc_triples.ecsv
        WDS_assoc_triples.vot

        For the nonassocited triple+ systems:
        
        WDS_nonassoc_triples.csv
        WDS_nonassoc_triples.ecsv
        WDS_nonassoc_triples.vot
            
    """

    # Prepare the indices for use in writing via a vector-like
    # process
    assoc_groupNum = assoc_TriplePlus[:,0]
    assoc_stars = assoc_TriplePlus[:,1]
     
    # Check that assoc_groupNum and assoc_stars have the same lengths
    if len(assoc_groupNum) == len(assoc_stars):
        pass
    else:
        # return 2 for this error
        return 2

    nonassoc_groupNum = nonassoc_TriplePlus[:,0]
    nonassoc_stars = nonassoc_TriplePlus[:,1]
    nonassoc_flags = nonassoc_TriplePlus[:,2]

    if len(nonassoc_groupNum) == len(nonassoc_stars):
        if len(nonassoc_groupNum) == len(nonassoc_flags):
            pass
    else:
        # return 3 for this error
        return 3

    
    # Define the output file names
    # associated filenams
    assoc_csv_filename = save_path+'/WDS_assoc_triples.csv'
    assoc_ecsv_filename = save_path+'/WDS_assoc_triples.ecsv'
    assoc_votable_filename = save_path+'/WDS_assoc_triples.vot'
    # nonassociated filenames
    nonassoc_csv_filename = save_path+'/WDS_nonassoc_triples.csv'
    nonassoc_ecsv_filename = save_path+'/WDS_nonassoc_triples.ecsv'
    nonassoc_votable_filename = save_path+'/WDS_nonassoc_triples.vot'
    
    # Test printing to screen the first 20 associated groupNum and star 
    # indices
    print(" The first 20 associated group numbers and stars are: \n")
    for i in range(0,21):
        print(WDS_grouped.groups[assoc_groupNum[i]][assoc_stars[i]]['WDS_Identifier','Discovr','Comp'])

    print(" \n")
    # Write out the results for the 3+ star WDS entries
    # Extract and write out the associated group data

    # Attempt at a vector method in next line produces repetitions of 
    # the first few entries. 
    # assoc_triple_data = WDS_grouped.groups[assoc_groupNum][assoc_stars]
    # nonassoc_triple_data = WDS_grouped.groups[nonassoc_groupNum][nonassoc_stars]
    
    # Try a loop over indices of groupNum and stars for associated stars
    # and: extract assoc groups and stars via a list. 
    assoc_triple_list = []

    # Loop over all groups and associated stars, append the information 
    # as a list. Each entry in the list is an astropy.table.row.Row object
    for j in range(0,len(assoc_groupNum)):
        assoc_triple_list.append(WDS_grouped.groups[assoc_groupNum[j]][assoc_stars[j]])
        if j%200==0:
            print(f"added row number {j} to assoc_triple_list \n")

    # Now convert the assoc_triple_data 'list' class into a Table class
    assoc_triple_data = Table(rows=assoc_triple_list, names=WDS_grouped.colnames)
    # Write out the triple data
    assoc_triple_data.write(assoc_csv_filename, format='ascii.csv', overwrite=True)
    assoc_triple_data.write(assoc_ecsv_filename, format='ascii.ecsv', overwrite=True)
    assoc_triple_data.write(assoc_votable_filename, format='votable', overwrite=True)    

    
    print("**********")
    print("Associated triple plus data written to the following files:")
    print(f"Files: \n {assoc_csv_filename},\n {assoc_ecsv_filename} and \
         \n {assoc_votable_filename} \n") 
    print("**********")

    # before moving to nonassociated groups delete the assoc_triple_list
    # del(assoc_triple_list)
    
    # Repeat list, loop task for nonassociated systems. The addition of  
    # the nonassociated flags needs to be considered
    nonassoc_triple_list = []
    
    # Loop over all groups and associated stars, append the information 
    # as a list. Each entry in the list is an astropy.table.row.Row object
    for k in range(0,len(nonassoc_groupNum)):
        row = WDS_grouped.groups[nonassoc_groupNum[k]][nonassoc_stars[k]]
        theFlag = nonassoc_flags[k]
        #convert Row class to a Table to access add column method
        aTempTable = Table(rows=[row])
        aTempTable.add_column(Column(name='nonassoc_flag',data=[theFlag]))
        # append the nonassoc data with the flag to the list
        # Zero index needed here in order for one row of 
        # many columns to become a list. Without the Zero 
        # index nonassoc_triple_list rows have a length of one.
        nonassoc_triple_list.append(aTempTable[0])
        if k%2000==0:
            print(f"added row number {k} to nonassoc_triple_list \n")

    # Now convert the assoc_triple_data 'list' class into a Table class
    # Recall we added nonassoc_flags to the output so a column name 
    # needs top be added. 
    nonassoc_colnames = WDS_grouped.colnames + ['nonassoc_flag']
    nonassoc_triple_data = Table(rows=nonassoc_triple_list, names=nonassoc_colnames)

    # Write out the nonassociated triple data
    nonassoc_triple_data.write(nonassoc_csv_filename, format='ascii.csv', overwrite=True)
    nonassoc_triple_data.write(nonassoc_ecsv_filename, format='ascii.ecsv',overwrite=True)
    nonassoc_triple_data.write(nonassoc_votable_filename, format='votable',overwrite=True)    
    
    # everything worked, return 0
    
    return 0

def makeIndices(WDS_grouped):

    """
    Creates the indices for the multiCompIndices and threePlusIndices 
    and saves them to pickle files for reuse. The pickle files saves 
    time as we avoid running this loop repetitively. If the user is happy 
    with the indices as they exist in the pickel file just read in the 
    pickel file and proceed!
    
    The index files are:
    
        multiCompIndices: indices identify all WDS entries with 
        two or more components
    
        threePlusIndices: indices identify all WDS entries with 
        three or more components
    
    INPUT:
    
        WDS_grouped: the grouped WDS catalog
    
    OUTPUT:
    
        indices.pkl: pickle file containing the two index files
    """
    
    print("You entered 'Yes'! Running makeIndices...")
    
    # A FOR loop to find WDS ID's with two or more components. 
    # Create a 1-d array to save these indices and define it as integer 
    # or the default is a float
    multiCompIndices = np.array([], dtype=np.int32)
    # The array multiCompNumber holds for each index in multiCompIndices
    # the number of components found for the corresponding WDS_Identifier
    multiCompNumber = np.array([], dtype=np.int32)
    multiCompCount = 0

    bigCountNumber = 2
    threePlusCountNumber = 3
    threePlusIndices = np.array([], dtype=np.int32)

    print("Organizing the indices for the WDS Components\n")
    for i in range(0,len(WDS_grouped.groups.keys)):
        numComps = len(WDS_grouped.groups[i]) 
        if i%5000==0:
            print('Component numbers checked for: ',i)   
        if  numComps >= bigCountNumber:
            # print(f"For index {i} {numComps} components found.")
            multiCompIndices = np.append(multiCompIndices, i)
            multiCompNumber  = np.append(multiCompNumber, numComps)
            multiCompCount +=1
            if numComps >= threePlusCountNumber:
                threePlusIndices = np.append(threePlusIndices, i) 
    
    print("\n")
    print(f"I found {multiCompCount} multi-component WDS systems.\n")

    # Write out indices to a file so we can read them in and avoid 
    # this function if desired
    indicesPickleFile = "indices.pkl"
    with open(indicesPickleFile, "wb") as binary_file:
        pickle.dump((threePlusIndices,multiCompIndices),binary_file)
    print(f"\n Indices saved to pickle file: {indicesPickleFile} \n")

    return 0

########################################################################

# For now, assume all work and files are in the current working directory
print("\n")
print("*************************************")
print("* Starting sort_query_results_22.py *")
print("*************************************\n")

# Start performance counter
startPerfTime = time.perf_counter()

save_path = os.getcwd()

# # Read in the query results table file

# qrt = query results table

print("Reading in the WDS Gaia matched and vetted query results\n")

qrt ='{0}/merge_sort_on_wds_gaia.ecsv'.format(save_path) 
query_results_table = QTable.read(qrt, header_start=0, data_start=1)

# Use group_by() method of Table so I can group all entries by the same WDS Identifier and 
# avoid some FOR loops 
# WDS_grouped = small_test_table.group_by('WDS_Identifier')
WDS_grouped = query_results_table.group_by('WDS_Identifier')
# Get indices of groups
group_indices = WDS_grouped.groups.indices

response = input("Do you want to run makeIndices? (yes/no): ").strip().lower()

if response == 'yes':
    makeIndices(WDS_grouped)
    with open("indices.pkl","rb") as binary_file:
        (threePlusIndices,multiCompIndices)=pickle.load(binary_file)
elif response == 'no':
    with open("indices.pkl","rb") as binary_file:
        (threePlusIndices,multiCompIndices)=pickle.load(binary_file)

# Attempt to create a table of primary star only entries via the unique function of astropy.table
# Note that WDS 00003-0654 LDS6079 is seelcted as by accident the second component was dropped off 
# in the 50 row selection
primaries_only = unique(WDS_grouped, keys=['WDS_Identifier'], keep='none')

# Select from the query_results_table the WDS stars that only have Gaia DR3 matches on the primaries
# report the length of the result
WDS_match_prim_gaia = unique(query_results_table, keys=['WDS_Identifier'], keep='none')
num_prim_only = len(WDS_match_prim_gaia)
print(f"There are {num_prim_only} WDS entries where only the primary has a Gaia DR3 match.")
print("Writing the primary star only results to their output files.\n")

# SAVE Primary only results TABLE AS ECSV and CSV and as a VOTable
# Two different ways of coding the ascii-write are 
# shown as an example for how to write these and one is a 
# little more concise.
# Start performance counter
startPerfTime = time.perf_counter()
# start CPU Processor counter
startProcTime = time.process_time()

ascii.write(WDS_match_prim_gaia, save_path+'/WDS_match_prim_gaia.ecsv', format='ecsv', overwrite=True)
stopEcsvPerfTime = time.perf_counter()
print('Primary only .ecsv file done.')
print(f"Performance time for ECSV writing=: {stopEcsvPerfTime - startPerfTime} seconds.\n")

WDS_match_prim_gaia.write(save_path+'/WDS_match_prim_gaia.csv', overwrite=True)
stopCsvPerfTime = time.perf_counter()
print('Primary only .csv file done.')
print(f"Performance time for CSV writing=: {stopCsvPerfTime - stopEcsvPerfTime} seconds.\n")

WDS_match_prim_gaia.write(save_path+'/WDS_match_prim_gaia.vot', format = 'votable', overwrite=True)
stopVotPerfTime = time.perf_counter()
print('Primary only .vot file done.')
print(f"Performance time for VOT writing=: {stopVotPerfTime - stopCsvPerfTime} seconds.\n")

stopPerfTime = time.perf_counter()
endProcTime = time.process_time()
elapProcTime = endProcTime - startProcTime
print(f"CPU Time for writing 3 files =: {elapProcTime} seconds.")
print(f"Total Performance time for Primary star only table writing=: \
{stopPerfTime - startPerfTime} seconds.\n")

# Before looping through multi-component systems add columns for testing if 
# multi-component systems are 
# associated or not. Using the small table as a test to begin with now.
tableLength = len(WDS_grouped)
WDS_grouped.add_column(np.zeros(tableLength)*u.arcsec, name='separation')
# save P.A. and maybe include it later
# WDS_grouped.add_column(np.zeros(tableLength)*u.deg, name='thetaGaia')
WDS_grouped.add_column(np.zeros(tableLength)*u.mas**2/u.yr**2, name='delta_mu_ra2')
WDS_grouped.add_column(np.zeros(tableLength)*u.mas**2/u.yr**2, name='delta_mu_dec2')
WDS_grouped.add_column(np.zeros(tableLength)*u.mas/u.yr, name='delta_mu')
WDS_grouped.add_column(np.zeros(tableLength)*u.mas/u.yr, name='sigma_delta_mu')
WDS_grouped.add_column(np.zeros(tableLength)*u.mas/u.yr, name='delta_mu_orbit')
WDS_grouped[0:3]

# Begin tests for associated and non-associated systems for pairs

non_assocCount = 0
maxSepCount = 0
assocCount = 0
non_assoc_indices = []
non_assoc_bySep_indices = []
assoc_indices = []
noncon_plaxCount = 0
nonassoc_byPlax_indices = []
pairCount = 0 # count number of pairs only
keplerCount = 0
noncon_pmKepCount = 0
nonassoc_pmKep_indices = []
unresolvedPairCount = 0

print("Checking WDS pairs first.")
print(f"There are {len(multiCompIndices)} entries to check. \n")

# Open the output ascii file
with open('pairsOutput.txt', 'w') as pairsOut:
    
# Looping through for systems with 2 components within the 
# systems identified via the indices recorded within multiCompIndices
    for anIndex in multiCompIndices:
        if anIndex%2000==0:
            print(f"\n Pair Association checked through multiCompIndex {anIndex}")
        if len(WDS_grouped.groups[anIndex]) == 2:
            pairCount+=1
            brightest = np.where(WDS_grouped.groups[anIndex]['phot_g_mean_mag'] == 
                np.min(WDS_grouped.groups[anIndex]['phot_g_mean_mag']))
            faintest  = np.where(WDS_grouped.groups[anIndex]['phot_g_mean_mag'] == 
                np.max(WDS_grouped.groups[anIndex]['phot_g_mean_mag']))
            parallax_to_use = WDS_grouped.groups[anIndex]['parallax'][brightest]
            ra_a, ra_b = WDS_grouped.groups[anIndex]['ra']
            dec_a, dec_b = WDS_grouped.groups[anIndex]['dec']
            par_a, par_b = WDS_grouped.groups[anIndex]['parallax']
            # use the ra and dec of components a and b to make SkyCoord objects
            coord_a = SkyCoord(ra_a, dec_a, frame='icrs')
            coord_b = SkyCoord(ra_b, dec_b, frame='icrs')
            # calculate the separation between the two objects
            WDS_grouped.groups[anIndex]['separation'][faintest] = coord_a.separation(coord_b).to(u.arcsec)
            # maximum separation test
            # Set pair_sep to the separation
            pair_sep = WDS_grouped.groups[anIndex]['separation'][faintest]
            # Max separation uses parallax in milli-arcseconds and returns arcsec
            # thus we make a unit conversion
            this_max_sep = 206.265*parallax_to_use*u.arcsec/u.mas
    
            # Set logical test equal to a result and use that result in the 'if'
            # statement. That is how to do logical tests with astropy units
            sep_test = pair_sep > this_max_sep
    
            if sep_test[0]:
                # system unbound due to Galactic tidal field: non-associated
                pairsOut.write(f"\n Pair unbound due to tidal field: {WDS_grouped.groups[anIndex]['WDS_Identifier'][0]} \n")
                pairsOut.write(f" Pair sep = {pair_sep} and max separation = {this_max_sep} \n")
                non_assocCount+=1
                maxSepCount+=1
                non_assoc_indices.append(anIndex)
                non_assoc_bySep_indices.append(anIndex)
                continue
            absParDiff = np.absolute(par_a - par_b)
            sig_par_a, sig_par_b = WDS_grouped.groups[anIndex]['parallax_error']
            plaxErrSum = sqrt(sig_par_a**2/(1. * u.mas)**2 \
                              + sig_par_b**2/(1. * u.mas)**2)
            # b_test = pair_sep > 4.0*u.arcsec
            # if b_test[0]: 
            # pdb.set_trace()
            if pair_sep.value[0] > 4.0:
                b = 3.0
            else: 
                b = 6.0
            # Test for consistency of parallaxes
            if absParDiff.value > b*plaxErrSum:
                # Parallaxes NOT consistent within the errors
                pairsOut.write(f"\n Parallaxes are inconsistent: {WDS_grouped.groups[anIndex]['WDS_Identifier'][0]} \n")
                pairsOut.write(f"Parallax difference = {absParDiff} and sigma check = {b*plaxErrSum} \n")
                non_assocCount+=1
                noncon_plaxCount+=1
                non_assoc_indices.append(anIndex)
                nonassoc_byPlax_indices.append(anIndex)
                continue
            # Start process for proper motions against a Keplerian orbit
            # print(f"I am in the Kepler velocity part for groupNum {groupNum}")
            # At the end of the loop: 
            #      keplerCount minus noncon_pmKepCount = assocCount
            keplerCount+=1
            # Check if the stars are unresolved
            if coord_a == coord_b:
                # unresolved, add them to the associated count
                pairsOut.write(f'\n *** In Pair stars = {anIndex} coord_a = coord_b ***')
                pairsOut.write(f'*** For Pair stars {anIndex} are unresolved and associated ***\n')
                pairsOut.write(f"{WDS_grouped.groups[anIndex]['WDS_Identifier','Discovr','Comp'][0]} \n")
                unresolvedPairCount+=1
                assocCount+=1
                assoc_indices.append(anIndex)
                continue
            pmra_a, pmra_b = WDS_grouped.groups[anIndex]['pmra']
            pmdec_a, pmdec_b = WDS_grouped.groups[anIndex]['pmdec']
            pmra_error_a, pmra_error_b = WDS_grouped.groups[anIndex]['pmra_error']
            pmdec_error_a, pmdec_error_b = WDS_grouped.groups[anIndex]['pmdec_error']
            mu_star_ra_a = pmra_a.value * cos(radians(dec_a.value))
            mu_star_ra_b = pmra_b.value * cos(radians(dec_b.value))
            delta_mu_ra2 = (mu_star_ra_a - mu_star_ra_b)**2
            delta_mu_dec2 = (pmdec_a.value - pmdec_b.value)**2
            delta_mu = (delta_mu_ra2 + delta_mu_dec2)**(1/2)
            sigma_term_1 = (pmra_error_a.value**2 + pmra_error_b.value**2) * delta_mu_ra2
            sigma_term_2 = (pmdec_error_a.value**2 + pmdec_error_b.value**2) * delta_mu_dec2
            sigma_delta_mu = (1/delta_mu) * (sigma_term_1 + sigma_term_2)**(1/2)
            # set a conversion factor to get correct units for delta_mu_orbit
            c_factor = (u.arcsec**(1/2))/(u.yr*u.mas**(1/2))
            delta_mu_orbit = 0.44*(parallax_to_use)**(3/2)*(pair_sep)**(-1/2)*c_factor
            # pdb.set_trace()
            if delta_mu > delta_mu_orbit.value + 2*sigma_delta_mu:
                pairsOut.write(f'\n *** For Pair stars {anIndex} vel > Keplerian estimate. ***\n')
                pairsOut.write(f"{WDS_grouped.groups[anIndex]['WDS_Identifier','Discovr','Comp'][0]}\n")  
                pairsOut.write(f'Delta_mu = {delta_mu} and is > {delta_mu_orbit.value + 2*sigma_delta_mu} \n')
                noncon_pmKepCount+=1
                non_assocCount+=1
                nonassoc_pmKep_indices.append(anIndex)
                non_assoc_indices.append(anIndex)
                continue
            # If we are here the pair is associated
            else:
                assocCount+=1
                assoc_indices.append(anIndex)

    print("The pair results are: ")
    print(f"\n Associated pairs found = {assocCount}")
    print(f"Non-associated pairs count = {non_assocCount}")
    print(f"Their sum should equal => {pairCount} \n")
    pairsOut.write(f"\n Associated pairs found = {assocCount}\n")
    pairsOut.write(f"Non-associated pairs count = {non_assocCount}\n")
    pairsOut.write(f"Their sum should equal => {pairCount} \n")
    
# Now, save associated and non-associated indices into 
# indices for the pairs ONLY. Why? Because we are treating 
# the pairs and 3-plus component groups separately. 
assoc_pair_indices = assoc_indices
non_assoc_pair_indices = non_assoc_indices
non_assoc_pair_bySep_indices = non_assoc_bySep_indices
non_assoc_pair_byPlax_indices = nonassoc_byPlax_indices
non_assoc_pair_pmKep_indices = nonassoc_pmKep_indices

# Write out the results for the pairs
print(f"A total of {pairCount} pairs were checked for associations.")
print(f"Writing out the associated results for {assocCount} pairs in three file formats. \n")

#ascii.write(WDS_grouped.groups[assoc_pair_indices], save_path+'/WDS_match_pairs_assoc.ecsv', format='ecsv', overwrite=True)
print('Associated matched pairs .ecsv file done.\n')

WDS_assoc_pair = WDS_grouped.groups[assoc_pair_indices]

#WDS_assoc_pair.write(save_path+'/WDS_match_pairs_assoc.csv', overwrite=True)
print('Associated matched pairs .csv file done.\n')

#WDS_assoc_pair.write(save_path+'/WDS_match_pairs_assoc.vot', format = 'votable', overwrite=True)
print('Associated matched pairs .vot file done.\n')

WDS_non_assoc_pair = WDS_grouped.groups[non_assoc_pair_indices]
# Add the non-associated flags to the WDS_non_assoc_pair table

non_assoc_flags = np.zeros(len(WDS_non_assoc_pair), dtype=int)
# Set flags as bySep  = 1
#              byPlax = 2
#              pmKep  = 3
# and convert the indices that are currently 
# class list -> np.array
flg_ind_1 = np.array(non_assoc_pair_bySep_indices)
flg_ind_2 = np.array(non_assoc_pair_byPlax_indices)
flg_ind_3 = np.array(non_assoc_pair_pmKep_indices)

# Get the WDS ID's from the main table of the of the WDS 
# entries that require the flags of 1, 2, or 3.
ids_to_flag_1 = WDS_grouped.groups[flg_ind_1]['WDS_Identifier']
ids_to_flag_2 = WDS_grouped.groups[flg_ind_2]['WDS_Identifier']
ids_to_flag_3 = WDS_grouped.groups[flg_ind_3]['WDS_Identifier']

# Set flags in non_assoc_table where 'WDS_Identifier' matches
non_assoc_flags[np.isin(WDS_non_assoc_pair['WDS_Identifier'], ids_to_flag_1)] = 1
non_assoc_flags[np.isin(WDS_non_assoc_pair['WDS_Identifier'], ids_to_flag_2)] = 2
non_assoc_flags[np.isin(WDS_non_assoc_pair['WDS_Identifier'], ids_to_flag_3)] = 3

# Add the new column to the WDS_non_assoc_pair table
WDS_non_assoc_pair.add_column(Column(data=non_assoc_flags, name='non_assoc_flag'))

# Write out the non-associated pair results

print(f"Writing out the non-associated results for {non_assocCount} pairs in three file formats. \n")

ascii.write(WDS_non_assoc_pair, save_path+'/WDS_non_assoc_pairs.ecsv', format='ecsv', overwrite=True)
print('Non-associated pairs .ecsv file done.\n')

WDS_non_assoc_pair.write(save_path+'/WDS_non_assoc_pairs.csv', overwrite=True)
print('Non-associated pairs .csv file done.\n')

WDS_non_assoc_pair.write(save_path+'/WDS_non_assoc_pairs.vot', format = 'votable', overwrite=True)
print('Non-associated pairs .vot file done.\n')

## End writing out results for pairs

# pdb.set_trace()

# Begin to work on results for WDS entries with 3 or more stars
print(f"Starting to check {len(threePlusIndices)} WDS entries with 3+ components.")
print(f"The maximum index number printed out will be {threePlusIndices[len(threePlusIndices)-1]} .\n")

# Create lists and counters for index number to write out for 
# 3+ WDS entries
threePlusIndexAssoc = []
assocThreePlusCount = 0
# count the primary stars found amongst the triples
# This must equal number of triples
triplePrimaryCount = 0

# Create a type None as a placeholder for the non-associated 
# array variable. This will be a type None until the first 
# non-associated component is identified in the 3-plus WDS entries
# Also one for the associated array and one to hold the primaries 

nonassoc_TriplePlus = None 
assoc_TriplePlus = None
primaries_TriplePlus = None

# Reset the non-associated index lists
non_assoc_bySep_indices = []
nonassoc_byPlax_indices = []
nonassoc_pmKep_indices = []

nonassocTriplesFile = "nonassocTriplesOut.txt"
nonassocTriplesOut = open(nonassocTriplesFile, 'w')
# Record systems with multiple primaries and count them
twoOrMorePrimaries = open("twoPlusPrimaries.txt", 'w')
manyPrimariesCount = 0

for anIndex in threePlusIndices:
    # create a list to hold the components to write out for 
    # the 3+ WDS entries that ARE associated. This is reset 
    # for each index number 

    # set groupNum equal to anIndex for now for reuse convenience
    # print(f"A start of threePlusIndices loop anIndex equals: {anIndex}")
    groupNum = anIndex
    if anIndex%2000==0:
        print(f"Association checked through threePlusIndices {anIndex}")
    # if WDS_grouped.groups[groupNum]['WDS_Identifier'][0] == '23599+6108':
      #   print(f" WDS ID =: {WDS_grouped.groups[groupNum]['WDS_Identifier'][0]} \n")
      #   pdb.set_trace()
    brightestMag = np.min(WDS_grouped.groups[groupNum]['phot_g_mean_mag'])
    brightest = np.where(WDS_grouped.groups[groupNum]['phot_g_mean_mag'] == brightestMag)[0]
    theComponents = WDS_grouped.groups[groupNum]['Comp']
    # Check if there are any masked entries in theComponents
    if theComponents.mask is not np.ma.nomask and np.any(theComponents.mask):
        # Allow for any masked entries mixed in with strings
        # Find indices of unmasked components
        theCompsUnmaskedIndex = np.where(~theComponents.mask)[0]
        # unmasked values of components
        theComponents = theComponents[~theComponents.mask]
        # But a gotcha is if all Components = '--' 
        if len(theComponents) == 0:
            theComponents = np.arange(0, len(WDS_grouped.groups[groupNum]['Comp']), 1)
        # Find the minimum value among unmasked elements
        min_value = min(theComponents)
        if min_value == 0:
            lowComponent = [0]
        else:
            orig_index_min_comp = theCompsUnmaskedIndex[theComponents == min_value]
            # lowComponent = np.where(theComponents == min(theComponents))[0]
            lowComponent = orig_index_min_comp

    else:
        lowComponent = np.where(theComponents == min(theComponents))[0]
    # 
    theThreePlusIndex = groupNum    
    commonBrightComponent = np.intersect1d(lowComponent,brightest)
    if len(commonBrightComponent) == 0:
        # Failed to find the commonBrightComponent with intersection
        # Use the Comp and find the brightest one and set to primary
        aPrimary = np.where(WDS_grouped.groups[groupNum]['phot_g_mean_mag'][lowComponent] 
            == min(WDS_grouped.groups[groupNum]['phot_g_mean_mag'][lowComponent]))[0]
        # Set the values for the primary star
        commonBrightComponent = aPrimary
        parallax_to_use = WDS_grouped.groups[theThreePlusIndex]['parallax'][commonBrightComponent]
        # Max separation uses parallax in milli-arcseconds and returns 
        # arcsec
        this_max_sep = 206.265*parallax_to_use*u.arcsec/u.mas
        ra_a = WDS_grouped.groups[theThreePlusIndex]['ra'][commonBrightComponent]
        dec_a = WDS_grouped.groups[theThreePlusIndex]['dec'][commonBrightComponent]
        coord_a = SkyCoord(ra_a, dec_a, frame='icrs')
        sig_par_a = WDS_grouped.groups[theThreePlusIndex]['parallax_error'][commonBrightComponent]
    elif len(commonBrightComponent) == 1:
        # Found a single bright component, set as primary
        parallax_to_use = WDS_grouped.groups[theThreePlusIndex]['parallax'][commonBrightComponent]        
        this_max_sep = 206.265*parallax_to_use*u.arcsec/u.mas
        ra_a = WDS_grouped.groups[theThreePlusIndex]['ra'][commonBrightComponent]
        dec_a = WDS_grouped.groups[theThreePlusIndex]['dec'][commonBrightComponent]
        coord_a = SkyCoord(ra_a, dec_a, frame='icrs')
        sig_par_a = WDS_grouped.groups[theThreePlusIndex]['parallax_error'][commonBrightComponent]
    else:
        # len(commonBrightComponent) >= 2
        manyPrimariesCount+=1        
        manyPrimariesToWrite = WDS_grouped.groups[110026]['WDS_Identifier'][0] \
        +' number of primaries = '+str(len(commonBrightComponent))
        twoOrMorePrimaries.write(manyPrimariesToWrite)
         
        
    # The commonBrightComponent i.e. the Primary star for this multi-component group 
    # is identified. Conduct the checks for associated/non-associated

    
    print(f"Entering stars loop for groupNum: {groupNum}")
    for stars in range(len(WDS_grouped.groups[theThreePlusIndex])):
        print(f"At top of for loop star = {stars}")
        print(f"{WDS_grouped.groups[groupNum]['WDS_Identifier','Discovr','Comp'][stars]}")
        #if groupNum == 11254:
        #    print(f"commonBrightComponent is: {commonBrightComponent}")
        if stars == commonBrightComponent:
            # do not evaluate primary against primary
            print(f"Primary star = {stars}, recording the primary.\n")
            triplePrimaryCount+=1
            if isinstance(primaries_TriplePlus,type(None)):
                primaries_TriplePlus = np.array([theThreePlusIndex,stars])
            else:
                primaries_TriplePlus = np.vstack((primaries_TriplePlus,[theThreePlusIndex,stars]))
            # The continue below loops back through the primary again, comment it out
            continue
        # else:
        # now loop primary over other components
        # And keep track of the indices and components
        group_number = theThreePlusIndex
        ra_b = WDS_grouped.groups[group_number]['ra'][stars]
        # print(f"ra_b returns {ra_b}")
        dec_b = WDS_grouped.groups[group_number]['dec'][stars]
        coord_b = SkyCoord(ra_b, dec_b, frame='icrs')
        sig_par_b = WDS_grouped.groups[group_number]['parallax_error'][stars]
        # calculate the separation between the two objects
        WDS_grouped.groups[group_number]['separation'][stars] = coord_a.separation(coord_b).to(u.arcsec)
        # Set pair_sep to a number without units of arcsec for comparison in the "If test"
        pair_sep = WDS_grouped.groups[group_number]['separation'][stars]
        # print(f"coord_a = {coord_a} and coord_b = {coord_b}")
        sep_test = pair_sep > this_max_sep
        
        if sep_test[0]:
            theTripleID = WDS_grouped.groups[group_number]['SOURCE_ID','WDS_Identifier','Discovr','Comp'][stars]
            nonassocTriplesOut.write(f"\n 3+ unbound due to tidal field: \n {theTripleID} \n")
            nonassocTriplesOut.write(f" 3+ sep = {pair_sep} and max separation = {this_max_sep} \n")
            print(f"Pair separation too large = {pair_sep} max_sep = {this_max_sep}")
            print("\n")
            non_assocCount+=1
            maxSepCount+=1
            flag=1
            if anIndex not in non_assoc_indices:
                non_assoc_indices.append(anIndex)
            if anIndex not in non_assoc_bySep_indices:
                non_assoc_bySep_indices.append(anIndex)
            if nonassoc_TriplePlus is None:
                nonassoc_TriplePlus = np.array([group_number,stars,flag])
            else:
                nonassoc_TriplePlus = extendFlagRows(nonassoc_TriplePlus,group_number,stars,flag)
            continue
        else:
            print(f"Pair sep = {pair_sep} max_sep = {this_max_sep} OK")
        par_b = WDS_grouped.groups[group_number]['parallax'][stars]
        absParDiff = np.absolute(parallax_to_use - par_b)
        sig_par_b = WDS_grouped.groups[group_number]['parallax_error'][stars]
        # sqrt needs to return a float, convert parallax errors to 
        # dimension-less quantities
        plaxErrSum = sqrt(sig_par_a**2/(1. * u.mas)**2 + \
                          sig_par_b**2/(1. * u.mas)**2)
        # b_test = pair_sep > 4.0*u.arcsec
        # pdb.set_trace()
        # if b_test:
        # pdb.set_trace()
        if pair_sep.value > 4.0:
            b = 3.0
        else: 
            b = 6.0
        # Test for consistency of parallaxes
        if absParDiff.value > b*plaxErrSum:
            # print(f"\n Parallaxes are inconsistent: {WDS_grouped.groups[group_number]['WDS_Identifier'][0]} ")
            # print(f"Parallax difference = {absParDiff} and sigma check = {b*plaxErrSum} \n")
            nonassocTriplesOut.write(f"\n Parallaxes are inconsistent: {WDS_grouped.groups[group_number]['WDS_Identifier'][0]} \n")
            nonassocTriplesOut.write(f"Parallax difference = {absParDiff} and sigma check = {b*plaxErrSum} \n")
            non_assocCount+=1
            noncon_plaxCount+=1
            flag=2
            if anIndex not in non_assoc_indices:
                non_assoc_indices.append(anIndex)
            if anIndex not in nonassoc_byPlax_indices:    
                nonassoc_byPlax_indices.append(anIndex)
            if nonassoc_TriplePlus is None:
                nonassoc_TriplePlus = np.array([group_number,stars,flag])
            else:
                nonassoc_TriplePlus = extendFlagRows(nonassoc_TriplePlus,group_number,stars,flag)  
            continue
        else:
            print(f"For star = {stars} parallaxes are consistent.")
        # Start process for proper motions against a Keplerian orbit
        print(f"I am in the Kepler velocity part for groupNum {group_number}")
        # keplerCount+=1
        # Check if the stars are unresolved
        if coord_a == coord_b:
            selected_columns = WDS_grouped.columns[:31]
            print(f"*** star = {stars} coord_a = coord_b ***")
            print(f"*** star {stars} are unresolved and associated ***""")
            assocThreePlusCount+=1
            print(f"assocThreePlusCount = {assocThreePlusCount}.\n")
            #if anIndex not in threePlusIndexAssoc:
            #    threePlusIndexAssoc.append(anIndex)
            threePlusIndexAssoc.append(anIndex)
            if assoc_TriplePlus is None:
                assoc_TriplePlus = np.array([group_number,stars])
            else:
                assoc_TriplePlus = extendAssocRows(assoc_TriplePlus,group_number,stars)
            continue
        pmra_a = WDS_grouped.groups[group_number]['pmra'][commonBrightComponent]
        pmdec_a = WDS_grouped.groups[group_number]['pmdec'][commonBrightComponent]
        pmra_b = WDS_grouped.groups[group_number]['pmra'][stars]
        pmdec_b = WDS_grouped.groups[group_number]['pmdec'][stars]
        pmra_error_a = WDS_grouped.groups[group_number]['pmra_error'][commonBrightComponent]
        pmdec_error_a = WDS_grouped.groups[group_number]['pmdec_error'][commonBrightComponent]
        pmra_error_b = WDS_grouped.groups[group_number]['pmra_error'][stars]
        pmdec_error_b = WDS_grouped.groups[group_number]['pmdec_error'][stars]
        mu_star_ra_a = pmra_a.value * cos(radians(dec_a.value))
        mu_star_ra_b = pmra_b.value * cos(radians(dec_b.value))
        delta_mu_ra2 = (mu_star_ra_a - mu_star_ra_b)**2
        delta_mu_dec2 = (pmdec_a.value - pmdec_b.value)**2
        delta_mu = (delta_mu_ra2 + delta_mu_dec2)**(1/2)
        sigma_term_1 = (pmra_error_a.value**2 + pmra_error_b.value**2) * delta_mu_ra2
        sigma_term_2 = (pmdec_error_a.value**2 + pmdec_error_b.value**2) * delta_mu_dec2
        sigma_delta_mu = (1/delta_mu) * (sigma_term_1 + sigma_term_2)**(1/2)
        c_factor = (u.arcsec**(1/2))/(u.yr*u.mas**(1/2))
        delta_mu_orbit = 0.44*(parallax_to_use)**(3/2)*(pair_sep)**(-1/2)*c_factor
        if delta_mu > delta_mu_orbit.value + 2*sigma_delta_mu:
            #print(f"*** For triple+ stars {anIndex} vel > Keplerian estimate. ***\n""")
            #print(f"{WDS_grouped.groups[group_number]['WDS_Identifier','Discovr','Comp'][0]}")  
            #print(f"Delta_mu = {delta_mu} and is > {delta_mu_orbit.value + 2*sigma_delta_mu} \n")
            nonassocTriplesOut.write(f"\n *** For 3+ stars {anIndex} vel > Keplerian estimate. ***\n""")
            nonassocTriplesOut.write(f"{WDS_grouped.groups[group_number]['WDS_Identifier','Discovr','Comp'][0]}\n")  
            nonassocTriplesOut.write(f"Delta_mu = {delta_mu} and is > {delta_mu_orbit.value + 2*sigma_delta_mu} \n")
            noncon_pmKepCount+=1
            non_assocCount+=1
            flag=3
            if anIndex not in nonassoc_pmKep_indices:
                nonassoc_pmKep_indices.append(anIndex)
            if anIndex not in non_assoc_indices:
                non_assoc_indices.append(anIndex)
            if nonassoc_TriplePlus is None:
                nonassoc_TriplePlus = np.array([group_number,stars,flag])
            else:
                nonassoc_TriplePlus = extendFlagRows(nonassoc_TriplePlus,group_number,stars,flag)      
            continue
        else:
            print(f"For star = {stars} PM good for bound orbit.")
            # We are here, so all checks are good and this pair 
            # component 'stars' for index = groupNum is associated
            assocThreePlusCount+=1
            print(f"assocThreePlusCount = {assocThreePlusCount}.\n")
            # if anIndex not in threePlusIndexAssoc:
            #    threePlusIndexAssoc.append(anIndex)
            threePlusIndexAssoc.append(anIndex)
            if assoc_TriplePlus is None:
                assoc_TriplePlus = np.array([group_number,stars])
            else:
                assoc_TriplePlus = extendAssocRows(assoc_TriplePlus,group_number,stars)

nonassocTriplesOut.close()
twoOrMorePrimaries.close()
print(f"\n I found {manyPrimariesCount} with multiple primaries.")
print(f" The systems, if any, are listed in the file twoPlusPrimaries.txt. \n")
# Stop the timer

# Get the truncated 31 column copy of WDS_grouped for ease of writing the 
# associated and non-associated triple results

writeIt = writeTriples(assoc_TriplePlus,nonassoc_TriplePlus,WDS_grouped,save_path)

if writeIt == 2:
    raise RuntimeError("writeTriples: assoc_groupNum and assoc_stars have unequal lengths")
elif writeIt == 3:
    raise RuntimeError("writeTriples: length problem with nonassoc file lengths ")
else:
    print(f"writeIt successfully completed! \n")

stopPerfTime = time.perf_counter()

print("\n Completed checking the 3+ component WDS catalog entries.\n")
# Inform user of time to run the code
print(f"Total Performance time for sort_query_results =: {(stopPerfTime - startPerfTime)/60.0} minutes.\n")