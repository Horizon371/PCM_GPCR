from get_chembl_data import *
import chembl_querries
import constants
from preprocess_data import preprocess_bioactivity_data

querry = chembl_querries.GET_UNIPROT_BIOACTIVITES_FOR_GPCR
bioactivities_csv = f"{constants.DATA_FOLDER}/chembl_gpcr_bioactivities.csv"

run_chembl_querry(querry, bioactivities_csv)
preprocess_bioactivity_data(bioactivities_csv)