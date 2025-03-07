import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import rdFingerprintGenerator
import constants

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

def one_hot_encoding(sequence):
    # Initialize a matrix for the sequence (rows = len(sequence), columns = 20)
    encoding_matrix = np.zeros((len(sequence), len(AA_ALPHABET)), dtype=int)
    
    # Encode each amino acid in the sequence
    for i, aa in enumerate(sequence):
        idx = AA_ALPHABET.index(aa)  # Find index of amino acid in the alphabet
        encoding_matrix[i, idx] = 1  # Set corresponding position to 1
    
    return encoding_matrix.tolist()


def preprocess_molecule_data(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    sanitized = sanitize(mol)
    return calculate_fingerprint(sanitized)


def sanitize(mol: Chem.Mol):
    remover = SaltRemover()
    stripped = remover.StripMol(mol)
    Chem.RemoveStereochemistry( stripped ) 
    return stripped


def calculate_fingerprint(mol: Chem.Mol):
    if mol is not None:
        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        fp = mfpgen.GetFingerprint(mol)
        return list(fp)
    else:
        return None  # In case of invalid SMILES

def preprocess_bioactivity_data(bioactivities_csv):

    activities = pl.read_csv(bioactivities_csv).select(["canonical_smiles", "sequence", "standard_value"])

    activities = activities.with_columns([
        pl.col(["canonical_smiles"]).map_elements(preprocess_molecule_data).alias("fingerprint"),
        pl.col(["sequence"]).map_elements(one_hot_encoding).alias("encoding")
    ])

    print(activities)

    activities = activities.filter(
        ~pl.col("canonical_smiles").is_in(["", None]) | ~pl.col("fingerprint").is_null()
    )

    activities.write_parquet(f"{constants.DATA_FOLDER}/chembl_gpcr.parquet")