import polars as pl
from rdkit import Chem
import numpy as np
from rdkit.Chem import rdFingerprintGenerator
import molvs

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

md = molvs.metal.MetalDisconnector()
lfc = molvs.fragment.LargestFragmentChooser()
uc = molvs.charge.Uncharger()


def preprocess_bioactivity_data(bioactivities_csv, processed_data_path):
    activities = pl.read_csv(bioactivities_csv).select(["canonical_smiles", "sequence", "standard_value", "pchembl_value"])

    activities = activities.with_columns([
        pl.col(["canonical_smiles"]).map_elements(standardize_to_fingerprint).alias("fingerprint"),
        pl.col(["sequence"]).map_elements(one_hot_encoding).alias("encoding")
    ])

    activities = activities.filter(
        ~pl.col("canonical_smiles").is_in(["", None]) | ~pl.col("fingerprint").is_null()
    )

    activities = activities.filter(pl.col("fingerprint").list.len() > 0)
    activities = activities.drop_nans()
    activities.write_parquet(processed_data_path)


def standardize_to_fingerprint(smiles: str):
    standardized = standardize_smiles(smiles)
    if standardized:
        return calculate_fingerprint(standardized)
    return None


def standardize_smiles(smiles):
    """
    Standardize a SMILES string by disconnecting metals, choosing the largest organic fragment, and uncharging the molecule.
    """
    stdsmiles = molvs.standardize.standardize_smiles(smiles)
    stdmol = Chem.MolFromSmiles(stdsmiles)
    stdmol = md.disconnect(stdmol)
    stdmol = lfc.choose(stdmol)
    stdmol = uc.uncharge(stdmol)
    if not molvs.Validator().validate(stdmol):
        return stdmol
    return None
    

def calculate_fingerprint(mol: Chem.Mol):
    if mol is not None:
        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        fp = mfpgen.GetFingerprint(mol)
        return list(fp)
    else:
        return None  # In case of invalid SMILES


def one_hot_encoding(sequence):
    # Initialize a matrix for the sequence (rows = len(sequence), columns = 20)
    encoding_matrix = np.zeros((len(sequence), len(AA_ALPHABET)), dtype=int)
    
    # Encode each amino acid in the sequence
    for i, aa in enumerate(sequence):
        idx = AA_ALPHABET.index(aa)  # Find index of amino acid in the alphabet
        encoding_matrix[i, idx] = 1  # Set corresponding position to 1
    
    return encoding_matrix.tolist()