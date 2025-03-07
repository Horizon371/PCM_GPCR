GET_UNIPROT_BIOACTIVITES_FOR_GPCR = """
SELECT DISTINCT 
    cs.molregno,
    cs.canonical_smiles,
    td.chembl_id AS chembl_target_id,
    cse.accession AS uniprot_id,
    cse.sequence,
    act.pchembl_value,
    act.standard_value
FROM target_dictionary td
  JOIN assays a ON td.tid = a.tid
  JOIN activities act ON a.assay_id = act.assay_id
  JOIN molecule_dictionary md ON md.molregno = act.molregno
  JOIN compound_structures cs ON md.molregno = cs.molregno
  JOIN compound_properties cp on md.molregno = cp.molregno
  JOIN target_components tc ON td.tid = tc.tid
  JOIN component_class cc ON cc.component_id = cse.component_id
  JOIN protein_classification pc ON pc.protein_class_id = cc.protein_class_id
  JOIN component_go cg ON cg.component_id = cse.component_id
  JOIN go_classification gc ON gc.go_id = cg.go_id
  JOIN component_sequences cse ON tc.component_id = cse.component_id
  WHERE (
    (
    gc.pref_name LIKE '%G protein%' OR gc.pref_name LIKE '%G-protein%' OR gc.pref_name LIKE '%GPCR%')
    OR (pc.pref_name LIKE '%G protein%' OR pc.pref_name LIKE '%G-protein%' OR pc.pref_name LIKE '%GPCR%')
    )
    AND act.standard_relation = "="
    AND a.confidence_score = 9
    AND cse.organism = "Homo sapiens"
    AND cp.heavy_atoms <= 50
    AND act.standard_units = 'nM'
    AND act.standard_type = "IC50"
"""