# -*- coding: utf-8 -*-

"""Utils to be used in notebooks"""

from collections import defaultdict
from functools import lru_cache
from typing import Set, Tuple, Any, Dict
import time
from pubchempy import BadRequestError, get_synonyms, PubChemHTTPError
import chembl_downloader
import sqlite3

import requests
import networkx as nx
import obonet
import pandas as pd
from tqdm import tqdm


@lru_cache(maxsize=None)
def get_ncbiontology():
    """Wrapper to get the ncbiontology."""
    return obonet.read_obo(
        'http://purl.obolibrary.org/obo/ncbitaxon.obo'
    )


def get_genus_and_family_info_for_plants(
    plant_ids: Set[str]
) -> Tuple[Dict[Any, list], Dict[Any, list]]:
    """Return dictionaries mapping plant species to their genus/families."""
    # Load ncbitaxon taxonomy. This might take a couple of mins
    ncbitaxon_ontology = get_ncbiontology()

    # Get the childs of Viridiplantae (all plants)
    plant_childs = nx.ancestors(ncbitaxon_ontology, 'NCBITaxon:33090')

    # Subset the graph to make it faster to the relevant part (plants only)
    ncbitaxon_ontology = ncbitaxon_ontology.subgraph(plant_childs).copy()

    # Relabel nodes to be consistent with the dataset
    ncbitaxon_ontology = nx.relabel_nodes(
        ncbitaxon_ontology,
        {
            node: node.replace('NCBITaxon:', 'ncbitaxon:').strip()
            for node in ncbitaxon_ontology.nodes()
        },
    )

    family_nodes = set()
    genus_nodes = set()

    # Identify the famlies and genus nodes
    for id_, data in ncbitaxon_ontology.nodes(data=True):

        if 'property_value' not in data:
            continue

        # Group edges based on the different taxonomic levels of interst
        for value in data['property_value']:
            if value == 'has_rank NCBITaxon:family':
                family_nodes.add(id_)
            elif value == 'has_rank NCBITaxon:genus':
                genus_nodes.add(id_)

    genus_to_species = defaultdict(list)
    family_to_species = defaultdict(list)

    # Build groups for each family
    for node in tqdm(family_nodes, desc='order family'):
        children = set(nx.ancestors(ncbitaxon_ontology, node))

        for child in children:

            # Check plant is one of the ones in the dataset
            if child not in plant_ids:
                continue

            family_to_species[node].append(child)

    # Build groups for each genus
    for node in tqdm(genus_nodes, desc='order genus'):
        children = set(nx.ancestors(ncbitaxon_ontology, node))

        for child in children:

            # Check plant is one of the ones in the dataset
            if child not in plant_ids:
                continue

            genus_to_species[node].append(child)

    return genus_to_species, family_to_species


def ncbitaxon_curies_to_names(
    plant_ids: Set[str]
) -> Tuple[Dict[Any, list], Dict[Any, list]]:
    """Return dictionaries mapping plant species to their names."""
    # Load ncbitaxon taxonomy. This might take a couple of mins
    ncbitaxon_ontology = get_ncbiontology()

    # Get the childs of Viridiplantae (all plants)
    plant_childs = nx.ancestors(ncbitaxon_ontology, 'NCBITaxon:33090')

    # Subset the graph to make it faster to the relevant part (plants only)
    ncbitaxon_ontology = ncbitaxon_ontology.subgraph(plant_childs).copy()

    return {
        node.replace('NCBITaxon:', 'ncbitaxon:').strip(): data.get('name')
        for node, data in ncbitaxon_ontology.nodes(data=True)
        if node.replace('NCBITaxon:', 'ncbitaxon:').strip() in plant_ids
    }


def create_np_classifier_vectors(
    df: pd.DataFrame
) -> pd.DataFrame:
    """Create NP-Classifier vectors for each plant."""

    all_classes = set(df['class'].unique())
    plants = set(df['plant_curie'].unique())

    vector_data = []

    for plant_name in tqdm(plants, desc='Generated NP-Classifier vectors'):
        phytochem_df = df[df['plant_curie'] == plant_name]

        t = {
            'plant_name': plant_name,
        }

        for class_name in all_classes:
            temp = phytochem_df[phytochem_df['class'] == class_name]
            t[class_name] = len(temp['chemical_curie'].unique())

        vector_data.append(t)

    return pd.DataFrame(vector_data)


def collapse_vector_to_family(
    plant_chemical_df: pd.DataFrame,
    df: pd.DataFrame,
    family_to_species: dict
) -> pd.DataFrame:
    """Collapsing the NP-Classifier vectors to family level."""
    data = []
    skipped_empty = 0
    skipped_med = 0
    skipped_non_med = 0

    plant_type_map = plant_chemical_df[
        ['plant_curie', 'plant_type']
    ].set_index('plant_curie').to_dict()['plant_type']
    df['plant_type'] = df['plant_name'].map(plant_type_map)

    for family_curie in tqdm(family_to_species):
        tmp_df = df[df['plant_name'].isin(family_to_species[family_curie])]

        if tmp_df.empty:
            skipped_empty += 1
            continue

        tmp_df = tmp_df.drop(columns=['plant_name'])

        med_df = tmp_df[tmp_df['plant_type'] == 'Medicinal']
        med_df = med_df.drop(columns=['plant_type'])

        # Remove family pairs with chemical class less than 10
        if med_df.sum().sum() < 10:
            skipped_med += 1
            continue

        med_dict = med_df.sum().to_dict()
        med_dict['family'] = family_curie
        med_dict['ftype'] = 'Medicinal'

        non_med_df = tmp_df[tmp_df['plant_type'] == 'Non-medicinal']
        non_med_df = non_med_df.drop(columns=['plant_type'])

        # Remove family pairs with chemical class less than 10
        if non_med_df.sum().sum() < 10:
            skipped_non_med += 1
            continue

        non_med_dict = non_med_df.sum().to_dict()
        non_med_dict['family'] = family_curie
        non_med_dict['ftype'] = 'Non-medicinal'

        data.append(med_dict)
        data.append(non_med_dict)

    print('Empty skipped -', skipped_empty)
    print('Medicinal skipped -', skipped_med)
    print('Non-medicinal skipped -', skipped_non_med)
    return pd.DataFrame(data)


def get_pubchem_assays(cids: list) -> pd.DataFrame:
    """Get active and inactive bioassays from PubChem"""
    
    assay_df = pd.DataFrame(columns=[
        'compound_id',
        'assay_ids',
        'assay_class',
    ])

    for cid in tqdm(cids, desc='Pinging PubChem for bioassays'):
        # Get assays where compound is active
        active_data = requests.get(
            f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/aids/JSON?aids_type=active'
        ).json()

        inactive_data = requests.get(
            f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/aids/JSON?aids_type=active'
        ).json()

        # Get hit in assay if information present
        if 'InformationList' in active_data:
            d = active_data['InformationList']['Information'][0]['AID']
            active_df = pd.DataFrame({
                'compound_id': cid,
                'assay_ids': d,
                'assay_class': 'active',
            }, index=[0])

            assay_df = pd.concat([assay_df, active_df], ignore_index=True)

        if 'InformationList' in inactive_data:
            d = inactive_data['InformationList']['Information'][0]['AID']
            inactive_df = pd.DataFrame({
                'compound_id': cid,
                'assay_ids': d,
                'assay_class': 'inactive',
            }, index=[0])

            assay_df = pd.concat([assay_df, inactive_df], ignore_index=True)


    return assay_df


def create_genus_compound_vectors(df: pd.DataFrame):
    """Creating the chemical-specie dataframe for statistical analysis."""
    plants = set(df['plant_curie'].unique())
    
    vector_data = []
    
    for plant_name in tqdm(plants, desc='Generated compound vectors'):
        phytochem_df = df[df['plant_curie'] == plant_name]
        total_chemicals = phytochem_df['chemical_curie'].unique()
        chemicals_only_in_plant = set()

        for chem in total_chemicals:
            num_plants_with_chem = df[df['chemical_curie'] == chem]['plant_curie'].nunique()
            if num_plants_with_chem < 2:
                chemicals_only_in_plant.add(chem)
            
        
        t = {
            'plant_name': plant_name,
            '# chemicals': len(total_chemicals),
            '# chemicals in plant only': len(chemicals_only_in_plant)
        }

        vector_data.append(t)
            
    return pd.DataFrame(vector_data)


def create_taxon_compound_vectors(level_mapper: dict, df: pd.DataFrame, med_plants: set):
    """Creating the chemical-family/ chemical-genus dataframe for statistical analysis."""
    
    rows = []
    
    for taxon_name, species in tqdm(level_mapper.items()):
        tmp_df = df[df['plant_curie'].isin(species)]

        total_chemicals = tmp_df.chemical_curie.unique()

        unique_chemicals = set()

        med_plants_in_level = med_plants.intersection(set(species))

        for chem in total_chemicals:
            chem_plants = set(df[df['chemical_curie'] == chem]['plant_curie'].unique())

            non_level_plants = chem_plants - set(species)

            if len(non_level_plants) < 1:
                unique_chemicals.add(chem)
        
        rows.append({
            'name': taxon_name,
            '# chemicals': len(total_chemicals),
            '# level specific chemicals': len(unique_chemicals),
            '# med plants': len(med_plants_in_level),
            '# plants in level': len(species),
            'plants in level': species
        })
        
    return pd.DataFrame(rows)


def get_chembl_id(pubchem_idx: str):
    """Map Pubchem CID to ChEMBL id"""
    try:
        other_idenfitiers = get_synonyms(identifier=pubchem_idx)
    except (PubChemHTTPError, BadRequestError):  # too many request
        time.sleep(3)
        try:
            other_idenfitiers = get_synonyms(identifier=pubchem_idx)
        except BadRequestError:  # incorrect pubchem id
            return None

    if len(other_idenfitiers) < 1:
        return None

    other_idenfitiers = other_idenfitiers[0]

    for idx in other_idenfitiers['Synonym']:
        if idx.startswith('CHEMBL'):
            return idx

    return None


def get_chembl_assays() -> pd.DataFrame:
    """Get active and inactive bioassays from ChEMBL"""
    
    assay_sql = """
    SELECT
        MOLECULE_DICTIONARY.molregno as chembl_id,
        ACTIVITIES.pchembl_value,
        ASSAYS.chembl_id as assay_id,
        COMPONENT_SEQUENCES.accession as target
    FROM MOLECULE_DICTIONARY
    JOIN ACTIVITIES ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
    JOIN ASSAYS ON ACTIVITIES.assay_id == ASSAYS.assay_id
    JOIN TARGET_DICTIONARY on ASSAYS.tid == TARGET_DICTIONARY.tid
    JOIN TARGET_COMPONENTS on TARGET_DICTIONARY.tid == TARGET_COMPONENTS.tid
    JOIN COMPONENT_SEQUENCES on TARGET_COMPONENTS.component_id == COMPONENT_SEQUENCES.component_id
    WHERE
        ASSAYS.assay_type in ('F', 'B')
        and ACTIVITIES.standard_value is not null
        and ACTIVITIES.standard_relation is not null
        and ASSAYS.assay_organism == 'Homo sapiens'
    """

    _, version = chembl_downloader.download_extract_sqlite(return_version=True)
    print(f'Working on CheMBL version {version}')
    assay_df = chembl_downloader.query(sql=assay_sql, version=version)
  
    return assay_df


def _get_assays_from_db(db_path: str, out_path: str):
    """Get assays from ChEMBL 32 database with SQLite."""
    conn = sqlite3.connect(db_path)

    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
    try:
        assert len(cursor.fetchall()) > 1
    except AssertionError:
        print('Incorrect database. Please download the database again.')

    assay_sql = """
    SELECT
        MOLECULE_DICTIONARY.molregno as chembl_id,
        COMPOUND_STRUCTURES.canonical_smiles as smiles,
        ACTIVITIES.pchembl_value as pchembl_value,
        ASSAYS.chembl_id as assay_id
    FROM MOLECULE_DICTIONARY
    JOIN ACTIVITIES ON MOLECULE_DICTIONARY.molregno == ACTIVITIES.molregno
    JOIN COMPOUND_STRUCTURES ON MOLECULE_DICTIONARY.molregno == COMPOUND_STRUCTURES.molregno
    JOIN ASSAYS ON ACTIVITIES.assay_id == ASSAYS.assay_id
    WHERE
        ASSAYS.assay_type in ('F', 'B')
        and ACTIVITIES.standard_value is not null
        and ACTIVITIES.pchembl_value is not null
        and ASSAYS.assay_organism == 'Homo sapiens'
    """

    assay_df = pd.read_sql(assay_sql, con=conn)
    assay_df.to_csv(out_path, sep='\t', index=False)


def get_bioassay_metadata(db_path: str, out_path: str):
    """Get time of registry for bioassays in ChEMBL32."""
    conn = sqlite3.connect(db_path)

    cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
    try:
        assert len(cursor.fetchall()) > 1
    except AssertionError:
        print('Incorrect database. Please download the database again.')

    _sql = """
    SELECT
        DOCS.doc_id as doc_id,
        DOCS.year as year,
        ASSAYS.chembl_id as assay_id
    FROM DOCS
    JOIN ASSAYS ON ASSAYS.doc_id = DOCS.doc_id
    WHERE
        ASSAYS.assay_type in ('F', 'B')
        and ASSAYS.assay_organism == 'Homo sapiens'
    """

    chembl_df = pd.read_sql(_sql, con=conn)
    chembl_df.shape
    chembl_df.to_csv(out_path, sep='\t', index=False)