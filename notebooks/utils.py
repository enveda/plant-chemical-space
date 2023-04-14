# -*- coding: utf-8 -*-

"""Utils to be used in notebooks"""

from collections import defaultdict
from functools import lru_cache
from typing import Set, Tuple, Any, Dict

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
