import MDAnalysis.analysis.distances
from MDAnalysis import transformations
import os

import numpy as np
import pandas as pd
from tqdm import tqdm

from MDAnalysis.analysis import leaflet


pep_dict = {
    'WF3': 'FLGALIKGAIHGGRFIHGMIQNHH',
    'WF4': 'GWGSIFKHGRHAAKHIGHAAVNHYL',
    'WF2': 'GWGSFFKKAAHVGKHVGKAALTHYL',
    'WF1': 'GKGRWLERIGKAGFIIIGGALDHL',
    'WF1a': 'WLRRIGKGVKIIGGAALDHL',
    'WF1a1': 'GRRKRKWLRRIGKGVKIIGGAALDHL'
}

def get_universe(peptide_path, traj_name, tpr_name, files=False):
    xtc_file_path = None
    peptide_name = os.path.basename(peptide_path)

    if os.path.isfile(f'{peptide_path}/{traj_name}.xtc'):
        xtc_file_path = f'{peptide_path}/{traj_name}.xtc'
        print('FOUND xtc')

    if not xtc_file_path:
        print(f'{peptide_path}/{traj_name}.xtc')
        print(f'Warning!! No xtc found for peptide {peptide_name}')

    tpr_file_path = None

    if os.path.isfile(f'{peptide_path}/{tpr_name}.tpr'):
        tpr_file_path = f'{peptide_path}/{tpr_name}.tpr'
    else:
        print(f'Warning!! No tpr found for peptide {peptide_name}')

    u = MDAnalysis.Universe(tpr_file_path, xtc_file_path)

    if files:
        return u, xtc_file_path, tpr_file_path
    else:
        return u
    
def get_peptide_range(peptides, res_start, pep_num):
    if len(peptides)>1 and ('pg' not in peptides and 'pepg' not in peptides):
        peptides_list = [peptides[0],peptides[0], peptides[0], peptides[0],
                        peptides[1],peptides[1],peptides[1],peptides[1]]
    else:
        peptides_list = [peptides[0],peptides[0], peptides[0], peptides[0]]
    if pep_num == 8:
        peptides_list = [peptides[0],peptides[0], peptides[0], peptides[0],
                        peptides[0],peptides[0],peptides[0],peptides[0]]

    pep_range = {}
    print('peptides_list is ', peptides_list)
    for i, peptide in enumerate(peptides_list):
        res_end = res_start+len(pep_dict[peptide])
        pep_range[f'pep{i+1}'] = list(range(res_start, res_end))
        res_start = res_end
    return pep_range

def get_aa_sequence(universe, peptides, pep_num):
    '''
    universe: MDAnalysis Universe object
    
    pep_num: int
        How many peptides the simulation contains
    '''
    protein_atoms = universe.select_atoms('protein')
    prot_residues = protein_atoms.residues
    res_names = prot_residues.resnames
    res_ids = prot_residues.residues.resids
    pep_num_dict = get_peptide_range(peptides, res_ids[0], pep_num)
    new_dict = dict(zip(res_ids, res_names))

    return pep_num_dict, new_dict

def get_combinations(u, peptide, pep_num, pairs=None):
    '''
    Get list of residue residue combinations

    Input: 
        peptide: list of strin
            list of peptides
        pepnum: int 
            number of peptides of each type
        pairs: 
    '''
    _, new_dict = get_aa_sequence(u, peptide, pep_num)
    peptide_ranges = get_peptide_range(peptide, list(new_dict.keys())[0], pep_num)
    comb_list = []
    for i in peptide_ranges.keys():
        for j in peptide_ranges.keys():
            comb_list.extend([(x, y) for x in peptide_ranges[i] for y in peptide_ranges[j] if (y, x) not in comb_list])
    return list(set(comb_list))


def distance_matrix_time(u, peptide, pep_num):

    step_size = 100
    start, stop, step = u.trajectory.check_slice_indices(None, None, None)

    frames = np.arange(start, stop, step_size)

    n_frames = frames.size

    combinations = get_combinations(u, peptide, pep_num)

    list_arrays = []
    for _, ts in tqdm(enumerate(u.trajectory[frames]), total=n_frames):
        min_dist_per_time = np.zeros((len(combinations), 4))

        #Loop through permutations, create 2D array with amino acid and minimum distance
        for i, c in enumerate(combinations):
            selection1 = u.select_atoms('resid %s'%c[0]).center_of_mass()
            selection2 = u.select_atoms('resid %s'%c[1]).center_of_mass()
            dist = MDAnalysis.analysis.distances.distance_array(selection1, selection2, 
                                                                box=u.dimensions, result=None, backend='serial')
            min_dist = dist
            l = [c[0], c[1], min_dist[0].astype(float)[0], int(ts.time/1000)]
            min_dist_per_time[i, :] = l
        list_arrays.append(min_dist_per_time)
    all_lists = np.concatenate(list_arrays)
    
    return pd.DataFrame(all_lists, columns=['Peptide1','Peptide2', 'mindist', 'Time(ns)'])


membrane = 'pg'
for pep in ['WF2']:
    p=f"../{pep}_pg_2"
    u = get_universe(p, 'md_0_1_combined', 'md_0_1')
    workflow = [transformations.unwrap(u.atoms)]
    u.trajectory.add_transformations(*workflow)
    df_list = distance_matrix_time(u, [l for l in pep.split('_')], 8)
    df_list.to_csv(f'distance_over_time_{pep}.csv')