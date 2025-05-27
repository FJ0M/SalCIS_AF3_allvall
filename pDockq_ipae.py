import argparse
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pdb
import json

parser = argparse.ArgumentParser(description='''Calculate a predicted DockQ score for a predicted structure.''')
parser.add_argument('--pdbfile', nargs=1, type=str, default=sys.stdin, help='Path to pdbfile to be scored. Note that this file needs to contain at least two chains. The B-factor column is assumed to contain the plDDT score from AlphaFold.')
#parser.add_argument('--rankingfile', nargs=1, type=str, default=sys.stdin, help='Path to ranking_debug.json file.')


#####################FUNCTIONS#########################
def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record


def read_pdb(pdbfile):
    '''Read a pdb file predicted with AF and rewritten to contain all chains
    '''

    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            # Get CB - CA for GLY
            if record['atm_name'] == 'CB' or (record['atm_name'] == 'CA' and record['res_name'] == 'GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'], record['y'], record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'], record['y'], record['z']]]
                    chain_plddt[record['chain']] = [record['B']]

    # Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt

def calc_pdockq(chain_coords, chain_plddt, t,pae_file):
    '''Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    '''

    # Get coords and plddt per chain
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

    # Calc 2-norm
    mat = np.append(coords1, coords2, axis=0)
    a_min_b = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1, l1:]  # upper triangular --> first dim = chain 1
    contacts = np.argwhere(contact_dists <= t)
    #print(contacts,l1)
    if contacts.shape[0] < 1:
        pdockq = 0
        ppv = 0
    else:
        # Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:, 0])], plddt2[np.unique(contacts[:, 1])]]))
        # Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt * np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018

        # PPV
        PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192,
                        0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
                        0.8919553, 0.88570037, 0.87822061, 0.87116417, 0.86040801,
                        0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
                        0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
                        0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
                        0.63555449, 0.55890174])

        pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
                                      0.60150931, 0.58313803, 0.5647381, 0.54122438, 0.52314392,
                                      0.49659878, 0.4774676, 0.44661346, 0.42628389, 0.39990988,
                                      0.38479715, 0.3649393, 0.34526004, 0.3262589, 0.31475668,
                                      0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
                                      0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
                                      0.06968505, 0.02946438])
        inds = np.argwhere(pdockq_thresholds >= pdockq)
        if len(inds) > 0:
            ppv = PPV[inds[-1]][0]
        else:
            ppv = PPV[0]

        meanpae = get_pae_region(contacts,l1,pae_file)
        
    return pdockq, PPV, meanpae


def get_max_iptm_ptm(rankingfile):
    # Check if the JSON file name is provided as an argument
    if len(sys.argv) < 2:
        print("Please provide the JSON file name as an argument.")
        sys.exit(1)

    # Extract the JSON file name from the command-line argument
    json_file = rankingfile

    # Read the JSON file
    with open(json_file) as file:
        data = json.load(file)

    # Get the "iptm+ptm" dictionary
    iptm_ptm_dict = data['iptm+ptm']

    # Find the maximum value in the dictionary
    max_value = max(iptm_ptm_dict.values())

    # Print the maximum value
    return max_value

def get_pae_region(contacts,l1,pae_file):
    import json
    json_file_path = pae_file
    
    # Open the JSON file and read its contents
    with open(json_file_path, "r") as json_file:
        data = json.load(json_file)

    predicted_aligned_error_sliced = pd.DataFrame(data[0]['predicted_aligned_error'])
    region_values=[]
    for pair in contacts:
    # Slice the heatmap data to get the values for the specified region
        pair_2 = [pair[0],pair[1]+l1]
        region_values.append(predicted_aligned_error_sliced.iloc[pair_2])
    return (np.mean(region_values))
    
#################MAIN####################

# Parse args
args = parser.parse_args()
# Read chains
chain_coords, chain_plddt = read_pdb(args.pdbfile[0])
# Check chains
#if len(chain_coords.keys()) < 2:
 #   print('Only one chain in pdbfile', args.pdbfile[0])
  #  sys.exit()

# Calculate pdockq
t = 8  # Distance threshold, set to 8 Å

#max_iptm_ptm = get_max_iptm_ptm(args.rankingfile[0])
pdb_dirname = os.path.dirname(args.pdbfile[0])

file_basename = os.path.basename(args.pdbfile[0])

# Split the basename using underscore as the delimiter and select the first part
pae_file = pdb_dirname + '/pae_'+('_'.join(file_basename.split('_')[1:])).split('.')[0]+'.json'

#print(pdb_dirname,pae_file)
pdockq, ppv, meanpae = calc_pdockq(chain_coords, chain_plddt, t, pae_file)

print(args.pdbfile[0],np.round(pdockq, 3), np.round(meanpae,3))#, np.round(max_iptm_ptm,3))


