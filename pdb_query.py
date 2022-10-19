
import requests
import json
import os
import pypdb
from datetime import date


# r<=3; x-ray; has ligand; has ki/kd; protein only#
def get_pdbids():
    json_ = '%7B%0A%20%20%22query%22%3A%20%7B%0A%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%22nodes%22%3A%20%5B%0A%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%20%20%20%20%22nodes%22%3A%20%5B%0A%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22nodes%22%3A%20%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22service%22%3A%20%22text%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22exptl.method%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22exact_match%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%20%22X-RAY%20DIFFRACTION%22%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22service%22%3A%20%22text%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22rcsb_entry_info.diffrn_resolution_high.value%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22less_or_equal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%203%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22logical_operator%22%3A%20%22and%22%0A%20%20%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22nodes%22%3A%20%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22nodes%22%3A%20%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22service%22%3A%20%22text%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22rcsb_binding_affinity.value%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22exists%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22service%22%3A%20%22text%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22rcsb_binding_affinity.type%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22exact_match%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%20%22Kd%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22label%22%3A%20%22nested-attribute%22%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22group%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22logical_operator%22%3A%20%22and%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22nodes%22%3A%20%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22service%22%3A%20%22text%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22rcsb_binding_affinity.value%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22exists%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22service%22%3A%20%22text%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22rcsb_binding_affinity.type%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22exact_match%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%20%22Ki%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22label%22%3A%20%22nested-attribute%22%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5D%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22logical_operator%22%3A%20%22or%22%0A%20%20%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22service%22%3A%20%22text%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22rcsb_entry_info.selected_polymer_entity_types%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22exact_match%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%20%22Protein%20(only)%22%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%5D%2C%0A%20%20%20%20%20%20%20%20%22label%22%3A%20%22text%22%0A%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%22type%22%3A%20%22terminal%22%2C%0A%20%20%20%20%20%20%20%20%22service%22%3A%20%22text_chem%22%2C%0A%20%20%20%20%20%20%20%20%22parameters%22%3A%20%7B%0A%20%20%20%20%20%20%20%20%20%20%22attribute%22%3A%20%22rcsb_chem_comp_container_identifiers.comp_id%22%2C%0A%20%20%20%20%20%20%20%20%20%20%22operator%22%3A%20%22exists%22%2C%0A%20%20%20%20%20%20%20%20%20%20%22negation%22%3A%20false%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%5D%0A%20%20%7D%2C%0A%20%20%22return_type%22%3A%20%22entry%22%2C%0A%20%20%22request_options%22%3A%20%7B%0A%20%20%20%20%22return_all_hits%22%3A%20true%0A%20%20%7D%0A%7D'
    url = f'https://search.rcsb.org/rcsbsearch/v2/query?json={json_}'
    r = requests.get(url).content
    my_json = r.decode('utf8').replace("'", '"')
    data = json.loads(my_json)['result_set']
    pdb_ids = [entry['identifier'] for entry in data]
    file_name = str(date.today()).replace('-', '') + '.txt'
    file_path = '/Users/xiaotongxu/structural_datasets/pdb'
    with open(os.path.join(file_path, file_name), 'w') as f:
        for id in pdb_ids:
            f.write("%s\n" % id)
    return pdb_ids

#cross check with PDBBind

def checking(pdb_ids):
    pdbbind_path = '/Users/xiaotongxu/structural_datasets/data_summary'
    pdbbind_file = 'INDEX_general_PL_data.2020'
    pdb_bind = []
    with open((os.path.join(pdbbind_path, pdbbind_file)),'r') as f:
        lines = f.readlines()
        for line in lines:
            if '#' not in line:
                pdb_id = line.split()[0].upper()
                pdb_bind.append(pdb_id)
    not_in = []
    for id in pdb_ids:
        if id not in pdb_bind:
            not_in.append(id)
    return not_in, pdb_bind


def get_activity(overlap):
    pdbbind_path = '/Users/xiaotongxu/structural_datasets/data_summary'
    pdbbind_file = 'INDEX_general_PL_data.2020'
    file_name = str(date.today()).replace('-', '') + '_activity.txt'
    file_path = '/Users/xiaotongxu/structural_datasets/pdb'
    f_w = open(os.path.join(file_path, file_name),"wt")
    with open((os.path.join(pdbbind_path, pdbbind_file)),'r') as f:
        lines = f.readlines()
        for line in lines:
            pdb_id = line.split()[0].upper()
            if pdb_id in overlap:
                f_w.write(str(line))

def get_activity(not_in):
    multi_acctivity= []
    string_list = []
    for pdb_id in not_in:
        describe = pypdb.get_info(pdb_id)
        resolution = describe['pdbx_vrpt_summary']['pdbresolution']
        affinity_info = describe['rcsb_binding_affinity']
        if len(affinity_info) == 1:
            affinity_info = affinity_info[0]
            affinity = affinity_info['value']
            type = affinity_info['type']
            unit = affinity_info['unit']
            comp_id = affinity_info['comp_id']
            string = f'{pdb_id} {resolution} {affinity} {type} {unit} {comp_id}'
            string_list.append(string)
        else:
            print(pdb_id)
            multi_acctivity.append(pdb_id)
    
    file_name = str(date.today()).replace('-', '') + '_clean_activity.txt'
    file_path = '/Users/xiaotongxu/structural_datasets/pdb'
    with open(os.path.join(file_path, file_name), 'w') as f:
        for string in string_list:
            f.write("%s\n" % string)
    return multi_acctivity






pdb_ids = get_pdbids()
not_in, pdb_bind = checking(pdb_ids)
overlap  = list(set(pdb_ids)^set(not_in))
get_activity(overlap)
multi_acctivity = get_activity(not_in)
#print(len(a))
#print(a[1:4])
#get_pdbids()