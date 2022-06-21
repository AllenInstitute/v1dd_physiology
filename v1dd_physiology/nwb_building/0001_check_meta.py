import os, json
import NeuroAnalysisTools.core.FileTools as ft

'''
the original meta json file received from team pika is saved 
in the "meta_lims/meta_pika_original" folder.

the "meta_pika_xxxxxx.json" files in "meta_lims" folder is 
the curated version of these files. 

the modificaitons are fixing typo, depth mismatch, removing 
problematic sessions, standardize session/experiment names, etc.
'''

script_folder = os.path.dirname(os.path.realpath(__file__))
meta_folder = os.path.join(script_folder, 'meta_lims')

depth_pairs = {
    "6" : 500,
    "7" : 525,
    "8" : 550,
    "9" : 575,
    "a" : 600,
    "b" : 625,
    "c" : 650,
    "d" : 675,
    "e" : 700,
    "f" : 725
}

meta_fns = ft.look_for_file_list(
    source=meta_folder,
    identifiers=['meta_pika_'],
    file_type='json',
    is_full_path=False)
meta_fns.sort()
print('\n'.join(meta_fns))
print()

for meta_fn in meta_fns:

    print(f'\nanalyzing {meta_fn} ...')

    mouse_id = meta_fn[10:16]
    meta_path = os.path.join(meta_folder, meta_fn)

    with open(meta_path, 'r') as f:
        meta = json.load(f)

    for exp_id, exp_dict in meta.items():

        assert(str(exp_dict['external_specimen_name']) == mouse_id)

        assert(exp_dict['oe'].split('_')[1] == mouse_id)

        assert(int(exp_dict['oe'].split('_')[-1]) == exp_dict['imaging_depth'])

        if mouse_id != "409828":
            assert(exp_dict['os'] == exp_dict['oe'][:25])
        
        col = exp_dict['oe'][19]
        vol = exp_dict['oe'][24]

        if vol in ["1", "2", "3", "4", "5"]:
            vol = int(vol)
            assert(exp_dict["imaging_depth"] > 50 + (vol-1) * 96 -25)
            assert(exp_dict["imaging_depth"] < 50 + vol * 96)
        else:
            assert(col == '1')
            assert(exp_dict["imaging_depth"] == depth_pairs[vol])
        