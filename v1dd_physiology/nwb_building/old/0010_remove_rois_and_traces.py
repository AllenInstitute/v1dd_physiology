'''
this script removed all the rois and traces 
from the old segmentation
'''

import os
import h5py
import NeuroAnalysisTools.core.FileTools as ft

data_folder = r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\nwbs'

curr_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(curr_folder)

nwb_fns = ft.look_for_file_list(
    source=data_folder, 
    identifiers=['.nwb'], 
    file_type='nwb', 
    is_full_path=False
    )

# print('\n'.join(nwb_fns))

for nwb_i, nwb_fn in enumerate(nwb_fns):
    
    print(f'processing {nwb_fn}, {nwb_i+1} / {len(nwb_fns)} ...')

    ff = h5py.File(os.path.join(data_folder, nwb_fn), 'r+')

    grps_to_del = []
    for key in ff['processing'].keys():
        if key.startswith('l0') or key.startswith('rois'):
            grps_to_del.append(f'/processing/{key}')

    for grp in grps_to_del:
        del ff[grp]

    del ff['analysis/lims_roi_ids'] 

    ff.close()

