'''
this script repacks all nwbs to save disk space 
'''

import os
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

for nwb_i, nwb_fn in enumerate(nwb_fns):
    
    print(f'processing {nwb_fn}, {nwb_i+1} / {len(nwb_fns)} ...')

    nwb_fn_new = f'{nwb_fn[:-8]}.nwb'

    os.system(f'h5repack {os.path.join(data_folder, nwb_fn)} {os.path.join(data_folder, nwb_fn_new)}')

    