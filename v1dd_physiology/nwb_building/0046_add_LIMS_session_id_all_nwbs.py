import os, h5py
import NeuroAnalysisTools.core.FileTools as ft
import utils as utils

data_folder = r'\\allen\programs\mindscope\workgroups' \
              r'\surround\v1dd_in_vivo_new_segmentation\data\nwbs'

fns = ft.look_for_file_list(
    source=data_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )
fns.sort()
print('\n'.join(fns))

for fn in fns:
    mouse_id = fn[1:7]
    sess_n = fn[8:10]

    meta = utils.get_all_sessions(mouse_id=mouse_id)
    sess_id = meta[sess_n]['session']
    print(f'{fn}: {sess_id}')

    ff = h5py.File(os.path.join(data_folder, fn), 'a')
    if 'session_id' in ff['general']:
        del ff['general/session_id']
    ff['general/session_id'] = sess_id
    ff.close()

