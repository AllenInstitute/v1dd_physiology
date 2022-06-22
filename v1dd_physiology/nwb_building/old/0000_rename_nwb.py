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

print('\n'.join(nwb_fns))

for nwb_fn in nwb_fns:
    os.rename(
        os.path.join(data_folder, nwb_fn), 
        os.path.join(data_folder, f'{nwb_fn[:-4]}_old.nwb'))

