import os
import pandas as pd
import NeuroAnalysisTools.core.FileTools as ft
import utils as utils

nwb_folder = r"\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\nwbs"

curr_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(curr_folder)

fns = ft.look_for_file_list(
    source=nwb_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

paths = []

for fn in fns:
    mid = fn.split('_')[0][1:]
    sess_id = fn[8:10]

    meta_sess = utils.get_all_sessions(mouse_id=mid)
    sess_dict = meta_sess[sess_id]
    specimen = 'specimen_{}'.format(sess_dict['specimen'])
    session = 'session_{}'.format(sess_dict['session'])
    sess_path = os.path.join(r'\\allen\programs\braintv\workgroups\v1column\rezaa\deepdive_volumes\v1DD_EM',
                             specimen, session)
    # print(sess_path)
    paths.append([fn[0:10], sess_path])

paths_df = pd.DataFrame(data=paths, columns=['sess_name', 'sess_path'])
print(paths_df)

paths_df.to_csv(os.path.join(r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\stimulus_tables', 
                             '3p_sess_paths.csv'))