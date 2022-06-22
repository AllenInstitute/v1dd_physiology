import os
import NeuroAnalysisTools.NwbTools as nt
import NeuroAnalysisTools.core.FileTools as ft
import utils as utils

data_folder = r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\nwbs'

fns = ft.look_for_file_list(
    source=data_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

# fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

for fn in fns:

    sess_folder = utils.get_lims_session_path_from_session_name(fn[0:10])

    # print(sess_folder)

    sync_path = ft.look_for_unique_file(sess_folder, identifiers=['_sync.h5'],
                                        file_type='h5', is_full_path=True)

    nwb_path = os.path.join(data_folder, fn)
    nwb_f = nt.RecordedFile(nwb_path)
    nwb_f.add_sync_data(f_path=sync_path, by_label=True)
    nwb_f.close()