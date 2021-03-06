import os
import shutil
import NeuroAnalysisTools.core.FileTools as ft

sess_to_move = [
    'M427836_1e', # no stimulus table
    'M427836_1f', # no stimulus table
    # 'M438833_1e', # not the right depth
    # 'M438833_1f', # not the right depth
]

nwb_folder = r"\\allen\programs\mindscope\workgroups\surround" \
             r"\v1dd_in_vivo_new_segmentation\data\nwbs"

save_folder = os.path.join(nwb_folder, 'problematic')

if not os.path.isdir(save_folder):
    os.makedirs(save_folder)

for sess in sess_to_move:

    fn = ft.look_for_unique_file(
        source=nwb_folder,
        identifiers=[sess],
        file_type='nwb',
        is_full_path=False)

    print(fn)

    if fn is not None: 

        file_path = os.path.join(nwb_folder, fn)
        save_path = os.path.join(save_folder, fn)

        if os.path.isfile(save_path):
            os.remove(save_path)

        shutil.move(file_path, save_path)