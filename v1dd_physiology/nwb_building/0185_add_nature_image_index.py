import os
import h5py
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft
import pandas as pd

nwb_folder = r'\\allen\programs\mindscope\workgroups\surround' \
             r'\v1dd_in_vivo_new_segmentation\data\nwbs'
index_fp = r"\\allen\programs\mindscope\workgroups\surround" \
           r"\v1dd_in_vivo_new_segmentation\data\stim_movies" \
           r"\natural_images_order.xlsx"

curr_folder = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_folder)

ni_12 = pd.read_excel(io=index_fp, sheet_name='natural_images_short')
ni_12_ind = np.array(ni_12['img_id_0_base'], dtype=np.uint)
ni_12_des = list(ni_12['description'])
ni_12_des = [n.encode('utf-8') for n in ni_12_des]

ni_236 = pd.read_excel(io=index_fp, sheet_name='natural_images_long')
ni_236_ind = np.array(ni_236['img_id_0_base'], dtype=np.uint)
ni_236_des = list(ni_236['description'])
ni_236_des = [n.encode('utf-8') for n in ni_236_des]

# print(ni_12)
# print(ni_236)

nwb_fns = ft.look_for_file_list(source=nwb_folder,
                                file_type='nwb',
                                identifiers=[''],
                                is_full_path=False)
nwb_fns.sort()


for nwb_i, nwb_fn in enumerate(nwb_fns):

    print(f'processing {nwb_fn}, {nwb_i+1} / {len(nwb_fns)} ...')
    nwb_f = h5py.File(os.path.join(nwb_folder, nwb_fn), 'a')
    ni12_grp = nwb_f['stimulus/presentation/natural_images_12']
    ni12_ind_dset = ni12_grp.create_dataset('image_index', data=ni_12_ind)
    ni12_ind_dset.attrs['neurodata_type'] = 'Custom'
    ni12_des_dset = ni12_grp.create_dataset('image_description', data=ni_12_des)
    ni12_des_dset.attrs['neurodata_type'] = 'Custom'

    ni236_grp = nwb_f['stimulus/presentation/natural_images']
    ni236_ind_dset = ni236_grp.create_dataset('image_index', data=ni_236_ind)
    ni236_ind_dset.attrs['neurodata_type'] = 'Custom'
    ni236_des_dset = ni236_grp.create_dataset('image_description', data=ni_236_des)
    ni236_des_dset.attrs['neurodata_type'] = 'Custom'

    nwb_f.close()