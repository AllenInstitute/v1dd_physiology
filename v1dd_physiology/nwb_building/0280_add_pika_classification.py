import os, h5py
import pandas as pd
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft

pred_path = r"\\allen\programs\mindscope\workgroups\surround" \
            r"\v1dd_in_vivo_new_segmentation\classifier_results" \
            r"\classifications.csv"

nwb_folder = r"\\allen\programs\mindscope\workgroups\surround" \
           r"\v1dd_in_vivo_new_segmentation\data\nwbs"    

curr_folder = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_folder)

df_pred = pd.read_csv(pred_path)
print(df_pred)        

nwb_paths = ft.look_for_file_list(
    source=nwb_folder,
    identifiers=[],
    file_type='nwb',
    is_full_path=True)
nwb_paths.sort()
# print('\n'.join(nwb_paths))

for nwb_i, nwb_p in enumerate(nwb_paths):

    _, nwb_fn = os.path.split(nwb_p)

    print(f'processing {nwb_fn}, {nwb_i + 1}/{len(nwb_paths)} ...')

    nwb_f = h5py.File(nwb_p, 'a')

    if 'roi_classification_pika' in nwb_f['analysis']:
        del nwb_f['analysis/roi_classification_pika']

    pika_grp = nwb_f['analysis'].create_group('roi_classification_pika')
    pika_grp.attrs['description'] = '''This group stores the results from 
    a neural network classifier trained by team Pika. The classifier predicts 
    if an roi is a valid cell body or not. The dataset 'score' is the 
    confidence score generated by the classifier and the dataset 'prediction' 
    is the boolean marker using 0.5 as threshold on the 'score'.
    '''

    plane_keys = [k for k in nwb_f['processing'].keys() if k.startswith('rois_and_traces_plane')]
    plane_ns = [f'plane{pi}' for pi in range(len(plane_keys))]

    for plane_n in plane_ns:

        plane_grp = nwb_f[f'processing/rois_and_traces_{plane_n}']
        exp_id = plane_grp['experiment_id'][()].decode()
        print(f'\t{plane_n}: {exp_id}')

        if 'roi_list' in plane_grp['ImageSegmentation/imaging_plane']:

            roi_names_p = plane_grp['ImageSegmentation/pipeline_roi_names']
            roi_names = plane_grp['ImageSegmentation/imaging_plane/roi_list']

            df_pred_plane = df_pred[df_pred['experiment_id'] == int(exp_id)].copy()
            df_pred_plane = df_pred_plane.sort_values(by='roi_id')

            assert(df_pred_plane.shape[0] == len(roi_names) == len(roi_names_p))

            pika_score = np.array(df_pred_plane['y_score'])
            pika_pred = np.array(df_pred_plane['y_pred'] == 'cell')
        
        else: # no roi

            roi_names_p = np.array([])
            roi_names = np.array([])
            pika_score = np.array([])
            pika_pred = np.array([])

        pika_plane_grp = pika_grp.create_group(plane_n)
        pika_plane_grp['pipeline_roi_names'] = roi_names_p
        pika_plane_grp['roi_names'] = roi_names
        pika_plane_grp['score'] = pika_score
        pika_plane_grp['prediction'] = pika_pred
        pika_plane_grp['prediction'].attrs['threshold'] = 0.5

