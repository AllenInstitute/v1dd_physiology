'''
run by 
/allen/aibs/mat/Peter/pet-python-36/bin/python 
'''

import sys
code_base_dir = '/allen/programs/braintv/workgroups/cortexmodels/peterl/forward_model' #IF YOU MOVE THE CODE, PLEASE EDIT THIS AS WELL!!!
sys.path.extend([code_base_dir])
import numpy as np
from scipy.signal import resample_poly
import h5py
import os
import sys
from joblib import Parallel, delayed
from importlib.machinery import SourceFileLoader
l0 = SourceFileLoader(
    'L0_analysis',
    '/allen/programs/braintv/workgroups'
    '/cortexmodels/peterl/forward_model/l0_analysis_deepscope.py').load_module()

output_dir = '/allen/programs/mindscope/workgroups/surround' \
             '/v1dd_in_vivo_new_segmentation/data/event_detection'

curr_folder = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_folder)

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

nwb_path = sys.argv[1]
print(f'processing {nwb_path} ...')
f = h5py.File(nwb_path, 'r')

# # example nwb_file
# nwb_path = '/allen/programs/mindscope/workgroups/surround' \
#            '/v1dd_in_vivo_new_segmentation/data/nwbs/M409828_1a_20180627.nwb'

nwb_folder, nwb_fn = os.path.split(os.path.abspath(nwb_path))

output_prefix = os.path.splitext(nwb_fn)[0]

#loop over all 6 planes. Upsample by a factor of 5 to get ~30Hz...
plane_keys = [p for p in f['processing'].keys() if p.startswith('rois_and_traces_plane')]
plane_num = len(plane_keys)

for i in range(plane_num):

    output_fn = f'{output_prefix}_{i}.npz'
    output_path = os.path.join(output_dir, output_fn)
    if os.path.isfile(output_path):
        print('file {} already exists. skip.'.format(output_fn))
    else:
        dff = f['processing'][f'rois_and_traces_plane{i}']['DfOverF']['dff_raw']['data'][()]
        ts = f['processing'][f'rois_and_traces_plane{i}']['DfOverF']['dff_raw']['timestamps'][()]

        if len(dff) == 1:
            np.savez_compressed(output_path, dff=np.array([[np.nan]]), 
                ts=np.array([np.nan]),  events=np.array([[np.nan]]), 
                noise_stds=np.array([np.nan]), lambdas=np.array([np.nan]))
        else:
            dff30Hz = resample_poly(dff, 5, 1, axis=1)
            ts30Hz = resample_poly(ts, 5, 1)
            l0a = l0.L0_analysis(dff30Hz)  # needs to be array of arrays
            events = l0a.get_events()
            np.savez_compressed(
                output_path, dff=dff30Hz, ts=ts30Hz,  events=events, 
                noise_stds=l0a._noise_stds, lambdas=l0a.lambdas)
