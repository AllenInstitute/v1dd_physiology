import os
import json
import datetime
import NeuroAnalysisTools.PreprocessPipeline as pp
import NeuroAnalysisTools.core.FileTools as ft
import utils as utils

data_folder = r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\nwbs'

curr_folder = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_folder)

psr = pp.Preprocessor()

meta_path = os.path.join(
    os.path.realpath(curr_folder), 
    'meta_lims', 
    'deepdive_EM_volumes_mouse_ids.json'
    )

with open(meta_path, 'r') as f:
    meta = json.load(f)

mouse_ids = list(meta.keys())
mouse_ids.remove('slc1') # this mouse is excluded

for mid in mouse_ids:
    
    print(f'\nprocessing mouse_id: {mid} ...')

    metam = utils.get_meta_mouse(mid)
    dob = metam['date_of_birth']
    prod = metam['prod']

    sesses = utils.get_all_sessions(metam['mouse_id'])
    sesses_2p = utils.pick_2p_sessions(sesses)

    print(sesses_2p)

    for sn, smeta in sesses_2p.items():
        print(f'\tadding session: {sn} ...')
        
        lims_path = utils.get_lims_session_path(
            sess=smeta, 
            prod=prod, 
            btv_path=r"\\allen\programs\braintv"
            )
        print(f'\tsess path: {lims_path}')

        date = utils.get_session_date(lims_path)
        age = (datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8])) -
               datetime.date(int(dob[0:4]), int(dob[4:6]), int(dob[6:8]))).days

        # plane_depth = smeta['depth']

        psr.generate_nwb_file(
            save_folder=data_folder,
            date=date,
            mouse_id=metam['mouse_id'],
            session_id=sn,
            experimenter='',
            genotype=metam['genotype'], 
            sex=metam['gender'],
            age=age,
            indicator='GCaMP6s',
            imaging_rate='6.6',
            imaging_depth='multi',
            imaging_location='primary visual cortex',
            imaging_device='two-photon scope',
            imaging_excitation_lambda='940 nm',
            is_date_back=True)
