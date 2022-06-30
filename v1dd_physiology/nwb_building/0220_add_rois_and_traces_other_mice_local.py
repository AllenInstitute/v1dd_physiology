import os
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.NwbTools as nt
import utils as utils

nwb_path = r"\\allen\programs\mindscope\workgroups\surround" \
           r"\v1dd_in_vivo_new_segmentation\data\nwbs" \
           r"\M416296_12_20181214.nwb"

curr_folder = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_folder)

nwb_folder, nwb_fn = os.path.split(nwb_path)
print(f'data path: {nwb_folder}')
print(f'adding rois and traces to {nwb_fn} ...')

sess_dict = utils.get_all_experiments_pika_meta(nwb_fn[:10])
print(sess_dict.keys())

# nwb_f = nt.RecordedFile(nwb_path)

for plane_n, plane_dict in sess_dict.items():

    print(f'\t{plane_n} ...')

    rois_traces_dict = utils.get_rois_and_traces_other_mice(
        exp_id=plane_dict['experiment'])

    plane_dict.update(rois_traces_dict)
    plane_dict.update({'plane_num' : len(sess_dict)})
    # print(plane_dict.keys())

    # utils.add_rois_and_traces_to_nwb(
    #     nwb_f=nwb_f, 
    #     plane_n=plane_n, 
    #     plane_dict=plane_dict
    #     )

# nwb_f.close()