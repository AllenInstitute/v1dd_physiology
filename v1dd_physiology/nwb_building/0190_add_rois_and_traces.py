import sys, os
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.NwbTools as nt
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)
import utils as utils

nwb_path = sys.argv[1]
nwb_folder, nwb_fn = os.path.split(nwb_path)
print(f'data path: {nwb_folder}')
print(f'adding rois and traces to {nwb_fn} ...')

sess_dict = utils.get_all_experiments_pika_meta(nwb_fn[:10])

# nwb_f = nt.RecordedFile(os.path.join(nwb_folder, fn))

for plane_n, plane_dict in sess_dict.items():

    rois_traces_dict = utils.get_rois_and_traces(
        exp_id=plane_dict['experiment'],
        folder_noisy='/allen/programs/mindscope/workgroups/surround/motion_correction_labeling_2022',
        folder_denoised='/allen/programs/mindscope/workgroups/surround/denoising_labeling_2022/denoised_movies',
        folder_roi='/allen/programs/mindscope/workgroups/surround/denoising_labeling_2022/segmentations',
        folder_noisy_trace_raw='/allen/programs/mindscope/workgroups/surround/trace_extraction_2022/traces_2022',
        folder_noisy_trace_demixed='/allen/programs/mindscope/workgroups/surround/trace_extraction_2022/demix_2022',
        folder_noisy_trace_subtracted='/allen/programs/mindscope/workgroups/surround/trace_extraction_2022/neuropil_2022',
        folder_noisy_trace_dff='/allen/programs/mindscope/workgroups/surround/trace_extraction_2022/dff_2022',
        )

    plane_dict.update(rois_traces_dict)
    print(plane_dict.keys())

    # utils.add_rois_and_traces_to_nwb(
    #     nwb_f=nwb_f, 
    #     plane_n=plane_n, 
    #     plane_dict=plane_dict
    #     )

# nwb_f.close()

