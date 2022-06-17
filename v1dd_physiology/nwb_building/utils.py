import os, h5py, json
import numpy as np
import tifffile as tf
import matplotlib.pyplot as plt
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.core.ImageAnalysis as ia


script_folder = os.path.dirname(os.path.realpath(__file__))


def get_meta_mouse(mouse_id):
    """
    return a dictionary of the meta data of the input mouse id 
    """

    path = os.path.join(
        script_folder, 
        'meta_lims', 
        'deepdive_EM_volumes_mouse_ids.json'
        )
    with open(path, 'r') as f:
        meta = json.load(f)

    _ = meta.pop('slc1') # this mouse is excluded

    if len(str(mouse_id)) < 6:
        metam = meta[mouse_id]
        metam['mouse_name'] = mouse_id
    else:
        metam = {}
        for key, value in meta.items():
            if value['mouse_id'] == str(mouse_id):
                metam = value
                metam['mouse_name'] = key
    return metam


def get_all_sessions(mouse_id):
    """
    return a dictionary of the meta data of all sessions from 
    the input mouse_id (six digit string)
    """

    json_folder = os.path.join(script_folder, 'meta_lims')

    json_path = ft.look_for_unique_file(
        source=json_folder, 
        identifiers=['deepdive_EM_volumes', mouse_id], 
        file_type='json', 
        print_prefix='', 
        is_full_path=True,
        is_verbose=True)

    # print(json_path)

    with open(json_path, 'r') as f:
        sesses = json.load(f)

    return sesses


def get_all_experiments_pika_meta(sess_n):
    """
    given a session name, e.g. M409828_11
    return the metadata of all experiments (planes) 
    in this session 

    This uses the modified pika meta file

    Parameters
    ----------
    sess_n : string
        length should be 10, e.g. M409828_11

    Returns
    -------
    sess_meta : dictionary
        {plane_n : 
            {
            'experiment' : experiment_id, string,
            'depth': int, plane depth in microns
            }
        }
    """
    json_folder = os.path.join(script_folder, 'meta_lims')

    json_path = ft.look_for_unique_file(
        source=json_folder, 
        identifiers=['meta_pika', sess_n[1:7]], 
        file_type='json', 
        print_prefix='', 
        is_full_path=True,
        is_verbose=True)

    with open(json_path, 'r') as f:
        sesses = json.load(f)

    sess_dict = {}
    depths = []
    for exp_i, exp_dict in sesses.items():
        oe_parts = exp_dict["oe"].split('_')
        col = oe_parts[2][3]
        vol = oe_parts[3][3]
        if (col == sess_n[8]) and (vol == sess_n[9]):
            sess_dict.update({exp_i: exp_dict})
            depths.append(int(oe_parts[-1]))

    depths.sort()
    
    sess_meta = {}
    for dep_i, depth in enumerate(depths):
        for exp_i, exp_dict in sess_dict.items():
            oe_parts = exp_dict["oe"].split('_')
            if int(oe_parts[-1]) == depth:
                sess_meta.update({f'plane{dep_i}': 
                                    {'experiment': exp_i,
                                     'depth': depth}})

    return sess_meta


def pick_3p_sessions(sesses):
    """
    from the dictonary of all session pick the three-photon sessions.
    """

    sesses_3p = {}
    for sess_n, sess in sesses.items():
        if len(sess_n) == 2 and \
        sess_n[0] == '1' and \
        sess_n[1] not in ['1', '2', '3', '4', '5']:
            
            sesses_3p.update({sess_n: sess})

    return sesses_3p


def get_lims_session_path(
    sess, prod, btv_path=r"\\allen\programs\braintv"
    ):
    """
    given session meta dictionary 
    return the lims path to the data of this session
    """

    specimen = f'specimen_{sess["specimen"]}'
    session = f'ophys_session_{sess["session"]}'
    prod = f'prod{prod}'

    sess_path = os.path.join(btv_path, 'production/neuralcoding', prod, specimen, session)

    return sess_path


def get_session_date(session_path):
    """
    given lims session_path as string, return a string of experiment date as a string "yyyymmdd"
    """
    platform_path = ft.look_for_unique_file(session_path, identifiers=['_platform.json'],
                                            file_type='json', is_full_path=True)
    with open(platform_path, 'r') as f:
        dict_platform = json.load(f)

    date = dict_platform['registration']['reticle_image']['acquired_at']
    date = date[0:4] + date[5:7] + date[8:10]
    return date


def get_lims_session_path_from_session_name(
    session_name,
    btv_path=r"\\allen\programs\braintv"
    ):
    """
    given session name (e.g. 'M409828_1a') return 
    lims path to the session
    """
    assert(len(session_name) == 10)

    mouse_id = session_name[1:7]
    key = session_name[-2:]

    meta_mouse = get_meta_mouse(mouse_id=mouse_id)
    meta_sess = get_all_sessions(mouse_id=mouse_id)
    sess_dict = meta_sess[key]
    specimen = 'specimen_{}'.format(sess_dict['specimen'])
    session = 'ophys_session_{}'.format(sess_dict['session'])
    prod = 'prod{}'.format(meta_mouse['prod'])

    sess_path = os.path.join(btv_path, 'production/neuralcoding', prod, specimen, session)

    return sess_path


def get_vasmap(
    session_name, 
    btv_path=r"\\allen\programs\braintv"
    ):

    sess_path = get_lims_session_path_from_session_name(
        session_name=session_name,
        btv_path=btv_path
        )

    vasmap_wf_path = ft.look_for_unique_file(source=sess_path,
                                             identifiers=['vasculature.tif'],
                                             file_type='tif',
                                             is_full_path=True)

    plane_fns = [f for f in os.listdir(sess_path) if f[0:17] == 'ophys_experiment_' and
                 os.path.isdir(os.path.join(sess_path, f))]
    plane_fns.sort()
    plane_path = os.path.join(sess_path, plane_fns[0])

    vasmap_2p_path = ft.look_for_unique_file(source=plane_path,
                                             identifiers=['averaged_surface.tif'],
                                             file_type='tif',
                                             is_full_path=True)

    vasmap_wf = tf.imread(vasmap_wf_path).squeeze()
    vasmap_2p = tf.imread(vasmap_2p_path)

    return vasmap_wf, vasmap_2p


def analyze_lsn_movie(mov, grid_size=9.3):
    """
    give the locally sparse noise movie for displaying,
    return, at each location of either sign, the frame indices
    that the probe is on. assuming the center of the frame is (0, 0)
    (altitude, azimuth), top: higher altitude, left: lower azimuth
    :param mov: 3d array, movie of sparse noise. each
                probe, is one pixel
    :param grid_size: float, grid size in degrees
    :return probe_list: list, each element represents a probe
                        at a give location and sign. This element
                        is a list it self:
                        (altitude,
                         azimuth,
                         sign,
                         list of frame indices that the probe is on)

    """

    azis = (np.arange(mov.shape[2]) - ((mov.shape[2] - 1) / 2)) * grid_size
    alts = (np.arange(mov.shape[1]) - ((mov.shape[1] - 1) / 2)) * grid_size
    alts = alts[::-1]

    # print(alts)
    # print(azis)

    probe_list = []

    for i in range(mov.shape[1]):
        for j in range(mov.shape[2]):

            azi = azis[j]
            alt = alts[i]

            pix = mov[:, i, j]

            on = np.where(pix == 1)
            probe_list.append([alt, azi, 1, list(on[0])])

            off = np.where(pix == -1)
            probe_list.append([alt, azi, -1, list(off[0])])

            # print(probe_list)

    return probe_list


def get_roi_mask_list_from_pipeline_json(json_path, img_h=512, img_w=512):
    """
    given the roi_dictionary from the "<experiment_id>_rois.json" file

    return roi_list: list
         each element is a NeuroAnalysisTools.core.ImageAnalysis.ROI object
    return roi_id_list: 1d array
         each element is a integer as roi_id
    """

    with open(json_path, 'r') as f:
        dict_all = json.load(f)

    roi_list = []
    roi_id_list = []
    curr_id = None

    for roi_i, roi_dict in enumerate(dict_all):

        # check if the id is incremental which is important to 
        # match rois with their traces
        if curr_id is None:
            curr_id = roi_dict['id']
        else:
            if roi_dict['id'] <= curr_id:
                raise ValueError('the roi ids in this json file is not incremental.')
            else:
                curr_id = roi_dict['id']
        
        # generate mask
        mask = np.zeros((img_h, img_w), dtype=np.uint8)
        mask_small = np.array(roi_dict['mask_matrix'], dtype=np.uint8)
        y = roi_dict['y']
        x = roi_dict['x']
        mask[y: y + mask_small.shape[0],
             x: x + mask_small.shape[1]] = mask_small

        roi_list.append(ia.ROI(mask))
        roi_id_list.append(roi_dict['id'])

    return roi_list, np.array(roi_id_list)


def get_traces_from_pipeline_h5(h5_path):
    """
    get traces from h5 files generated in pipeline (from team Pika)

    :return roi_ids: 1d array, ids for each roi, should be incremental
    :return traces: 2d array, roi x time
        the order of traces should match the roi_ids

    """
    ff = h5py.File(h5_path, 'r')
    roi_ids = ff['roi_names'][()]
    traces = ff['data'][()]
    roi_ids = np.array([int(r.decode()) for r in roi_ids])  
    
    # sort roi ids and traces
    sort_ind = np.argsort(roi_ids)
    roi_ids = roi_ids[sort_ind]
    traces = traces[sort_ind, :]

    return roi_ids, traces


def get_neuropil_subtracted_traces_from_pipeline_h5(h5_path):

    ff = h5py.File(h5_path, 'r')
    roi_ids = ff['roi_names'][()]
    traces = ff['FC'][()]
    rs = ff['r'][()]
    rmses = ff['RMSE'][()]
    roi_ids = np.array([int(r.decode()) for r in roi_ids]) 

    # sort roi ids, traces, rs and rmses
    sort_ind = np.argsort(roi_ids)
    roi_ids = roi_ids[sort_ind]
    traces = traces[sort_ind, :]
    rs = rs[sort_ind]
    rmses = rmses[sort_ind]

    return roi_ids, traces, rs, rmses


def get_dff_traces_from_pipeline_h5(h5_path):

    ff = h5py.File(h5_path, 'r')
    roi_ids = ff['roi_names'][()]
    traces = ff['data'][()]
    num_frame_small_baseline = ff['num_small_baseline_frames'][()]
    dff_sigmas = ff['sigma_dff'][()]
    roi_ids = np.array([int(r.decode()) for r in roi_ids])  

    # sort roi ids, traces, num_frame_small_baseline and sigmas_dff
    sort_ind = np.argsort(roi_ids)
    roi_ids = roi_ids[sort_ind]
    traces = traces[sort_ind, :]
    num_frame_small_baseline = num_frame_small_baseline[sort_ind]
    dff_sigmas = dff_sigmas[sort_ind]

    return roi_ids, traces, num_frame_small_baseline, dff_sigmas


def get_rois_and_traces(
    exp_id,
    folder_noisy=r'\\allen\programs\mindscope\workgroups\surround\motion_correction_labeling_2022',
    folder_denoised=r'\\allen\programs\mindscope\workgroups\surround\denoising_labeling_2022\denoised_movies',
    folder_roi=r'\\allen\programs\mindscope\workgroups\surround\denoising_labeling_2022\segmentations',
    folder_noisy_trace_raw=r'\\allen\programs\mindscope\workgroups\surround\trace_extraction_2022\traces_2022',
    folder_noisy_trace_demixed=r'\\allen\programs\mindscope\workgroups\surround\trace_extraction_2022\demix_2022',
    folder_noisy_trace_subtracted=r'\\allen\programs\mindscope\workgroups\surround\trace_extraction_2022\neuropil_2022',
    folder_noisy_trace_dff=r'\\allen\programs\mindscope\workgroups\surround\trace_extraction_2022\dff_2022',
    ):
    
    """
    Collect rois and traces information
    """

    plane_dict = {}

    # ================= noisy projection ========================
    proj_max_noisy_path = os.path.join(
        folder_noisy, exp_id, f'{exp_id}_filtered_max.h5')
    proj_max_noisy = h5py.File(proj_max_noisy_path, 'r')['max_projection'][()]

    proj_avg_noisy_path = os.path.join(
        folder_noisy, exp_id, f'{exp_id}_avg.h5')
    proj_avg_noisy = h5py.File(proj_avg_noisy_path, 'r')['avg_projection'][()]
    
    plane_dict.update({
        'projection_max_noisy' : proj_max_noisy,
        'projection_avg_noisy' : proj_avg_noisy
        })
    # print('aaaa')
    # ================= noisy projection ========================


    # ================= denoised projection ========================
    proj_max_denoi_path = os.path.join(
        folder_denoised, exp_id, f'{exp_id}_max.h5')
    proj_max_denoi = h5py.File(proj_max_denoi_path, 'r')['max_projeciton'][()] # matching the typo

    proj_avg_denoi_path = os.path.join(
        folder_denoised, exp_id, f'{exp_id}_avg.h5')
    proj_avg_denoi = h5py.File(proj_avg_denoi_path, 'r')['avg_projection'][()]

    proj_cor_denoi_path = os.path.join(
        folder_denoised, exp_id, f'{exp_id}_corr.h5')
    proj_cor_denoi = h5py.File(proj_cor_denoi_path, 'r')['corr_projection'][()]

    plane_dict.update({
        'projection_max_denoised' : proj_max_denoi,
        'projection_avg_denoised' : proj_avg_denoi,
        'projection_cor_denoised' : proj_cor_denoi
        })
    # print('bbbb')
    # ================= denoised projection ========================

    # f, axs = plt.subplots(ncols=3, nrows=2, figsize=(10, 6))
    # f.suptitle(f'experiment: {exp_id}')
    # axs[0, 0].imshow(proj_avg_noisy)
    # axs[0, 0].set_title('noisy mean')
    # axs[0, 1].imshow(proj_max_noisy)
    # axs[0, 1].set_title('noisy max')
    # axs[1, 0].imshow(proj_avg_denoi)
    # axs[1, 0].set_title('denoised mean')
    # axs[1, 1].imshow(proj_max_denoi)
    # axs[1, 1].set_title('denoised max')
    # axs[1, 2].imshow(proj_cor_denoi)
    # axs[1, 2].set_title('denoised corr')
    # for ax in axs.flat:
    #     ax.set_axis_off()
    # plt.tight_layout()
    # plt.show()

    # ================= get roi list ========================
    roi_json_path = ft.look_for_unique_file(
        source=os.path.join(folder_roi, exp_id),
        identifiers=[exp_id],
        file_type='json',
        is_full_path=True)
    roi_list, roi_id_list = get_roi_mask_list_from_pipeline_json(roi_json_path)
    plane_dict.update({
        'roi_list' : roi_list,
        'roi_id_list': [f'{exp_id}_{r:04d}' for r in roi_id_list]
        })
    # print('cccc')
    # ================= get roi list ========================


    # ================= noise traces ========================
    # the true roi ids in the right order to be checked
    roi_ids_target = np.array([int(r[-4:]) for r in plane_dict['roi_id_list']])
    
    # get raw noisy traces
    trace_noisy_raw_path = os.path.join(
        folder_noisy_trace_raw, exp_id, 'roi_traces.h5')
    roi_ids, traces_noisy_raw = get_traces_from_pipeline_h5(trace_noisy_raw_path)
    assert(np.array_equal(roi_ids, roi_ids_target))
    # print('dddd')

    # get neuropil noisy traces
    trace_noisy_neuropil_path = os.path.join(
        folder_noisy_trace_raw, exp_id, 'neuropil_traces.h5')
    roi_ids, traces_noisy_neuropil = get_traces_from_pipeline_h5(trace_noisy_neuropil_path)
    assert(np.array_equal(roi_ids, roi_ids_target))
    # print('eeee')

    # get demixed noisy traces
    trace_noisy_demixed_path = os.path.join(
        folder_noisy_trace_demixed, exp_id, f'{exp_id}_demixed_traces.h5')
    roi_ids, traces_noisy_demixed = get_traces_from_pipeline_h5(trace_noisy_demixed_path)
    assert(np.array_equal(roi_ids, roi_ids_target))
    # print('ffff')

    # get subtract noisy traces
    trace_noisy_subtracted_path = os.path.join(
        folder_noisy_trace_subtracted, exp_id, 'neuropil_correction.h5')
    roi_ids, traces_noisy_subtracted, rs, rmses \
        = get_neuropil_subtracted_traces_from_pipeline_h5(trace_noisy_subtracted_path)
    assert(np.array_equal(roi_ids, roi_ids_target))
    # print('gggg')

    # get dff traces
    trace_noisy_dff_path = os.path.join(
        folder_noisy_trace_dff, exp_id, f'{exp_id}_dff.h5')
    roi_ids, traces_noisy_dff, num_frame_small_baseline, dff_sigmas \
        = get_dff_traces_from_pipeline_h5(trace_noisy_dff_path)
    assert(np.array_equal(roi_ids, roi_ids_target))
    # print('hhhh')

    plane_dict.update({
        'traces_noisy_raw': traces_noisy_raw,
        'traces_noisy_neuropil': traces_noisy_neuropil,
        'traces_noisy_demixed': traces_noisy_demixed,
        'traces_noisy_subtracted': traces_noisy_subtracted,
        'traces_noisy_dff': traces_noisy_dff,
        'subtraction_r': rs,
        'subtraction_rmse': rmses,
        'dff_frame_num_small_baseline': num_frame_small_baseline,
        'dff_sigma': dff_sigmas
        })
    # ================= noise traces ========================

    return plane_dict


def add_rois_and_traces_to_nwb(nwb_f=nwb_f, plane_n=plane_n, 
    plane_dict=plane_dict):

    pass

