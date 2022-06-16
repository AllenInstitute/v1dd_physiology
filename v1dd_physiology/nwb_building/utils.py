import os
import json
import numpy as np
import tifffile as tf
import NeuroAnalysisTools.core.FileTools as ft

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