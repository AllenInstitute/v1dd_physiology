import os, h5py
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.SingleCellAnalysis as sca


db_path = r"\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data"


# ========================== FETCHING DATA FROM NWB FILES ================================

def get_all_sessions(database_path=db_path):
    """
    return all session ids in the database

    Parameters
    ----------
    database_path: string
        path to the database base folder

    Returns
    -------
    sessions: list of strings
        list of all session ids in the database
    """

    nwb_fns = ft.look_for_file_list(
        source=os.path.join(database_path, 'nwbs'),
        identifiers=[],
        file_type='nwb',
        is_full_path=False)
    
    sessions = [n[0:10] for n in nwb_fns]
    return sessions


def get_nwb_path(session_id, database_path=db_path):
    """
    Given session id, return the path to the corresponding nwb path

    Parameters
    ----------
    session_id: string
        M<mouse_id>_<column_id><volume_id>, e.g. 'M409828_11'
    database_path: string
        path to the database base folder

    Returns
    -------
    nwb_path: string
        path to the corresponding nwb file
    """

    nwb_path = ft.look_for_unique_file(
        source=os.path.join(database_path, 'nwbs'),
        identifiers=[session_id],
        file_type='nwb',
        is_full_path=True)

    return nwb_path


def check_nwb_integrity(nwb_f, verbose=True):
    """
    Check the data integrity of an nwb file, and 
    print results. Will not raise Exceptions. 

    Please note this check is not exhaustive.

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    verbose : bool

    Returns
    -------
    msg : string
        results of this check
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    msg = f'\nCheck integrity of {os.path.split(nwb_f.filename)[1]} ...'
    
    # check vasculature map
    for vm in ['vasmap_mp', 'vasmap_mp_rotated', 'vasmap_wf', 'vasmap_wf_rotated']:
        if vm not in nwb_f['acquisition/images']:
            msg += f'\n\tCannot find "{vm}" in "/acquisition/images".'.expandtabs(4)

    # check raw time sync data
    for sync in ['digital_2p_vsync',
                 'digital_cam1_exposure',
                 'digital_cam2_exposure',
                 'digital_stim_photodiode',
                 'digital_stim_vsync']:
        if (f'{sync}_fall' not in nwb_f['acquisition/timeseries']) or \
        (f'{sync}_rise' not in nwb_f['acquisition/timeseries']):

            msg += f'\n\tCannot find "{sync}" signal in "/acquisition/timeseries".'.expandtabs(4)

    # check stimulus
    for stim in ['drifting_gratings_full',
                 'drifting_gratings_windowed',
                 'locally_sparse_noise',
                 'natural_images',
                 'natural_images_12',
                 'natural_movie',
                 'spontaneous']:
        if stim not in nwb_f['stimulus/presentation']:
            msg += f'\n\tCannot find "{stim}" in "/stimulus/presentation".'.expandtabs(4)

    # check stimulus onset times
    for stim_on_t in ['onsets_drifting_gratings_full',
                      'onsets_drifting_gratings_windowed',
                      'probe_onsets_lsn']:
        if stim_on_t not in nwb_f['analysis']:
            msg += f'\n\tCannot find "{stim_on_t}" in "/analysis".'.expandtabs(4)

    # check eye tracking
    if 'eye_tracking_right' not in nwb_f['processing']:
        msg += '\n\tCannot find "eye_tracking_right" in "/processing".'.expandtabs(4)

    # check locomotion
    if 'locomotion' not in nwb_f['processing']:
        msg += '\n\tCannot find "locomotion" in "/processing".'.expandtabs(4)

    # check imaging planes
    rtkeys = [k for k in nwb_f['processing'] if k.startswith('rois_and_traces_plane')]
    msg += f'\n\tTotal number of imaging plane(s): {len(rtkeys)}.'.expandtabs(4)

    for plane_i in range(len(rtkeys)):
        msg += f'\n\tChecking plane{plane_i} ...'.expandtabs(4)
        rtkey = f'rois_and_traces_plane{plane_i}'
        evkey = f'l0_events_plane{plane_i}'

        if not rtkey in nwb_f['processing']:
            msg += f'\n\t\tCannot find "{rtkey}" in "/processing".'.expandtabs(4)
        elif not evkey in nwb_f['processing']:
            msg += f'\n\t\tCannot find "{evkey}" in "/processing".'.expandtabs(4)
        else:
            seg_path = f'ImageSegmentation/imaging_plane'
            if not seg_path in nwb_f[f'processing/{rtkey}']:
                msg += f'\n\t\tCannot find segmentation results ("{seg_path}") in ' \
                       f'"/processing/{rtkey}".'.expandtabs(4)
            else:

                # check traces
                trace_shapes = {}
                for ftkey in ['f_raw',
                              'f_raw_demixed',
                              'f_raw_neuropil',
                              'f_raw_subtracted']:
                    if not f'Fluorescence/{ftkey}' in nwb_f[f'processing/{rtkey}']:
                        msg += f'\n\t\tCannot find extracted traces ("Fluorescence/{ftkey}") in ' \
                               f'"/processing/{rtkey}"'.expandtabs(4)
                    else:
                        trace_shapes.update({ftkey : nwb_f[f'processing/{rtkey}/Fluorescence/{ftkey}/data'].shape})

                # check dF/F
                if not 'DfOverF/dff_raw' in nwb_f[f'processing/{rtkey}']:
                    msg += f'\n\t\tCannot find extracted dF/F traces ("DfOverF/dff_raw") in ' \
                           f'"/processing/{rtkey}"'.expandtabs(4)
                else:
                    trace_shapes.update({'dff_raw' : nwb_f[f'processing/{rtkey}/DfOverF/dff_raw/data'].shape})

                # check events
                if not 'DfOverF/l0_events' in nwb_f[f'processing/{evkey}']:
                    msg += f'\n\t\tCannot find extracted events ("DfOverF/l0_events") in ' \
                           f'/processing/{evkey}'.expandtabs(4)
                else:
                    trace_shapes.update({'l0_events' : nwb_f[f'processing/{evkey}/DfOverF/l0_events/data'].shape})

                # check projection images
                for ri in ['correlation_projection_denoised',
                           'max_projection_denoised',
                           'max_projection_raw',
                           'mean_projection_denoised',
                           'mean_projection_raw']:
                    if f'reference_images/{ri}' not in nwb_f[f'processing/{rtkey}/{seg_path}']:
                        msg += f'\n\t\tCannot find projection imagings ("reference_images/{ri}")' \
                               f'in "/processing/{rtkey}/{seg_path}"'.expandtabs(4)

                roi_ns = [rk for rk in nwb_f[f'processing/{rtkey}/{seg_path}'].keys() if 
                          rk.startswith('roi_') and (not rk.endswith('_list'))]
                roi_ns.sort()
                roi_num = len(roi_ns)

                msg += f'\n\t\tThis plane has {len(roi_ns)} rois.'.expandtabs(4)

                # print(trace_shapes)

                # check trace shape
                if roi_num == 0: # zero ori
                    for t_n, t_shape in trace_shapes.items():
                        if t_shape != (1, 1):
                            msg += f'\n\t\t{t_n} does not have shape: (1, 1).'
                else:
                    # check roi names
                    roi_list = nwb_f[f'processing/{rtkey}/{seg_path}/roi_list'][()]
                    roi_list = [r.decode() for r in roi_list]
                    if not np.array_equal(np.array(roi_ns), np.array(roi_list)):
                        msg += f'\n\t\troi group names and roi_list do not match.'.expandtabs(4)

                    num_time_sample = []
                    for t_n, t_shape in trace_shapes.items():

                        if t_shape[0] != roi_num:
                            msg += f'\n\t\tThe first dimension of {t_n} ({t_shape[0]}) does not ' \
                                   f'match the number of rois ({roi_num}).'.expandtabs(4)

                        num_time_sample.append(t_shape[1])

                    if len(set(num_time_sample)) != 1:
                        msg += f'\n\t\tNumber of timestamps among different traces are not the ' \
                               f'same ({set(num_time_sample)}).'.expandtabs(4)

                # check prediction shape
                cla_grp_key = f'/analysis/roi_classification_pika/plane{plane_i}'
                if not cla_grp_key in nwb_f:
                    msg += f'\n\t\tCannot find pika classification results: "{cla_grp_key}"'.expandtabs(4)
                else:
                    cla_grp = nwb_f[cla_grp_key]
                    for cla_key in ['pipeline_roi_names',
                                    'prediction',
                                    'roi_names',
                                    'score']:
                        if cla_key not in cla_grp:
                            msg += f'\n\t\tCannot find {cla_key} in {cla_grp_key}.'.expandtabs(4)
                        else:
                            len_pred = nwb_f[f'{cla_grp_key}/{cla_key}'].shape[0]
                            if len_pred != roi_num:
                                msg += f'\n\t\tThe length of prediction results "{cla_grp_key}/{cla_key}" ({len_pred}) ' \
                                      f'does not match the number of rois ({roi_num}).'.expandtabs(4)

    msg += '\nDone.'

    if verbose:
        print(msg)
    
    return msg


def get_scope_type(nwb_f):
    """
    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    scope_type : string
        if two-photon return "2p"
        if three-photon return "3p"
    """
    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    if nwb_f['general/optophysiology/imaging_plane_0/device'][()] == b'two-photon scope':
        return "2p"
    elif nwb_f['general/optophysiology/imaging_plane_0/device'][()] == b'three-photon scope':
        return "3p"
    else:
        raise LookupError("Do not understand imaging device.")


def get_vasculature_map(nwb_f,  type='wf', is_standard='False'):
    """
    get vasculature map of a particular session

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode
    type : str, 
        'wf' for wide field, 'mp' for multi-photon (2p or 3p)
    is_standard : bool
        if False, original image acquired;
        if True, rotated to match standard orientation 
        (up: anterior, left: lateral).
    
    Returns
    ------- 
    vasmap: 2d array
        vasculature_map
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    vasmap_path = f'acquisition/images/vasmap_{type}'

    if is_standard:
        vasmap_path = f'{vasmap_path}_rotated'

    vasmap = nwb_f[vasmap_path][()]

    return vasmap


def get_session_id(nwb_f):
    """
    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode


    Returns
    -------     
    session_id: string
        10-character session id 'M<mouse_id>_<column_id><volume_id>'.
        for example: M409828_11


    """
    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    nwb_fn = os.path.split(nwb_f.filename)[1]
    return nwb_fn[0:10]


def get_lims_session_id(nwb_f):
    """
    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    sess_id : str
        LIMS ophys session id, unquie to every nwb file
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    return nwb_f['general/session_id'][()].decode()


def get_windowed_grating_location(nwb_f):
    """
    return altitude and azimuth location of the center of windowed gratings

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    alt : float
        altitude of the windowed drifting grating circle in degrees
    azi : float
        azimuth of the windowed drifting grating circle in degrees
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    alt = nwb_f['stimulus/presentation/drifting_gratings_windowed/center_alt'][()]
    azi = nwb_f['stimulus/presentation/drifting_gratings_windowed/center_azi'][()]
    
    return alt, azi


def get_windowed_grating_diameter(nwb_f):
    """
    return windowed drifting grating diameter in degrees 

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    dgcw_dia : float
        diameter of the windowed drifing grating circle in degrees
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    dgcw_dia = nwb_f['stimulus/presentation/drifting_gratings_windowed/diameter_deg'][()]

    return dgcw_dia


def get_plane_names(nwb_f):
    """
    return plane names in a session 

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    plane_ns : list of string
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    plane_ns = [k[-6:] for k in nwb_f['processing'].keys() if k.startswith('rois_and_traces_plane')]
    
    return plane_ns


def get_lims_experiment_id(nwb_f, plane_n):
    """
    get the LIMS experiment id for a imaging plane

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    Returns
    -------
    exp_id : str
        LIMS experiment id for this imaging plane
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    exp_id = nwb_f[f'processing/rois_and_traces_{plane_n}'
                   f'/experiment_id'][()].decode()

    return exp_id


def get_plane_depth(nwb_f, plane_n):
    """
    get the depth in microns for a imaging plane

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    Returns
    -------
    depth : int
        micorns under pia
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')
    
    depth = nwb_f[f'processing/rois_and_traces_{plane_n}/imaging_depth_micron'][()]

    return depth


def get_plane_projections(nwb_f, plane_n):
    """
    get the depth in microns for a imaging plane

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    Returns
    -------
    proj_raw_mean : 2d array
    proj_raw_max : 2d array
    proj_denoised_mean : 2d array
    proj_denoised_max : 2d array
    proj_denoised_corr : 2d array
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    img_grp = nwb_f[f'processing/rois_and_traces_{plane_n}/ImageSegmentation'
                    f'/imaging_plane/reference_images']

    proj_raw_mean = img_grp['mean_projection_raw/data'][()]
    proj_raw_max = img_grp['max_projection_raw/data'][()]
    proj_denoised_mean = img_grp['mean_projection_denoised/data'][()]
    proj_denoised_max = img_grp['max_projection_denoised/data'][()]
    proj_denoised_corr = img_grp['correlation_projection_denoised/data'][()]

    return proj_raw_mean, proj_raw_max, proj_denoised_mean, \
        proj_denoised_max, proj_denoised_corr


def get_roi_ns(nwb_f, plane_n):
    """
    get the roi names for a imaging plane

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    Returns
    -------
    roi_ns : list of string
        roi counting should be continuous, with the last 4 digit as 
        index, which can be used to retrieve saved traces.
        ['roi_0000', 'roi_0001', ....]
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    # roi_ns = nwb_f[f'processing/rois_and_traces_{plane_n}'
    #                f'/ImageSegmentation/imaging_plane/roi_list'][()]
    # roi_ns = [r.decode() for r in roi_ns]
    # roi_ns.sort()

    roi_ns = nwb_f[f'processing/rois_and_traces_{plane_n}'
                   f'/ImageSegmentation/imaging_plane'].keys()
    roi_ns = [r for r in roi_ns if r.startswith('roi_')]
    roi_ns = [r for r in roi_ns if not r.endswith('_list')]
    roi_ns.sort()
    
    return roi_ns


def get_pika_roi_ids(nwb_f, plane_n):
    """
    get the ori_ids from pika for a imaging plane

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    Returns
    -------
    roi_ids : list of string
        the format is <session_id>_<roi_id>, for example:
        ['795018590_0000', '795018590_0001', '795018590_0002', ...]
        Note: roi counting may not be continuous
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')
    
    roi_grp = nwb_f[f'processing/rois_and_traces_{plane_n}/ImageSegmentation'
                    f'/imaging_plane']
    roi_ns = get_roi_ns(nwb_f=nwb_f, plane_n=plane_n)
    roi_ids = [roi_grp[f'{r}/roi_description'][()].decode() for r in roi_ns]
    return roi_ids


def get_plane_traces(nwb_f, plane_n, trace_type):
    """
    get the activity traces for rois in a imaging plane.

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    trace_type : string
        type of trace to be extracted. Should be one of the 
        following:
            'raw'
            'demixed'
            'neuropil'
            'subtracted'
            'dff'
            'events'

    Returns
    -------
    traces : 2d array
        activity traces of all rois in this imaging plane with 
        specified tracetype. (roi x time)
    trace_ts : 1d array
        timestamps for the activity trace in seconds, 
        len(trace_ts) == traces.shape[1]
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    if trace_type == 'raw':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw']
    elif trace_type == 'demixed':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw_demixed']
    elif trace_type == 'neuropil':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw_neuropil']
    elif trace_type == 'subtracted':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw_subtracted']
    elif trace_type == 'dff':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/DfOverF/dff_raw']
    elif trace_type == 'events':
        trace_grp = nwb_f[f'processing/l0_events_{plane_n}' \
                          f'/DfOverF/l0_events']
    else:
        raise LookupError(f'Do not understand "trace_type", should be '
            f'one of the following ["raw", "demixed", "neuropil", "sutracted", '
            f'"dff", "events"]. Got "{trace_type}".')

    traces = trace_grp['data'][()]
    trace_ts = trace_grp['timestamps'][()]
    
    assert(traces.shape[1] == trace_ts.shape[0])

    return traces, trace_ts


def get_pika_roi_id(nwb_f, plane_n, roi_n):
    """
    get the pika roi id of an roi

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    roi_n : string
        name of the roi, e.g. 'roi_0005'

    Returns
    -------
    pika_roi_id : str
        the format is <session_id>_<roi_id>, for example:
        '795018590_0000'. Note: roi counting may not be continuous
    """

    roi_id = nwb_f[f'processing/rois_and_traces_{plane_n}/ImageSegmentation'
                   f'/imaging_plane/{roi_n}/roi_description'][()].decode()
    return roi_id


def get_pika_classifier_score(nwb_f, plane_n, roi_n):
    """
    get the pika classifier score of an roi

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    roi_n : string
        name of the roi, e.g. 'roi_0005'

    Returns
    -------
    score : float
        score should be in the range from 0 to 1. 
        The larger the more likely for this roi to be a cell soma.
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    roi_ind = int(roi_n[-4:])
    return nwb_f[f'analysis/roi_classification_pika/{plane_n}/score'][roi_ind]


def get_roi_mask(nwb_f, plane_n, roi_n):
    """
    get the binary roi mask of an roi

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    roi_n : string
        name of the roi, e.g. 'roi_0005'

    Returns
    -------
    roi_mask : 2d array
        binary mask of the roi
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    plane_grp = nwb_f[f'processing/rois_and_traces_{plane_n}/ImageSegmentation']

    width = plane_grp['img_width'][()]
    height = plane_grp['img_height'][()]

    mask = np.zeros((height, width), dtype=np.uint8)
    
    roi_grp = plane_grp[f'imaging_plane/{roi_n}']
    pixels = roi_grp['pix_mask'][()]
    mask[pixels[1, :], pixels[0, :]] = 1

    return mask


def get_single_trace(nwb_f, plane_n, roi_n, trace_type):
    """
    get the activity trace for an roi.

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    plane_n : string
        name of the plane, should be 'plane0', 'plane1', ...

    roi_n : string
        name of the roi, e.g. 'roi_0005'

    trace_type : string
        type of trace to be extracted. Should be one of the 
        following:
            'raw'
            'demixed'
            'neuropil'
            'subtracted'
            'dff'
            'events'

    Returns
    -------
    trace : 1d array
        activity trace of specified tracetype
    trace_ts : 1d array
        timestamps for the activity trace in seconds, 
        should be the same shape as trace
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    roi_ind = int(roi_n[-4:])

    if trace_type == 'raw':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw']
    elif trace_type == 'demixed':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw_demixed']
    elif trace_type == 'neuropil':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw_neuropil']
    elif trace_type == 'subtracted':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/Flurescence/f_raw_subtracted']
    elif trace_type == 'dff':
        trace_grp = nwb_f[f'processing/rois_and_traces_{plane_n}' \
                          f'/DfOverF/dff_raw']
    elif trace_type == 'events':
        trace_grp = nwb_f[f'processing/l0_events_{plane_n}' \
                          f'/DfOverF/l0_events']
    else:
        raise LookupError(f'Do not understand "trace_type", should be '
            f'one of the following ["raw", "demixed", "neuropil", "sutracted", '
            f'"dff", "events"]. Got "{trace_type}".')

    trace = trace_grp['data'][roi_ind, :]
    trace_ts = trace_grp['timestamps'][()]

    assert(trace.shape == trace_ts.shape)

    return trace, trace_ts


def get_stim_list(nwb_f):
    """
    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    stims : list of strings
        list of all displayed stimuli
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    stims = list(nwb_f['stimulus/presentation'].keys())

    return stims


def get_dgc_onset_times(nwb_f, dgc_type):
    """
    get onset times of each drifting grating condition

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    dgc_type: string
        should be 'windowed' or 'full'

    Returns
    -------
    dgc_onset_times : dictionary
        {condition_name : onset_timestamps in seconds (1d array)}
    """
    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    dgc_grp = nwb_f[f'analysis/onsets_drifting_gratings_{dgc_type}']

    dgc_onset_times = {}
    for dgc_condi in sorted(dgc_grp):
        dgc_onset_times.update(
            {dgc_condi : dgc_grp[dgc_condi]['onset_ts_sec'][()]}
            )

    return dgc_onset_times


def get_lsn_onset_times(nwb_f):
    """
    get onset times of each locally sparse noise square.

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    lsn_onset_times : dictionary
        {square_name : onset_timestamps in seconds (1d array)}
    """
    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    lsn_grp = nwb_f['analysis/probe_onsets_lsn']

    lsn_onset_times = {}
    for probe in sorted(lsn_grp):
        lsn_onset_times.update(
            {probe : lsn_grp[probe][()]}
            )

    return lsn_onset_times

    
# ========================== FETCHING DATA FROM NWB FILES ================================



# ==================== FETCHING DATA FROM RESPONSE METRICS FILES =========================

def get_rm_path(nwb_f):
    """
    given the nwb file object, return the path to the 
    corresponding response matrix file path 

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode

    Returns
    -------
    rm_path : string
        path to the corresponding response matrix file
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    nwb_folder, nwb_fn = os.path.split(nwb_f.filename)
    base_folder = os.path.split(nwb_folder)[0]

    rm_path = os.path.join(
        base_folder, 
        'response_matrices', 
        f'response_matrix_{nwb_fn[0:10]}.hdf5')
    
    return rm_path


def get_session_id_from_rmf(rm_f):
    """
    return 10-character session id from response matrix file

    Parameters
    ----------
    rm_f : hdf5 File object
        response matrix file, should be in read-only mode

    Returns
    -------
    sess_id : 10-character string
        e.g., 'M409828_11'
    """
    if rm_f.mode != 'r':
        raise OSError('The response matrix file should be opened in read-only mode.')

    rm_fn = os.path.split(rm_f.filename)[1]
    return rm_fn[16:26]


def get_strf(rm_f, plane_n, roi_n, trace_type='events'):
    """
    get spatial temporal receptive field of a given roi defined by 
    response matrix file, plane_n, and roi_n.

    Parameters
    ----------
    rm_f : hdf5 File object
        response matrix file, should be in read-only mode

    plane_n : 6-character string
        plane name, e.g., 'plane0', 'plane1', ...

    roi_n : 8-character string
        roi name, e.g., 'roi_0000', 'roi_0100', ....

    trace_type : string
        type of response to be extracted, should be either 'dff' or 'events'

    Returns
    -------
    strf : NeuroAnalysisTools.SingleCellAnalysis.SpatialTemporalReceptiveField
    object
        spatial temporal receptive field of the defined roi.
    """

    if rm_f.mode != 'r':
        raise OSError('The response matrix file should be opened in read-only mode.')

    roi_ind = int(roi_n[-4:])
    sess_id = get_session_id_from_rmf(rm_f)

    h5_grp = rm_f[f'strf/{plane_n}']
    
    if trace_type in ['dff', 'events']:
    
        sta_ts = h5_grp.attrs['sta_timestamps']

        probe_ns = list(h5_grp.keys())
        probe_ns.sort()

        locations = []
        signs = []
        traces = []
        trigger_ts = []

        for probe_i, probe_n in enumerate(probe_ns):

            locations.append([float(probe_n[3:9]), float(probe_n[13:19])])
            signs.append(int(probe_n[24:26]))

            traces.append(h5_grp[f'{probe_n}/sta_{trace_type}'][roi_ind, :, :])
            trigger_ts.append(h5_grp[f'{probe_n}'].attrs['global_trigger_timestamps'])

        return sca.SpatialTemporalReceptiveField(
            locations=locations, signs=signs, traces=traces, time=sta_ts,
            trigger_ts=trigger_ts, name=f'{sess_id}_{plane_n}_{roi_n}',
            locationUnit='degree', trace_data_type=f'sta_{trace_type}')
    
    else:
        raise LookupError('Do not under stand "trace_type", should be "dff" of "events".')


def get_greedy_rf(rm_f, plane_n, roi_n, trace_type='events', alpha='0.100'):
    """
    get greedy pixel-wise receptive field of a given roi defined by 
    response matrix file, plane_n, and roi_n.

    Parameters
    ----------
    rm_f : hdf5 File object
        response matrix file, should be in read-only mode

    plane_n : 6-character string
        plane name, e.g., 'plane0', 'plane1', ...

    roi_n : 8-character string
        roi name, e.g., 'roi_0000', 'roi_0100', ....

    trace_type : string
        type of response to be extracted, should be either 'dff' or 'events'

    alpha : 5-character string
        significance level, should be '0.100' or '0.050'

    Returns
    -------
    greedy_rf_on : NeuroAnalysisTools.SingleCellAnalysis.SpatialReceptiveField
    object
    
    greedy_rf_off : NeuroAnalysisTools.SingleCellAnalysis.SpatialReceptiveField
    object
    """
    
    if rm_f.mode != 'r':
        raise OSError('The response matrix file should be opened in read-only mode.')
    
    if alpha not in ['0.100', '0.050']:
        raise LookupError(f"cannot understand parameter 'alpha' "
                          f"({alpha}), should be '0.100' or '0.050'.")

    if trace_type == 'dff':
        h5_grp_key = f'greedy_rfs_1000ms_alpha{alpha}_{trace_type}/{plane_n}'
    elif trace_type == 'events':
        h5_grp_key = f'greedy_rfs_0500ms_alpha{alpha}_{trace_type}/{plane_n}'
    else:
        raise LookupError(f"cannot understand parameter 'trace_type' "
                          f"({trace_type}), should be 'dff' or 'events'.")

    # print(h5_grp_key)
    h5_grp = rm_f[h5_grp_key]
    roi_ind = int(roi_n[-4:])

    alts_on = h5_grp['greedy_pixelwise_rfs_on'].attrs['alt_positions']
    azis_on = h5_grp['greedy_pixelwise_rfs_on'].attrs['azi_positions']
    map_on = h5_grp['greedy_pixelwise_rfs_on'][roi_ind, :, :]
    assert (len(alts_on) == map_on.shape[0])
    assert (len(azis_on) == map_on.shape[1])
    rf_on = sca.SpatialReceptiveField(mask=map_on, altPos=alts_on, aziPos=azis_on,
                                      sign='ON', dataType=trace_type)

    alts_off = h5_grp['greedy_pixelwise_rfs_off'].attrs['alt_positions']
    azis_off = h5_grp['greedy_pixelwise_rfs_off'].attrs['azi_positions']
    map_off =h5_grp['greedy_pixelwise_rfs_off'][roi_ind, :, :]
    assert (len(alts_off) == map_off.shape[0])
    assert (len(azis_off) == map_off.shape[1])
    rf_off = sca.SpatialReceptiveField(mask=map_off, altPos=alts_off, aziPos=azis_off,
                                       sign='OFF', dataType=trace_type)

    return rf_on, rf_off


def get_dgcrm(rm_f, plane_n, roi_n, trace_type='events', dgc_type='full'):
    """
    get drifting grating response matrix of a given roi defined by 
    response matrix file, plane_n, and roi_n.

    Parameters
    ----------
    rm_f : hdf5 File object
        response matrix file, should be in read-only mode

    plane_n : 6-character string
        plane name, e.g., 'plane0', 'plane1', ...

    roi_n : 8-character string
        roi name, e.g., 'roi_0000', 'roi_0100', ....

    trace_type : string
        type of response to be extracted, should be either 'dff' or 'events'

    dgc_type : string
        type of drifting grating circle, should be 'full' of 'windowed'

    Returns
    -------
    dgcrm : NeuroAnalysisTools.SingleCellAnalysis.DriftingGratingResponseMatrix
    object
        drifting gration response matrix of the defined roi.
    """

    if rm_f.mode != 'r':
        raise OSError('The response matrix file should be opened in read-only mode.')
    
    roi_ind = int(roi_n[-4:])

    if dgc_type in ['full', 'windowed']:
            grp_n = f'responses_drifting_gratings_{dgc_type}/{plane_n}'
    else:
        raise LookupError(f"cannot understand parameter 'dgc_type' "
                          f"({dgc_type}), should be 'full' or 'windowed'.")

    h5_grp = rm_f[grp_n]

    if trace_type in ['dff', 'events']:
        dgcrm = sca.get_dgc_response_matrix_from_nwb(
            h5_grp=h5_grp, roi_ind=roi_ind, trace_type=f'sta_{trace_type}')
    else:
        raise LookupError(f"cannot understand parameter 'trace_type' "
                          f"({trace_type}), should be 'events' or 'dff'.")

    return dgcrm


# ==================== FETCHING DATA FROM RESPONSE METRICS FILES =========================


if __name__ == "__main__":

    
    # # ========================================================================
    # nwb_path = r"\\allen\programs\mindscope\workgroups\surround" \
    #            r"\v1dd_in_vivo_new_segmentation\data\nwbs" \
    #            r"\M427836_25_20190426.nwb"

    # nwb_f = h5py.File(nwb_path, 'r')
    # check_nwb_integrity(nwb_f=nwb_f, verbose=True)
    # # ========================================================================


    # # ========================================================================
    # import matplotlib.pyplot as plt
    # nwb_path = r"\\allen\programs\mindscope\workgroups\surround" \
    #            r"\v1dd_in_vivo_new_segmentation\data\nwbs" \
    #            r"\M427836_25_20190426.nwb"
    # plane_n = 'plane3'
    # roi_n = 'roi_0034'

    # nwb_f = h5py.File(nwb_path, 'r')
    # rm_path = get_rm_path(nwb_f=nwb_f)
    # rm_f = h5py.File(rm_path, 'r')

    # strf = get_strf(rm_f=rm_f, plane_n=plane_n, roi_n=roi_n, 
    #     trace_type='dff')
    # print(strf.name)
    # strf.plot_traces(yRange=[0, 0.02])
    # plt.show()

    # dgcrm = get_dgcrm(rm_f=rm_f, plane_n=plane_n, roi_n=roi_n, 
    #     trace_type='dff', dgc_type='windowed')
    # # ========================================================================


    # ========================================================================
    import matplotlib.pyplot as plt
    nwb_path = r"\\allen\programs\mindscope\workgroups\surround" \
               r"\v1dd_in_vivo_new_segmentation\data\nwbs" \
               r"\M409828_13_20181213.nwb"
    plane_n = 'plane0'

    nwb_f = h5py.File(nwb_path, 'r')
    rm_path = get_rm_path(nwb_f=nwb_f)
    rm_f = h5py.File(rm_path, 'r')

    roi_n = 'roi_0086'

    rf_on, rf_off = get_greedy_rf(rm_f=rm_f, plane_n=plane_n, roi_n=roi_n,
        trace_type='events', alpha='0.100')
    f, axs = plt.subplots(ncols=1, nrows=2, figsize=(6,8))
    f.suptitle(roi_n)
    rf_on.plot_rf2(plot_axis=axs[0], cmap='Reds')
    rf_off.plot_rf2(plot_axis=axs[1], cmap='Blues')
    plt.tight_layout()
    plt.show()

    # for i in range(100):

    #     roi_n = f'roi_{i:04d}'

    #     rf_on, rf_off = get_greedy_rf(rm_f=rm_f, plane_n=plane_n, roi_n=roi_n,
    #         trace_type='events', alpha='0.100')
    #     f, axs = plt.subplots(ncols=1, nrows=2, figsize=(6,8))
    #     f.suptitle(roi_n)
    #     rf_on.plot_rf2(plot_axis=axs[0], cmap='Reds')
    #     rf_off.plot_rf2(plot_axis=axs[1], cmap='Blues')
    #     plt.tight_layout()
    #     plt.show()
    # ========================================================================