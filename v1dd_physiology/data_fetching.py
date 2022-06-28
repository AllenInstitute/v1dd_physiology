import os, h5py
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft

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


def get_scope_type(nwb_f):
    """
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


def get_projection_images(nwb_f, plane_n, is_plot=False):
    """
    return mean projection and max projection of a given imaging plane

    Parameters
    ----------
    nwb_f : hdf5 File object
        should be in read-only mode
    plane_n : string
        "plane0", "plane1", ...
    is_plot : bool

    Returns
    -------
    mean_proj : 2d array
        mean projection 
    max_proj : 2d array
        max projection
    """

    if nwb_f.mode != 'r':
        raise OSError('The nwb file should be opened in read-only mode.')

    mean_proj = nwb_f['/processing/rois_and_traces_{}/ImageSegmentation' \
                      '/imaging_plane/reference_images/mean_projection' \
                      '/data'.format(plane_n)][()]
    max_proj = nwb_f['/processing/rois_and_traces_{}/ImageSegmentation' \
                      '/imaging_plane/reference_images/max_projection' \
                      '/data'.format(plane_n)][()]
    nwb_f.close()

    if is_plot:
        sl = 2 # plot saturation level [0-100]
        f = plt.figure(figsize=(10, 5))
        f.suptitle(plane_n)
        ax0 = f.add_subplot(121)
        ax0.set_title('mean_projection')
        ax0.set_axis_off()
        ax0.imshow(
            mean_proj, 
            vmin=np.percentile(mean_proj[:], sl),
            vmax=np.percentile(mean_proj[:], 100-sl), 
            cmap='gray', 
            interpolation='nearest'
            )

        ax1 = f.add_subplot(122)
        ax1.set_title('max_projection')
        ax1.set_axis_off()
        ax1.imshow(
            max_proj, 
            vmin=np.percentile(max_proj[:], sl),
            vmax=np.percentile(max_proj[:], 100-sl), 
            cmap='gray', 
            interpolation='nearest'
            )

    return mean_proj, max_proj


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

    roi_ns = nwb_f[f'processing/rois_and_traces_{plane_n}'
                   f'/ImageSegmentation/imaging_plane/roi_list'][()]
    roi_ns = [r.decode() for r in roi_ns]
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

# ========================== FETCHING DATA FROM NWB FILES ================================



# ==================== FETCHING DATA FROM RESPONSE METRICS FILES =========================


# ==================== FETCHING DATA FROM RESPONSE METRICS FILES =========================