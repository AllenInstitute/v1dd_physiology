import h5py


# ========================== FETCHING DATA FROM NWB FILES ================================

def get_projection_images(nwb_path, plane_n, is_plot=False):
    """
    return mean projection and max projection of a given imaging plane
    """

    nwb_f = h5py.File(nwb_path, 'r')

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


def get_vasculature_map(nwb_path,  type='wf', is_standard='False'):
    """
    get vasculature map of a particular session

    :param nwb_f: hdf5 File object
    :param type: str, 'wf' for wide field or '2p' for 2p photon
    :param is_standard: bool, if False, original image acquired;
                              if True, rotated to match standard orientation
                                       up: anterior, left: lateral
    :return vasmap: 2d array
    """

    nwb_f = h5py.File(nwb_path, 'r')

    vasmap_path = f'acquisition/images/vasmap_{type}'

    if is_standard:
        vasmap_path = f'{vasmap_path}_rotated'

    vasmap = nwb_f[vasmap_path][()]

    nwb_f.close()

    return vasmap


def get_windowed_grating_location(nwb_path):
    """
    return altitude and azimuth location of the center of windowed gratings
    """
    nwb_f = h5py.File(nwb_path, 'r')
    alt = nwb_f['stimulus/presentation/drifting_gratings_windowed/center_alt'][()]
    azi = nwb_f['stimulus/presentation/drifting_gratings_windowed/center_azi'][()]
    nwb_f.close()
    return alt, azi


def get_windowed_grating_diameter(nwb_path):
    """
    return windowed drifting grating diameter in degrees 
    """
    nwb_f = h5py.File(nwb_path, 'r')
    dgcw_dia = nwb_f['stimulus/presentation/drifting_gratings_windowed/diameter_deg'][()]
    nwb_f.close()
    return dgcw_dia



# ========================== FETCHING DATA FROM NWB FILES ================================



# ==================== FETCHING DATA FROM RESPONSE METRICS FILES =========================


# ==================== FETCHING DATA FROM RESPONSE METRICS FILES =========================