import os
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.NwbTools as nt
import utils as utils

nwb_folder = r"\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\nwbs"

curr_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(curr_folder)

fns = ft.look_for_file_list(
    source=nwb_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

for fi, fn in enumerate(fns):

    print('processing {}, {} / {}'.format(fn, fi + 1, len(fns)))

    sess_path = utils.get_lims_session_path_from_session_name(
        session_name=fn[0:10])

    log_path = ft.look_for_unique_file(source=sess_path,
                                       identifiers=['_stim.pkl'],
                                       file_type='pkl',
                                       is_full_path=True)
    # print(log_path)
    log = ft.loadFile(log_path)

    dx = log[b'items'][b'foraging'][b'encoders'][0][b'dx']

    # '''
    # dx the movement of yoked stimulus controlled by the mouse,
    #
    # the exact value is how many pixels the stimulus had moved
    # between two visual frames, but the conversion from this
    # pixel to running distance (cm) is not clear.
    # My guess is that mouse movements are equal to the
    # psysical size of pixels.
    #
    # so:
    # the monior is 52 x 32.5 cm and pixel resolution is 1920 x 1200.
    # so each pixel is 0.027 cm
    # '''
    # dis = np.cumsum(dx) * 0.027

    '''
    the below conversion is from Saskia
    '''
    dxcm = ((dx / 360) * 5.5036 * np.pi * 2)
    dis = np.cumsum(dxcm)

    nwb_f = nt.RecordedFile(os.path.join(nwb_folder, fn))

    description = '''
        1d locomotion position calculated from the "/items/foraging/encoders/0/dx"
        position = np.cumsum(dx) * 0.027. 0.027 is the pixel size in cm.
                  '''
    lm_series = nwb_f.create_timeseries('SpatialSeries', name='distance', modality='other',)
    lm_series.set_data(dis, unit='cm', conversion=np.nan, resolution=np.nan)
    lm_series.set_time_as_link('/acquisition/timeseries/digital_stim_vsync_rise')
    lm_series.set_comments('')
    lm_series.set_description(description)
    lm_series.set_value('num_samples', dis.shape[0])
    lm_series.set_value('reference_frame', '')
    lm_series.set_source('BrainObservatory behavior stage with CamStim software')

    lm_mod = nwb_f.create_module('locomotion')
    lm_interf = lm_mod.create_interface("Position")
    lm_interf.add_timeseries(lm_series)
    lm_series.finalize()
    lm_interf.finalize()
    lm_mod.finalize()

    # print('for debug ...')

    nwb_f.close()