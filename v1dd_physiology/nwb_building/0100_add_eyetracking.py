import os, h5py, cv2
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.NwbTools as nt

data_folder = r'\\allen\programs\mindscope\workgroups\surround' \
              r'\v1dd_in_vivo_new_segmentation\data\nwbs'
eye_folder = r'\\allen\programs\mindscope\workgroups\surround' \
             r'\v1dd_in_vivo_new_segmentation\data\eye_tracking_movies'

diagonal_length = 9.0 # mm
side = 'right'
nasal_dir = 'right'
eyetracking_ts_name = 'digital_cam2_exposure_rise'

curr_folder = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_folder)

fns = ft.look_for_file_list(
    source=data_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

# fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

for fi, fn in enumerate(fns):

    print(f'processing {fn}, {fi + 1}/{len(fns)} ...')

    nwb_path = os.path.join(data_folder, fn)
    eye_path = os.path.join(eye_folder, f'{fn[0:10]}_eye_ellipse.hdf5')

    eye_f = h5py.File(eye_path, 'r')
    ell_data = eye_f['ellipse'][()]
    eye_f.close()

    mov_fn = f'{fn[0:10]}_eye_ellipse.avi'
    mov_path = os.path.join(eye_folder, mov_fn)
    mov = cv2.VideoCapture(mov_path)
    mov_shape = (int(mov.get(cv2.CAP_PROP_FRAME_HEIGHT)),
                 int(mov.get(cv2.CAP_PROP_FRAME_WIDTH)))
    mov.release()

    description = '''This is a PupilTracking timeseries. The eyetracking movie is recorded by 
                the pipeline stage and feature points were extracted by DeepLabCut using pipeline eyetracking
                model. For each movie frame, each circumference of corneal reflection, eye lid and pupil were 
                labeled by 12 points using DeepLabCut and the model. Then an ellipse was fit to each of the object. 
                Each ellipse is represent by 5 numbers (see below). The tracking data is saved in the "data" 
                field. Data should be an array with shape (n, 15). n is the eyetracking movie frame number. 
                data[:, 0:5] is the fitted ellipses for the corneal reflection, data[:, 5:10] is the fitted ellipse 
                for the eye, and data[:, 10:15] is the fitted ellipse for the pupil. Each five number of an ellipse 
                represent:
                center_elevation: in mm, small means ventral large means dorsal, this is relative to the movie FOV
                center_azimuth: in mm, small means nasal, large means temporal, this is realtive to the movie FOV
                long_axis_length: in mm
                short_axis_lengt: in mm
                angle_long_axis: in degrees, counterclockwise, 0 is the temporal.
                '''

    pixel_size = diagonal_length / np.sqrt(mov_shape[0] ** 2 + mov_shape[1] ** 2)
    if nasal_dir == 'right':
        pass
    elif nasal_dir == 'left':
        ell_data[:, 1] = mov_shape[1] - ell_data[:, 1]
        ell_data[:, 6] = mov_shape[1] - ell_data[:, 6]
        ell_data[:, 11] = mov_shape[1] - ell_data[:, 11]
        ell_data[:, 4] = (180 - ell_data[:, 4]) % 360
        ell_data[:, 9] = (180 - ell_data[:, 9]) % 360
        ell_data[:, 14] = (180 - ell_data[:, 14]) % 360
    else:
        raise ValueError('\tDo not understand "nasal_dir" ({}). '
                         'Should be "left" or "right"'.format(nasal_dir))

    ell_data[0:4] = ell_data[0:4] * pixel_size
    ell_data[5:9] = ell_data[5:9] * pixel_size
    ell_data[10:14] = ell_data[10:14] * pixel_size

    nwb_f = nt.RecordedFile(nwb_path)
    nwb_f.add_eyetracking_general(ts_path=eyetracking_ts_name, data=ell_data,
                                  module_name='eye_tracking', side=side, comments='',
                                  description=description, source='')
    nwb_f.close()








