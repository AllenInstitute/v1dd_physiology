import os
import cv2
import h5py
import numpy as np
import NeuroAnalysisTools.PreprocessPipeline as pp
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.DeepLabCutTools as dlct

data_folder = r'\\allen\programs\mindscope\workgroups\surround' \
              r'\v1dd_in_vivo_new_segmentation\data\eye_tracking_movies\3p'
et_side = 'right'
et_confidence_thr = 0.7
et_point_num_thr = 11
et_ellipse_fit_function = cv2.fitEllipse
et_is_generate_labeled_movie = True
et_diagonal_length = 9.0
et_nasal_dir = 'right'

curr_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(curr_folder)

dlc_result_paths = ft.look_for_file_list(
    source=data_folder,
    identifiers=['DLC_resnet50_universal_eye_tracking', 'M438833_1b'],
    # identifiers=['DLC_resnet50_universal_eye_tracking'],
    file_type='h5',
    is_full_path=True)
dlc_result_paths.sort()
print('\n'.join(dlc_result_paths))

for dlcr_i, dlcr_p in enumerate(dlc_result_paths):

    dlcr_fn = os.path.split(dlcr_p)[1]

    print('\nprocessing {}, {} / {}'.format(dlcr_fn, dlcr_i + 1, len(dlc_result_paths)))
    save_path = os.path.join(data_folder, '{}_ellipse.hdf5'.format(dlcr_fn[0:14]))

    if os.path.isfile(save_path):
        print('\tEllipse fit results already exists. skip.'
              'Path: \n\t{}'.format(os.path.realpath(save_path)))
    else:
        print('\tfitting ellipse ...')
        df_pts = dlct.read_data_file(dlcr_p, is_verbose=False)
        df_ell = dlct.get_all_ellipse(df_pts=df_pts, lev_thr=et_confidence_thr,
                                      num_thr=et_point_num_thr, fit_func=et_ellipse_fit_function,
                                      is_verbose=False)

        save_f = h5py.File(save_path, 'x')
        dset = save_f.create_dataset('ellipse', data=np.array(df_ell))
        dset.attrs['columns'] = list(df_ell.columns)
        save_f.close()

        # print('\tgenerating labeled movie ...')
        mov_path_raw = os.path.join(data_folder, '{}.avi'.format(dlcr_fn[0:14]))
        mov_path_lab = os.path.join(data_folder, '{}_ellipse.avi'.format(dlcr_fn[0:14]))
        dlct.generate_labeled_movie(mov_path_raw=mov_path_raw,
                                    mov_path_lab=mov_path_lab,
                                    df_ell=df_ell,
                                    fourcc='XVID',
                                    is_verbose=True,
                                    fps=30.)