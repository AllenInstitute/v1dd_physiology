import os
# import shutil
import utils as utils
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.core.ImageAnalysis as ia
import cv2
import numpy as np
import matplotlib.pyplot as plt

data_folder = r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\nwbs'
save_folder = r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\eye_tracking_movies\3p'

fns = ft.look_for_file_list(
    source=data_folder, 
    # identifiers=[''],
    identifiers=['M438833_1b'],
    file_type='nwb',
    is_full_path=False
    )

fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

pairs = []

for fn in fns:

    sess_folder = utils.get_lims_session_path_from_session_name(fn[0:10])
    # print(sess_folder)

    video_path = ft.look_for_unique_file(
        source=sess_folder,
        identifiers=['_eye'],
        file_type='avi',
        is_full_path=True)

    save_path = os.path.join(save_folder, f'{fn[0:10]}_eye.avi')
    pairs.append((save_path, video_path))

print('\n'.join([str(p) for p in pairs]))

for pair_i, pair in enumerate(pairs):

    src = pair[1]
    trg = pair[0]

    print(f'processing {pair[0]}, {pair_i+1}/{len(pairs)} ...')

    # rotate the move 90 degree, only for 3p
    video_capture = cv2.VideoCapture(src)
    frame_num = int(video_capture.get(cv2.CAP_PROP_FRAME_COUNT))
    frame_h = int(video_capture.get(cv2.CAP_PROP_FRAME_HEIGHT))
    frame_w = int(video_capture.get(cv2.CAP_PROP_FRAME_WIDTH))
    fps = video_capture.get(cv2.CAP_PROP_FPS)


    mov_r = cv2.VideoWriter(
        filename=trg, 
        fourcc=cv2.VideoWriter_fourcc(*'XVID'), 
        fps=fps, 
        frameSize=(frame_h, frame_w),
        isColor=False)

    for frame_i in range(frame_num):
        _, frame = video_capture.read()
        # print(frame.shape)
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        # frame_r = frame.transpose(1, 0, 2)[::-1, :, :]
        frame_r = frame.transpose()[::-1, :]
        frame_r = cv2.equalizeHist(frame_r)

        # f, axs = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
        # # axs[0].imshow(np.mean(frame, axis=2))
        # # axs[1].imshow(np.mean(frame_r, axis=2))
        # axs[0].imshow(frame)
        # axs[1].imshow(frame_r)
        # plt.show()

        mov_r.write(frame_r)

    mov_r.release()
