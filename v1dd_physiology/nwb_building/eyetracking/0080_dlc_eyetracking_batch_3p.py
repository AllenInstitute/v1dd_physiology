"""
this is a script version of the jupyter notebook for Allen Institute pipeline eyetracking analysis
https://github.com/AllenInstitute/dlc-eye-tracking/blob/master/DLC%20Eye%20Tracking%20and%20Ellipse%20Fitting.ipynb

packages:
https://github.com/AlexEMG/DeepLabCut

https://github.com/AllenInstitute/dlc-eye-tracking/
by Peter Ledochowitsch, MAT, Allen Institute

This notebook demonstrates the necessary steps to use DeepLabCut for Eye Tracking.
It assumes that DeepLabCut is installed.

This notebook illustrates how to:

do eye and pupil tracking
fit ellipses
generate labeled video
This uses a pretrained model! Make sure that the GPU is available or the kernel might die!


***
this script is not part of the NeuroAnalysisTools analysis scripts
please use a environment with DeepLabCut and TensorFlow installed to run
"""

import os
import tensorflow as tf
import deeplabcut

print('tensorflow version: {}'.format(tf.__version__))
print('deeplabcut version: {}'.format(deeplabcut.__version__))

# ======================================== windows ==================================================
# this path point to the backed up model save in Jun's online storage disks
# path_config = r"Z:\deeplabcut_models\fromPeter\universal_eye_tracking-peterl-2019-07-10\config.yaml"

# this path point to the model save on the \\allen server
path_config = r"\\allen\programs\braintv\workgroups\cortexmodels\peterl\visual_behavior" \
              r"\DLC_models\universal_eye_tracking-peterl-2019-07-10\config.yaml"
# ======================================== windows ==================================================


# ========================================= ubuntu ==================================================
# this path point to the backed up model save in Jun's online storage disks
# from Jun's ubuntu machine
# path_config = '/media/data4/deeplabcut_models/fromPeter/universal_eye_tracking-peterl-2019-07-10' \
#               '/config_backup.yaml'

# this path point to the model save on the //allen server from Jun's ubuntu machine
# path_config = '/media/cortexmodels/peterl/visual_behavior/DLC_models' \
#               '/universal_eye_tracking-peterl-2019-07-10/config.yaml'
# ========================================= ubuntu ==================================================


mov_folder = r'\\allen\programs\mindscope\workgroups\surround\v1dd_in_vivo_new_segmentation\data\eye_tracking_movies\3p'
mov_fns = os.listdir(mov_folder)
mov_fns = [f for f in mov_fns if f[-4:] == '.avi']
# mov_fns = [f for f in mov_fns if (f[-4:] == '.avi') and ('M438833_1b' in f)]
mov_paths = [os.path.join(mov_folder, fn) for fn in mov_fns]
print('number of movies: {}'.format(len(mov_paths)))
print('\n'.join(mov_paths))

# deeplabcut.analyze_videos(path_config,mov_paths)  #can accept list of video files
deeplabcut.create_labeled_video(path_config,mov_paths) #can accept list of video files
