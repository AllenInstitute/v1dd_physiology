import os, h5py
import numpy as np
import tifffile as tf
import NeuroAnalysisTools.core.FileTools as ft
import utils as utils

nwb_folder = r"\\allen\programs\mindscope\workgroups\surround" \
             r"\v1dd_in_vivo_new_segmentation\data\nwbs"

stim_mov_path = r"\\allen\programs\mindscope\workgroups\surround" \
                r"\v1dd_in_vivo_new_segmentation\data\stim_movies" \
                r"\stim_locally_sparse_nois_16x28_displayed.tif"

grid_size = 9.3 # degree

curr_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(curr_folder)

## analyze locally sparse noise movie
## ==========================================================================
lsn_mov = tf.imread(stim_mov_path)
probe_list = utils.analyze_lsn_movie(mov=lsn_mov, grid_size=grid_size)
probe_ns = ['alt{:06.1f}_azi{:06.1f}_sign{:02.0f}'.format(p[0], p[1], p[2])
            for p in probe_list]

probe_dict = {}
for probe_i, probe_n in enumerate(probe_ns):
    probe_dict.update({probe_n:probe_list[probe_i][3]})

# print(probe_dict)
probe_ns.sort()
## ==========================================================================


fns = ft.look_for_file_list(
    source=nwb_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

# fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

for fi, fn in enumerate(fns):

    print('processing {}, {} / {}'.format(fn, fi + 1, len(fns)))

    nwb_f = h5py.File(os.path.join(nwb_folder, fn), 'a')

    stim_array = nwb_f['/stimulus/presentation/locally_sparse_noise/data'][()]

    # del nwb_f['analysis/probe_onsets_lsn']
    save_grp = nwb_f['analysis'].create_group('probe_onsets_lsn')

    for probe_n in probe_ns:
        frame_list = probe_dict[probe_n]
        onset_ts = []
        for frame_i in frame_list:
            onsets = stim_array[:, 0][stim_array[:, 2] == frame_i]
            onset_ts = onset_ts + list(onsets)

        # print(onset_ts)
        save_grp.create_dataset(probe_n, data=np.array(onset_ts))

    nwb_f.close()

