import sys, os, h5py
import numpy as np
import pandas as pd
import tifffile as tf
import matplotlib.pyplot as plt
import NeuroAnalysisTools.core.FileTools as ft
import greedy_pixelwise_rf as gprf
import v1dd_physiology.data_fetching as daf

nwb_path = sys.argv[1]
nwb_folder, nwb_fn = os.path.split(nwb_path)
print(f'processing {nwb_fn} ...')

# trace_type = 'dff' # 'dff' or 'events'
trace_type = sys.argv[2]
if trace_type == 'dff':
    time_offsets = [0., 1.0]
elif trace_type == 'events':
    time_offsets = [0., 0.5]

alphas = [0.05, 0.1]
grid_size = 9.3
gray_level = 0

# template_path = r"\\allen\programs\mindscope\workgroups\surround" \
#                 r"\v1dd_in_vivo_new_segmentation\data\stim_movies" \
#                 r"\stim_locally_sparse_nois_16x28_displayed.tif"

template_path = "/allen/programs/mindscope/workgroups/surround" \
                "/v1dd_in_vivo_new_segmentation/data/stim_movies" \
                "/stim_locally_sparse_nois_16x28_displayed.tif"

template = tf.imread(template_path)
azis = (np.arange(template.shape[2]) - ((template.shape[2] - 1) / 2)) * grid_size
alts = (np.arange(template.shape[1]) - ((template.shape[1] - 1) / 2)) * grid_size
alts = alts[::-1]

tw_len = time_offsets[1] - time_offsets[0]
tw_str = f'{tw_len * 1000:04.0f}ms'

nwb_f = h5py.File(nwb_path, 'r')
rm_path = daf.get_rm_path(nwb_f=nwb_f)
rm_f = h5py.File(rm_path, 'a')
plane_ns = daf.get_plane_names(nwb_f=nwb_f)

for alpha in alphas:

    print(f'\tsignificance level: {alpha}:')

    alpha_grp_n = f'greedy_rfs_{tw_str}_alpha{alpha:5.3f}_{trace_type}'
    if alpha_grp_n in rm_f:
        del rm_f[alpha_grp_n]
    alpha_grp = rm_f.create_group(alpha_grp_n)

    for plane_i, plane_n in enumerate(plane_ns):

        print(f'\t\t{plane_n}, {plane_i + 1} / {len(plane_ns)}  ...')

        stim_table = nwb_f['stimulus/presentation/locally_sparse_noise/data'][()]
        stim_table = pd.DataFrame(data=stim_table, columns=('start_ts', 'end_ts', 'frame'))

        roi_ns = daf.get_roi_ns(nwb_f=nwb_f, plane_n=plane_n)

        if len(roi_ns) > 0:

            if plane_n in alpha_grp:
                del alpha_grp[plane_n]
                
            plane_grp = alpha_grp.create_group(plane_n)

            _, ts = daf.get_single_trace(nwb_f=nwb_f, plane_n=plane_n, 
                roi_n='roi_0000', trace_type=trace_type)

            # convert timestamps to frame numbers
            imaging_frames = np.zeros((len(stim_table), 2), dtype=np.int64)
            for i_sweep, start_ts in enumerate(stim_table['start_ts'].values):
                end_ts = stim_table['end_ts'].values[i_sweep]
                if start_ts < ts[0]:
                    print('Invalid timestamp: too early')
                elif end_ts > ts[-1]:
                    print('Invalid timestamp: too late')
                else:
                    imaging_frames[i_sweep, 0] = np.argwhere(ts >= (start_ts + time_offsets[0]))[0, 0]
                    imaging_frames[i_sweep, 1] = np.argwhere(ts <= (start_ts + time_offsets[1]))[-1, 0]
            stim_table['start'] = imaging_frames[:, 0]
            stim_table['end'] = imaging_frames[:, 1]
            # convert timestamps to frame numbers finished

            rfs_on = []
            rfs_off = []

            for roi_i, roi_n in enumerate(roi_ns):
                
                if roi_i % 50 == 0:
                    print(f'\t\t\t{roi_n}, {roi_i + 1} / {len(roi_ns)}.')

                trace, _ = daf.get_single_trace(nwb_f=nwb_f, plane_n=plane_n, 
                    roi_n=roi_n, trace_type=trace_type)
                trace = np.where(np.isfinite(trace), trace, 0)
                rf_on, rf_off = gprf.get_receptive_field_greedy(L0_events=trace,
                                                                stimulus_table=stim_table,
                                                                LSN_template=template,
                                                                GRAY_VALUE=gray_level,
                                                                alpha=alpha,
                                                                sweep_response_type='mean')

                assert (rf_on.shape[0] == len(alts))
                assert (rf_on.shape[1] == len(azis))
                assert (rf_off.shape[0] == len(alts))
                assert (rf_off.shape[1] == len(azis))

                rfs_on.append(rf_on)
                rfs_off.append(rf_off)

                # plt.figure()
                # ax1 = plt.subplot(211)
                # ax1.imshow(rf_on, cmap='Reds', interpolation='none', origin='lower')
                # ax2 = plt.subplot(212)
                # ax2.imshow(rf_off, cmap='Blues', interpolation='none', origin='lower')
                # plt.show()

            rfs_on = np.array(rfs_on)
            rfs_off = np.array(rfs_off)

            rf_on_dset = plane_grp.create_dataset(
                'greedy_pixelwise_rfs_on', data=rfs_on, compression='lzf')
            rf_on_dset.attrs['alt_positions'] = alts
            rf_on_dset.attrs['azi_positions'] = azis

            rf_off_dset = plane_grp.create_dataset(
                'greedy_pixelwise_rfs_off', data=rfs_off, compression='lzf')
            rf_off_dset.attrs['alt_positions'] = alts
            rf_off_dset.attrs['azi_positions'] = azis

nwb_f.close()
rm_f.close()
