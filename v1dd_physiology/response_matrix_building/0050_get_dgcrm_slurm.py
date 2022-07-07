'''
this script extracts response matrices of every roi to
drifting gratings and save them into
response_matrices files
'''


import sys, os, h5py
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.core.TimingAnalysis as ta
import v1dd_physiology.data_fetching as daf

# stim_type = 'full' # 'windowed' or 'full'

nwb_path = sys.argv[1]
stim_type = sys.argv[2]

t_win = [-1., 3.]


def get_sta(arr, arr_ts, trigger_ts, frame_start, frame_end):
    sta_arr = []

    for trig in trigger_ts:
        trig_ind = ta.find_nearest(arr_ts, trig)

        if trig_ind + frame_end < arr.shape[1]:
            curr_sta = arr[:, (trig_ind + frame_start): (trig_ind + frame_end)]
            # print(curr_sta.shape)
            sta_arr.append(curr_sta.reshape((curr_sta.shape[0], 1, curr_sta.shape[1])))

    sta_arr = np.concatenate(sta_arr, axis=1)
    return sta_arr


nwb_folder, nwb_fn = os.path.split(nwb_path)
print(f'processing {nwb_fn}.')

nwb_f = h5py.File(os.path.join(nwb_folder, nwb_fn), 'r')
rm_path = daf.get_rm_path(nwb_f=nwb_f)
# print(rm_path)
rm_f = h5py.File(rm_path, 'a')

onset_grp = nwb_f[f'analysis/onsets_drifting_gratings_{stim_type}']
grating_ns = list(onset_grp.keys())
grating_ns.sort()
# print(grating_ns)

plane_ns = daf.get_plane_names(nwb_f=nwb_f)
plane_ns.sort()

res_grp_key = f'responses_drifting_gratings_{stim_type}'
if res_grp_key in rm_f:
    del rm_f[res_grp_key]
res_grp = rm_f.create_group(res_grp_key)

for plane_n in plane_ns:

    traces_dff, ts = daf.get_plane_traces(
        nwb_f=nwb_f, plane_n=plane_n, trace_type='dff')

    traces_eve, _ = daf.get_plane_traces(
        nwb_f=nwb_f, plane_n=plane_n, trace_type='events')

    if len(traces_dff) > 1:

        res_grp_plane = res_grp.create_group(plane_n)

        frame_dur = np.mean(np.diff(ts))
        frame_start = int(np.floor(t_win[0] / frame_dur))
        frame_end = int(np.ceil(t_win[1] / frame_dur))

        t_axis = np.arange(frame_end - frame_start) * frame_dur + (frame_start * frame_dur)
        res_grp_plane.attrs['sta_timestamps'] = t_axis

        for grating_n in grating_ns:

            grating_onsets = onset_grp[f'{grating_n}/onset_ts_sec'][()]

            curr_grating_grp = res_grp_plane.create_group(grating_n)
            curr_grating_grp.attrs['global_trigger_timestamps'] = grating_onsets
            curr_grating_grp.attrs['sta_traces_dimenstion'] = 'roi x trial x timepoint'

            sta_dff = get_sta(arr=traces_dff, arr_ts=ts, trigger_ts=grating_onsets,
                              frame_start=frame_start, frame_end=frame_end)

            sta_eve = get_sta(arr=traces_eve, arr_ts=ts, trigger_ts=grating_onsets,
                              frame_start=frame_start, frame_end=frame_end)

            curr_grating_grp.create_dataset('sta_dff', data=sta_dff, compression='lzf')

            curr_grating_grp.create_dataset('sta_events', data=sta_eve, compression='lzf')

rm_f.close()
nwb_f.close()            