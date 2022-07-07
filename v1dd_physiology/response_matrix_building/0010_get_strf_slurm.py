import os, sys, h5py
import numpy as np
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.core.TimingAnalysis as ta
import v1dd_physiology.data_fetching as daf

time_win = [-0.5, 1.5]


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


nwb_path = sys.argv[1]
nwb_folder, nwb_fn = os.path.split(nwb_path)

print(f'processing {nwb_fn} ...')

nwb_f = h5py.File(os.path.join(nwb_folder, nwb_fn), 'r')
rm_path = daf.get_rm_path(nwb_f=nwb_f)
rm_f = h5py.File(rm_path, 'a')

onset_grp = nwb_f['analysis/probe_onsets_lsn']
probe_ns = list(onset_grp.keys())
probe_ns.sort()

plane_ns = daf.get_plane_names(nwb_f)
plane_ns.sort()

if 'strf' in rm_f:
    del rm_f['strf']
strf_grp = rm_f.create_group(f'strf')

for plane_n in plane_ns:

    traces_dff, trace_ts = daf.get_plane_traces(
        nwb_f=nwb_f, plane_n=plane_n, trace_type='dff')

    traces_eve, trace_ts = daf.get_plane_traces(
        nwb_f=nwb_f, plane_n=plane_n, trace_type='events')

    if len(traces_dff) > 1: # there are rois in the plane

        strf_grp_plane = strf_grp.create_group(plane_n)

        frame_dur = np.mean(np.diff(trace_ts))
        frame_start = int(np.floor(time_win[0] / frame_dur))
        frame_end = int(np.ceil(time_win[1] / frame_dur))
        t_axis = np.arange(frame_end - frame_start) * frame_dur + (frame_start * frame_dur)

        strf_grp_plane.attrs['sta_timestamps'] = t_axis

        for probe_i, probe_n in enumerate(probe_ns):

            probe_onsets = onset_grp[probe_n][()]

            curr_probe_grp = strf_grp_plane.create_group(probe_n)
            curr_probe_grp.attrs['global_trigger_timestamps'] = probe_onsets
            curr_probe_grp.attrs['sta_traces_dimenstion'] = 'roi x trial x timepoint'

            sta_dff = get_sta(
                arr=traces_dff, arr_ts=trace_ts, trigger_ts=probe_onsets, 
                frame_start=frame_start, frame_end=frame_end)

            sta_eve = get_sta(
                arr=traces_eve, arr_ts=trace_ts, trigger_ts=probe_onsets, 
                frame_start=frame_start, frame_end=frame_end)

            curr_probe_grp.create_dataset(f'sta_dff', data=sta_dff, compression='lzf')
            curr_probe_grp.create_dataset(f'sta_events', data=sta_eve, compression='lzf')

rm_f.close()
nwb_f.close()