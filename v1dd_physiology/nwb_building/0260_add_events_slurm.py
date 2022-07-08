import os, sys
import numpy as np
import NeuroAnalysisTools.NwbTools as nt
import NeuroAnalysisTools.core.FileTools as ft

event_folder = "//allen/programs/mindscope/workgroups/surround" \
               "/v1dd_in_vivo_new_segmentation/data/event_detection"
downsample_r = 5

nwb_path = sys.argv[1]
nwb_folder, nwb_fn = os.path.split(nwb_path)

nwb_f = nt.RecordedFile(os.path.join(nwb_folder, nwb_fn))

plane_keys = [k for k in nwb_f.file_pointer['processing'].keys() \
              if k.startswith('rois_and_traces_plane')]
plane_num = len(plane_keys)

sess_id = os.path.splitext(nwb_fn)[0]

for plane_i in range(plane_num):

    if f'l0_events_plane{plane_i}' not in list(nwb_f.file_pointer['processing'].keys()):

        print(f'\tadding events of plane{plane_i} in {nwb_fn}')
        event_path = os.path.join(event_folder, f'{sess_id}_{plane_i}.npz')

        if not os.path.isfile(event_path):
            raise LookupError

        event_f = np.load(event_path)
        events = event_f['events']
        noise_stds = event_f['noise_stds']
        lambdas = event_f['lambdas']

        if len(events) == 1: # zero-roi imaging plane

            e_mo = nwb_f.create_module(f'l0_events_plane{plane_i}')
            e_if = e_mo.create_interface('DfOverF')
            e_ts = nwb_f.create_timeseries('RoiResponseSeries', 'l0_events')
            e_ts.set_data(np.array([[np.nan]]), unit='au', conversion=np.nan, resolution=np.nan)
            e_ts.set_value('data_format', 'roi (row) x time (column)')
            e_ts.set_description(
                'events detected by L0 event detection algorithm from DfOverF traces. '
                'Code by Peter Ledochowitsch.')
            e_ts.set_time_as_link(f'processing/rois_and_traces_plane{plane_i}/DfOverF/dff_raw')
            e_ts.set_value_as_link('segmentation_interface',
                                   f'processing/rois_and_traces_plane{plane_i}/ImageSegmentation')
            e_ts.set_value_as_link('roi_names',
                                   f'processing/rois_and_traces_plane{plane_i}/DfOverF/dff_raw/roi_names')
            e_ts.set_value('num_samples', 0)
            e_ts.set_value('noise_stds', np.array([np.nan]))
            e_ts.set_value('lambdas', np.array([np.nan]))
            e_if.add_timeseries(e_ts)
            e_ts.finalize()
            e_if.finalize()
            e_mo.finalize()

        else:

            events_d = []
            for i in range(downsample_r):
                events_d.append(events[:, i::downsample_r])
            events_d = np.mean(np.array(events_d), axis=0)

            # check data integrity
            roi_num_tot = len(nwb_f.file_pointer[f'processing/rois_and_traces_plane{plane_i}/'
                                                 f'ImageSegmentation/pipeline_roi_names'])
            ts_num = len(nwb_f.file_pointer[f'processing/rois_and_traces_plane{plane_i}/'
                                            f'DfOverF/dff_raw/timestamps'])
            if events_d.shape[0] != roi_num_tot:
                raise ValueError(f'number of evert traces ({events_d.shape[0]}) '
                                 f'does not match number of rois ({roi_num_tot}).')
            if events_d.shape[1] != ts_num:
                raise ValueError(f'number of downsampled event timestamps '
                                 f'({events_d.shape[1]}) does not match number of '
                                 f'df/f timestamps ({ts_num}).')

            e_mo = nwb_f.create_module(f'l0_events_plane{plane_i}')
            e_if = e_mo.create_interface('DfOverF')
            e_ts = nwb_f.create_timeseries('RoiResponseSeries', 'l0_events')
            e_ts.set_data(events_d, unit='au', conversion=np.nan, resolution=np.nan)
            e_ts.set_value('data_format', 'roi (row) x time (column)')
            e_ts.set_description(
                'events detected by L0 event detection algorithm from DfOverF traces. '
                'Code by Peter Ledochowitsch.')
            e_ts.set_time_as_link(f'processing/rois_and_traces_plane{plane_i}/DfOverF/dff_raw')
            e_ts.set_value_as_link('segmentation_interface',
                                   f'processing/rois_and_traces_plane{plane_i}/ImageSegmentation')
            e_ts.set_value_as_link('roi_names',
                                   f'processing/rois_and_traces_plane{plane_i}/DfOverF/dff_raw/roi_names')
            e_ts.set_value('num_samples', events_d.shape[1])
            e_ts.set_value('noise_stds', noise_stds)
            e_ts.set_value('lambdas', lambdas)
            e_if.add_timeseries(e_ts)
            e_ts.finalize()
            e_if.finalize()
            e_mo.finalize()
        
nwb_f.close()