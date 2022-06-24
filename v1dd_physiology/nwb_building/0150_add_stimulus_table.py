import os, h5py
import numpy as np
import pandas as pd
import NeuroAnalysisTools.core.FileTools as ft
import NeuroAnalysisTools.NwbTools as nt
import utils as utils

nwb_folder = r"\\allen\programs\mindscope\workgroups\surround" \
             r"\v1dd_in_vivo_new_segmentation\data\nwbs"
             
stim_table_p = r"\\allen\programs\mindscope\workgroups\surround" \
               r"\v1dd_in_vivo_new_segmentation\data\stimulus_tables" \
               r"\stim_table.hdf5"

curr_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(curr_folder)

st_f = h5py.File(stim_table_p, 'r')

fns = ft.look_for_file_list(
    source=nwb_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

# fns = [fn for fn in fns if fn[9] in '12345']
fns.sort()
print('\n'.join(fns))

for fi, fn in enumerate(fns):

    print('processing {}, {} / {}'.format(fn, fi + 1, len(fns)))

    sess_n = fn[0:10]

    #spontaneous
    spt = pd.read_hdf(stim_table_p,
                      key=f'{sess_n}/spontaneous',
                      mode='r')

    dgf = pd.read_hdf(stim_table_p,
                      key=f'{sess_n}/drifting_gratings_full',
                      mode='r')

    dgw = pd.read_hdf(stim_table_p,
                      key=f'{sess_n}/drifting_gratings_windowed',
                      mode='r')
    dgw_alt = st_f[f'{sess_n}/drifting_gratings_windowed'].attrs['center_alt']
    dgw_azi = st_f[f'{sess_n}/drifting_gratings_windowed'].attrs['center_azi']
    dgw_dia = st_f[f'{sess_n}/drifting_gratings_windowed'].attrs['diameter_deg']

    lsn = pd.read_hdf(stim_table_p,
                      key=f'{sess_n}/locally_sparse_noise',
                      mode='r')

    ni = pd.read_hdf(stim_table_p,
                     key=f'{sess_n}/natural_images',
                     mode='r')

    ni12 = pd.read_hdf(stim_table_p,
                       key=f'{sess_n}/natural_images_12',
                       mode='r')

    nm = pd.read_hdf(stim_table_p,
                     key=f'{sess_n}/natural_movie',
                     mode='r')

    # print(nm)

    nwb_f = nt.RecordedFile(os.path.join(nwb_folder, fn))

    spt_ts = nwb_f.create_timeseries('TimeSeries', 'spontaneous', 'stimulus')
    spt_ts.set_time(spt['Start'])
    spt_ts.set_data(spt.to_numpy(), unit='', conversion=np.nan, resolution=np.nan)
    spt_ts.set_value('column_names', [n.encode('utf-8') for n in spt.columns])
    spt_ts.set_value('duration_sec', 300.)
    spt_ts.set_description('mean gray for spontaneous activities.')
    spt_ts.finalize()

    dgf_ts = nwb_f.create_timeseries('TimeSeries', 'drifting_gratings_full', 'stimulus')
    dgf_ts.set_time(dgf['Start'])
    dgf_ts.set_data(dgf.to_numpy(), unit='', conversion=np.nan, resolution=np.nan)
    dgf_ts.set_value('column_names', [n.encode('utf-8') for n in dgf.columns])
    dgf_ts.set_value('TF_Hz', 1.)
    dgf_ts.set_value('contrast', 0.8)
    dgf_ts.set_value('duration_sec', 2)
    dgf_ts.set_value('iteartion', 8)
    dgf_ts.set_description('Full field drifting grating, 2 seconds on, 1 second gap in between')
    dgf_ts.finalize()

    dgw_ts = nwb_f.create_timeseries('TimeSeries', 'drifting_gratings_windowed', 'stimulus')
    dgw_ts.set_time(dgw['Start'])
    dgw_ts.set_data(dgw.to_numpy(), unit='', conversion=np.nan, resolution=np.nan)
    dgw_ts.set_value('column_names', [n.encode('utf-8') for n in dgw.columns])
    dgw_ts.set_value('TF_Hz', 1.)
    dgw_ts.set_value('contrast', 0.8)
    dgw_ts.set_value('duration_sec', 2)
    dgw_ts.set_value('iteartion', 8)
    dgw_ts.set_value('center_alt', dgw_alt)
    dgw_ts.set_value('center_azi', dgw_azi)
    dgw_ts.set_value('diameter_deg', dgw_dia)
    dgw_ts.set_description('Windowed drifting grating, 2 seconds on, 1 second gap in between')
    dgw_ts.finalize()

    lsn_ts = nwb_f.create_timeseries('TimeSeries', 'locally_sparse_noise', 'stimulus')
    lsn_ts.set_time(lsn['Start'])
    lsn_ts.set_data(lsn.to_numpy(), unit='', conversion=np.nan, resolution=np.nan)
    lsn_ts.set_value('column_names', [n.encode('utf-8') for n in lsn.columns])
    lsn_ts.set_value('size_deg', 9.)
    lsn_ts.set_value('duration_sec', 0.33)
    lsn_ts.set_description('Locally sparse noise, frame is the frame id of the '
                           'presaved tiff file. Frames in this file was displayed '
                           'at 3 Hz.')
    lsn_ts.finalize()

    ni_ts = nwb_f.create_timeseries('TimeSeries', 'natural_images', 'stimulus')
    ni_ts.set_time(ni['Start'])
    ni_ts.set_data(ni.to_numpy(), unit='', conversion=np.nan, resolution=np.nan)
    ni_ts.set_value('column_names', [n.encode('utf-8') for n in ni.columns])
    ni_ts.set_value('iteration', 8)
    ni_ts.set_value('duration_sec', 0.33)
    ni_ts.set_description('Natural images with 118 images presented in a random order, '
                          'but fixed with two different seeds. So for each iteration '
                          'the frame indices range from 0-235. Each of the two seeded '
                          'sequences was presented 4 times. The images were presented '
                          'at 3 Hz.')
    ni_ts.finalize()

    ni12_ts = nwb_f.create_timeseries('TimeSeries', 'natural_images_12', 'stimulus')
    ni12_ts.set_time(ni12['Start'])
    ni12_ts.set_data(ni12.to_numpy(), unit='', conversion=np.nan, resolution=np.nan)
    ni12_ts.set_value('column_names', [n.encode('utf-8') for n in ni12.columns])
    ni12_ts.set_value('iteration', 40)
    ni12_ts.set_value('duration_sec', 0.33)
    ni12_ts.set_description('Natural images with 12 images presented in a fixed order, '
                            'The whole sequence was presented 40 times and the images '
                            'were presented at 3 Hz.')
    ni12_ts.finalize()

    nm_ts = nwb_f.create_timeseries('TimeSeries', 'natural_movie', 'stimulus')
    nm_ts.set_time(nm['Start'])
    nm_ts.set_data(nm.to_numpy(), unit='', conversion=np.nan, resolution=np.nan)
    nm_ts.set_value('column_names', [n.encode('utf-8') for n in nm.columns])
    nm_ts.set_value('frame_rate_Hz', 30.)
    nm_ts.set_description('Natural movie. 2 minutes clip 8 repeats or 0.5 minutes'
                          'clip 30 repeats. Frames of the movie were displayed at '
                          '30 Hz.')
    nm_ts.finalize()

    nwb_f.close()
