"""
Reza's stimulus table does not have position of
windowed drifting gratings. This script will, for
each session, go to original dislplay log .pkl
file and find the location of the windowed drifting
gratings and add it to the converted stimulus table
hdf5 file
"""

import os, h5py
import NeuroAnalysisTools.core.FileTools as ft
import utils as utils

st_path = r"\\allen\programs\mindscope\workgroups\surround" \
          r"\v1dd_in_vivo_new_segmentation\data\stimulus_tables" \
          r"\3p_stim_table.hdf5"

curr_folder = os.path.dirname(os.path.realpath(__file__))
os.chdir(curr_folder)

st_f = h5py.File(st_path, 'a')
sess_ns = list(st_f.keys())
sess_ns.sort()

for sess_i, sess_n in enumerate(sess_ns):

    print('processing {}, {} / {} ...'.format(sess_n, sess_i + 1, len(sess_ns)))

    sess_path = utils.get_lims_session_path_from_session_name(
        session_name=sess_n
        )

    log_path = ft.look_for_unique_file(source=sess_path,
                                       identifiers=['_stim.pkl'],
                                       file_type='pkl',
                                       is_full_path=True)

    # print(log_path)
    log = ft.loadFile(log_path)
    stim_text = log[b'stimuli'][0][b'stim'].decode('utf-8')
    start_ind = stim_text.index('pos=array') + 11
    stim_text = stim_text[start_ind:]
    end_ind = stim_text.index(']), ')
    stim_text = stim_text[0: end_ind]
    mid_ind = stim_text.index(', ')
    azi = float(stim_text[0: mid_ind])
    alt = float(stim_text[mid_ind+2:])
    print((azi, alt))

    dgw_grp = st_f[sess_n]['drifting_gratings_windowed']
    dgw_grp.attrs['center_azi'] = azi
    dgw_grp.attrs['center_alt'] = alt
    dgw_grp.attrs['diameter_deg'] = 30

st_f.close()