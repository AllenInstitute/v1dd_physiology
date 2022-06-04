import os, h5py
import numpy as np
import pandas as pd
import NeuroAnalysisTools.core.FileTools as ft

nwb_folder = r"\\allen\programs\mindscope\workgroups\surround" \
             r"\v1dd_in_vivo_new_segmentation\data\nwbs"

stim_grp_n = 'drifting_gratings_windowed'

curr_folder = os.path.dirname(os.path.abspath(__file__))
os.chdir(curr_folder)

fns = ft.look_for_file_list(
    source=nwb_folder, 
    identifiers=[''], 
    file_type='nwb',
    is_full_path=False
    )

fns = [fn for fn in fns if fn[9] not in '12345']
fns.sort()
print('\n'.join(fns))

for fi, fn in enumerate(fns):

    print('processing {}, {} / {}'.format(fn, fi + 1, len(fns)))

    nwb_f = h5py.File(os.path.join(nwb_folder, fn), 'a')
    dgc_grp = nwb_f['analysis'].create_group(f'onsets_{stim_grp_n}')

    stim_df = pd.DataFrame(data=nwb_f[f'stimulus/presentation/{stim_grp_n}/data'][()],
                           columns=['Start', 'End', 'TF', 'SF', 'Ori'])

    dgc_alt = nwb_f[f'stimulus/presentation/{stim_grp_n}/center_alt'][()]
    dgc_azi = nwb_f[f'stimulus/presentation/{stim_grp_n}/center_azi'][()]
    dgc_rad = nwb_f[f'stimulus/presentation/{stim_grp_n}/diameter_deg'][()] / 2.

    dgc_ns = []
    for stim_i, stim_row in stim_df.iterrows():
        if np.isnan(stim_row['TF']): # blank trial
            dgc_ns.append(f'alt{dgc_alt:06.1f}'
                          f'_azi{dgc_azi:06.1f}'
                          f'_sf0.00'
                          f'_tf00.0'
                          f'_dire000'
                          f'_con0.80'
                          f'_rad{dgc_rad:03.0f}')
        else:
            dgc_ns.append(f'alt{dgc_alt:06.1f}'
                          f'_azi{dgc_azi:06.1f}'
                          f'_sf{stim_row["SF"]:04.2f}'
                          f'_tf{stim_row["TF"]:04.1f}'
                          f'_dire{stim_row["Ori"]:03.0f}'
                          f'_con0.80'
                          f'_rad{dgc_rad:03.0f}')

    stim_df['dgc_n'] = dgc_ns

    dgc_ns = list(set(dgc_ns))
    dgc_ns.sort()

    # print(stim_df)
    # _ = [print(n) for n in dgc_ns]

    iteration = nwb_f[f'stimulus/presentation/{stim_grp_n}/iteartion'][()]
    for dgc_n in dgc_ns:
        starts = np.array(stim_df[stim_df['dgc_n'] == dgc_n]['Start'])

        # the number of times display for each condition does not
        # always equal iteration !!!
        # if len(starts) != iteration:
        #     print(f'{dgc_n}, num of ts: {len(starts)}, iteration: {iteration}')

        condi_grp = dgc_grp.create_group(dgc_n)
        condi_grp.create_dataset('onset_ts_sec', data=starts)

    nwb_f.close()