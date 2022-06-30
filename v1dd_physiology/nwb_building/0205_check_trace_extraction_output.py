import os
import utils
import NeuroAnalysisTools.core.FileTools as ft

base_folder = r"\\allen\programs\mindscope\workgroups\surround\trace_extraction_2022"

db_path = r"\\allen\programs\mindscope\workgroups\surround" \
          r"\v1dd_in_vivo_new_segmentation\data"

nwb_fns = ft.look_for_file_list(
    source=os.path.join(db_path, 'nwbs'),
    identifiers=[],
    file_type='nwb',
    is_full_path=False)

sessions = [n[0:10] for n in nwb_fns if '409828' not in n]

exps = []
for sess in sessions:
    curr_exps = utils.get_all_experiments_pika_meta(sess)
    for _, curr_exp in curr_exps.items():
        exps.append(curr_exp['experiment'])

# print(exps)

for exp in exps:

    flag = 0

    folder_traces = os.path.join(base_folder, 'traces_2022', exp)
    if not os.path.isdir(folder_traces):
        if flag == 0:
            print(f'\nexperiment: {exp}:')
            flag = 1
        print('\tDid not find traces folder.')
    else:
        traces = [n for n in os.listdir(folder_traces) if n == 'roi_traces.h5']
        if len(traces) == 0:
            if flag == 0:
                print(f'\nexperiment: {exp}:')
            flag = 1
            print('\tDid not find roi traces.')

        neuropil = [n for n in os.listdir(folder_traces) if n == 'neuropil_traces.h5']
        if len(neuropil) == 0:
            if flag == 0:
                print(f'\nexperiment: {exp}:')
            flag = 1
            print('\tDid not find neuropil traces.')


    folder_demixed = os.path.join(base_folder, 'demix_2022', exp)
    if not os.path.isdir(folder_demixed):
        if flag == 0:
            print(f'\nexperiment: {exp}:')
        flag = 1
        print('\tDid not find demixed folder.')
    else:
        demixed = [n for n in os.listdir(folder_demixed) if '_demixed_traces.h5' in n]
        if len(demixed) == 0:
            if flag == 0:
                print(f'\nexperiment: {exp}:')
            flag = 1
            print('\tDid not find demixed traces.')

    folder_subtracted = os.path.join(base_folder, 'neuropil_2022', exp)
    if not os.path.isdir(folder_subtracted):
        if flag == 0:
            print(f'\nexperiment: {exp}:')
        flag = 1
        print('\tDid not find neuropil folder.')
    else:
        subtracted = [n for n in os.listdir(folder_subtracted) if n == 'neuropil_correction.h5']
        if len(subtracted) == 0:
            if flag == 0:
                print(f'\nexperiment: {exp}:')
            flag = 1
            print('\tDid not find neuropil subtracted traces.')

    folder_dff = os.path.join(base_folder, 'dff_2022', exp)
    if not os.path.isdir(folder_dff):
        if flag == 0:
            print(f'\nexperiment: {exp}:')
        flag = 1
        print('\tDid not find dff folder.')
    else:
        dff = [n for n in os.listdir(folder_dff) if '_dff.h5' in n]
        if len(dff) == 0:
            if flag == 0:
                print(f'\nexperiment: {exp}:')
            flag = 1
            print('\tDid not find dff traces.')





