import os
import h5py
import pickle
import pandas as pd

"""
the stim table from reza are .pkl files generated by python2
which cannot be read by python3. So this script uses python2
read this files and save them into a single .hdf5 file that 
is readable by python3.

for this a "2p_temp" environment was created.
"""

sess_p_path = r"\\allen\programs\mindscope\workgroups\surround" \
              r"\v1dd_in_vivo_new_segmentation\data\stimulus_tables\sess_paths.csv"

save_path = r"\\allen\programs\mindscope\workgroups\surround" \
            r"\v1dd_in_vivo_new_segmentation\data\stimulus_tables"

sess_p_df = pd.read_csv(sess_p_path, index_col=0)
# print(sess_p_df)

save_fp = os.path.join(save_path, 'stim_table.hdf5')
# save_f = h5py.File(save_fp, 'x')

for sess_i, sess_row in sess_p_df.iterrows():

    sess_n = sess_row['sess_name']
    print('processing {}, {} / {}'.format(sess_n, sess_i + 1, len(sess_p_df)))

    st_path = os.path.join(sess_row['sess_path'], 'stim_table.pkl')

    with open(st_path, 'rb') as file:
        st_dict = pickle.load(file)

    # save_grp = save_f.create_group(sess_n)
    for stim_n, stim_df in st_dict.items():
        # stim_grp = save_grp.create_group(stim_n)
        stim_df.to_hdf(save_fp, key='{}/{}'.format(sess_n, stim_n))