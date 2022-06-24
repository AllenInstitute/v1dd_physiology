
This folder contains code to build nwbs database.

It heavily relied on the [NeuroAnalysisTools](https://github.com/zhuangjun1981/NeuroAnalysisTools) package. The most frequently used class is the 

`NeuroAnalysisTools.NwbTools.RecordedFile` class which is a subclass of `NeuroAnalysisTools.external.nwb.NWB` class which follows NWB 1.0 standard.

Once nwb/respone_table files are built, the user does not need to interact with these codes and they codes should be just viewed as record keeping.

Importantly, the `meta_lims` folder contains the meta data of all mice which were used to query data for nwb building.

The workflow is the run the list of scripts in order to build the database step by step. The job of each script is described in its title. Some odd thing are listed below.

 * the `eyetracking` folder contains the script to preprocess the eyetracking movies
 * `0045_remove_problematic_nwbs_3p.py` remove "M427836\_1e" and "M427836_1f" sessions, since their stimulus table can not be found from Reza Abbasi Asl's preprocessing folder.
 * `0120_get_stim_table_path.py` get the path to the stimulus table from Reza Abbasi Asl's preprocessing folder for each session.
 * `0130_conver_stimulus_table.py` convert the stimulus table (in the format of python2 pickle file) into hdf5 format.
 * the `*add_rois_and_traces*` and `*event_detection*` scripts can be run on slurm server by running the correspondent .sh shell script on "hpc-login".
 * the `utils` file contains many functions which are called by the scripts in this folder.