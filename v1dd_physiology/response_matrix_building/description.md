## Description
The scripts in this folder calculate response matrix for each roi from 
nwb files. 

In most cases, the response matrix means stimulus triggered average, and 
it is saved as 3d array ("roi x trial x timepoint"). The global trigger 
time for each specific stimulus can be found in the group attributes.

Response matrices are only extracted from dF/F and events traces.

## Workflow  

  1. run `0010_get_strf_slurm.py` on slurm cluster
  2. run `0060_get_dgcrm_w_slurm.sh` on slurm cluster
  3. run `0070_get_dgcrm_f_slurm.sh` on slurm cluster
  4. run `0090_get_greedy_rf_dff_slurm.sh` on slurm cluster
  5. run `0100_get_greedy_rf_events_slurm.sh` on slurm cluster