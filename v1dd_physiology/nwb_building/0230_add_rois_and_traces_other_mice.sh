#!/bin/bash
#SBATCH -J v1dd_event    # Job name
#SBATCH -N 1   # number of nodes
#SBATCH -c 35  # number of cores (per node?)
#SBATCH --mem=16gb                         # Job memory request (per node)
#SBATCH --time=10:00:00                     # Time limit hrs:min:sec
#SBATCH --array=0-35 

#SBATCH --partition braintv                 # Partition used for processing
#SBATCH --mail-user=junz@alleninstitute.org # Where to send mail  
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)

# %A" is replaced by the job ID and "%a" with the array index
#SBATCH -o /allen/programs/mindscope/workgroups/surround/jun_testing/slurm_output/add_roi&trace_%A_%a.out
#SBATCH -e /allen/programs/mindscope/workgroups/surround/jun_testing/slurm_output/add_roi&trace_%A_%a.err
 
pwd; hostname; date

search_dir=/allen/programs/mindscope/workgroups/surround/v1dd_in_vivo_new_segmentation/data/nwbs

# nwb_list=("$search_dir"/M416296*.nwb)
# nwb_list=("$search_dir"/M427836*.nwb)
nwb_list=("$search_dir"/M438833*.nwb)

# echo "${nwb_list[@]}"

# nwb_list=("${nwb_list[@]:1:1}")
# echo "${nwb_list[@]}"

a1=${nwb_list[`expr $SLURM_ARRAY_TASK_ID`]}

python_path="/home/junz/anaconda3/envs/analysis/bin/python"
script_path="/allen/programs/mindscope/workgroups/surround/v1dd_in_vivo_new_segmentation/v1dd_physiology/v1dd_physiology/nwb_building/0210_add_rois_and_traces_other_mice_slurm.py"
echo "adding rois and traces"
"$python_path" "$script_path" $a1

date