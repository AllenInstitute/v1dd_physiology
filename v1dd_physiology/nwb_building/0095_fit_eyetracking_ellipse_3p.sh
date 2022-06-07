#!/bin/bash
#SBATCH -J v1dd_event    # Job name
#SBATCH -N 1   # number of nodes
#SBATCH -c 32  # number of cores (per node?)
#SBATCH --mem=16gb                         # Job memory request (per node)
#SBATCH --time=10:00:00                     # Time limit hrs:min:sec
#SBATCH --array=0-32 # 0-32 

#SBATCH --partition braintv                 # Partition used for processing
#SBATCH --mail-user=junz@alleninstitute.org # Where to send mail  
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)

# %A" is replaced by the job ID and "%a" with the array index
#SBATCH -o /allen/programs/mindscope/workgroups/surround/jun_testing/eyetracking_output/ellipse_%A_%a.out
#SBATCH -e /allen/programs/mindscope/workgroups/surround/jun_testing/eyetracking_output/ellipse_%A_%a.err
 
pwd; hostname; date

search_dir=/allen/programs/mindscope/workgroups/surround/v1dd_in_vivo_new_segmentation/data/eye_tracking_movies/3p

h5_list=("$search_dir"/*.h5)
echo "${h5_list[@]}"

# h5_list=("${h5_list[@]:1:1}")
# echo "${h5_list[@]}"

a1=${h5_list[`expr $SLURM_ARRAY_TASK_ID`]}

# source activate analysis
echo "Run eyetracking ellipse fitting on multiple nodes"
/home/junz/anaconda3/envs/analysis/bin/python /allen/programs/mindscope/workgroups/surround/v1dd_in_vivo_new_segmentation/v1dd_physiology/v1dd_physiology/nwb_building/0095_fit_eyetracking_ellipse_3p_slurm.py $a1

date