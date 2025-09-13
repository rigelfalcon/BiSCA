#!/bin/bash -l
#SBATCH --job-name=BxaBsTri15k24030116
#SBATCH --account=def-pvaldes # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=0-3:00:00         # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20      # adjust this if you are using parallel commands
#SBATCH --mem=80G             # adjust this according to the memory requirement per node you need
#SBATCH --mail-user=yingwangrigel@gmail.com # adjust this to match your email address
#SBATCH --mail-type=ALL


# Choose a version of MATLAB by loading a module:
module load StdEnv/2023
module load matlab/2023b
matlab -nodisplay -r "run('setup.m'); addpath('./example/iEEG/'); Ntask=1;itask=1;npertask=1772;taskname='BxaBsTri15k24030116';step3_hoxialpha_test_tri_3model_taskmng"




