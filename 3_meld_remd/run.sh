#!/bin/bash
#SBATCH --job-name=cpu-MELD.job
#SBATCH --error=equilibrate.err
#SBATCH --output=equilibrate.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.desilva@ufl.edu
#SBATCH --time=5-00:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=hpg-turin
##SBATCH --gpus=1
#SBATCH --account=alberto.perezant
#SBATCH --qos=alberto.perezant
##SBATCH --distribution=cyclic:cyclic
date;hostname;pwd

cd $SLURM_SUBMIT_DIR

echo "Nodes     :   $SLURM_JOB_NUM_NODES"
echo "MPI ranks :   $SLURM_NTASKS"
#echo "GPUs      :   $SLURM_GPUS"

module purge

# -- local amber --
ml gcc/14.2.0 openmpi/5.0.7 #cuda/12.8.1
source /orange/alberto.perezant/t.desilva/amber_setup/pmemd24/amber.sh

#export UCX_NET_DEVICES="mlx5_4:1,mlx5_7:1,mlx5_8:1,mlx5_9:1,mlx5_10:1,mlx5_13:1,mlx5_14:1,mlx5_15:1"

#echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
#echo "SLURM_JOB_GPUS=$SLURM_JOB_GPUS"

#export UCX_TLS=rc,cuda_copy,cuda_ipc
#export OMPI_MCA_coll=^hcoll

#echo "Starting 2 replicas on 1 GPUs"

#mpirun -np 2 $PMEMDHOME/bin/pmemd.MPI -ng 2 -groupfile equilibrate.groupfile
#srun --mpi=${HPC_PMIX} $PMEMDHOME/bin/pmemd.MPI -groupfile equilibrate.groupfile
srun --mpi=pmix_v5 $PMEMDHOME/bin/pmemd.MPI -O -i meld.mdin -o meld.mdout -c ../gpu_eq/equilibrate.rst.001 -r meld.rst -x meld.nc -inf meld.mdinfo -p ../1_system_setup/3gb1.prmtop

