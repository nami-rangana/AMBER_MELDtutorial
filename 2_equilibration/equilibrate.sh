#!/bin/bash
#SBATCH --job-name=eqMELD
#SBATCH --error=equilibrate.err
#SBATCH --output=equilibrate.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.desilva@ufl.edu
#SBATCH --time=45:00 
#SBATCH --nodes=4
#SBATCH --ntasks=12 
#SBATCH --ntasks-per-node=3 
#SBATCH --gpus=4 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --partition=hpg-turin
#SBATCH --distribution=block:block
#SBATCH --account=alberto.perezant
#SBATCH --qos=alberto.perezant

date;hostname;pwd

cd $SLURM_SUBMIT_DIR

echo "========== JOB INFO =========="
echo "Job ID       : $SLURM_JOB_ID"
echo "Nodes        : $SLURM_JOB_NUM_NODES"
echo "MPI ranks    : $SLURM_NTASKS"
echo "Node list    : $SLURM_JOB_NODELIST"
echo "=============================="

module purge

# -- local amber --
#ml cuda/12.8.1 gcc/14.2.0 openmpi/5.0.10
ml cuda/13.2.1  gcc/14.2.0 openmpi/5.0.10
source /orange/alberto.perezant/t.desilva/amber/install/amber.sh

#export OMPI_MCA_coll_han_priority=0
#export UCX_TLS=tcp,sm,self

#mpirun -np $SLURM_NTASKS $AMBERHOME/bin/pmemd.cuda.MPI -ng 12 -groupfile equilibrate.groupfile
srun --mpi=pmix_v5 $AMBERHOME/bin/pmemd.cuda.MPI -ng 12 -groupfile equilibrate.groupfile
#
#
