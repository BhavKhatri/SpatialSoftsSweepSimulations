#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:ompthreads=32:mem=64g
#PBS -N SummaryStats_32_threads
## stdout/stderr (*.o<jobid>/*.e<jobid>) are written here instead of $PBS_O_WORKDIR.
## NB: this folder must already exist at submission time (PBS does not create it),
## so run `mkdir -p JOB_OUTPUTS` once in $PBS_O_WORKDIR before qsub-ing this script.
#PBS -o JOB_OUTPUTS/
#PBS -e JOB_OUTPUTS/

## CONTSIM, MU and M are captured at submission time, e.g.:
##   qsub -v CONTSIM="",MU="1.2e-7",M="256" ScriptSummaryStats32.pbs.sh
## Left unset, all three default to "" (empty), which SampleSummaryStats_Parallel_Script.jl
## treats as "use default behaviour" for ContinueSim, the Species-default μ, and the M
## hardcoded in that script (handy for quick local testing with a small M that never
## needs to reach the cluster's submitted jobs).

# cd SampleSummaryStats
# mkdir SummaryStats_N=1e4_1e11_D=0.01_1_s=1e-2_1_Area=1e4_1e7_Nsamp=5000_M=1000
# # cd SummaryStas_N=1e4_1e11_D=0.01_1_s=1e-2_1_Area=1e4_1e7_Nsamp=5000_M=1000
# cd ..

## NB values for ncpus and mem are allocated
## to each node (specified by select=N)
##
## to direct output to cwd, use $PBS_O_WORKDIR:
## specify LOGFILE found in ~/ during execution then moved to cwd on job completion
##
cd $PBS_O_WORKDIR
JOBNUM=`echo $PBS_JOBID | sed 's/\..*//'`
LOGFILE=${PBS_JOBNAME}.o${JOBNUM}

#########################################
##                                     ##
## Output some useful job information. ##
##                                     ##
#########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: job number is $JOBNUM
echo PBS: logfile is $LOGFILE
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
#echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

lscpu

echo ------------------------------------------------------
export PATH=$PATH:~/.juliaup/bin
export JULIA_NUM_THREADS=32

# cd RunOutputs

julia SampleSummaryStats_Parallel_Script.jl "$CONTSIM" "$MU" "$M"
# ~/julia-1.10.5/bin/julia < SampleSummaryStats_Parallel_Script.jl
# Median_Over_ReplicatesABC_parallelthreads_script.jl
