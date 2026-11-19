#!/bin/bash
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=12:00:00
#PBS -J 1-110


## Load python module:
module load tools/prod
module load MATLAB/2023a_Update_3

## Copy files and parameters:
##pwd

cd ${PBS_O_WORKDIR}

cp -r src/ $TMPDIR
cp -r resources $TMPDIR
cp process_results_HPCv4.py $TMPDIR
cd $TMPDIR


pwd

##cp -r $PBS_O_WORKDIR/src .
##cp -r $PBS_O_WORKDIR/resources .
mkdir -p outputs

cd src

## Write the ${PBS_ARRAY_INDEX} (which becomes the country ISO) here:
echo ${PBS_ARRAY_INDEX} > countries_to_run.txt

ls -l


##echo "Hello world" > outputs/temp1.txt

matlab -nosplash -nodisplay -nojvm -batch main_script -logfile myTestRun_output.log


# Make directory if it does not exist:
mkdir -p ${PBS_O_WORKDIR}/outputs
cd ..
#python3 process_results_HPCv2.py outputs/${PBS_ARRAY_INDEX}
python3 process_results_HPCv4.py outputs/ ${PBS_ARRAY_INDEX}

## Move the output files to my home dir:
##mv outputs/${PBS_ARRAY_INDEX}/summary_outputs*.csv ${PBS_O_WORKDIR}/outputs/
mv outputs/summary_outputs*.csv ${PBS_O_WORKDIR}/outputs/
mv outputs/elim_outputs*.csv ${PBS_O_WORKDIR}/outputs/
#mv outputs/BD_coverage*.csv ${PBS_O_WORKDIR}/outputs/
#mv outputs/HepB3_coverage*.csv ${PBS_O_WORKDIR}/outputs/
#mv outputs/PAP_coverage*.csv ${PBS_O_WORKDIR}/outputs/

now=$(date +"%d_%m_%Y")
mkdir -p $EPHEMERAL/outputs\_$now/${PBS_ARRAY_INDEX}
mv outputs/* $EPHEMERAL/outputs\_$now/${PBS_ARRAY_INDEX}

