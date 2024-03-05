#!/bin/bash
#----------------------------------------------------------------------------
#- Laminar vortex shedding across a 2D cylinder at Re=100
#- Source      : Martin Einarsve, https://github.com/meinarsve/CFDwOpenFoam
#- Adapted by  : Nishant Kumar
#- Date        : 11/05/2022
#----------------------------------------------------------------------------
#SBATCH --job-name=laminarVortexShedding
#SBATCH --output=log.solve
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=0-12:00:00
#SBATCH --partition=imb-resources
#SBATCH --mem=93G
#SBATCH --ear=off
##SBATCH --mem-per-cpu=100

echo "#############################"
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "SLURM_TASKS_PER_NODE:" $SLURM_TASKS_PER_NODE
echo "#############################"

#- Ensure only owner can read the output
umask 0077

export SLURM_COMP_VERBOSE=3
export SLURM_LOADER_VERBOSE=3

#- Number of cores to perform POD
np=12
#- Mesh size identifier
ms="m0.5"

#- Run locally or on scratch (1: enabled, 0: disabled)
moveToScratch=0
#- If run on scratch, move everything to ${SLURM_SUBMIT_DIR} on finish (1: enabled, 0: disabled)
moveFromScratch=0
#- Cronjob to erase coordinates in postProcessing/internalField/*/cloud_*.xy (1: enabled, 0: disabled)
cronEraseFlag=0

#- Erase coordinates: Specify last time directory in postProcessing to ignore 
#-      'none'  : Run for all time steps
#-      'auto'  : Automatically detect last time directory and replace with a (int/float) value
#-      <value> : (int/float) Specify valid last time step from postProcessing/internalField
#- NOTE: Only applicable for cronEraseFlag=0.
lastTimeErase='none'

#- Source the bash profile and then call the appropriate OpenFOAM version function
#- so that all the modules and environment variables get set
source /gpfs/home/nkumar001/.bash_profile
#- Load modules
module purge
export MODULEPATH=$HOME/centos_7/easybuild/modules/all/Core:$MODULEPATH
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 OpenFOAM/2.4.0
source ${EBROOTOPENFOAM}/OpenFOAM-${EBVERSIONOPENFOAM}/etc/bashrc

dt=$(date '+%Y%m%d') # %H:%M:%S

#- Set directory where the job will run
if [ $moveToScratch -eq 1 ]; then
    #- Place the job in the scratch directory
    SCRDIR=/scratch/$USER/$(realpath ${SLURM_SUBMIT_DIR} | cut -d'/' -f5-).$ms.run${dt}
    mkdir -p $SCRDIR
    # 
    #- Copy the input files to ${TMP}
    echo "Copying from ${SLURM_SUBMIT_DIR}/ to ${SCRDIR}/"
    /usr/bin/rsync -arvx "${SLURM_SUBMIT_DIR}/" ${SCRDIR}/
    #
    cd $SCRDIR
# elif [ $moveToScratch -eq 0 ]; then
#     #- Run from current directory
#     cd ${0%/*} || exit 1
fi
parentDIR=$(pwd)

#- Automatically detect last time directory for coordinate erase
#- NOTE: Replace variable $lastTimeErase with the time value instead of creating a new variable
timeEraseFile=$parentDIR/postProcessing/internalField/listDif.tmp
if [ "$lastTimeErase" == "auto" -a -f "$timeEraseFile" ]; then
    lastTimeErase=$(tac $timeEraseFile |egrep -m 1 .)
fi

#- Scheduler for modifying output files
endFlag=0
stage="runTime"
eraseCmd="$parentDIR/system/sampling/eraseCoordinatesAuto.sh \
    $parentDIR $endFlag $lastTimeErase $cronEraseFlag \
    > $parentDIR/log.eraseCoordinates.$stage 2>&1"
if [ $cronEraseFlag -eq 1 ]; then
    cronJob="*/1 * * * * $eraseCmd"
    ( crontab -l | grep -v -F "$eraseCmd" ; echo "$cronJob" ) | crontab -
fi

echo "Starting OpenFOAM job at: " $(date) " using " $nprocs " cores"
start=`date +%s.%N`

#- Copy relevant pointCloud file
cp ${parentDIR}/system/sampling/pointCloud.${ms}.dat ${parentDIR}/system/sampling/pointCloud.dat

#- Launch using primitive commands
#--------------------------------------------------
#- Convert gmsh to OpenFOAM
#- NOTE: Modification in constant/polyMesh/boundary: 
#-      frontAndBack
#-          type            empty;
#-           // physicalType    patch;
boundaryFile="constant/polyMesh/boundary"
if [ -f "$boundaryFile" ]; then
    mv $boundaryFile ${boundaryFile}.bak
else
    echo ERROR: Failed to find ${boundaryFile}. Run gmshToFoam.
    exit 1
fi
gmshToFoam cylinder.${ms}.msh > log.gmshToFoam 2>&1

#- Restore constant/polyMesh/boundary
if [ -f $boundaryFile.bak ]; then
    mv ${boundaryFile}.bak $boundaryFile
fi

#- Decompose mesh
checkMesh > log.checkMesh.1 2>&1
renumberMesh -overwrite > log.renumberMesh 2>&1
decomposePar > log.decomposePar 2>&1
checkMesh > log.checkMesh.2 2>&1

#- Solve
mpirun -np $np icoFoam -parallel < /dev/null > log.icoFoam 2>&1
#--------------------------------------------------

end=`date +%s.%N`
echo "Ending OpenFOAM job at: " $(date)
echo "Runtime:" $( echo "$end - $start" | bc -l )

#- Erase coordinates
if [ $cronEraseFlag -eq 0 ]; then
    start=`date +%s.%N`
    $eraseCmd
    end=`date +%s.%N`
    echo "Time elapsed for erasing coordinates:" $( echo "$end - $start" | bc -l )
fi

#- Modify output files left out by scheduler and stop
if [ $cronEraseFlag -eq 1 ]; then
    endFlag=1
    stage="endTime"
    $eraseCmd

    #- Remove crontab entries
    # /usr/bin/crontab -r
    ( crontab -l | grep -v -F "$eraseCmd" ) | crontab -
fi

#- Move results
if [ $moveFromScratch -eq 1 ]; then
    #- Job done, copy everything back
    echo "Copying from ${SCRDIR}/ to ${SLURM_SUBMIT_DIR}/"
    /usr/bin/rsync -avx --exclude "*_0.gz" --exclude "phi*.gz" --exclude "ddt*.gz" ${SCRDIR}/ "${SLURM_SUBMIT_DIR}/"
      
    #- Delete temporary files
    [ $? -eq 0 ] && /bin/rm -rf ${SCRDIR}
fi

# ----------------------------------------------------------------- end-of-file