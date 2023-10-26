#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################

#@BATCH_JOBNAME@RUN_N
#@BATCH_TIME@RUN_T
#@RUN_P
#@RUN_Q
#@BATCH_GROUP
#@BATCH_JOINOUTERR

umask 022

limit stacksize unlimited
limit coredumpsize 0

#######################################################################
#                     Set Model Run Parameters
#######################################################################

set NH    = @USE_NONHYDRO
set FV_NX = @FV_NX
set FV_NY = @FV_NY
set NX = $FV_NX
set NY = $FV_NY
@  MODEL_NPES = $NX * $NY

set RUN_DT = @DT

set AGCM_IM = @AGCM_IM
set AGCM_JM = `expr $AGCM_IM \* 6`
set AGCM_LM = @AGCM_LM

# Set OpenMP Threads
set N_OMP = 1

# Set number of tracers
set N_TRACERS = 2
# Q is always included, so pass N-1 tracers to the layout
@ N_M1_TRACERS = $N_TRACERS - 1

set FV3_NPX = `expr $AGCM_IM \+ 1`
set FV3_NPY = `expr $AGCM_IM \+ 1`
set FV3_NPZ = ${AGCM_LM}

set HIST_IM = @HIST_IM
set HIST_JM = @HIST_JM
set OUTPUT_GRID = "PC${HIST_IM}x${HIST_JM}-DC"

set BEG_DATE = '18910301 000000'
set END_DATE = '29990311 000000'
set JOB_SGMT = '00000001 000000'

set USE_SHMEM    = @USE_SHMEM
set USE_IOSERVER = @USE_IOSERVER

if ($USE_IOSERVER == 0) then
   set IOS_NODES = 0
else
   set IOS_NODES = @IOS_NDS
endif

# Check for Over-Specification of CPU Resources
# ---------------------------------------------
if ($?SLURM_NTASKS) then
   set  NCPUS = $SLURM_NTASKS
else if ($?PBS_NODEFILE) then
   set  NCPUS = `cat $PBS_NODEFILE | wc -l`
else
   set  NCPUS = NULL
endif

if ( $NCPUS != NULL ) then

   if ( $USE_IOSERVER == 1 ) then

      set NCPUS_PER_NODE = @NCPUS_PER_NODE

      @ NODES  = `echo "( ($MODEL_NPES + $NCPUS_PER_NODE) + ($IOS_NODES * $NCPUS_PER_NODE) - 1)/$NCPUS_PER_NODE" | bc`
      @ NPES   = $NODES * $NCPUS_PER_NODE

      if( $NPES > $NCPUS ) then
         echo "CPU Resources are Over-Specified"
         echo "--------------------------------"
         echo "Allotted  NCPUs: $NCPUS"
         echo "Requested NCPUs: $NPES"
         echo ""
         echo "Specified NX: $NX"
         echo "Specified NY: $NY"
         echo ""
         echo "Specified IOSERVER_NODES: $IOS_NODES"
         echo "Specified cores per node: $NCPUS_PER_NODE"
         exit
      endif

   else

      @ NPES = $MODEL_NPES

      if( $NPES > $NCPUS ) then
         echo "CPU Resources are Over-Specified"
         echo "--------------------------------"
         echo "Allotted  NCPUs: $NCPUS"
         echo "Requested NCPUs: $NPES"
         echo ""
         echo "Specified NX: $NX"
         echo "Specified NY: $NY"
         exit
      endif

   endif

else
   # This is for the desktop path

   @ NPES = $MODEL_NPES

endif

#######################################################################
#             Experiment Specific Environment Variables
#######################################################################

setenv GEOSDIR @GEOSDIR
setenv GEOSBIN @GEOSBIN
setenv GEOSETC @GEOSETC

set TAG = `cat $GEOSETC/.FV3_VERSION`
set RUN_CMD = "$GEOSBIN/esma_mpirun -np "

setenv ARCH `uname`
source $GEOSBIN/g5_modules
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BASEDIR}/${ARCH}/lib:${GEOSDIR}/lib
module list

#########################################
# Set MPI Environment Variables
#########################################

@SETENVS

         setenv EXETAG "$TAG"
         setenv FV3PRC '@FV_PRECISION'
         setenv EXPID   "c${AGCM_IM}_L${AGCM_LM}_T${N_TRACERS}_${NX}x${NY}_${N_OMP}threads"
         setenv EXPDSC  "c${AGCM_IM}_L${AGCM_LM}_T${N_TRACERS}_${NX}x${NY}_${N_OMP}threads"
         setenv EXPDIR  @EXPDIR
         setenv SCRDIR  $EXPDIR/scratch_${EXPID}_${EXETAG}-${FV3PRC}
if ($NH) setenv SCRDIR  ${SCRDIR}_NH.$$

#######################################################################
#                 Create Experiment Scratch-Directory
#######################################################################

if (! -e $SCRDIR            ) mkdir -p $SCRDIR
cd $SCRDIR
/bin/rm -rf *

#######################################################################
#          Create LIVE RC Files from Templates for Current Run
#######################################################################

/bin/rm -f ioserver.nml
cat >   ioserver.nml << EOF
&ioserver
/
EOF

/bin/rm -f ExtData.rc
cat >      ExtData.rc << EOF
USE_EXTDATA: .false.
EOF

/bin/rm -f logging.yaml
cat >      logging.yaml << EOF
schema_version: 1
locks:
   mpi: {class: MpiLock, comm: MPI_COMM_WORLD}
formatters:
   plain: {class: Formatter, format: '%(message)a'}
   basic: {class: Formatter, format: '%(short_name)a15~: %(level_name)a~: %(message)a'}
   mpi: {class: MpiFormatter, format: '%(mpi_rank)i4.4~: %(name)~: %(level_name)a~: %(message)a', comm: MPI_COMM_WORLD}
handlers:
   console: {class: streamhandler, formatter: basic, unit: OUTPUT_UNIT, level: INFO}
   console_plain: {class: streamhandler, formatter: plain, unit: OUTPUT_UNIT, level: INFO}
   warnings: {class: FileHandler, filename: warnings_and_errors.log, lock: mpi, level: WARNING, formatter: basic}
   errors: {class: StreamHandler, formatter: basic, unit: ERROR_UNIT, level: ERROR}
   mpi_shared: {class: FileHandler, filename: allPEs.log, formatter: mpi, comm: MPI_COMM_WORLD, lock: mpi, rank_keyword: rank, level: DEBUG}
root:
   handlers: [warnings, errors, console]
   level: WARNING
   root_level: WARNING
loggers:
   errors: {handlers: [errors], level: ERROR}
   CAP: {level: WARNING, root_level: INFO}
   MAPL: {handlers: [mpi_shared], level: WARNING, root_level: INFO}
   MAPL.profiler: {handlers: [console_plain], propagate: FALSE, level: WARNING, root_level: INFO}
EOF

/bin/rm -f CAP.rc
cat >      CAP.rc << EOF
MAPLROOT_COMPNAME: DYN
        ROOT_NAME: DYN
ROOT_CF: AGCM.rc
HIST_CF: HISTORY.rc
LATLON: 0
BEG_DATE:     ${BEG_DATE}
END_DATE:     ${END_DATE}
JOB_SGMT:     ${JOB_SGMT}
NUM_SGMT:     1
HEARTBEAT_DT: ${RUN_DT}
USE_SHMEM: ${USE_SHMEM}
MAPL_ENABLE_TIMERS: YES
MAPL_ENABLE_MEMUTILS: NO
PRINTSPEC: 0  # (0: OFF, 1: IMPORT & EXPORT, 2: IMPORT, 3: EXPORT)
EOF

/bin/rm -f AGCM.rc
cat >      AGCM.rc << EOF
# Model Resolution and Timestep Parameters
# ----------------------------------------
            NX: ${NX}
            NY: ${NY}
            IM: ${AGCM_IM}
            JM: ${AGCM_JM}
            LM: ${AGCM_LM}
       AGCM_IM: ${AGCM_IM}
       AGCM_JM: ${AGCM_JM}
       AGCM_LM: ${AGCM_LM}
IOSERVER_NODES: @IOS_NDS
      GRIDNAME: PE${AGCM_IM}x${AGCM_JM}-CF
 DYN.GRID_TYPE: Cubed-Sphere
  DYN.GRIDNAME: PE${AGCM_IM}x${AGCM_JM}-CF
        DYN.NF: 6
        DYN.LM: ${AGCM_LM}
  DYN.IM_WORLD: ${AGCM_IM}
FV3_STANDALONE: 1
        DYCORE: FV3
     COLDSTART: 1
       CASE_ID: 1
  CASE_TRACERS: ${N_M1_TRACERS}
  HEARTBEAT_DT: ${RUN_DT}
     ADIABATIC: 1
        FV_OFF: 0
      fix_mass: 0
  RECORD_FINAL: NO
# Set the number of parallel I/O processes to use when
# RESTART_TYPE and or CHECKPOINT_TYPE are set to pbinary or pnc4
#---------------------------------------------------------------
NUM_READERS: @NUM_READERS
NUM_WRITERS: @NUM_WRITERS
# AGCM Model Restart Files
# ------------------------
DYN_INTERNAL_CHECKPOINT_FILE: fvcore_internal_checkpoint
DYN_INTERNAL_CHECKPOINT_TYPE: pnc4
DYN_INTERNAL_HEADER:          1
EOF

/bin/rm -f HISTORY.rc
cat >      HISTORY.rc << EOF
VERSION: 1
EXPID:  ${EXPID}
EXPDSC: ${EXPDSC}

COLLECTIONS:  'inst3_3d_diag'
              'inst1_2d_diag'
              ::

GRID_LABELS: PC@HIST_IMx@HIST_JM-DC
::

PC@HIST_IMx@HIST_JM-DC.GRID_TYPE: LatLon
PC@HIST_IMx@HIST_JM-DC.IM_WORLD: @HIST_IM
PC@HIST_IMx@HIST_JM-DC.JM_WORLD: @HIST_JM
PC@HIST_IMx@HIST_JM-DC.POLE: PC
PC@HIST_IMx@HIST_JM-DC.DATELINE: DC
PC@HIST_IMx@HIST_JM-DC.LM: @AGCM_LM


  inst1_2d_diag.format:    'CFIO',
  inst1_2d_diag.template:  '%y4%m2%d2_%h2%n2z.nc4',
  inst1_2d_diag.mode:      'instantaneous',
  inst1_2d_diag.grid_label: PC@HIST_IMx@HIST_JM-DC,
  inst1_2d_diag.frequency:  010000,
  inst1_2d_diag.duration:   010000,
  inst1_2d_diag.fields:     'PS'       , 'DYN'           ,
                            'PHIS'     , 'DYN'           ,
                            'DXC'      , 'DYN'           ,
                            'DYC'      , 'DYN'           ,
                            'AREA'     , 'DYN'           ,
                            ::

  inst3_3d_diag.format:    'CFIO',
  inst3_3d_diag.template:  '%y4%m2%d2_%h2%n2z.nc4',
  inst3_3d_diag.mode:      'instantaneous',
  inst3_3d_diag.frequency:  030000,
  inst3_3d_diag.duration:   030000,
  inst3_3d_diag.vscale:     100.0,
  inst3_3d_diag.vunit:      'hPa',
  inst3_3d_diag.vvars:      'log(PLE)' , 'DYN'          ,
  inst3_3d_diag.levels:      1000 975 950 925 900 875 850 825 800 775 750 725 700 650 600 550 500 450 400 350 300 250 200 150 100,
  inst3_3d_diag.grid_label: PC@HIST_IMx@HIST_JM-DC,
  inst3_3d_diag.fields:     'PS'       , 'DYN'           ,
                            'PHIS'     , 'DYN'           ,
                            'U'        , 'DYN'           ,
                            'V'        , 'DYN'           ,
                            'T'        , 'DYN'           ,
                            'AREA'     , 'DYN'           ,
                            ::
EOF

         set hydrostatic='.true.'
if ($NH) set hydrostatic='.false.'

set GRID_INPUT = "'INLINE'"

/bin/rm -f input.nml
cat >      input.nml << EOF
&fv_core_nml
       npx = ${FV3_NPX}
       npy = ${FV3_NPX}
       npz = ${FV3_NPZ}
       adiabatic = .true.
       hydrostatic = ${hydrostatic}
       make_nh = .T.
       fv_debug = .F.
       fv_sg_adj = -1
       n_sponge = -1
       n_zfilter = 0
       compute_coords_locally = .false.
/

&fv_grid_nml
/
#       grid_file = $GRID_INPUT

&main_nml
/

&test_case_nml
       test_case = 5
/

&fms_io_nml
/

&fms_nml
        print_memory_usage=.false.
        domains_stack_size = 24000000
/
EOF

setenv OMP_NUM_THREADS $N_OMP
echo "OMP_NUM_THREADS $OMP_NUM_THREADS"
if ($N_OMP > 1) then
   setenv I_MPI_PIN_DOMAIN omp
   setenv KMP_AFFINITY compact
   setenv KMP_STACKSIZE 16m
endif

@SINGULARITY_BUILD  #######################################################################
@SINGULARITY_BUILD  #          Settings for Singularity - EXPERIMENTAL
@SINGULARITY_BUILD  #######################################################################
@SINGULARITY_BUILD
@SINGULARITY_BUILD  # Detect if StandAlone_FV3_Dycore.x is in the current directory
@SINGULARITY_BUILD  # -------------------------------------------------------------
@SINGULARITY_BUILD
@SINGULARITY_BUILD  # If you are using singularity, set the path to the singularity sandbox here
@SINGULARITY_BUILD  setenv SINGULARITY_SANDBOX @SINGULARITY_SANDBOX
@SINGULARITY_BUILD
@SINGULARITY_BUILD  echo "We are running under Singularity"
@SINGULARITY_BUILD  echo ""
@SINGULARITY_BUILD
@SINGULARITY_BUILD  # Detect if StandAlone_FV3_Dycore.x is in the experiment directory
@SINGULARITY_BUILD  if (-e $EXPDIR/StandAlone_FV3_Dycore.x) then
@SINGULARITY_BUILD     echo "Found StandAlone_FV3_Dycore.x in $EXPDIR"
@SINGULARITY_BUILD
@SINGULARITY_BUILD     echo "NOTE: Testing has shown Singularity only works when running with"
@SINGULARITY_BUILD     echo "      the StandAlone_FV3_Dycore.x executable directly from the installation bin directory"
@SINGULARITY_BUILD     echo ""
@SINGULARITY_BUILD     echo "      So, we will *ignore* the local StandAlone_FV3_Dycore.x and "
@SINGULARITY_BUILD     echo "      instead use $GEOSBIN/StandAlone_FV3_Dycore.x"
@SINGULARITY_BUILD     echo ""
@SINGULARITY_BUILD  endif
@SINGULARITY_BUILD
@SINGULARITY_BUILD  echo "Using StandAlone_FV3_Dycore.x from $GEOSBIN"
@SINGULARITY_BUILD  setenv FV3EXE $GEOSBIN/StandAlone_FV3_Dycore.x
@SINGULARITY_BUILD
@SINGULARITY_BUILD  # Load the Singularity module
@SINGULARITY_BUILD  module load singularity
@SINGULARITY_BUILD
@SINGULARITY_BUILD  # Set Singularity Bind Paths. Note: These are dependent on where you are running.
@SINGULARITY_BUILD  # By default, we'll assume you are running this script from NOBACKUP
@SINGULARITY_BUILD  setenv REAL_BIND_PATH @REAL_BIND_PATH
@SINGULARITY_BUILD  setenv SINGULARITY_BIND_PATH "-B ${NOBACKUP}:${NOBACKUP},${REAL_BIND_PATH}:${REAL_BIND_PATH}"
@SINGULARITY_BUILD
@SINGULARITY_BUILD  # Set a variable to encapsulate all Singularity details
@SINGULARITY_BUILD  setenv SINGULARITY_RUN "singularity exec $SINGULARITY_BIND_PATH $SINGULARITY_SANDBOX"

@NATIVE_BUILD  # Detect if StandAlone_FV3_Dycore.x is in the experiment directory
@NATIVE_BUILD  if (-e $EXPDIR/StandAlone_FV3_Dycore.x) then
@NATIVE_BUILD     echo "Found StandAlone_FV3_Dycore.x in $EXPDIR"
@NATIVE_BUILD
@NATIVE_BUILD     echo "Copying $EXPDIR/StandAlone_FV3_Dycore.x to $SCRDIR"
@NATIVE_BUILD     /bin/cp $EXPDIR/StandAlone_FV3_Dycore.x $SCRDIR/StandAlone_FV3_Dycore.x
@NATIVE_BUILD
@NATIVE_BUILD     setenv FV3EXE $SCRDIR/StandAlone_FV3_Dycore.x
@NATIVE_BUILD  else
@NATIVE_BUILD     echo "Using StandAlone_FV3_Dycore.x from $GEOSBIN"
@NATIVE_BUILD
@NATIVE_BUILD     setenv FV3EXE $GEOSBIN/StandAlone_FV3_Dycore.x
@NATIVE_BUILD  endif

#######################################################################
#                          Run the Model
#######################################################################
echo "  "
echo "***** USING **** $FV3EXE *********************"

if( $USE_SHMEM == 1 ) $GEOSBIN/RmShmKeys_sshmpi.csh >& /dev/null

if( $USE_IOSERVER == 1) then
   set IOSERVER_OPTIONS = "--npes_model $MODEL_NPES --nodes_output_server $IOS_NODES"
else
   set IOSERVER_OPTIONS = ""
endif

@SINGULARITY_BUILD  $RUN_CMD $NPES $SINGULARITY_RUN $FV3EXE $IOSERVER_OPTIONS --logging_config logging.yaml |& tee ${SCRDIR}.log
@NATIVE_BUILD  $RUN_CMD $NPES $FV3EXE $IOSERVER_OPTIONS --logging_config logging.yaml |& tee ${SCRDIR}.log

if( $USE_SHMEM == 1 ) $GEOSBIN/RmShmKeys_sshmpi.csh >& /dev/null

set rc =  $status
echo       Status = $rc
if ($rc != 0) then
  echo This Job Failed
  RmShmKeys >& /dev/null
  exit 0
endif

exit 0
