### Notes 8/28/15

To build executable
make ifort CORE=atmosphere

To clean executable
make clean CORE=atmosphere

To run 4 tasks w/ 16 threads per task
cd 120km
bsub < run_4task_16tpt

where run_4task_16tpt
```
#!/bin/csh
#
# LSF batch script to run an MPI application
#
#BSUB -P SCIS0005            # project code
#BSUB -W 0:15                # wall-clock time (hrs:mins)
#BSUB -n 4                   # number of tasks in job         
#BSUB -R "span[ptile=1]"     # run 1 MPI tasks per node
#BSUB -J run_4task_16tpt     # job name
#BSUB -o run.%J.out          # output file name in which %J is replaced by the job ID
#BSUB -e run.%J.err          # error file name in which %J is replaced by the job ID
#BSUB -q regular             # queue

limit stacksize unlimited
setenv OMP_NUM_THREADS 16
setenv MP_TASK_AFFINITY core:${OMP_NUM_THREADS}
setenv LD_LIBRARY_PATH ${NETCDF}/lib:$LD_LIBRARY_PATH

mpirun.lsf ../atmosphere_model
```

##### Only the following files contain OpenMP directives

core_atmosphere/mpas_atm_core.F
framework/mpas_threading.F
core_atmosphere/dynamics/mpas_atm_time_integration.F
core_atmosphere/physics/mpas_atmphys_driver.F
core_atmosphere/physics/mpas_atmphys_driver_microphysics.F

##### Thread loop bounds calculated in mpas_threading_init (framework/mpas_threading.F)

Call tree:
mpas() in driver/mpas.F
mpas_init() in driver/mpas_subdriver.F
mpas_threading_init() in framework/mpas_threading.F

```
#ifdef _OPENMP
!$OMP PARALLEL PRIVATE(threadid)
            threadid = OMP_get_thread_num()

            block % cellThreadStart(threadid+1) = (threadid * nCells / block % nThreads) + 1
            block % cellThreadEnd(threadid+1)   = ((threadid+1) * nCells / block % nThreads)
            block % cellSolveThreadStart(threadid+1) = (threadid * nCellsSolve / block % nThreads) + 1
            block % cellSolveThreadEnd(threadid+1)   = ((threadid+1) * nCellsSolve / block % nThreads)
            block % edgeThreadStart(threadid+1) = (threadid * nEdges / block % nThreads) + 1
            block % edgeThreadEnd(threadid+1)   = ((threadid+1) * nEdges / block % nThreads)
            block % edgeSolveThreadStart(threadid+1) = (threadid * nEdgesSolve / block % nThreads) + 1
            block % edgeSolveThreadEnd(threadid+1)   = ((threadid+1) * nEdgesSolve / block % nThreads)
            block % vertexThreadStart(threadid+1) = (threadid * nVertices / block % nThreads) + 1
            block % vertexThreadEnd(threadid+1)   = ((threadid+1) * nVertices / block % nThreads)
            block % vertexSolveThreadStart(threadid+1) = (threadid * nVerticesSolve / block % nThreads) + 1
            block % vertexSolveThreadEnd(threadid+1)   = ((threadid+1) * nVerticesSolve / block % nThreads)
!$OMP END PARALLEL
```
##### Look into how to get timings by process and thread

### Notes 9/16/15

##### Adding timing

Added per-thread timings to following files:
core_atmosphere/mpas_atm_core.F
core_atmosphere/dynamics/mpas_atm_time_integration.F
core_atmosphere/physics/mpas_atmphys_driver.F

Not needed for these files
core_atmosphere/physics/mpas_atmphys_driver_microphysics.F
framework/mpas_threading.F

