### Notes 8/28/15

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
