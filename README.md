This repository contains scripts for launching and post-processing the simulations in the paper `The importance of geometry in the corneal micropocket assay.`


# Setting Up

Requirements to run simulations:

 * Microvessel Chaste Python package (and its dependencies)
 
 Requirements to do postprocessing:

 * Microvessel Chaste Python package (and its dependencies)
 * VTK Python package
 
 
# Running Simulations
 
 The simulation launcher is in `test/python/batch_simulation.py`. To run simulations do:
 
 
 ```python
 mpirun -np $NUM_CPUS python test/python/batch_simulation.py
 ```
 
 The collection of different simulations is in `src/python/cornea/simulations`. To change the simulation run by the launcher modify the line:
 
  ```python
from cornea.simulations.fixed_case5 import study, master_work_dir, study_data
 ```
 
 in `test/python/batch_simulation.py` to the desired simulation. Results will be in the directory specified by the `$CHASTE_TEST_OUTPUT` environment variable.
 
 
 # Postprocessing Simulations
 
 To postprocess simulations do:
 
  ```python
 mpirun -np $NUM_CPUS python test/python/batch_postprocess.py
 ```
 
 Change the line `work_dir = "Python/Cornea/Submission/Fixed_Case5"` to the desired output directory. Some postprocessing operations will work in parallel, but most only work in serial (requiring `mpirun -np 2`).
 
After postprocessing the plots used in the paper will be included in the simulation folders. Screenshots have been generated in Paraview using the `vtp`, `vti` and `vtu` files produced during the simulations.