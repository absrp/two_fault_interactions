### Project description 
This repository contains a set of Jupyter Notebooks to analyze the outputs of boundary element code Motorcycle. These scripts can be used to reproduce the analysis in Rodriguez Padilla et al. (in review). 

### Motorcycle input files
The models analyzed by scripts in this repository are generated using Sylvain Barbot's spectral boundary element code Motorcycle. 
[Motorcycle read the docs - Barbot et al.](https://motorcycle.readthedocs.io/en/latest/)

The analysis in this project spans the following model parameter suite:

Model suite
| D/L | Ru | DX (m) | N2 |
|---|---|---|---|
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 6.7 | 20 | 4096 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 8.3 | 20 | 4096 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 10.6 | 20 | 4096 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 13.3 | 10 | 8192 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 17.1 | 10 | 8192 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 22.2 | 5 | 16384 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 33.3 | 5 | 16384 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 47.6 | 2.5 | 32768 |
| 0.04, 0.12, 0.2, 0.34, 0.6, 1, 2, 4, 6, 12, 24 | 66.7 | 2.5 | 32768 |

where D/L is the scaled distance between the faults (distance/length of rate weakening patch), Ru is the instability ratio, DX is the cell size, and N2 is the number of cells on the fault. 

The input files to run the Motorcycle models associated with this project are:

Single fault benchmarks:
- single_fault_ref_antiplane.sh
- single_fault_ref_inplane.sh
  
Two-fault models:
-two_faults_antiplane.sh
-two_faults_inplane.sh

### Model postprocessing and analysis scripts
1. utils.py
    Functions required to run the Jupyter Notebooks.
2. model_visualizer.ipynb
   Visualize model outputs and calculate alignment and correlation coefficients. 
3. single_fault_benchmark.ipynb
   Visualie and analyze single fault benchmarks.
4. single_fault_benchmark_times.ipynb
   Measure interevent times for a single fault benchmark. 
5. two_fault_times.ipynb
   Measure interevent times for a two fault model. 
6. map_observations.ipynb
    Make map of sites with observations of synchronization and alternation from   paleoseismic chronologies and incremental slip rates. Make illustration of sync and alternation behaviors.
7. comparison_inter_event_times.ipynb
   Compare interevent times for each two fault model and the corresponding single fault benchmark for the same value of Ru.
8. evolution_interaction_metrics_with_time.ipynb
   Plot evolution of the alignment and correlation coefficients for a model with its corresponding catalog.
9. measure_alignment_corr_coeff_data.ipynb
    Measure alignment and correlation coefficients for geologic datasets.
10. corr_metrics_stress_analysis.ipynb
    Plot complete model phase space, correlation and alignment coefficients for geologic datasets, and Coulomb stress analysis.
   
### Contact
Please direct any questions, suggestions or concerns to: 

Email - alba.rodriguez@usu.edu, alba@caltech.edu
Website - [absrp@github.io](absrp@github.io)

### Acknowledgements
We thank Sylvain Barbot for responding to questions about Motorcycle

### Manuscript link
Stay tuned!
