# rec_shape_defect
Code with the article Reconstruction of smooth shape defects in waveguides using locally resonant frequencies

Main programs: 
  - visualization.m
Programm to visualize the propagation of a wavefield in a perturbed waveguide at resonant frequencies. It is used to produce Figure 4-5-6. 

- rec_shape_increasing.m 
Program to recover increasing width defects in waveguides in 2D. This program generates the data using FEM/PML methods and use the reconstruction method presented in the article to provide an approximation of the width defect. It is used to produce Figure 15. 

- rec_shape_general_case.m 
Program to recover general width defects in waveguides in 2D. This program generates the data using FEM/PML methods and use the reconstruction method presented in the article to provide an approximation of the width defect. It is used to produce Figure 16. 
