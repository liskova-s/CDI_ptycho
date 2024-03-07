# Lensless imaging
This is a repository for a research project on Lensless imaging methods, drawned up (Nov 2023 - Sep 2024).  
Currently undergoing modifications and major fixes.  
#### _Main goals of the project:_ 
- _implement various numerical solutions to the diffraction problem_
- _research lensless imaging methods (CDI + ptychography)_
- _implement both imaging methods_
#### _Current state_
_The Fraunhofer, Fresnel and Rayleigh-Sommmerfeld freespace propagation models have been researched and the first implementations have been made. Currently they are being tested and tuned with respect to computing time vs. result accuracy._
#### _Outline of the following steps_
- _starting with building and testing freespace propagators_
- _using propagators to simulate a diffraction pattern of a given object_
- _implementing iterative CDI algoritmhs for retrieval of the object_
- _adding experimental parameters (such as noise and dynamic range of the camera) to the CDI simulation_
- _implementing ptychography based algoritmhs for object retrieval_
