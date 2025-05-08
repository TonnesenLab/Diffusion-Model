# Diffusion Flux Model by Paula Giménez Mínguez

This app allows the user to model diffusion of particles in a field of view provided as a tiff image. 
Diffusion depends on the pixel intensity so that pixels with high intensity will allow particles to 
diffuse freely and pixels with low intensity representy structures will hinder diffusion until prevent it. 

## Objective

The model is meant to be used with SUPER RESOLUTION SHADOW IMAGES and it has been applied to 
study diffusion of transmitters in the brain extracellular space. Althought it can be applicable 
to study diffusion processes in different areas and at several spatio-temporal scales

## Respository Estructure

DiffusionFluxModel
- src/ # MATLAB source code
- data/ # Entry data
- results/ # generated results
- README.md # This file
- LICENSE.txt # BSD 3 clause license

## Rquisites

- MATLAB R2019b or newer
- Toolboxes:
  - Curve Fitting Toolbox 3.5.10
  - Image Processing Toolbox 11.0
  - Statistics and Machine Learning Toolbox 11.6

## Execution of the code

-execute DiffusionModelApp.mlapp 

## LICENSE
This project holds under the BSD 3-Clause license. For more information check 
the LICENSE.txt file

## Recomended Citing
Si usas este código, por favor cita el siguiente artículo:

scss
Copiar
Editar
Apellido, N. (202X). Título del artículo. *Nombre de la revista*, volumen(número), páginas. DOI.
? Enlace al código
Repositorio en GitHub: https://github.com/usuario/repositorio

DOI Zenodo (si aplica): https://doi.org/xxx
