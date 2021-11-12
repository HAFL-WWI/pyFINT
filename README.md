# pyFINT: Python implementation of FINT

## Overview
pyFINT is a Python implementation of the FINT (Find INdividual Trees) tool from ecorisQ for detecting individual trees in a normalized crown height model. The program has been created as part of the swiss FINT-CH research project. 

## TL;DR
`docker build --rm -t pyfint .`
`docker run --rm -v $(pwd)/output:/pyfint/output -ti pyfint`

## Usage
The program has been written for Python 3.x. See requirements.txt for dependencies. The classes of the python version were kept close to the original C++ implementation. Their use is therefore not necessarily very pythonic. Details about the FINT algorithm and the program outputs can be found in the original tool's user manual: https://www.ecorisq.org/docs/FINT_manual_EN.pdf. Compared to the original version, pyFINT only supports normalized surface models in the GeoTIFF or Esri ASCII format as input and it doesn't support the calculation of a normalized surface model form a surface and terrain model. Additionally, pyFINT supports optional resizing and gauss filtering of the surface model, not yet available in the original C++ program. The example_main.py shows examples for the most common use cases of pyFINT.

Inputs:
- Normalized surface model as GeoTIFF or ASCII
- Optional: digital terrain model as GeoTIFF or ASCII
Outputs:
- *.TXT and *.CSV files with the detected individual trees
- schema.ini with an explanation of the columns in the *.TXT and *.CSV
- Optional: If a resizing or filtering option is configured, the modified normalized surface model is output as a GeoTIFF 

## Links
- FINT-CH Project Report (in German): https://www.ecorisq.org/relevant-literature/fint/58-final-report-fint-ch-project-in-german/file
- Download page ecorisQ: https://www.ecorisq.org/ecorisq-tools
