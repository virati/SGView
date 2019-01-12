# Spectrogram Viewer
A MATLAB based spectrogram viewer that imports generic timeseries data. Beyond time-frequency displays and interactive interval-based PSD estimation, this tool does some basic gain compression assessment


## Installation
The folder is self-contained and includes related libraries.


### Dependencies
* Filtfilthd
* 

## Frontend
The frontend is focused on time-frequency domain displays of data. UI is available for more advanced functionality. The primary functionality of interest for contemporary (2013-2018) clinical electrophysiology is likely

![Frontend of SGView](imgs/SG_View_frontend.png)

## Disclaimer
The viewer is provided *as is* and without guarantees or liability.

## Notes
My original goal for this was overly-ambitious. I wanted it to include everything analysis-wise that existed in 2013 to address frequency-domain electrophysiology. It took a while but I realized that the field of electrophysiology was nowhere near the level that neuroimaging was and there really was no way to do a simple frontend for electrophysiology analysis. There are simply way too many knobs. You can see this attempt in my earlier commits.

