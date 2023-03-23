# epilepticActivity_IEEG (Drs. Maya Geva-Sagiv and Shdema Epstein)

This toolbox was designed for detecting pathological activity (interictal spikes) in iEEG data recorded from epileptic patients. 
Several methods for spikes detection are implemented in the class, the commonly used is - detectTimes. 
Its basic idea is to find points in the data which pass one or more thresholds of -
(i) amplitude following highpass; (ii) unfiltered amplitude (iii) signal gradient.  
See spikeDetector_example.m for an example to how to use the class for IED detection

Other methods which were implemented - 
(1) spike detection using wavelet analysis (based on West et al 2003). 
(2) spike detection using Taeger energy (based on Zaveri et al 2014).
these methods were not validated on data, so should be used with caution.

IISView is a validation gui that allows loading detections, mark and delete manually.

