# epilepticActivity_IEEG

This toolbox was designed for detecting pathological activity (interictal spikes) in iEEG data recorded from epileptic patients. 
Several methods for spikes detection are implemented in the class, the commonly used one that - detectTimes. 
Its basic idea is to find points in the data which satisfy one of several conditions regarding values of -
-- amplitude after bandpass, 
-- unfiltered amplitude 
-- gradient  


Other methods which were implemented - 
-- spike detection using wavelet analysis (based on West et al 2003) 
-- using Taeger energy (based on Zaveri et al 2014), 
these methods were not validated on data, so should be used with caution.

IISView is a validation gui that allows loading detections, mark and delete manually.

