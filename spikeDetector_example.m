IIS_det = SpikeWaveDetector; % see SpikeWaveDetector.docx for details on functions and default parameters

dataFolder = header.processed_MACRO; % extracted iEEG files
% IIS_det.samplingRate default is 1kHz - update if different

channels = [1, 9, 16, 24, 31, 38, 45, 53, 60, 67]; % deepest contact for each electrode

% Extract interictal activities from all channels
clear peakTimes mlink passedConditions
for ii = 1:length(channels)
    filename = fullfile(dataFolder, sprintf('CSC%d.mat',channels(ii)));
    mlink{ii} = matfile(filename);
    [peakTimes{ii}, passedConditions{ii}]= detectTimes(IIS_det, mlink{ii}.data,true);
end

% Plot IED detected on given channels
plotSpikeWaves(IIS_det, mlink{ii}.data, peakTimes{ii},1, passedConditions{ii});
