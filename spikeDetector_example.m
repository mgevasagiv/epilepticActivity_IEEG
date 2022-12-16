IIS_det = SpikeWaveDetector;
load('E:\Data_p\p487\EXP3\p487_EXP3_dataset.mat');
load('E:\Data_p\MACRO_MONTAGE\p487\EXP3\MacroMontage.mat')
dataFolder = header.processed_MACRO;
channels = [1, 9, 16, 24, 31, 38, 45, 53, 60, 67]; % deepest contact for each electrode

clear peakTimes mlink passedConditions
for ii = 1:length(channels)
    filename = fullfile(dataFolder, sprintf('CSC%d.mat',channels(ii)));
    mlink{ii} = matfile(filename);
    [peakTimes{ii}, passedConditions{ii}]= detectTimes(IIS_det, mlink{ii}.data,true);
end

plotSpikeWaves(IIS_det, mlink{ii}.data, peakTimes{ii},1, passedConditions{ii});
