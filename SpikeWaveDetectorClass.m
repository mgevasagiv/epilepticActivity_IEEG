classdef SpikeWaveDetectorClass < handle
    
    properties
        
        %general constants
        samplingRate = 1000;
        minDistSpikes = 50; % minimal distance for 'different' spikes - in miliseconds
        
        %plotting constants
        plotBeforeAfter = 5000; %constant for the plotting method - how much to plot before and after the peak, in ms
        blockSizePlot = 1;
        
        % constants for detection based on frequency analysis - the
        % thresholds for the standard deviation are based on the papers
        % Andrillon et al (envelope condition) and Starestina et al (other
         %conditions)
        SDthresholdEnv = 5; %threshold in standard deviations for the envelope after bandpass (HP)
        SDthresholdAmp = 20; %threshold in standard deviations for the amplitude
        SDthresholdGrad = 20; %threshold in standard deviations for the gradient
        SDthresholdConjAmp = 5; %threshold in standard deviations for the amplitude for the conjunction of amp&grad condition
        SDthresholdConjGrad = 5; %threshold in standard deviations for the gradient for the conjunction of amp&grad condition
        SDthresholdConjEnv = 5; %threshold in standard deviations for the HP for the conjunction of amp&HP condition
        useEnv = true;
        useAmp = false;
        useGrad = false;
        useConjAmpGrad = true;
        useConjAmpEnv = false;
        isDisjunction = false;
        blockSizeSec = 30; % filter and find peaks at blocks of X seconds - based on Andrillon et al
        
        %         %the bandpass range is based on Andrillon et al
        %         lowCut = 50; %low bound for band pass
        %         highCut = 150; %high bound for band pass
        %         minLengthSpike = 5; %a spike is detected if there are points for X ms passing the threshold - in ms, based on Andrillon et al

        %the highpass range is based on Staresina et al 
        lowCut = 250;  % low bound for high pass (Staresina)
        highCut = inf; % high bound for band pass
        minLengthSpike = 5;
        maxLengthSpike = 70; % mSec % Andrillon
        
        conditionsArrayTrueIfAny = false;
        percentageOfNansAllowedArounsSpike = 0.1; % how many NaNs are allowed in the vicinity of the spike (vicinity = minDistSpikes/2 before and after)
        HighBandPassScore = 11; % this is for debugging - finding especially high STDs for HP
        
        
        %constants from original version of frequency analysis detection
        nPointsBlockSizeFiltfilt = 2*10^6; % This constant is used in the original version of the code
        
        %constants for bandpass
        defaultFilterOrder = 1;
        nanWarning = 0.01;
        
        %constants for detections based on wavelets - taken from the paper
        %West et al
        scaleForCWS = [1:30];
        thresholdForScale1 = 400;
        thresholdForScale2 = 400;
        scale1 = 3;
        scale2 = 7;
        scaleSlow = 28;
        c1 = 1;
        c2 = 1;
        tau = 0.125; %in seconds
        
        %constants for detection based on Taegar energy - taken from the
        %paper Zaveri et al
        minFreqT = 10;
        maxFreqT = 70;
        derBefore = 8; %ms, based on the paper
        derAfter = 12; %ms
        threshMah = 200;
        nanThresh = 0.5;
 
        % sleep scoring
        NREM = 1;
        
    end
    
    methods
        
        %% frequency analysis detection
        
        function [peakTimes, peakStats] = detectTimes(obj, data, returnPeakStats, sleepScoringVec)
            
            %Frequency analysis - based on Selective neuronal lapses precede human cognitive lapses
            %following sleep deprivation, Nir et al, 2017
            %peakStats - a struct with the fields:
            %passedConditions - an array of size N*D, where N is the number
            %of detected peaks and D is the number of conditions checked in
            %the code (a detection will occur if the disjunction of the
            %conditions is true). At index i,j - the array is true if the
            %peak at index i in peakTimes passed condition j. The order of
            %conditions: 1. Envelope of the signal after bandpass is above
            %a set threshold, 2. Amplitude of the signal is above a set
            %threshold, 3. Gradient (between consecutive time points) of
            %the signal is above a set threshold, 4. The amplitude & the
            %gradient are both above a set threshold which is smaller than
            %the thresholds used in conditions 2 & 3. 5. The amplitude and
            %the signal after bandpass are both above a threshold smaller
            %than 1 & 2
            %indsPerPeak - a cell the length of the number of peaks, each
            %element contains all the indices of that peak
            %zscoresPerPeaksEnv - a cell the length of the number of peaks, each
            %element contains all the HP zscores at the indices of that peak
            %zscoresPerPeaksAmp - a cell the length of the number of peaks, each
            %element contains all the amplitude zscores at the indices of that peak
            %zscoresPerPeaksGrad - a cell the length of the number of peaks, each
            %element contains all the gradient zscores at the indices of that peak
            %zscoresPerPeaksMax - an array N*3, N the number of peaks. For
            %each peak stores the maximal zscores of HP (index 1),
            %amplitude (index 2), and gradient (index 3) for that peak.
            
             if nargin < 4
                useSleepScoring = 0;
             else
                 useSleepScoring = 1;
             end
            
            if nargin < 3
                returnPeakStats = false;
            end
            
            peakTimes = [];
            if returnPeakStats
                passedConditions = [];
                zscoresPerPeaksMax = [];
                zscoresPerPeaksEnv = {};
                zscoresPerPeaksAmp = {};
                zscoresPerPeaksGrad = {};
                indsPerPeak = {};
            end
            %replace nans by zeros
            originalData = data;
            data(isnan(data)) = 0;
            
            % Detect IIS on NREM sleep
            if useSleepScoring
                sleepScoring = sleepScoringVec;
                data(sleepScoring ~= obj.NREM) = 0;
            end
            zsAmp_all = zscore(data); % zscore over the entire NREM session   
            
            pointsInBlock = obj.blockSizeSec*obj.samplingRate;
            nBlocks = floor(length(data)/pointsInBlock);
            ind = 1;
            for iBlock = 1:nBlocks
                %use 3 conditions: absolute amplitude above a threshold,
                %gradient above threshold, and envelope of the signal after
                %a bandpass above a threshold
                
                currBlock = data((iBlock-1)*pointsInBlock+1:iBlock*pointsInBlock);
                nCurrBlock = length(currBlock);
                
                % amplitude
                if obj.useAmp || obj.useConjAmpGrad || obj.useConjAmpEnv
                    zsAmp = zsAmp_all((iBlock-1)*pointsInBlock+1:iBlock*pointsInBlock); 
                                % Nov 20, changing from block-based z-scoring 
                                %         (that detected many spindles) to NREM-vec zscoring
                                % Dec 21, bug fix - checking abs values
                    pointsPassedThreshAmplitude = abs(zsAmp) > obj.SDthresholdAmp;
                    pointsPassedThreshAmplitudeLowThresh = abs(zsAmp) > obj.SDthresholdConjAmp;
                else
                    if obj.isDisjunction
                        pointsPassedThreshAmplitude = false(1,nCurrBlock);
                    else
                        pointsPassedThreshAmplitude = true(1,nCurrBlock);
                    end
                end
                
                % gradient
                if obj.useGrad || obj.useConjAmpGrad
                    dataGradient = [0 diff(currBlock)];
                    zsGrad = zscore(dataGradient);
                    pointsPassedThreshGradient = zsGrad > obj.SDthresholdGrad;
                    pointsPassedThreshGradientLowThresh = zsGrad > obj.SDthresholdConjGrad;
                else
                    if obj.isDisjunction
                        pointsPassedThreshGradient = false(1,nCurrBlock);
                    else
                        pointsPassedThreshGradient = true(1,nCurrBlock);
                    end
                end
                
                
                %bandpass and envelope
                if obj.useEnv || obj.useConjAmpEnv
                    %first perform bandpass filtering
                    filteredBlock = obj.bandpass(currBlock, obj.samplingRate, obj.lowCut, obj.highCut);
                    %find envelope
                    envBlock = abs(hilbert(filteredBlock));
                    %find points which pass the threshold as set by number of
                    %SDs as compared to the current block
                    zsEnv = zscore(envBlock);
                    pointsPassedThreshEnv = zsEnv > obj.SDthresholdEnv;
                    pointsPassedThreshEnvLowThresh = zsEnv > obj.SDthresholdConjEnv;
                else
                    if obj.isDisjunction
                        pointsPassedThreshEnv = false(1,nCurrBlock);
                    else
                        pointsPassedThreshEnv = true(1,nCurrBlock);
                    end
                end
                
                % conjunction of amplitude & gradient with lower thresholds
                if obj.useConjAmpGrad
                    pointsPassedThreshAmpGradLowThresh = pointsPassedThreshGradientLowThresh & pointsPassedThreshAmplitudeLowThresh;
                else
                    pointsPassedThreshAmpGradLowThresh = false(1,nCurrBlock);
                end
                
                % conjunction of amplitude & HP with lower thresholds
                if obj.useConjAmpEnv
                    pointsPassedThreshAmpEnvLowThresh = pointsPassedThreshEnvLowThresh & pointsPassedThreshAmplitudeLowThresh;
                else
                    pointsPassedThreshAmpEnvLowThresh = false(1,nCurrBlock);
                end
                
                if obj.isDisjunction
                    %if isDisjuction is true - points are detected as threshold if any of the conditions
                    %is met
                    pointsPassedThresh = pointsPassedThreshEnv | pointsPassedThreshGradient | pointsPassedThreshAmplitude | ...
                        pointsPassedThreshAmpGradLowThresh | pointsPassedThreshAmpEnvLowThresh;
                else
                    %if isDisjuction is false - points are detected as threshold if all of the conditions
                    %are met
                    pointsPassedThresh = pointsPassedThreshEnv & pointsPassedThreshGradient & pointsPassedThreshAmplitude;
                end
                               
                %findSequence is called in order to check whether there is
                %a sequence of points lasting X ms to pass the threshold
                %(in paper - 5 ms, due to adapatations of code - 1 ms, i.e. no minimal length), and separates between
                %different spikes (merges points which are close
                %together to one spike)
                if sum(pointsPassedThresh)>0
                    if returnPeakStats
                        [currPeaks,allPeakInds] = obj.findSequences(currBlock, pointsPassedThresh);
                        nCurrPeaks = length(currPeaks);
                        currPassedConditions = [];
                        currZscoresPerPeaksMax = [];
                        currIndsPerPeak = cell(1,nCurrPeaks);
                        currZscoresPerPeaksEnv = cell(1,nCurrPeaks);
                        currZscoresPerPeaksAmp = cell(1,nCurrPeaks);
                        currZscoresPerPeaksGrad = cell(1,nCurrPeaks);
                        for iPeak = 1:nCurrPeaks

                            if obj.conditionsArrayTrueIfAny
                                currPassedConditions = [currPassedConditions; any(pointsPassedThreshEnv(allPeakInds{iPeak})) any(pointsPassedThreshAmplitude(allPeakInds{iPeak})) ...
                                    any(pointsPassedThreshGradient(allPeakInds{iPeak})) any(pointsPassedThreshAmpGradLowThresh(allPeakInds{iPeak})) ...
                                    any(pointsPassedThreshAmpEnvLowThresh(allPeakInds{iPeak}))];
                            else
                                seqToFind = ones(1,obj.minLengthSpike);
                                currPassedConditions = [currPassedConditions; ~isempty(strfind(pointsPassedThreshEnv(allPeakInds{iPeak}),seqToFind)) ...
                                    ~isempty(strfind(pointsPassedThreshAmplitude(allPeakInds{iPeak}),seqToFind)) ~isempty(strfind(pointsPassedThreshGradient(allPeakInds{iPeak}),seqToFind)) ...
                                    ~isempty(strfind(pointsPassedThreshAmpGradLowThresh(allPeakInds{iPeak}),seqToFind)) ~isempty(strfind(pointsPassedThreshAmpEnvLowThresh(allPeakInds{iPeak}),seqToFind))];
                            end
                            currZscoresPerPeaksMax = [currZscoresPerPeaksMax; max(zsEnv(allPeakInds{iPeak})) max(zsGrad(allPeakInds{iPeak})) max(zsAmp(allPeakInds{iPeak}))];
                            currIndsPerPeak{iPeak} = allPeakInds{iPeak};
                            currZscoresPerPeaksEnv{iPeak} = zsEnv(allPeakInds{iPeak});
                            currZscoresPerPeaksAmp{iPeak} = zsAmp(allPeakInds{iPeak});
                            currZscoresPerPeaksGrad{iPeak} = zsGrad(allPeakInds{iPeak});
                        end
                    else
                        currPeaks = obj.findSequences(currBlock, pointsPassedThresh);
                    end
                    
                    %remove spikes with too many NaNs around them
                    isnanAtPeak = false(1,length(currPeaks));
                    spikeVicinity = round((obj.minDistSpikes/1000)*obj.samplingRate/2); % convert from ms
                    for iPeak = 1:length(currPeaks)
                        origPeakInd = currPeaks(iPeak)+(iBlock-1)*pointsInBlock;
                        %handle indices at beginning or end of the original data
                        if origPeakInd <= spikeVicinity
                            startPoint = 1;
                        else
                            startPoint = origPeakInd-spikeVicinity;
                        end
                        if origPeakInd+spikeVicinity > length(originalData)
                            endPoint = length(originalData);
                        else
                            endPoint = origPeakInd+spikeVicinity;
                        end
                        dataAroundSpike = originalData(startPoint:endPoint);
                        isnanAtPeak(iPeak) = isnan(originalData(origPeakInd)) | ...
                            sum(isnan(dataAroundSpike))/length(dataAroundSpike) >= obj.percentageOfNansAllowedArounsSpike;
                        %sometimes the data has zeros instead of NaN - next
                        %code is to deal with it
                        if ~isnanAtPeak(iPeak)
                            isnanAtPeak(iPeak) = sum(dataAroundSpike==0)/length(dataAroundSpike) >= obj.percentageOfNansAllowedArounsSpike;
                        end
                    end
                    currPeaks = currPeaks(~isnanAtPeak);
                    
                    %remove the nan peaks also from the statistics
                    if returnPeakStats
                        passedConditions = [passedConditions; currPassedConditions(~isnanAtPeak,:)];
                        zscoresPerPeaksMax = [zscoresPerPeaksMax; currZscoresPerPeaksMax(~isnanAtPeak,:)];
                        
                        indsPerPeak = [indsPerPeak currIndsPerPeak(~isnanAtPeak)];
                        zscoresPerPeaksEnv = [zscoresPerPeaksEnv currZscoresPerPeaksEnv(~isnanAtPeak)];
                        zscoresPerPeaksAmp = [zscoresPerPeaksAmp currZscoresPerPeaksAmp(~isnanAtPeak)];
                        zscoresPerPeaksGrad = [zscoresPerPeaksGrad currZscoresPerPeaksGrad(~isnanAtPeak)];
                        
                    end
                else
                    currPeaks = [];
                end
                
                peakTimes = [peakTimes currPeaks+(iBlock-1)*pointsInBlock];   
            end
            
            if returnPeakStats
                peakStats.zscoresPerPeaksMax = zscoresPerPeaksMax;
                peakStats.zscoresPerPeaksEnv = zscoresPerPeaksEnv;
                peakStats.zscoresPerPeaksAmp = zscoresPerPeaksAmp;
                peakStats.zscoresPerPeaksGrad = zscoresPerPeaksGrad;
                peakStats.indsPerPeak = indsPerPeak;
                peakStats.passedConditions = passedConditions;
            end
            
        end
        
        function scores = getScoresTimesPointsFreq(obj, data, timePoints)
            %The method receives the data and a set of time points and
            %checks for each time point whether the detection using
            %frequency analysis recognizes a peak around that point. scores
            %is an array the size of timePoints where 0 means no peak was
            %recognized and 1 means a peak was recognized.
            
            nTP = length(timePoints);
            scores = zeros(1,nTP);
            pointsInBlock = obj.blockSizeSec*obj.samplingRate;
            %for recognition - use a block in which the point is at its
            %center
            pointsBeforeAfter = round(pointsInBlock/2);
            
            for iTP = 1:nTP
                %next code is for the edges - in case it's not possible to
                %build a block in which the point is exactly at the center
                %because it's at the start/end of the data
                if timePoints(iTP)>pointsBeforeAfter
                    minPoint = timePoints(iTP)-pointsBeforeAfter;
                    peakPoint = pointsBeforeAfter;
                else
                    minPoint = 1;
                    peakPoint = timePoints(iTP);
                end
                
                if timePoints(iTP)+pointsBeforeAfter>length(data)
                    maxPoint = length(data);
                else
                    maxPoint = timePoints(iTP)+pointsBeforeAfter;
                end
                currBlock = data(minPoint:maxPoint);
                
                %the rest of the code is the frequency detection algorithm
                %as already implemented - band pass, envelope, thresholds
                %by SD distance from mean of block, and a sequence finding
                filteredBlock = obj.bandpass(currBlock, obj.samplingRate, obj.lowCut, obj.highCut);
                %find envelope
                envBlock = abs(hilbert(filteredBlock));
                pointsPassedThresh = zscore(envBlock(~isnan(envBlock))) > obj.SDthreshold;
                if sum(pointsPassedThresh)>0
                    currPeaks = obj.findSequences(currBlock, pointsPassedThresh);
                else
                    currPeaks = [];
                    scores(iTP) = 0;
                end
                
                %obj.minDistSpikes is in ms - translate to number of data
                %points
                distSpikePoints = round(obj.minDistSpikes*obj.samplingRate/1000);
                if ~isempty(currPeaks)
                    %any peak distanced minDistSpikes from the time point
                    %counts as a recognition of spike
                    if any(abs(currPeaks-peakPoint) <= distSpikePoints)
                        scores(iTP) = 1;
                    end
                end
            end
        end
        
        function peakStats = getStats(obj, data, peakTimes)
            
            %receives precalculated peakTimes and data and returns statistics. 
            %The indices of each peak are the given peak time +-
            %obj.minDistSpike ms.
            %peakStats - a struct with the fields:
            %passedConditions - an array of size N*D, where N is the number
            %of detected peaks and D is the number of conditions checked in
            %the code (a detection will occur if the disjunction of the
            %conditions is true). At index i,j - the array is true if the
            %peak at index i in peakTimes passed condition j. The order of
            %conditions: 1. Envelope of the signal after bandpass is above
            %a set threshold, 2. Amplitude of the signal is above a set
            %threshold, 3. Gradient (between consecutive time points) of
            %the signal is above a set threshold, 4. The amplitude & the
            %gradient are both above a set threshold which is smaller than
            %the thresholds used in conditions 2 & 3. 5. The amplitude and
            %the signal after bandpass are both above a threshold smaller
            %than 1 & 2
            %indsPerPeak - a cell the length of the number of peaks, each
            %element contains all the indices of that peak
            %zscoresPerPeaksEnv - a cell the length of the number of peaks, each
            %element contains all the HP zscores at the indices of that peak
            %zscoresPerPeaksAmp - a cell the length of the number of peaks, each
            %element contains all the amplitude zscores at the indices of that peak
            %zscoresPerPeaksGrad - a cell the length of the number of peaks, each
            %element contains all the gradient zscores at the indices of that peak
            %zscoresPerPeaksMax - an array N*3, N the number of peaks. For
            %each peak stores the maximal zscores of HP (index 1),
            %amplitude (index 2), and gradient (index 3) for that peak.
            
            %peakTimes = [];
            passedConditions = [];
            zscoresPerPeaksMax = [];
            zscoresPerPeaksEnv = {};
            zscoresPerPeaksAmp = {};
            zscoresPerPeaksGrad = {};
            indsPerPeak = {};
            
            %translate peakTimes to indices
            peakTimes = round(peakTimes*obj.samplingRate/1000);
            
            %replace nans by zeros
            originalData = data;
            data(isnan(data)) = 0;
            
            pointsInBlock = obj.blockSizeSec*obj.samplingRate;
            nBlocks = floor(length(data)/pointsInBlock);
            ind = 1;
            distSpikePoints = round(obj.minDistSpikes*obj.samplingRate/1000);
            
            for iBlock = 1:nBlocks
                %use 3 conditions: absolute amplitude above a threshold,
                %gradient above threshold, and envelope of the signal after
                %a bandpass above a threshold
                               
                currStartInd = (iBlock-1)*pointsInBlock+1;
                currEndInd = iBlock*pointsInBlock;
                currBlock = data(currStartInd:currEndInd);
                nCurrBlock = length(currBlock);
                currPeaks = peakTimes(peakTimes>=currStartInd & peakTimes<=currEndInd);
                currPeaks = currPeaks-currStartInd+1;
                nCurrPeaks = length(currPeaks);
                
                % amplitude
                zsAmp = zscore(currBlock);
                pointsPassedThreshAmplitude = zsAmp > obj.SDthresholdAmp;
                pointsPassedThreshAmplitudeLowThresh = zsAmp > obj.SDthresholdConjAmp;
                
                % gradient
                dataGradient = [0 diff(currBlock)];
                zsGrad = zscore(dataGradient);
                pointsPassedThreshGradient = zsGrad > obj.SDthresholdGrad;
                pointsPassedThreshGradientLowThresh = zsGrad > obj.SDthresholdConjGrad;
                
                
                %bandpass and envelope
                
                %first perform bandpass filtering
                filteredBlock = obj.bandpass(currBlock, obj.samplingRate, obj.lowCut, obj.highCut);
                %find envelope
                envBlock = abs(hilbert(filteredBlock));
                %find points which pass the threshold as set by number of
                %SDs as compared to the current block
                zsEnv = zscore(envBlock);
                pointsPassedThreshEnv = zsEnv > obj.SDthresholdEnv;
                pointsPassedThreshEnvLowThresh = zsEnv > obj.SDthresholdConjEnv;
                
                % conjunction of amplitude & gradient with lower thresholds
                
                pointsPassedThreshAmpGradLowThresh = pointsPassedThreshGradientLowThresh & pointsPassedThreshAmplitudeLowThresh;
                pointsPassedThreshAmpEnvLowThresh = pointsPassedThreshEnvLowThresh & pointsPassedThreshAmplitudeLowThresh;
                
                currPassedConditions = [];
                currZscoresPerPeaksMax = [];
                currZscoresPerPeaksEnv = cell(1,nCurrPeaks);
                currZscoresPerPeaksAmp = cell(1,nCurrPeaks);
                currZscoresPerPeaksGrad = cell(1,nCurrPeaks);
                currIndsPerPeak = cell(1,nCurrPeaks);
                               
                for iPeak = 1:nCurrPeaks    
                    allPeakInds = max(currPeaks(iPeak)-distSpikePoints,1):min(currPeaks(iPeak)+distSpikePoints,nCurrBlock);
                    if obj.conditionsArrayTrueIfAny
                        currPassedConditions = [currPassedConditions; any(pointsPassedThreshEnv(allPeakInds)) any(pointsPassedThreshAmplitude(allPeakInds)) ...
                            any(pointsPassedThreshGradient(allPeakInds)) any(pointsPassedThreshAmpGradLowThresh(allPeakInds)) ...
                            any(pointsPassedThreshAmpEnvLowThresh(allPeakInds))];
                    else
                        seqToFind = ones(1,obj.minLengthSpike);
                        currPassedConditions = [currPassedConditions; ~isempty(strfind(pointsPassedThreshEnv(allPeakInds),seqToFind)) ...
                            ~isempty(strfind(pointsPassedThreshAmplitude(allPeakInds),seqToFind)) ~isempty(strfind(pointsPassedThreshGradient(allPeakInds),seqToFind)) ...
                            ~isempty(strfind(pointsPassedThreshAmpGradLowThresh(allPeakInds),seqToFind)) ~isempty(strfind(pointsPassedThreshAmpEnvLowThresh(allPeakInds),seqToFind))];
                    end
                    currZscoresPerPeaksMax = [currZscoresPerPeaksMax; max(zsEnv(allPeakInds)) max(zsGrad(allPeakInds)) max(zsAmp(allPeakInds))];
                    currIndsPerPeak{iPeak} = allPeakInds;
                    currZscoresPerPeaksEnv{iPeak} = zsEnv(allPeakInds);
                    currZscoresPerPeaksAmp{iPeak} = zsAmp(allPeakInds);
                    currZscoresPerPeaksGrad{iPeak} = zsGrad(allPeakInds);
                end
                 
                passedConditions = [passedConditions; currPassedConditions];
                zscoresPerPeaksMax = [zscoresPerPeaksMax; currZscoresPerPeaksMax];
                
                indsPerPeak = [indsPerPeak currIndsPerPeak];
                zscoresPerPeaksEnv = [zscoresPerPeaksEnv currZscoresPerPeaksEnv];
                zscoresPerPeaksAmp = [zscoresPerPeaksAmp currZscoresPerPeaksAmp];
                zscoresPerPeaksGrad = [zscoresPerPeaksGrad currZscoresPerPeaksGrad];
                
            end
                
            peakStats.zscoresPerPeaksMax = zscoresPerPeaksMax;
            peakStats.zscoresPerPeaksEnv = zscoresPerPeaksEnv;
            peakStats.zscoresPerPeaksAmp = zscoresPerPeaksAmp;
            peakStats.zscoresPerPeaksGrad = zscoresPerPeaksGrad;
            peakStats.indsPerPeak = indsPerPeak;
            peakStats.passedConditions = passedConditions;
        end
        
        %code from Maya
        function [peakTimes] = detectSAWsegments(obj, EEG_filename, header ,channel)
            
            %This is the original code from Maya which I rewrote according
            %to the paper
            
            % For example
            % [peakTimes] = detect_SAW_segments('REC1M_EEG.mat');
            %
            % PeakTimes is a two-column matrix where each row is a detected spike-wave complex;
            % the first column is the segment number, the second column is the position (ms) within that segment
            % Based on YN code, modified MGS Feb 2017
            
            %             if (nargin < 2)
            %                 SDthreshold = 8;
            %             end
            
            data = load(EEG_filename,'data');
            data = data.data;
            peakTimes = [];
            thresholds = [];
            hpSignal = [];
            
            if obj.nPointsBlockSizeFiltfilt > length(data)
                obj.nPointsBlockSizeFiltfilt = length(data);
            end
            for ii_block = 1 : floor( length(data) / obj.nPointsBlockSizeFiltfilt ),
                
                idx_block = (1:obj.nPointsBlockSizeFiltfilt)+(ii_block-1)*obj.nPointsBlockSizeFiltfilt ;
                filteredData = obj.bandpass(data(idx_block), obj.samplingRate, obj.lowCut, obj.highCut);
                hpSignal = [hpSignal(:)', filteredData(:)'];
                
            end % The "problematic discontinuities" in data_filt_ripples - on the borders of 2 blocks - are SMALL in size
            
            if length( data( idx_block(end)+1 : end ) ) > 1
                filteredData = obj.bandpass(data(idx_block(end)+1:end), obj.samplingRate, obj.lowCut, obj.highCut);
                hpSignal = [hpSignal(:)', filteredData(:)'];
            end
            
            if length(hpSignal) ~= length(data)
                error('wrong hp calc')
            end
            
            peakTimesWithinSegment = find(zscore(hpSignal(~isnan(hpSignal))) > obj.SDthreshold);
            
            if (~isempty(peakTimesWithinSegment) > 0)
                
                % Here we need to reduce number since it detects a lot of neighboring points and this makes subsequent computations slow+
                diffBetweenDetectedPoints = diff(peakTimesWithinSegment);
                RowsWithBigEnoughDifferences = find(diffBetweenDetectedPoints > obj.minDistSpikes);
                finalRowsToTake = RowsWithBigEnoughDifferences+1;
                finalPeakTimesWithinSegment = peakTimesWithinSegment(finalRowsToTake);
                
                peakTimes = finalPeakTimesWithinSegment;
            end
            
            %             fileIdentifier = sprintf('%sE%dC%03d',header.id,header.experimentNum,channel);
            %             saveName = sprintf('SpikeAndWaveTimesFor_%s_threshold=%d.mat',fileIdentifier,obj.SDthreshold);
            %
            %             %%
            %             fprintf(1, 'Detected %d spike and waves for %s\n', length(peakTimes), fileIdentifier);
            %             save(saveName, 'peakTimes', 'thresholds');
            
        end
        
        
        %% wavelet analsys detection - this code was not validated or debugged, so if utilized should be used with caution
        
        function peakTimes = detectTimesWavelet(obj,data)
            
            %based on the paper Wavelet analysis of epileptic spikes -
            %Latka, Was, Kozik, West (2003). The heart of the algorithm is
            %implemented in the method detectTimeWavelets1Block
            
            peakTimes = [];
            
            %calculate block size - we find the peaks seperately for each
            %block
            pointsInBlock = obj.blockSizeSec*obj.samplingRate;
            nBlocks = floor(length(data)/pointsInBlock);
            
            for iBlock = 1:nBlocks
                currBlock = data((iBlock-1)*pointsInBlock+1:iBlock*pointsInBlock);
                currPeaks = obj.detectTimeWavelets1Block(currBlock);
                peakTimes = [peakTimes currPeaks+(iBlock-1)*pointsInBlock];
            end
        end
        
        function peakTimes = detectTimeWavelets1Block(obj, dataBlock, toPlot)
            %The method finds peaks in one block using the wavelet method according to the Latka,
            %West paper
            
            if nargin < 3
                toPlot = false;
            end
            
            %performs wavelet analysis on the block with the mexican hat
            %wavelet
            wt = cwt(dataBlock,obj.scaleForCWS,'mexh');
            
            delaySize = round(obj.tau*obj.samplingRate);
            lengthBlock = length(dataBlock)-delaySize+1;
            %according to the paper - two functions are calculated for two
            %scales (by default scale 3 and scale 7). The low scales (3,7)
            %represent the fast spike element, the high scale (by default
            %30) represent the slow wave which follows. There is a time
            %delay between them - be default 0.125 seconds. The function
            %calculated will receive a high value if there is a fast
            %element followed by a slow element. C1 and C2 should be
            %manipulated according the relative contribution of the
            %elements given the expected shape of the spike-wave
            
            funcAtScale1 = obj.c1*wt(obj.scale1,1:lengthBlock)+obj.c1*wt(obj.scaleSlow,delaySize:end);
            funcAtScale2 = obj.c2*wt(obj.scale2,1:lengthBlock)+obj.c2*wt(obj.scaleSlow,delaySize:end);
            
            funcAtScale1 = funcAtScale1.^2./var(funcAtScale1);
            funcAtScale2 = funcAtScale2.^2./var(funcAtScale2);
            
            %peaks are recognized where both functions are above threshold
            allPeaks = find(funcAtScale1>obj.thresholdForScale1 & funcAtScale2>obj.thresholdForScale2);
            
            %obj.minDistSpikes is in ms - translate to number of data
            %points
            distSpikePoints = round(obj.minDistSpikes*obj.samplingRate/1000);
            %next sections is for merging peaks which are close together to
            %one spike
            if ~isempty(allPeaks)
                seqDists = diff(allPeaks);
                %two peaks which are in a distance of minDistSpikes + tau
                %are considered the same spike
                seqBeginnings = [allPeaks(1) allPeaks(find(seqDists>distSpikePoints+obj.tau*1000)+1)];
                nSeqs = length(seqBeginnings);
                peakTimes = zeros(1,nSeqs);
                for iPeak = 1:nSeqs
                    %set the time point of the peak to be the largest
                    %absolute value (i.e. positive or negative large peak) within the range of minDistSpike+tau
                    %(the additional tau element is because it defines the
                    %shift of the slow wave from the fast spike, but both
                    %are considered the same event). Of course if only
                    %positive or negative peaks are preferred as reference
                    %point it can be easily changed
                    if seqBeginnings(iPeak)+obj.tau*1000+obj.minDistSpikes < length(dataBlock)
                        [~,maxPoint] = max(abs(dataBlock(seqBeginnings(iPeak):seqBeginnings(iPeak)+obj.tau*1000+obj.minDistSpikes)));
                    else
                        [~,maxPoint] = max(abs(dataBlock(seqBeginnings(iPeak):end)));
                    end
                    peakTimes(iPeak) = maxPoint+seqBeginnings(iPeak);
                end
            else
                peakTimes = [];
            end
            
            %plot of the different elements accordign to which peaks will
            %be recognized
            if toPlot
                nPoints = length(dataBlock);
                nPlots = 3;
                subplot(nPlots,1,1);
                plot(dataBlock);
                title('data');
                xlim([1 nPoints]);
                subplot(nPlots,1,2);
                plot(funcAtScale1);
                title('scale 3 sum');
                xlim([1 nPoints]);
                subplot(nPlots,1,3);
                plot(funcAtScale2);
                title('scale 7 sum');
                xlim([1 nPoints]);
            end
        end
        
        function scores = getScoresTimesPointsWavelet(obj, data, timePoints)
            %The method receives the data and a set of time points and
            %checks for each time point what is the maximal value of the functions calculated
            %in the wavelet recognition process around the time point. i.e.,
            %given a threshold for wavelet detection T, for each time
            %point at index i, if scores(i)>=T the points will be detected
            %as a peak and vice versa
            %This code is actually not general because there are two threshold  -
            %for each scale, and they can be different
            
            nTP = length(timePoints);
            scores = zeros(1,nTP);
            pointsInBlock = obj.blockSizeSec*obj.samplingRate;
            
            %for recognition - use a block in which the point is at its
            %center
            pointsBeforeAfter = round(pointsInBlock/2);
            
            for iTP = 1:nTP
                %next code is for the edges - in case it's not possible to
                %build a block in which the point is exactly at the center
                %because it's at the start/end of the data
                if timePoints(iTP)>pointsBeforeAfter
                    minPoint = timePoints(iTP)-pointsBeforeAfter;
                    peakPoint = pointsBeforeAfter;
                else
                    minPoint = 1;
                    peakPoint = timePoints(iTP);
                end
                
                if timePoints(iTP)+pointsBeforeAfter>length(data)
                    maxPoint = length(data);
                else
                    maxPoint = timePoints(iTP)+pointsBeforeAfter;
                end
                currBlock = data(minPoint:maxPoint);
                scores(iTP) = obj.getScoreWavelets1Block(currBlock,peakPoint);
            end
        end
        
        function scoreBlock = getScoreWavelets1Block(obj, dataBlock, eventInd)
            
            %the next code is same as wavelet detection (see documentation
            %inside detectTimeWavelets1Block
            wt = cwt(dataBlock,obj.scaleForCWS,'mexh');
            wt = wt.^2;
            wt = wt/var(dataBlock);
            delaySize = round(obj.tau*obj.samplingRate);
            lengthBlock = length(dataBlock)-delaySize+1;
            funcAtScale1 = obj.c1*wt(obj.scale1,1:lengthBlock)+obj.c1*wt(obj.scaleSlow,delaySize:end);
            funcAtScale2 = obj.c2*wt(obj.scale2,1:lengthBlock)+obj.c2*wt(obj.scaleSlow,delaySize:end);
            
            %looks for the maximal scores around the area of the time point
            eventInds = [eventInd - obj.minDistSpikes - obj.tau*1000:eventInd + obj.minDistSpikes];
            
            %maximal score of the minimimum of the two functions (as the
            %confition to detect a peak is that both functions will be
            %above the threshold).
            scoreBlock = max(min(funcAtScale1(eventInds),funcAtScale2(eventInds)));
            
        end
        
        
        %% Taeger energy detection - this code was not validated or debugged, so if utilized should be used with caution
        
        function peakTimes = detectTimesTaeger(obj,data, toPlot)
            
            %based on the paper - Automatic detection of prominent
            %interictal spikes in intracranial EEG,
            % Gaspard, Alawadri, Zaveri, 2014
            % The heart of the algorithm is implemented in the method detectTimeTeager1Block
            
            if nargin < 3
                toPlot = false;
            end
            
            peakTimes = [];
            
            pointsInBlock = obj.blockSizeSec*obj.samplingRate;
            nBlocks = floor(length(data)/pointsInBlock);
            
            for iBlock = 1:nBlocks
                currBlock = data((iBlock-1)*pointsInBlock+1:iBlock*pointsInBlock);
                %do not find peak on blocks with too many NaNs
                if sum(isnan(currBlock))/length(currBlock)>obj.nanThresh
                    continue;
                end
                currPeaks = obj.detectTimeTeager1Block(currBlock, toPlot);
                peakTimes = [peakTimes currPeaks+(iBlock-1)*pointsInBlock];
                if toPlot
                    pause;
                end
            end
        end
        
        function peakTimes = detectTimeTeager1Block(obj, dataBlock, toPlot)
            
            %Detects peaks for one block according to the paper by Gaspard,
            %Zaveri
            
            if nargin < 3
                toPlot = false;
            end
            
            %bandpass
            filteredBlock = obj.bandpass(dataBlock, obj.samplingRate, obj.minFreqT, obj.maxFreqT);
            nPoints = length(filteredBlock);
            
            %calculate Taeger energy
            energSig=energyop(filteredBlock);
            %calculate the first derivative of the block
            funDer = ppval(fnder(spline([1:nPoints],filteredBlock)),[1:nPoints]);
            %build the triples of the derivative before and after the peak
            %and the energy at the peak
            dataTriplets = [funDer(1:nPoints-obj.derBefore-obj.derAfter)' energSig(obj.derBefore+1:nPoints-obj.derAfter) funDer(obj.derBefore+obj.derAfter+1:nPoints)'];
            %calculate the Mahalanobis of each point from the general
            %distribution of the block - in order to detect the outliers
            %which are recognized as peaks
            mahDist = mahal(dataTriplets,dataTriplets);
            
            %peaks are recognized as points which pass a threshold for the
            %distance
            allPeaks = find(mahDist>obj.threshMah)';
            
            %obj.minDistSpikes is in ms - translate to number of data
            %points
            distSpikePoints = round(obj.minDistSpikes*obj.samplingRate/1000);
            
            %the next code merges together peaks which are close together
            %as one spike
            if ~isempty(allPeaks)
                seqDists = diff(allPeaks);
                seqBeginnings = [allPeaks(1) allPeaks(find(seqDists>distSpikePoints)+1)];
                nSeqs = length(seqBeginnings);
                peakTimes = zeros(1,nSeqs);
                for iPeak = 1:nSeqs
                    lowerLimit = seqBeginnings(iPeak);
                    if length(seqBeginnings)>iPeak
                        upperLimit = seqBeginnings(iPeak+1)-1;
                    else
                        upperLimit = length(dataBlock);
                    end
                    currInds = allPeaks(allPeaks >= lowerLimit & allPeaks <= upperLimit);
                    %the point which defines a peaks is chosen to be the
                    %point with the maximal data value among the points
                    %which passed the threshold
                    [~,maxPoint] = max(dataBlock(currInds));
                    peakTimes(iPeak) = currInds(maxPoint);
                end
            else
                peakTimes = [];
            end
            
            if toPlot
                hold off;
                nPoints = length(dataBlock);
                nPlots = 2;
                subplot(nPlots,1,1);
                plot(dataBlock);
                title('data');
                xlim([1 nPoints]);
                hold off;
                subplot(nPlots,1,2);
                plot(mahDist);
                title('mah distance');
                xlim([1 nPoints]);
                hold off;
            end
            
        end
        
        %% help functions
        
        function BP = bandpass(obj, timecourse, SamplingRate, low_cut, high_cut, filterOrder)
            
            %bandpass code - from Maya
            if (nargin < 6)
                filterOrder = obj.defaultFilterOrder;
            end
            
            % Maya GS - handle NAN values
            indices = find(isnan(timecourse));
            if length(indices) > obj.nanWarning*length(timecourse)
                warning('many NaN values in filtered signal')
            end
            timecourse(indices) = 0;
            %
            if high_cut == inf
                [b, a] = butter(filterOrder,(low_cut/SamplingRate)*2,'high');
            else
                [b, a] = butter(filterOrder, [(low_cut/SamplingRate)*2 (high_cut/SamplingRate)*2]);
            end
            BP = filtfilt(b, a, timecourse );
            BP(indices) = NaN;
        end
        
        function plotSpikeWaves(obj, data, peakTimes, blockSizeToPlot, peakStats, plotZScores)
            
            %plots the peak times
            %receives the data, the timing of the spikes, and how many
            %blocks (spikes) to plot in each subplot (recommended is 5,
            %default is 1)
            %peakStats - the output of detectTimes with statistics on
            %peaks. If provided, for each peak the title will include
            %information on which detection condition it passed and the
            %maximal zscore for each parameter (envelope of HP, amplitude,
            %gradient)
            
            if nargin < 4 || isempty(blockSizeToPlot)
                blockSizeToPlot = obj.blockSizePlot;
            end
            
            if nargin < 5
                plotConditionsData = false;
            else
                plotConditionsData = true;
            end
            
            if nargin < 6 || ~plotConditionsData
                plotZScores = false;
            end
                 
            nPeaks = length(peakTimes);
            indBlock = 1;
            
            blockSizeData = obj.blockSizeSec*obj.samplingRate;
            
            for iPeak = 1:nPeaks
                if plotZScores
                    subplot(4,blockSizeToPlot,indBlock);
                else
                    subplot(1,blockSizeToPlot,indBlock);
                end
                
                %plot the block in which the peak was detected
                blockNum = floor(peakTimes(iPeak)/blockSizeData)+1;
                minPoint = (blockNum-1)*blockSizeData+1;
                maxPoint = blockNum*blockSizeData;
                peakPoint = mod(peakTimes(iPeak),blockSizeData);
                
                currData = data(minPoint:maxPoint);
                hold off;
                plot(currData);
                hold all;
                if ~plotConditionsData
                    plot(peakPoint,min(currData)*2,'marker','*','color','k');
                else
                    currInds = peakStats.indsPerPeak{iPeak};
                    plot(currInds,min(currData)*2,'marker','*','color','k');
%                     [~,maxIndHP] = max(peakStats.zscoresPerPeaksEnv{iPeak});
%                     [~,maxIndAmp] = max(peakStats.zscoresPerPeaksAmp{iPeak});
%                     [~,maxIndGrad] = max(peakStats.zscoresPerPeaksGrad{iPeak});
%                     plot(currInds(maxIndHP),min(currData)*2.1,'marker','*','color','k');
%                     plot(currInds(maxIndAmp),min(currData)*2.2,'marker','*','color','b');
%                     plot(currInds(maxIndGrad),min(currData)*2.3,'marker','*','color','g');
%                     
%                     indsPassedThreshHP = peakStats.zscoresPerPeaksEnv{iPeak}>=obj.SDthresholdEnv;
%                     indsPassedThreshAmp = peakStats.zscoresPerPeaksAmp{iPeak}>=obj.SDthresholdAmp;
%                     indsPassedThreshGrad = peakStats.zscoresPerPeaksGrad{iPeak}>=obj.SDthresholdGrad;
%                     if sum(indsPassedThreshHP)>0
%                         plot(peakStats.indsPerPeak{iPeak}(indsPassedThreshHP),min(currData)*2.1,'marker','*','color','r');
%                     end
%                     if sum(indsPassedThreshAmp)>0
%                         plot(peakStats.indsPerPeak{iPeak}(indsPassedThreshAmp),min(currData)*2.2,'marker','*','color','b');
%                     end
%                     if sum(indsPassedThreshGrad)>0
%                         plot(peakStats.indsPerPeak{iPeak}(indsPassedThreshGrad),min(currData)*2.3,'marker','*','color','g');
%                     end
                end
                
                if plotConditionsData
                    title(['spike #', num2str(iPeak),' passed conditions ',num2str(peakStats.passedConditions(iPeak,:)),' max zscores: HP(red) = ', num2str(peakStats.zscoresPerPeaksMax(iPeak,1)), ' Amp(blue) = ', ...
                        num2str(peakStats.zscoresPerPeaksMax(iPeak,3)), ' Grad(green) = ', num2str(peakStats.zscoresPerPeaksMax(iPeak,2))]);
                end
                hold off;

                if plotZScores
                    currData = data(minPoint:maxPoint);
                    
                    zsAmp = zscore(currData);
                    
                    dataGradient = [0 diff(currData)];
                    zsGrad = zscore(dataGradient);
                    
                    filteredBlock = obj.bandpass(currData, obj.samplingRate, obj.lowCut, obj.highCut);
                    %find envelope
                    envBlock = abs(hilbert(filteredBlock));
                    zsHP = zscore(envBlock);
                    
                    subplot(4,blockSizeToPlot,indBlock+blockSizeToPlot);
                    plot(zsHP,'color','r');
                    title('HP zscore');
                    hold all;
                    line([0,blockSizeData], [obj.SDthresholdEnv,obj.SDthresholdEnv],'color','k');
                    line([0,blockSizeData], [obj.SDthresholdConjEnv,obj.SDthresholdConjEnv],'color','k');
                    hold off;
                    
                    subplot(4,blockSizeToPlot,indBlock+blockSizeToPlot*2);
                    plot(zsAmp,'color','b');
                    title('Amplitude zscore');
                    hold all;
                    line([0,blockSizeData], [obj.SDthresholdAmp,obj.SDthresholdAmp],'color','k');
                    line([0,blockSizeData], [obj.SDthresholdConjAmp,obj.SDthresholdConjAmp],'color','k');
                    hold off;
                    
                    subplot(4,blockSizeToPlot,indBlock+blockSizeToPlot*3);
                    plot(zsGrad,'color','g');
                    title('Gradient zscore');
                    hold all;
                    line([0,blockSizeData], [obj.SDthresholdGrad,obj.SDthresholdGrad],'color','k');
                    line([0,blockSizeData], [obj.SDthresholdConjGrad,obj.SDthresholdConjGrad],'color','k');
                    hold off;
                end
                
                indBlock = indBlock+1;
                if indBlock > blockSizeToPlot || iPeak == nPeaks
                    indBlock = 1;
                    pause;
                end
            end
        end
        
        
        
        
        
        
        
    end
    
    methods (Access = private)
        function [seqPeaks,peakAllInds] = findSequences(obj,data, dataThresh)
            %help function for frequency analysis detection -
            %receives the vector of binary values and finds sequence of ones in length obj.minLengthSpike, which are
            %separated from each other by at least obj.minDistSpikes, returns peak
            %times where each. Remove sequences that are above
            %obj.maxLengthSpike.
            
            %obj.minLengthSpike and minDistSpikes are in ms - translate to number of data
            %points
            numConsSpikes = round(obj.minLengthSpike*obj.samplingRate/1000);
            distSpikePoints = round(obj.minDistSpikes*obj.samplingRate/1000);
            maxLength = round(obj.maxLengthSpike*obj.samplingRate/1000);
            
            %finding sequences of 1's in length obj.minLengthSpike
            accDataBlock = dataThresh;
            threshInds = find(dataThresh);
            for iSeqLength = 1:numConsSpikes-1
                accDataBlock = accDataBlock(1:end-1)+dataThresh(iSeqLength+1:end);
            end
            seqInds = find(accDataBlock>=numConsSpikes);
            seqPeaks = [];
            peakAllInds = {};
            rmvDet = [];
            
            %merge peaks which are close together to one spike
            if ~isempty(seqInds)
                seqDists = diff(seqInds);
                distinctSeq = find(seqDists>distSpikePoints)+1;
                seqBeginnings = [seqInds(1) seqInds(distinctSeq)];
                nSeqs = length(seqBeginnings);
                seqPeaks = zeros(1,nSeqs);
                for iPeak = 1:nSeqs
                    lowerLimit = seqBeginnings(iPeak);
                    if length(seqBeginnings)>iPeak
                        upperLimit = seqInds(distinctSeq(iPeak)-1)+numConsSpikes-1;
                    else
                        upperLimit = seqInds(end)+numConsSpikes-1;
                    end
                    currInds = lowerLimit:upperLimit;
                    
                    if length(currInds) > maxLength
                        rmvDet(end+1) = iPeak;
                    end
                    
                    peakAllInds{iPeak} = currInds;
                    %the point which defines a peaks is choden to be the
                    %point with the maximal data value among the points
                    %which passed the threshold
                    [~,maxPoint] = max(data(currInds));
                    seqPeaks(iPeak) = currInds(maxPoint);

                end
                seqPeaks(rmvDet) = [];
                peakAllInds(rmvDet) = [];
            end
        end
    end
end