function SCORED_TRIALS = scoreTrials(DATA,varargin)
%SCORED_TRIALS = scoreTrials(DATA) 
%
%Scores network trials, where DATA is a 'simulationResults.mat' file
%
% Input DATA is the structure loading of the 'simulationDataX.mat' file
% (for X some integer).
%     i.e. DATA = load(sprintf('simulationData%i.mat',number));
%          SCORED_TRIALS = scoreTrials(DATA);
%
% Output SCORED_TRIALS is a column vector where 1 means correct trial, 0
% means incorrect trial, and NaN means omitted trial.
%
% Optional 'Name',Value pair arguments can be provided after DATA:
%     'binSize': width of bin used for calculating firing rates, [s]
%     'stateThreshType': if 'relative', the threshold value argument is
%                            treated as a fraction of the cluster's maximum
%                            firing rate
%                        if 'absolute', the threshold value argument is a
%                            firing rate
%     'stateThresh': the threshold value argument (see above for
%                        interpretation)
%     'stateMinTime': firing rates must be at least at threshold value for
%                         at least this long to be considered online
%     'time1': time after stimulus when we start checking for clusters to
%                  be online, [s]
%     'time2': time after stimulus when we stop checking for clusters to be
%                  online, [s]
%
% Trials were already scored during the simulation, but this function
% theoretically allows for more flexibility (i.e. you could change the
% scoring rule after having already run the simulation)
%
% As written here, the function calculates the firing rates of the left and
% right action clusters, thresholds them to binarize them, and then infers
% choices from whether or not each cluster came online during a specific
% window relative to stimulus onset (only left -> chose left, only right ->
% chose right, else -> omitted)
%
% -LL
%

% parse input
if isempty(varargin)
    % default parameters
    binSize = 0.05;
    stateThreshType = 'relative';
    stateThresh = 0.4;
    stateMinTime = 0.05;
    time1 = 0.5;
    time2 = 3.0;
else
    % name, value pair arguments
    ind = find(cellfun(@(x)strcmpi(x,'binSize'),varargin),1);
    if ~isempty(ind), binSize = varargin{ind+1}; else, binSize = 0.05; end
    ind = find(cellfun(@(x)strcmpi(x,'stateThreshType'),varargin),1);
    if ~isempty(ind), stateThreshType = varargin{ind+1}; else, stateThreshType = 'relative'; end
    ind = find(cellfun(@(x)strcmpi(x,'stateThresh'),varargin),1);
    if ~isempty(ind), stateThresh = varargin{ind+1}; else, stateThresh = 0.4; end
    ind = find(cellfun(@(x)strcmpi(x,'stateMinTime'),varargin),1);
    if ~isempty(ind), stateMinTime = varargin{ind+1}; else, stateMinTime = 0.05; end
    ind = find(cellfun(@(x)strcmpi(x,'time1'),varargin),1);
    if ~isempty(ind), time1 = varargin{ind+1}; else, time1 = 0.05; end
    ind = find(cellfun(@(x)strcmpi(x,'time2'),varargin),1);
    if ~isempty(ind), time2 = varargin{ind+1}; else, time2 = 0.05; end
end
% calculate action cluster firing rates
binEdges = -DATA.warmup/1000:binSize:DATA.totalT/1000;
trialSequence = DATA.stimuli;
nTrials = length(trialSequence);
firingRates = cell(nTrials,4);
allFirings = DATA.firings_all;
lActInd = DATA.inds.actLE;
rActInd = DATA.inds.actRE;
for trial = 1:nTrials
    firings = allFirings{trial};
    lActSpikes = firings(ismember(firings(:,2),lActInd),:);
    rActSpikes = firings(ismember(firings(:,2),rActInd),:);
    lActFR = zeros(1,length(binEdges)-1);
    rActFR = zeros(1,length(binEdges)-1);
    for i = 1:length(binEdges)-1
        lActFR(i) = sum(lActSpikes(:,1)/1000>=binEdges(i)&lActSpikes(:,1)/1000<binEdges(i+1))/length(lActInd)/binSize;
        rActFR(i) = sum(rActSpikes(:,1)/1000>=binEdges(i)&rActSpikes(:,1)/1000<binEdges(i+1))/length(rActInd)/binSize;
    end
    firingRates{trial,1} = lActFR;
    firingRates{trial,2} = rActFR;
end
% threshold action cluster firing rates to determine 'on' or 'off'
for trial = 1:nTrials
    if strcmp(stateThreshType,'relative')
        lActFR_max = max(cellfun(@max,firingRates(:,1)));
        rActFR_max = max(cellfun(@max,firingRates(:,2)));
        lActState = firingRates{trial,1}>stateThresh*lActFR_max;
        rActState = firingRates{trial,2}>stateThresh*rActFR_max;
    elseif strcmp(stateThreshType,'absolute')
        lActState = firingRates{trial,1}>stateThresh;
        rActState = firingRates{trial,2}>stateThresh;
    else
        error('Invalid value for optional argument ''stateThreshType''. Valid values are ''relative'' and ''absolute''.');
    end
    % remove states that are too short
    temp_onsets = find(diff([0 lActState 0])==1); temp_offsets = find(diff([0 lActState 0])==-1);
    temp_inds = find(binSize*(temp_offsets-temp_onsets)<stateMinTime);
    if ~isempty(temp_inds)
        for ind = temp_inds
            lActState(temp_inds(ind):min(temp_offsets(ind)-1,length(lActState))) = 0;
        end
    end
    temp_onsets = find(diff([0 rActState 0])==1); temp_offsets = find(diff([0 rActState 0])==-1);
    temp_inds = find(binSize*(temp_offsets-temp_onsets)<stateMinTime);
    if ~isempty(temp_inds)
        for ind = temp_inds
            rActState(temp_inds(ind):min(temp_offsets(ind)-1,length(rActState))) = 0;
        end
    end
    firingRates{trial,3} = lActState;
    firingRates{trial,4} = rActState;
end
% apply rule for classifying correct trials
corrChoices = reshape(double(cellfun(@(x)ismember(x,{'Maltose','Octaacetate'}),trialSequence)),[],1);
choices = NaN(nTrials,1);
SCORED_TRIALS = NaN(nTrials,1);
%b1 = ceil((DATA.stim.tStart+time1)/binSize);
%b2 = floor((DATA.stim.tStart+time2)/binSize);
b1 = floor((DATA.stim.tStart+time1)/binSize)+1;
b2 = ceil((DATA.stim.tStart+time2)/binSize);
choices(arrayfun(@(x)any(firingRates{x,3}(b1:b2))&~any(firingRates{x,4}(b1:b2)),(1:nTrials)')) = 0;  
choices(arrayfun(@(x)~any(firingRates{x,3}(b1:b2))&any(firingRates{x,4}(b1:b2)),(1:nTrials)')) = 1;
SCORED_TRIALS(choices==corrChoices) = 1; 
SCORED_TRIALS(choices~=corrChoices & ~isnan(choices)) = 0; % exclude no-decision trials (this was done in experimental data)
end