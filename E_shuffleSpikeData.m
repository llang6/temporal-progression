%% shuffleSpikeData
%
% Create two surrogate datasets from pre-packaged data for use in HMM
% control analyses
%
% We have pre-packaged data in 'spikes' and 'win_train' for each
% experimental and simulation session.
%
% This script will load the pre-packaged data and form two new shuffled versions of
% 'spikes' -- one based on a "circular" shuffle and one based on a "swap"
% shuffle.
% 
% The circular shuffle pushes each spike train ahead by a random amount,
% and the portion of the train that "goes off the edge" comes back to the
% front as if the front and end were connected. This shuffle preserves
% auto-correlations (and firing rates) and disrupts cross-correlations.
%
% The swap shuffle bins ensemble spike trains and randomly permutes them
% across time. This shuffle preserves cross-correlations (and firing rates)
% and disrupts auto-correlations.
%
% You can simply use the already-formed files in /ProcessedData for the next step.
%
% Refer to the definition of the main (local) function as well as the
% optional output plots to see how the shuffling works.
%
% -LL
%

%% dependencies
% requires access to data:
% /ProcessedData/Experiment/spikes_expX.mat for X in 1:21
% /ProcessedData/Experiment/win_train_expX.mat for X in 1:21
% /ProcessedData/Simulation/spikes_simX.mat for X in 245:254
% /ProcessedData/Simulation/win_train_simX.mat for X in 245:254
%
% requires access to functions:
% loadVar (in +fun folder)

%% parameters
sessionsExp = 1;
sessionsSim = 245;
binSize = 0.005; % bin size for swap shuffle, [s]
isPlot = true;

%% setup 
addpath(pwd); % gives access to helper functions in /+fun

%% shuffle all the data
% experiment
for session = sessionsExp
    spikes = fun.loadVar(sprintf('%s/ProcessedData/Experiment/spikes_exp%i.mat',pwd,session));
    win_train = fun.loadVar(sprintf('%s/ProcessedData/Experiment/win_train_exp%i.mat',pwd,session));
    [spikes_shuff_circ,spikes_shuff_swap,circMap,swapMap] = main(spikes,win_train,binSize,isPlot);

    %save(sprintf('%s/ProcessedData/Experiment/spikes_exp%i_shuff_circ.mat',pwd,session),spikes_shuff_circ);
    %save(sprintf('%s/ProcessedData/Experiment/spikes_exp%i_shuff_swap.mat',pwd,session),spikes_shuff_swap);

end
% simulation
for session = sessionsSim
    spikes = fun.loadVar(sprintf('%s/ProcessedData/Simulation/spikes_sim%i.mat',pwd,session));
    win_train = fun.loadVar(sprintf('%s/ProcessedData/Simulation/win_train_sim%i.mat',pwd,session));
    [spikes_shuff_circ,spikes_shuff_swap,circMap,swapMap] = main(spikes,win_train,binSize,isPlot);
    
    %save(sprintf('%s/ProcessedData/Simulation/spikes_sim%i_shuff_circ.mat',pwd,session),spikes_shuff_circ);
    %save(sprintf('%s/ProcessedData/Simulation/spikes_sim%i_shuff_swap.mat',pwd,session),spikes_shuff_swap);

end

%% plot results of shuffling an example trial
if isPlot
    trial = 20;
    figure(1); clf;
    subplot(2,1,1);
    [numTrial,numNeurons] = size(spikes);
    for i = 1:numNeurons
        if isempty(spikes(trial,i).spk), continue; end
        x = [spikes(trial,i).spk'; spikes(trial,i).spk'];
        y = [(i-0.4)*ones(1,length(spikes(trial,i))); (i+0.4)*ones(1,length(spikes(trial,i)))];
        plot(x,y,'k'); hold on;
        firstSpike = spikes(trial,i).spk(1);
        plot([firstSpike; firstSpike],[(i-0.4); (i+0.4)],'r'); hold on;
    end
    xlim(win_train(trial,:)); ylim([0,numNeurons+1]); yticks(1:numNeurons);
    xlabel('Time [s]','fontsize',22); ylabel('Neurons','fontsize',22);
    set(gca,'TickDir','out','color','none','box','off');
    title(sprintf('Trial %i spike trains',trial),'fontsize',22);
    subplot(2,1,2);
    for i = 1:numNeurons
        if isempty(spikes_shuff_circ(trial,i).spk), continue; end
        x = [spikes_shuff_circ(trial,i).spk'; spikes_shuff_circ(trial,i).spk'];
        y = [(i-0.4)*ones(1,length(spikes_shuff_circ(trial,i))); ...
             (i+0.4)*ones(1,length(spikes_shuff_circ(trial,i)))];
        plot(x,y,'k'); hold on;
        firstSpike = circMap(i,trial);
        plot([firstSpike; firstSpike],[(i-0.4); (i+0.4)],'r'); hold on;
    end
    xlim(win_train(trial,:)); ylim([0,numNeurons+1]); yticks(1:numNeurons);
    xlabel('Time [s]','fontsize',22); ylabel('Neurons','fontsize',22);
    set(gca,'TickDir','out','color','none','box','off');
    title(sprintf('Trial %i spike trains (circular shuffle)',trial),'fontsize',22);

    figure(2); clf;
    subplot(2,1,1);
    for i = 1:numNeurons
        if isempty(spikes(trial,i).spk), continue; end
        x = [spikes(trial,i).spk'; spikes(trial,i).spk'];
        y = [(i-0.4)*ones(1,length(spikes(trial,i))); ...
             (i+0.4)*ones(1,length(spikes(trial,i)))];
        plot(x,y,'k'); hold on;
    end
    numBins = ceil((win_train(trial,2)-win_train(trial,1))/binSize);
    xline(win_train(trial,1)+binSize*(swapMap(1,trial)-1),'r'); hold on; 
    xline(win_train(trial,1)+binSize*swapMap(1,trial),'r'); hold on;
    xline(win_train(trial,1)+binSize*(swapMap(2,trial)-1),'b'); hold on; 
    xline(win_train(trial,1)+binSize*swapMap(2,trial),'b'); hold on;
    xline(win_train(trial,1)+binSize*(swapMap(3,trial)-1),'g'); hold on; 
    xline(win_train(trial,1)+binSize*swapMap(3,trial),'g'); hold on;
    xlim([win_train(trial,1)-binSize,win_train(trial,2)+binSize]); xticks(win_train(trial,:));
    ylim([0,numNeurons+1]); yticks(1:numNeurons);
    xlabel('Time [s]','fontsize',22); ylabel('Neurons','fontsize',22);
    set(gca,'TickDir','out','color','none','box','off');
    title(sprintf('Trial %i spike trains',trial),'fontsize',22);
    subplot(2,1,2);
    for i = 1:numNeurons
        if isempty(spikes_shuff_swap(trial,i).spk), continue; end
        x = [spikes_shuff_swap(trial,i).spk'; spikes_shuff_swap(trial,i).spk'];
        y = [(i-0.4)*ones(1,length(spikes_shuff_swap(trial,i))); ...
             (i+0.4)*ones(1,length(spikes_shuff_swap(trial,i)))];
        plot(x,y,'k'); hold on;
    end
    xline(win_train(trial,1)+binSize*(swapMap(4,trial)-1),'r'); hold on; 
    xline(win_train(trial,1)+binSize*swapMap(4,trial),'r'); hold on;
    xline(win_train(trial,1)+binSize*(swapMap(5,trial)-1),'b'); hold on; 
    xline(win_train(trial,1)+binSize*swapMap(5,trial),'b'); hold on;
    xline(win_train(trial,1)+binSize*(swapMap(6,trial)-1),'g'); hold on; 
    xline(win_train(trial,1)+binSize*swapMap(6,trial),'g'); hold on;
    xlim([win_train(trial,1)-binSize,win_train(trial,2)+binSize]); xticks(win_train(trial,:));
    ylim([0,numNeurons+1]); yticks(1:numNeurons);
    xlabel('Time [s]','fontsize',22); ylabel('Neurons','fontsize',22);
    set(gca,'TickDir','out','color','none','box','off');
    title(sprintf('Trial %i spike trains (swap shuffle)',trial),'fontsize',22);
end

%% main shuffling function definition
function [SPIKES_SHUFF_CIRC,SPIKES_SHUFF_SWAP,varargout] = main(SPIKES,WIN_TRAIN,BIN_SIZE,IS_PLOT)
% Takes in SPIKES and WIN_TRAIN and outputs two shuffled versions of
%     SPIKES: a circular-shuffled version and a swap-shuffled version.
% A parameter BIN_SIZE is necessary for the swap-shuffling.
% There is also an option controlled by IS_PLOT that will return two
%     additional matrices (first one tells you where each neuron's first
%     spike ends up after the circular shuffle and second one tells you
%     where the first, middle, and last bin end up after the swap shuffle)
%     for use in the plotting portion of the script.
% The plots should help visualize how the shuffling procedures work.

% extract data
[numTrials, numNeurons] = size(SPIKES);
% initialize
SPIKES_SHUFF_CIRC = struct();
SPIKES_SHUFF_SWAP = struct();
if IS_PLOT, circMap = NaN(numNeurons,numTrials); swapMap = NaN(3,numTrials); end
% loop over trials
for trial = 1:numTrials
    % independent circular shuffle
    for i = 1:numNeurons
        randShift = (WIN_TRAIN(trial,2) - WIN_TRAIN(trial,1))*rand;
        division = WIN_TRAIN(trial,2) - randShift;
        train = SPIKES(trial,i).spk';
        newTrain = [train(train >= division) - (division - WIN_TRAIN(trial,1)), ...
                    train(train < division) + randShift]';
        SPIKES_SHUFF_CIRC(trial,i).spk = newTrain;
        if IS_PLOT
            if ~isempty(train)
                if train(1) >= division
                    circMap(i,trial) = train(1) - (division - WIN_TRAIN(trial,1));
                else
                    circMap(i,trial) = train(1) + randShift;
                end
            end
        end
    end
    % column swap shuffle
    numBins = ceil((WIN_TRAIN(trial,2)-WIN_TRAIN(trial,1))/BIN_SIZE);
    newBins = randperm(numBins);
    for i = 1:numNeurons
        train = SPIKES(trial,i).spk;
        newTrain = zeros(length(train),1);
        for j = 1:length(train)
            oldBin = min(floor((train(j)-WIN_TRAIN(trial,1))/BIN_SIZE) + 1,numBins);
            newTrain(j) = train(j) + (newBins(oldBin) - oldBin)*BIN_SIZE;
        end
        SPIKES_SHUFF_SWAP(trial,i).spk = newTrain;
    end
    if IS_PLOT
       % find some good example bins
       allSpikes = []; for n = 1:size(SPIKES,2), allSpikes = [allSpikes;SPIKES(trial,n).spk]; end
       numSpikes = arrayfun(@(i)sum((allSpikes>WIN_TRAIN(trial,1)+(i-1)*BIN_SIZE)&(allSpikes<WIN_TRAIN(trial,1)+i*BIN_SIZE)),1:numBins);
       [~,ind] = sort(numSpikes,'descend');
       swapMap(1,trial) = ind(1);
       swapMap(2,trial) = ind(2);
       swapMap(3,trial) = ind(3);
       swapMap(4,trial) = newBins(ind(1));
       swapMap(5,trial) = newBins(ind(2));
       swapMap(6,trial) = newBins(ind(3));
    end
end
if IS_PLOT
    varargout{1} = circMap;
    varargout{2} = swapMap;
end
end
