%% fitHMMs
%
% Fit Hidden Markov Models to spiking data from each experimental session
% or from each simulated session.
%
% This script loads the 'spikes' structure and 'win_train' matrix for each
% experimental or simulated session and proceeds to fit an HMM to each one.
%
% Actually fitting the HMMs requires two package folders ('+hmm' and '+aux').
% These package folders contain code written by Luca Mazzucato:
% https://github.com/mazzulab/contamineuro_2019_spiking_net
%
% You can simply use the already-fitted HMMs in /HMMData for the next step. 
%
% Note that if you fit a new HMM using this script it may not be exactly
% the same as the already-fitted one because the algorithm uses random
% initial conditions and the random number generator seed is not set prior
% to fitting.
%
% -LL
%

%% dependencies
% requires access to data:
% /ProcessedData/Experiment/spikes_expX.mat for X in 1:21
% /ProcessedData/Experiment/spikes_expX_shuff_circ.mat for X in 1:21
% /ProcessedData/Experiment/spikes_expX_shuff_colswap for X in 1:21
% /ProcessedData/Experiment/win_train_expX.mat for X in 1:21
% /ProcessedData/Simulation/spikes_simX.mat for X in 245:254
% /ProcessedData/Simulation/spikes_simX_shuff_circ.mat for X in 245:254
% /ProcessedData/Simulation/spikes_simX_shuff_colswap.mat for X in 245:254
% /ProcessedData/Simulation/win_train_simX.mat for X in 245:254
%
% requires access to functions:
% loadVar (in +fun)
% funHMM (main function in +hmm, which requires several other functions in +aux and +hmm folders)

%% parameters
expOrSim = 'simulation'; % 'experiment', 'simulation'
data = 'unshuffled'; % 'unshuffled', 'shuffled_circular', 'shuffled_swap'
stimType = 'strong'; % 'original', 'long', 'strong' (only matters if expOrSim is 'simulation')
autoSave = false;

%% setup
addpath(pwd); % get access to package functions
inputFiles1 = {}; inputFiles2 = {};
if strcmp(data,'shuffled_circular')
    append = '_shuff_circ';
elseif strcmp(data,'shuffled_swap')
    append = '_shuff_swap';
else
    append = '';
end
if strcmp(expOrSim,'experiment')
    for i = 1:21
        inputFiles1 = [inputFiles1,sprintf('%s/ProcessedData/Experiment/spikes_exp%i%s.mat',pwd,i,append)];
        inputFiles2 = [inputFiles2,sprintf('%s/ProcessedData/Experiment/win_train_exp%i.mat',pwd,i)];
    end
elseif strcmp(expOrSim,'simulation')
    switch stimType
        case 'original'
            append2 = '';
        case 'long'
            append2 = '_longStim';
        case 'strong'
            append2 = '_strongStim';
    end
    for i = 245:254
        inputFiles1 = [inputFiles1,sprintf('%s/ProcessedData/Simulation/spikes_sim%i%s%s.mat',pwd,i,append,append2)];
        inputFiles2 = [inputFiles2,sprintf('%s/ProcessedData/Simulation/win_train_sim%i%s.mat',pwd,i,append2)];
    end
end

%% HMM analysis
for i = 1:length(inputFiles1)
    % load data 
    spikes = fun.loadVar(inputFiles1{i});
    win_train = fun.loadVar(inputFiles2{i});
    [ntrials, gnunits] = size(spikes);
    % HMM parameters
    MODELSEL = 'BIC'; % 'XVAL', 'AIC'
    DATAIN = struct('spikes',spikes,'win',win_train,'METHOD',MODELSEL); 
    % pass all to main HMM function
    res = hmm.funHMM(DATAIN);
    
    % save file
    if autoSave
        if strcmp(expOrSim,'experiment')
            save(sprintf('%s/HMMData/Experiment/HMM_exp%i%s.mat',pwd,i,append),'res','spikes','win_train');
        elseif strcmp(expOrSim,'simulation')
            save(sprintf('%s/HMMData/Simulation/HMM_sim%i%s%s.mat',pwd,i+244,append,append2),'res','spikes','win_train');
        end
    end
    
end
