%% fitHMMs
%
% Fit Hidden Markov Models to spiking data from each experimental session
% or from each simulated session.
%
% This script loads the 'spikes' structure and 'win_train' matrix for each
% experimental or simulated session and proceeds to fit an HMM to each one.
%
% To actually fit the HMMs, you will need to download two additional
% package folders ('+hmm' and '+aux') and save them in your AnalysisDemo
% directory. The package folders contain code written by Luca Mazzucato.
%
% You can save the outputs (it will take a very long time to run) or
% download the already-fitted HMMs and use them for the next step. 
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
% requires access to functions:
% loadVar (in +fun)
% funHMM (main function in +hmm, which requires several other functions in +aux and +hmm folders)

%% parameters
expOrSim = 'experiment'; % 'experiment', 'simulation'
data = 'original'; % 'original', 'shuffled_circular', 'shuffled_swap'
autoSaveOutput = false;

%% setup
homeDir = pwd; addpath(homeDir); % get access to package functions
inputFiles1 = {}; inputFiles2 = {};
if strcmp(data,'shuffled_circular')
    append = '_shuff_circ';
elseif strcmp(data,'shuffled_swap')
    append = '_shuff_swap';
else
    append = '';
end
if strcmp(expOrSim,'experiment')
    cd('ProcessedData'); cd('Experiment'); loadDir = pwd; cd(homeDir);
    if autoSaveOutput
        cd('HMMData'); cd('Experiment'); saveDir = pwd; cd(homeDir);
    end
    for i = 1:21
        inputFiles1 = [inputFiles1,sprintf('spikes_exp%i%s.mat',i,append)];
        inputFiles2 = [inputFiles2,sprintf('win_train_exp%i.mat',i)];
    end
elseif strcmp(expOrSim,'simulation')
    cd('ProcessedData'); cd('Simulation'); loadDir = pwd; cd(homeDir);
    if autoSaveOutput
        cd('HMMData'); cd('Simulation'); saveDir = pwd; cd(homeDir);
    end
    for i = 245:254
        inputFiles1 = [inputFiles1,sprintf('spikes_sim%i%s.mat',i,append)];
        inputFiles2 = [inputFiles2,sprintf('win_train_sim%i.mat',i)];
    end
end

%% HMM analysis
for i = 1:length(inputFiles1)
    % load data 
    cd(loadDir);
    spikes = fun.loadVar(inputFiles1{i});
    win_train = fun.loadVar(inputFiles2{i});
    cd(homeDir);
    [ntrials, gnunits] = size(spikes);
    % HMM parameters
    MODELSEL = 'BIC'; % 'XVAL', 'AIC'
    DATAIN = struct('spikes',spikes,'win',win_train,'METHOD',MODELSEL); 
    % pass all to main HMM function
    res = hmm.funHMM(DATAIN);
    % save file
    if autoSaveOutput
        cd(saveDir);
        if strcmp(expOrSim,'experiment')
            save(sprintf('HMM_exp%i%s.mat',i,append));
        elseif strcmp(expOrSim,'simulation')
            save(sprintf('HMM_sim%i%s.mat',i+244,append));
        end
        cd(homeDir);
    end
end
