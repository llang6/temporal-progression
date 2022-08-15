%% prePackageExpDataForHMM
%
% Pre-package the raw experimental data for input to HMM analysis
%
% All raw experimental data (spiking data, event timestamps, and trial
% types) are contained in the 'SessionInfo' structure.
%
% The HMM script requires two inputs: 'spikes' and 'win_train'.
%
% This script loads 'SessionInfo' and converts the data inside into
% 'spikes' and 'win_train'.
%
% You can choose to save the outputs or download the already-formed files
% for use in the next step.
%
% -LL 
%

%% dependencies
% requires access to data:
% /RawData/Experiment/SessionInfo.mat

%% parameters
sessions = 1; % any value or range of values from 1 to 21
preSampTime = 0.1; % spike trains clipped before preSampTime [s] before taste event
postDecTime = 0.1; % spike trains clipped after postDecTime [s] after decision event
minFR = 2; % neurons firing at < minFR [Hz] will be removed
autoSaveOutput = false;

%% setup and load data
homeDir = pwd; addpath(homeDir);
if ~exist('SessionInfo','var')
    cd('RawData'); cd('Experiment');
    load('SessionInfo.mat'); 
    cd(homeDir);
end
if autoSaveOutput
    cd('ProcessedData'); cd('Experiment');
    saveDir = pwd;
    cd(homeDir);
end

%% loop over sessions
for session = sessions
    fprintf('\nPre-packaging experiment session %i...',session);
    % check parameters for errors
    if ~ismember(session,1:21), error('Session %i does not exist',session); end
    if SessionInfo(session).TrialEvents.SampleTimes(1) < preSampTime % first case
        error('''preSampTime'' too large for trial 1');
    end
    for i = 1:1:length(SessionInfo(session).TrialEvents.DecisionTimes)-1 % middle cases
        if (SessionInfo(session).TrialEvents.DecisionTimes(i)+postDecTime) > (SessionInfo(session).TrialEvents.SampleTimes(i+1)-preSampTime)
            error('Window overlap between trials %3i and %3i',i,i+1);
        end
    end
    if (SessionInfo(session).TrialEvents.DecisionTimes(end)+postDecTime) > SessionInfo(session).TimeLength % last case
        error('''postDecTime'' too large for last trial');
    end
    % construct 'spikes' struct and 'win_train' matrix
    spikes = []; 
    spikes(SessionInfo(session).NumTrials,SessionInfo(session).NumNeurons).spk = []; 
    win_train = NaN(SessionInfo(session).NumTrials,2); 
    for trial = 1:SessionInfo(session).NumTrials 
        for neuron = 1:SessionInfo(session).NumNeurons 
            % everything is in terms of alignment to decision event
            spikes_align_dec = SessionInfo(session).Neuron(neuron).SpikeTrain - SessionInfo(session).TrialEvents.DecisionTimes(trial);
            pre_time = SessionInfo(session).TrialEvents.SampleTimes(trial) - SessionInfo(session).TrialEvents.DecisionTimes(trial) - preSampTime;
            post_time = postDecTime;
            spikes(trial,neuron).spk = (spikes_align_dec(spikes_align_dec>=pre_time & spikes_align_dec<=post_time))';
            win_train(trial,:) = [pre_time,post_time];
        end
    end
    % remove quiet neurons
    quietNeurons = [];
    for i = 1:1:length(SessionInfo(session).Neuron)
        aveFR = sum(arrayfun(@(trial)length(spikes(trial,i).spk),1:size(spikes,1)))/sum(diff(win_train')); 
        if aveFR < minFR
            quietNeurons = [quietNeurons,i];
        end
    end
    spikes_active = spikes;
    spikes_active(:,quietNeurons) = [];
    % save output
    if autoSaveOutput
        cd(saveDir);
        save(sprintf('spikes_exp%i.mat',session),'spikes_active');
        save(sprintf('win_train_exp%i.mat',session),'win_train');
        cd(homeDir);
    end
    fprintf('Done.\n');
end