%% prePackageSimDataForHMM
%
% Pre-package the raw simulation data for input to HMM analysis
%
% All raw simulation data (spiking data, event timestamps, and trial
% types) are contained in the session's 'simulationResults' file.
%
% The HMM script requires two inputs: 'spikes' and 'win_train'.
%
% This script loads 'simulationData' and converts the data inside into
% 'spikes' and 'win_train'. There is random sampling involved -- the HMM
% will not be fit to the spike trains of every neuron in the simulated
% network.
%
% You can simply use the already-formed files in /ProcessedData for the next step.
%
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Simulation/simulationDataX.mat for X in 245:254

%% parameters
sessions = 245; % 245-254;
preStim = 0; % 0;
postStim = 3; % 3;
nSamp = 1; % 1; % sample this many neurons from each E cluster
downsample = 14; % 14; % final size of sample (if < nSamp*(# clusters), we are 'downsampling')

%% setup
homeDir = pwd; addpath(homeDir);
cd('RawData'); cd('Simulation'); 
loadDir = pwd; cd(homeDir);

%% loop over sessions
for session = sessions
    fprintf('\nPre-packaging simulation session %i',session);
    cd(loadDir); % could also add load/save dirs to path in setup, but this avoids confusion
    data = load(sprintf('simulationData%i.mat',session));
    cd(homeDir);
    firings_all = data.firings_all;
    ntrials = numel(firings_all);
    stimStart = data.stim.tStart; trialEnd = data.totalT/1000;
    if stimStart - preStim < 0, preStim = stimStart; end
    if stimStart + postStim > trialEnd, postStim = trialEnd - stimStart; end
    leftBound = (stimStart - preStim)*1000; rightBound = (stimStart + postStim)*1000;
    win_train = [-1*preStim*ones(ntrials,1), postStim*ones(ntrials,1)];
    sampNeurons = NaN(nSamp,data.net.Q); % indices of sampled neurons
    nCount = 0;
    for cluster = [data.clustersWithRoles,setdiff(1:data.net.Q,data.clustersWithRoles)]
        nCount = nCount + 1;
        indices = (data.net.indE(cluster)):(data.net.indE(cluster+1)-1);
        FRs = NaN(length(indices),1);
        for i = 1:length(indices)
            frates = zeros(1,ntrials);
            for j = 1:ntrials
                myArray = firings_all{j};
                mySpikes = myArray(myArray(:,2)==indices(i),1);
                frates(j) = sum(mySpikes>leftBound & mySpikes<rightBound)/(preStim+postStim);
            end
            FRs(i) = mean(frates);
        end
        if sum(FRs>=2)<2 
            fprintf('\nWARNING! Cluster %i has fewer than 2 neurons with firing rates at least 2 Hz. Randomly choosing %i and moving on...\n',cluster,nSamp);  
        else
            indices = indices(FRs>=2);
        end
        rN = indices(randperm(length(indices))); % random sampling
        sampNeurons(:,nCount) = reshape(rN(1:nSamp),[],1);
    end
    sampNeurons = reshape(sampNeurons,1,[]);
    sampNeurons = sampNeurons(sort(randperm(length(sampNeurons),downsample))); % random sub-sampling
    spikes = struct();
    for trial = 1:ntrials
        for neuron = 1:length(sampNeurons)
            f = firings_all{trial};
            f = f(f(:,2)==sampNeurons(neuron),1)/1000-stimStart;
            spikes(trial,neuron).spk = f(f>=win_train(trial,1)&f<=win_train(trial,2));
        end
    end
    
    %save(sprintf('spikes_sim%i.mat',session),'spikes');
    %save(sprintf('win_train_sim%i.mat',session),'win_train');

    fprintf('Done.\n');
end
