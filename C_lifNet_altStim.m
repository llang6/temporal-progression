%% lifNet_altStim
%
% Reload the original 10 networks from lifNet and run new simulations with
% modified stimuli
%
% -LL
%

%% dependencies
% requires access to data:
% /RawData/Simulation/simulationDataX.mat for X in 245:254
%
% requires access to functions:
% getActionGate (in +fun)
% simulation (in +fun)
% rasterPlot (in +fun)
% spikes2FR (in +fun)
% threshold (in +fun)

%% options and parameters
%close all; clearvars; clc; % start fresh
%rng(9995); % make simulation results repeatable
isPlot = true;
% --- stimulus properties ---
stim_tauFall = 0.160; % decay time constant, [s]
stim_Gain = 2.00; % strength/amplitude, [*100% increase in external current]

%% simulations
for network_session = 1:10
    
    clearvars -except stim_tauFall stim_Gain isPlot network_session;

    load(sprintf('%s/RawData/Simulation/simulationData%i.mat',pwd,244+network_session));

    % ==== fixed parameters ====
    % ---- time ----
    totalT = 3000; % trial time, [ms]

    % ==== setup ====
    timerStart_total = tic;
    rngAtStart = rng;
    addpath(pwd);

    % ==== stimulus ====
    % construct time vector
    timev = 0:dt:(totalT+warmup); 
    % construct the stimulus
    stim.tauFall = stim_tauFall;
    stim.Gain = stim_Gain;
    stim.Curve = zeros(1,length(timev));
    if strcmp(stim.Profile,'square')     
        stim.Curve(timev>=(warmup+1000*stim.tStart) & timev<=(warmup+1000*stim.tStart+1000*stim.Duration)) = stim.Gain;
    elseif strcmp(stim.Profile,'alpha')  
        stim.Curve(timev>=warmup+1000*stim.tStart) = (1/((stim.tauFall*1000)-(stim.tauRise*1000)))*(exp(-(0:dt:(totalT-1000*stim.tStart))/(stim.tauFall*1000))-exp(-(0:dt:(totalT-1000*stim.tStart))/(stim.tauRise*1000)));
        stim.Curve = stim.Curve/max(abs(stim.Curve))*stim.Gain;
    end

    % ==== main simulations: baseline ====
    timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nSimulating baseline trial...');

    gate.Curve = fun.getActionGate(gate.StartMean,gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
    [spkdata,~] = fun.simulation(net,neu,stim,inds,'',gate.Switch,gate.Curve,S,isModifyWeights,warmup,BinW,timev,'',false,'');

    fprintf('Done. ');
    toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % extract spike data
    firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end

    % raster plot: baseline
    if isPlot
        figure(1); clf; hold all;
        fun.rasterPlot(firings,inds,indexReMap_inv,clustersWithRoles,totalT,net);
        title('Baseline trial (no stimulus)');
        drawnow;
    end

    % ==== main simulations: trials ====
    firings_all = cell(numTrials,1);
    inputs_all = cell(numTrials,1);
    firingRates = cell(numTrials,6);
    stimuli = [repmat({'Sucrose'},1,ceil(numTrials/4)),... 
               repmat({'Maltose'},1,ceil(numTrials/4)),...
               repmat({'Quinine'},1,ceil(numTrials/4)),...
               repmat({'Octaacetate'},1,ceil(numTrials/4))];

    for iter = 1:numTrials

        timerStart_trial = tic; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\nSimulating trial %i...',iter);

        gate.Curve = fun.getActionGate(gate.StartTimes(iter),gate.Duration,gate.Floor,gate.Ceiling,warmup,totalT,dt);
        [spkdata,~] = fun.simulation(net,neu,stim,inds,stimuli{iter},gate.Switch,gate.Curve,S,isModifyWeights,...
            warmup,BinW,timev,'',false,'');

        fprintf('Done. ');
        toc(timerStart_trial); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % extract spikes and calculate firing rates
        firings = []; for i = 1:length(spkdata), firings = [firings ; spkdata{i}]; end
        firings_all{iter} = firings;
        bins = -warmup/1000:binSize:totalT/1000;
        [~,frLStim] = fun.spikes2FR(firings,bins,[inds.sucE,inds.quiE]); firingRates{iter,1} = frLStim;
        [~,frRStim] = fun.spikes2FR(firings,bins,[inds.malE,inds.octE]); firingRates{iter,2} = frRStim;
        [~,frLCue] = fun.spikes2FR(firings,bins,inds.cueLE); firingRates{iter,3} = frLCue;
        [~,frRCue] = fun.spikes2FR(firings,bins,inds.cueRE); firingRates{iter,4} = frRCue;
        [~,frLAct] = fun.spikes2FR(firings,bins,inds.actLE); firingRates{iter,5} = frLAct;
        [~,frRAct] = fun.spikes2FR(firings,bins,inds.actRE); firingRates{iter,6} = frRAct;

        % raster plot: trial
        if isPlot
            % new figure for each trial if there won't be too many
            figure(2); clf; hold all;
            fun.rasterPlot(firings,inds,indexReMap_inv,clustersWithRoles,totalT,net);
            title(sprintf('Trial: %i, Stimulus: %s',iter,stimuli{iter}));
            drawnow;
        end

    end

    % ==== threshold firing rates and score performance ====
    corrChoices = reshape(double(cellfun(@(x)ismember(x,{'Maltose','Octaacetate'}),stimuli)),[],1);
    firingRates_thresh = cell(numTrials,6);
    frLAct_max = max(cellfun(@max,firingRates(:,5)));
    frRAct_max = max(cellfun(@max,firingRates(:,6)));
    frLCue_max = max(cellfun(@max,firingRates(:,3)));
    frRCue_max = max(cellfun(@max,firingRates(:,4)));
    for TR = 1:numTrials    
        firingRates_thresh{TR,5} = fun.threshold(firingRates{TR,5},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
        firingRates_thresh{TR,6} = fun.threshold(firingRates{TR,6},stateThresh*max(frLAct_max,frRAct_max),stateMinTime,binSize);
        firingRates_thresh{TR,3} = fun.threshold(firingRates{TR,3},stateThresh*max(frLCue_max,frRCue_max),stateMinTime,binSize);
        firingRates_thresh{TR,4} = fun.threshold(firingRates{TR,4},stateThresh*max(frLCue_max,frRCue_max),stateMinTime,binSize);
    end

    choices = NaN(numTrials,1);
    bin1 = floor((warmup/1000+stim.tStart+bound1)/binSize)+1; 
    bin2 = ceil((warmup/1000+stim.tStart+bound2)/binSize);
    for trial = 1:numTrials
        if any(firingRates_thresh{trial,5}(bin1:bin2)) && ~any(firingRates_thresh{trial,6}(bin1:bin2))
            choices(trial) = 0;
        elseif ~any(firingRates_thresh{trial,5}(bin1:bin2)) && any(firingRates_thresh{trial,6}(bin1:bin2))
            choices(trial) = 1;
        end
    end

    fprintf('\nPerformance by action cluster activation from %.1f to %.1f s (correct/performed/total): %i/%i/%i\n\n',...
            bound1,bound2,sum(choices==corrChoices),sum(~isnan(choices)),numTrials);

    toc(timerStart_total);

    % save(sprintf('%s/RawData/Simulation/simulationData%i_tau%i_gain%i',pwd,244+network_session,round(stim_tauFall*1000),round(stim_Gain*100)));
    % for tau = 705 ms, gain = 60%, file was renamed 'simulationData%i_longStim.mat'
    % for tau = 160 ms, gain = 200%, file was renamed 'simulationData%i_strongStim.mat'

end
