function fun_PlotTrial(data,firings,Params,trialSequence,iTrial)

data_i = data{iTrial};
firings_i = firings{iTrial};
stimuliPresent = trialSequence(iTrial,:);

% default plots:
plot_raster = 1; % raster plot
plot_currents = 1; % postsynaptic current for representative E and I neurons
plot_rates = 1; % average cluster firing rates
plot_stimulus = 1; % stimuli time course

% load simulation results
Params.Sim.plot_length = Params.Sim.plot_length;

% save plots
HOME = Params.HOME;
savedir = Params.savedir;
if ~exist(savedir,'dir'); mkdir(savedir); end
filesave = fullfile(savedir,sprintf('[Jp%0.03g]_[p%d]',Params.Jplus,Params.p));

Ne_plot = data_i.Ne_plot; % number of exc neuron to plot (typically all)
Ni_plot = data_i.Ni_plot; % number of inh neurons to plot  (typically all)
N_e = Ne_plot;
N_i = Ni_plot;
N_decision = Params.N_decision;
ind_plot = data_i.ind_plot; % [5; N_e+5]; % indices of neurons to plot
vplot = data_i.vplot; % store membrane potential for plots; rows=neurons, cols=time steps;
iEplot = data_i.iEplot; % store EPSC for plots; rows=neurons, cols=time steps;
iExtplot = data_i.iExtplot; % store IPSC for plots; rows=neurons, cols=time steps;
iIplot = data_i.iIplot; % store IPSC for plots; rows=neurons, cols=time steps;
VTh = data_i.VTh;
Jee = Params.Jee;
tau = data_i.tau;
p = data_i.p;

Ext = Params.Ext;
Network = Params.Network;
Sim = Params.Sim;
if any(strcmp(fieldnames(Params),'NcE'))
    NcE = Params.NcE;
elseif any(strcmp(fieldnames(Params),'popsize'))
    NcE = Params.popsize;
end

% keep only pops with at least 1 neuron
NcE = NcE(NcE>0);
clustermatrix = Params.clustermatrix;
clustermatrix = clustermatrix(1:numel(NcE),:);
cusumNcE = [0 cumsum(NcE)'];
fprintf('only pops sharing up to %d clusters are populated.\n',sum(clustermatrix(end,:)));

XLIM = [Sim.t_Start Sim.t_End];
if Sim.t_Start < -5
    XLIM(1) = -5;
end
if Sim.t_End > 5
    XLIM(2) = 5;
end

fprintf('--- Plots saved in %s...\n',filesave);

%--------------------
% FIRST PLOT: rasters
%--------------------
numfig = 1;
if plot_raster
    if Params.plotAllTrials
        if any(cellfun(@(x) ~isempty(x), firings))
            option = 'poporder';
            fieldNames = {'fieldNames','N_e','Ne_plot','Ni_plot','NcE','N_decision','XLIM','option','p','clustermatrix'};
            IN = aux.v2struct(fieldNames);
            figure(numfig);
            for tr = 1:size(trialSequence,1)
                clf;
                PlotAllNeurons(Ext,firings{tr},IN,trialSequence(tr,:),tr);
                if Params.autoSave
                    cd(savedir);
                    saveas(gcf,[filesave '_rasterTrial' num2str(tr) '.jpg']);
                    cd(HOME);
                end
            end
            numfig = numfig + 1;
        end
    else
        if any(any(firings_i))
            figure(numfig); clf;
            numfig = numfig + 1;
            % A. Plot population rasters
            option = 'poporder';
            fieldNames = {'fieldNames','N_e','Ne_plot','Ni_plot','NcE','N_decision','XLIM','option','p','clustermatrix'};
            IN = aux.v2struct(fieldNames);
            PlotAllNeurons(Ext,firings_i,IN,stimuliPresent,iTrial);
            % B. Ensemble rasterplot (10 units)
            if Params.autoSave
                cd(savedir);
                saveas(gcf,[filesave '_rasterTrial' num2str(iTrial) '.jpg']);
                cd(HOME);
            end
        end
    end
end

%--------------------
% SECOND PLOT: PSC
%--------------------
if plot_currents
    figure(numfig); clf;
    set(gcf, 'Position',  [100, 100, 800, 800])
    numfig = numfig+1;
    nplot = min(size(iIplot,1),2);
    np = nplot+1;
    % A. Plot membrane potentials
%         VTh = VTh*15 - 65;
%         vplot = vplot*15 - 65*ones(size(vplot)); % trasform to mV
    subplot(np,1,1);
    PlotMembrane(Sim,nplot,vplot/Jee-100,firings_i,ind_plot,VTh/Jee-100);
    xlim(XLIM);
    % D. Plot currents
    % bandpass IPSC
    hicutoff = 200; % hi-edge (Hz)
    srate = 10000; % sampling rate (Hz)
    PSCPlot(Sim,Ext,iIplot*Jee,iEplot*Jee,iExtplot*Jee,hicutoff,srate,VTh*Jee./(tau),ind_plot);
    xlim(XLIM);
    % SAVE PLOT
    if Params.autoSave
        cd(savedir);
        saveas(gcf,[filesave '_PSCTrial' num2str(iTrial) '.jpg']);
        cd(HOME);
    end
end

%----------------------------
% THIRD PLOT: summary stats
%----------------------------
if plot_rates
    if p > 1 && any(any(firings_i))
        
        figure(numfig); clf; 
        BinSize = 0.025;
        bins = Sim.t_End-Sim.plot_length:BinSize:Sim.t_End;
        
        % A. Cluster PSTH
        subplot(2,1,1)
        [rate, x_bins, numfigFlag] = PlotPopRates(bins,firings_i,N_e,N_i,NcE,N_decision,numfig);
        xlim(XLIM);
        % B. Single units PSTH from ensemble
        subplot(2,1,2)
        PlotPopRaster(rate,x_bins,bins);
        xlim(XLIM);
        numfig = numfig + 1 + numfigFlag;
        % SAVE PLOT
        if Params.autoSave
            cd(savedir);
            saveas(gcf,[filesave '_ratesTrial' num2str(iTrial) '.jpg']);
            cd(HOME);
        end
    end
end

if plot_stimulus
    if any(strcmp(fieldnames(Ext),'stim'))
        Tseq = Sim.t_Start:Sim.dt_step:Sim.t_End;
        nstim = numel(Ext.stim);
        Leg = cell(1,nstim);
        figure(numfig); clf; h = [];  numfig = numfig + 1; hold on;
        for n = 1:nstim
            subplot(nstim,1,n);
            if any(strcmp(Ext.stim(n).name,{'Sucrose','Quinine','Maltose','Octaacetate'}))
                y_inc = Ext.stim(n).profile_inc;
                y_dec = Ext.stim(n).profile_dec;
                indstim = Tseq>Ext.stim(n).interval(1) & Tseq<Ext.stim(n).interval(2);
                yplot_inc = zeros(1,numel(Tseq));
                yplot_dec = yplot_inc;
                yplot_inc(indstim) = 100*y_inc(Tseq(indstim));
                yplot_dec(indstim) = 100*y_dec(Tseq(indstim));
                mods = unique(Ext.stim(n).DiffGain_Inc);
                for i = 1:numel(mods)
                    plot(Tseq(2:end-1),mods(i)*yplot_inc(2:end-1),'r','linewidth',2); hold on;
                end
                mods = unique(Ext.stim(n).DiffGain_Dec);
                for i = 1:numel(mods)
                    plot(Tseq(2:end-1),mods(i)*yplot_dec(2:end-1),'b','linewidth',2);
                end
                Leg{n} = Ext.stim(n).name;
                xlab = ''; ylab = '';
                tt = sprintf('%s [',Ext.stim(n).name);
                for c = sort(Ext.stim(n).StimClust_Inc')
                    tt = [tt ' ' num2str(c)];
                end
                tt = [tt ' ;'];
                for c = sort(Ext.stim(n).StimClust_Dec')
                    tt = [tt ' ' num2str(c)];
                end
                tt = [tt ']-clusters'];
                if n == 1
                    ylab = ('Stim. gain (%)');
                    xlab = ('Time [s]');
                end
                fntsz = 15;
                aux.figset(gca,xlab,ylab,tt,fntsz);
            else
                y = Ext.stim(n).profile;
                indstim = Tseq>Ext.stim(n).interval(1) & Tseq<Ext.stim(n).interval(2);
                yplot = zeros(1,numel(Tseq));
                yplot(indstim) = 100*y(Tseq(indstim));
                plot(Tseq(2:end-1),yplot(2:end-1),'r','linewidth',2);
                Leg{n} = Ext.stim(n).name;
                xlab = ''; ylab = '';
                tt = sprintf('%s [',Ext.stim(n).name);
                for c = find(Ext.stim(n).StimClust)
                    tt = [tt ' ' num2str(c)];
                end
                tt = [tt ']-clusters'];
                if n == 1
                    ylab = ('Stim. gain (%)');
                    xlab = ('Time [s]');
                end
                fntsz = 15;
                aux.figset(gca,xlab,ylab,tt,fntsz);
            end
        end
        xlim(XLIM);
        if Params.autoSave
            cd(savedir);
            saveas(gcf,[filesave '_Stim.jpg']);
            cd(HOME);
        end
        % if gaussian stimulus, plot distribution
        % index of gaussian stimuli
        ind_gauss = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'gauss')), {Ext.stim(:).name},'uniformoutput',false)));
        if any(strcmp(Params.stimuli,'CSgauss'))
            figure(numfig); clf; h = [];  numfig = numfig+1; hold on;
            for i_g = 1:numel(ind_gauss)
                subplot(numel(ind_gauss),1,i_g);
                distr = Ext.stim(ind_gauss(i_g)).gauss;
                %                 [f,x]=hist(distr);
                %                 plot(x,f,'linewidth',3);
                OPTIONS = struct('colors','m','yvalues','count','XLim',[-1 1]);
                h = aux.PlotHistStairs(distr,30,OPTIONS);
                xlab = ''; ylab = ''; tt = ['distr. of ' Ext.stim(ind_gauss(i_g)).name '; std=' num2str(std(distr))];
                if i_g == 1
                    xlab = '\delta\nu_{ext}'; ylab = 'occurrence';
                end
                xlim([-1 1]);
                %                 ylim([0 min(200,numel(distr)/4)]);
                aux.figset(gca,xlab,ylab,tt,fntsz);
            end
            if Params.autoSave
                cd(savedir);
                saveas(gcf,[filesave '_StimGauss.jpg']);
                cd(HOME);
            end
        end
    end
end

%-----------------
% AUX FUN
%-----------------

% 1A. Total raster plot

function PlotAllNeurons(Ext,firings,IN,stimuliPresent,trialNum)

    hold on;
    % unpack
    aux.v2struct(IN);
    %
    cusumNcE = [0 cumsum(NcE)'];
    % Plot all excitatory neuron spikes in black
    plot(firings(firings(:,2)<=Ne_plot,1),firings(firings(:,2)<=Ne_plot,2),'k.','markersize',1);
    % Plot all inhibitory neuron spikes in red
    plot(firings(firings(:,2)>N_e & firings(:,2)<=N_e+Ni_plot,1),firings(firings(:,2)>N_e & firings(:,2)<=N_e+Ni_plot,2)-N_e+Ne_plot,'r.','markersize',1);
    % Plot left decision cluster spikes in green
    plot(firings(firings(:,2)>N_e+N_i & firings(:,2)<=N_e+N_i+N_decision,1),firings(firings(:,2)>N_e+N_i & firings(:,2)<=N_e+N_i+N_decision,2)-N_e-N_i-2*N_decision,'g.','markersize',1);
    % Plot right decision cluster spikes in cyan
    plot(firings(firings(:,2)>N_e+N_i+N_decision & firings(:,2)<=N_e+N_i+2*N_decision,1),firings(firings(:,2)>N_e+N_i+N_decision & firings(:,2)<=N_e+N_i+2*N_decision,2)-N_e-N_i-2*N_decision,'c.','markersize',1);
    if any(strcmp(fieldnames(Ext),'stim'))
        for n = 1:numel(Ext.stim)
            time_bin = [Ext.stim(n).interval(1) Ext.stim(n).interval(2)];
            if ismember(Ext.stim(n).name,stimuliPresent) && strcmp(Ext.stim(n).name,'US') && ~isempty(Ext.stim(n).ind) 
                % find pops with target neurons
                pop_ind = zeros(1,numel(cusumNcE)-1);
                for k = 1:numel(cusumNcE)-1
                    pop_ind(k) = any(Ext.stim(n).ind>=cusumNcE(k)+1 & Ext.stim(n).ind<=cusumNcE(k+1));
                end
                pop_ind = find(pop_ind);
                for k = 1:numel(pop_ind)
                    lower = (cusumNcE(pop_ind(k)));
                    upper = cusumNcE(pop_ind(k)+1);
                    indspikes = firings(:,2)>=lower & firings(:,2)<=upper;
                    indstim = firings(:,1)>=Ext.stim(n).interval(1) & firings(:,1)<Ext.stim(n).interval(2);
                    indspikes = indspikes & indstim;
                    plot(firings(indspikes,1),firings(indspikes,2),'m.','markersize',1);
                    hold on;
                end
            end
            % PLOT NEURONS THAT RESPOND TO STIMULUS
            if ismember(Ext.stim(n).name,stimuliPresent) && any(strcmp(Ext.stim(n).name,{'Sucrose','Quinine','Maltose','Octaacetate'})) && (~isempty(Ext.stim(n).ind_inc) || ~isempty(Ext.stim(n).ind_dec)) 
                % override plot neurons that explicitly respond to stimulus
                plot(firings(ismember(firings(:,2),Ext.stim(n).ind_inc),1),firings(ismember(firings(:,2),Ext.stim(n).ind_inc),2)-N_e+Ne_plot,'m.','markersize',1); hold on;
                plot(firings(ismember(firings(:,2),Ext.stim(n).ind_dec),1),firings(ismember(firings(:,2),Ext.stim(n).ind_dec),2)-N_e+Ne_plot,'b.','markersize',1); hold on;
                % find pops with target neurons
%                 pop_ind_inc = zeros(1,numel(cusumNcE)-1);
%                 pop_ind_dec = pop_ind_inc;
%                 for k = 1:numel(cusumNcE)-1
%                     pop_ind_inc(k) = any(Ext.stim(n).ind_inc>=cusumNcE(k)+1 & Ext.stim(n).ind_inc<=cusumNcE(k+1));
%                     pop_ind_dec(k) = any(Ext.stim(n).ind_dec>=cusumNcE(k)+1 & Ext.stim(n).ind_dec<=cusumNcE(k+1));
%                 end
%                 pop_ind_inc = find(pop_ind_inc);
%                 pop_ind_dec = find(pop_ind_dec);
%                 if ~isempty(intersect(pop_ind_inc,pop_ind_dec)), error('Clusters are both excited and inhibited by stimulus'); end
%                 for k = 1:numel(pop_ind_inc)
%                     lower = cusumNcE(pop_ind_inc(k));
%                     upper = cusumNcE(pop_ind_inc(k)+1);
%                     indspikes = firings(:,2)>=lower & firings(:,2)<=upper;
%                     indstim = firings(:,1)>=Ext.stim(n).interval(1) & firings(:,1)<Ext.stim(n).interval(2);
%                     indspikes = indspikes & indstim;
%                     plot(firings(indspikes,1),firings(indspikes,2),'m.','markersize',1);
%                     hold on;
%                 end
%                 for k = 1:numel(pop_ind_dec)
%                     lower = (cusumNcE(pop_ind_dec(k)));
%                     upper = cusumNcE(pop_ind_dec(k)+1);
%                     indspikes = firings(:,2)>=lower & firings(:,2)<=upper;
%                     indstim = firings(:,1)>=Ext.stim(n).interval(1)-0.5 & firings(:,1)<Ext.stim(n).interval(2)-0.5;
%                     indspikes = indspikes & indstim;
%                     plot(firings(indspikes,1),firings(indspikes,2),'b.','markersize',1);
%                     hold on;
%                 end
            end
            % PLOT NEURONS THAT SYNAPSE ONTO DECISION CLUSTERS
            % Option 1: Neuron-based wiring
            %
%             if any(strcmp(Ext.stim(n).name,{'Sucrose','Quinine','Maltose','Octaacetate'})) && (~isempty(Ext.stim(n).ind_inc) || ~isempty(Ext.stim(n).ind_dec))
%                 selectiveSNs = Ext.stim(1).ind_inc;
%                 selectiveQNs = Ext.stim(2).ind_inc;
%                 selectiveMNs = Ext.stim(3).ind_inc;
%                 selectiveONs = Ext.stim(4).ind_inc;
%                 selectiveLeftNeurons = [];
%                 selectiveRightNeurons = [];
%                 for neu = 1:1400
%                     if (ismember(neu,selectiveSNs) || ismember(neu,selectiveQNs)) && ~(ismember(neu,selectiveMNs) || ismember(neu,selectiveONs))
%                         selectiveLeftNeurons = [selectiveLeftNeurons neu];
%                     end
%                     if (ismember(neu,selectiveMNs) || ismember(neu,selectiveONs)) && ~(ismember(neu,selectiveSNs) || ismember(neu,selectiveQNs))
%                         selectiveRightNeurons = [selectiveRightNeurons neu];
%                     end
%                 end   
%                 plot(firings(ismember(firings(:,2),selectiveLeftNeurons),1),firings(ismember(firings(:,2),selectiveLeftNeurons),2),'g.','markersize',1); hold on;
%                 plot(firings(ismember(firings(:,2),selectiveRightNeurons),1),firings(ismember(firings(:,2),selectiveRightNeurons),2),'c.','markersize',1); hold on;
%             end
            %
            % Option 2: Cluster-based wiring
            %
            if any(strcmp(Ext.stim(n).name,{'Sucrose','Quinine','Maltose','Octaacetate'})) && (~isempty(Ext.stim(n).ind_inc) || ~isempty(Ext.stim(n).ind_dec))
                if any(strcmp(Ext.stim(n).name,{'Sucrose','Quinine'}))
                    for clu = Ext.stim(n).StimClust_Inc'
                        lower = cusumNcE(clu) + 1;
                        upper = cusumNcE(clu+1);
                        plot(firings(firings(:,2) >= lower & firings(:,2) <= upper,1),firings(firings(:,2) >= lower & firings(:,2) <= upper,2),'g.','markersize',1); hold on;
                    end
                elseif any(strcmp(Ext.stim(n).name,{'Maltose','Octaacetate'}))
                    for clu = Ext.stim(n).StimClust_Inc'
                        lower = cusumNcE(clu) + 1;
                        upper = cusumNcE(clu+1);
                        plot(firings(firings(:,2) >= lower & firings(:,2) <= upper,1),firings(firings(:,2) >= lower & firings(:,2) <= upper,2),'c.','markersize',1); hold on;
                    end
                end
            end
            stimInd = find(arrayfun(@(x) ismember(x.name,stimuliPresent), Ext.stim),1);
            if ~isempty(stimInd)
                if any(strcmp(Ext.stim(stimInd).name,{'Sucrose','Quinine','Maltose','Octaacetate'})) && (~isempty(Ext.stim(stimInd).ind_inc) || ~isempty(Ext.stim(stimInd).ind_dec)) 
                    plot(firings(ismember(firings(:,2),Ext.stim(stimInd).ind_inc),1),firings(ismember(firings(:,2),Ext.stim(stimInd).ind_inc),2)-N_e+Ne_plot,'m.','markersize',1); hold on;
                    plot(firings(ismember(firings(:,2),Ext.stim(stimInd).ind_dec),1),firings(ismember(firings(:,2),Ext.stim(stimInd).ind_dec),2)-N_e+Ne_plot,'b.','markersize',1); hold on;
                end
            end
            %
        end
    end
    line([0 0],[0 Ne_plot+Ni_plot],'color','b');
    hold on;
    % Cluster boundaries
    yline(-99.5,'k'); hold on;
    yline(0.5,'k'); hold on;
    for clust = Ext.stim(1).StimClust_Inc'
        yline(cusumNcE(clust)+0.5,'k'); hold on;
        yline(cusumNcE(clust+1)+0.5,'k'); hold on;
        text(3.05,(cusumNcE(clust)+cusumNcE(clust+1)+1)/2,'S');
    end
    for clust = Ext.stim(2).StimClust_Inc'
        yline(cusumNcE(clust)+0.5,'k'); hold on;
        yline(cusumNcE(clust+1)+0.5,'k'); hold on;
        text(3.05,(cusumNcE(clust)+cusumNcE(clust+1)+1)/2,'Q');
    end
    for clust = Ext.stim(3).StimClust_Inc'
        yline(cusumNcE(clust)+0.5,'k'); hold on;
        yline(cusumNcE(clust+1)+0.5,'k'); hold on;
        text(3.05,(cusumNcE(clust)+cusumNcE(clust+1)+1)/2,'M');
    end
    for clust = Ext.stim(4).StimClust_Inc'
        yline(cusumNcE(clust)+0.5,'k'); hold on;
        yline(cusumNcE(clust+1)+0.5,'k'); hold on;
        text(3.05,(cusumNcE(clust)+cusumNcE(clust+1)+1)/2,'O');
    end
    if isfield(Params,'dClusts')
        text(3.05,(cusumNcE(Params.dClusts(1))+cusumNcE(Params.dClusts(1)+1)+1)/2,'L');
        text(3.05,(cusumNcE(Params.dClusts(2))+cusumNcE(Params.dClusts(2)+1)+1)/2,'R');
    end
    ylim([-2*N_decision Ne_plot+Ni_plot]);
    xlab = 'Time [s]';
    ylab = 'Neurons';
    tt = 'All neurons ordered by population';
    tastes = {'Sucrose','Quinine','Maltose','Octaacetate'};
    tastesPresent = ismember(tastes,stimuliPresent);
    if sum(tastesPresent) > 0
        tt = {tt ; ['Trial ' num2str(trialNum) ' (' tastes{tastesPresent} ')']};
    else
        tt = {tt ; ['Trial ' num2str(trialNum)]};
    end
    aux.figset(gca,xlab,ylab,tt,15);

    xlim(XLIM);
    hold off

end

%---------------------------------------------------

% 2. PSC Plot

function PSCPlot(Sim,Ext,iIplot,iEplot,iExtplot,hicutoff,srate,VTh,ind_plot)
    stimcolors = {'m',[0.5,0.5,0.5], 'b', 'r', 'c'};
    nplot = min(2,size(iIplot,1))+1;
    dt = Sim.dt_step;
    hE = []; hI = []; hext = []; htot = []; hth = [];
    t_plot = max([Sim.t_End-Sim.plot_length,Sim.t_Start])+dt:dt:Sim.t_End;
    ExcInh = {'Exc','Inh'};
    Legplot = {'i_E','i_I','i_{ext}','i_{tot}','V_{th}/\tau'};
    for j = 1:nplot-1
        subplot(nplot,1,j+1);  hstim = [];
        hold on;
        [iISmooth] = aux.eegfilt(iIplot(j,:),srate,0,hicutoff);
        [iESmooth] = aux.eegfilt(iEplot(j,:),srate,0,hicutoff);
        [iExtSmooth] = aux.eegfilt(iExtplot(j,:),srate,0,hicutoff);
        hE(j) = plot(t_plot, iESmooth,'color','b','linewidth',1); % EPSC
        hI(j) = plot(t_plot, iISmooth,'color','r','linewidth',1); % IPSC
        hext(j) = plot(t_plot, iExtplot(j,:),'color','g','linewidth',1); % EPSC
        itot = iEplot(j,:) + iIplot(j,:) + iExtplot(j,:);
%         hth(j) = plot([t_plot(1) t_plot(end)],ones(1,2)*VTh(ind_plot(j))/tau(ind_plot(j)),'color','b','linewidth',1,'linestyle','-.');
        hth(j) = plot([t_plot(1) t_plot(end)],ones(1,2)*VTh(ind_plot(j)),'color','b','linewidth',1,'linestyle','-.');
        [itotSmooth] = aux.eegfilt(itot,srate,0,hicutoff);
        htot(j) = plot(t_plot,itotSmooth,'color','k','linewidth',1); % membrane potentials
        grid on
        if any(strcmp(fieldnames(Ext),'stim'))
            ii = 0;
            for n = 1:numel(Ext.stim)
                if ~ismember(Ext.stim(n).name,stimuliPresent), continue; end
                ii = ii + 1;
                time_bin = [Ext.stim(n).interval(1) Ext.stim(n).interval(2)];
                lower = min(min(iISmooth))*ones(1,numel(time_bin));
                %                     upper=max(max(iESmooth))*ones(1,numel(time_bin));
                %                     [~,~]=aux.jbfill(time_bin,lower,upper,'m',0,0,0.2);
                hstim(ii) = line(time_bin([1,end]),(1+0.1*ii)*lower([1,end]),'color',stimcolors{ii},'linewidth',5);
                Legplot = [Legplot,Ext.stim(n).name];
            end
        end
        xlab = 'Time [s]';
        ylab = 'mV/s';
        tt = sprintf('PSC to %s unit',ExcInh{j});
        fntsz = 15;
        aux.figset(gca,xlab,ylab,tt,fntsz);
        ind_lim = (t_plot>Sim.t_Start+0.1 & t_plot<Sim.t_End-0.1);
        legend([hE(1) hI(1) hext(1) htot(1) hth(1), hstim],Legplot);
        xlim([t_plot(1) t_plot(end)]);
        ylim([1.5*min(min(iISmooth(:,ind_lim)))...
            1.5*max([max(max(iESmooth(:,ind_lim))) max(max(iExtSmooth(:,ind_lim)))])]);
    end

end
%-------------------------------------------------

% 1A. Exc and Inh membrane potentials

function PlotMembrane(Sim,nplot,vplot,firings,ind_plot,VTh)
    ThStyle = {':','--'};
    dt = Sim.dt_step;
    colors = ['b'; 'r'];
    hh = []; cnt = 0;
    bins = min([Sim.t_End-Sim.plot_length,Sim.t_Start])+dt:dt:Sim.t_End;
    for j = 1:nplot
        cnt = cnt + 1;
        hh(cnt) = plot(bins, vplot(j,1:numel(bins)),'color',colors(j,:)); % membrane potentials
        hold on
        if any(any(firings))
            ind = find(firings(:,2)==ind_plot(j));
            if ~isempty(ind) % plot spike overshoots
                h = line([firings(ind,1), firings(ind,1)]',...
                    [VTh(ind_plot(j))*ones(size(ind)) VTh(ind_plot(j))*ones(size(ind))*2]');
                set(h,'color',colors(j,:));
                hold on
            end
        end
        cnt = cnt + 1;
        hh(cnt) = plot([bins(1) bins(end)],[VTh(ind_plot(j)) VTh(ind_plot(j))],'color',colors(j,:),'linestyle',ThStyle{j});
        hold on
    end

    legend(hh,'E','V_{Th,E}','I','V_{Th,I}');
    ylim([1.2*min(min(vplot(:,0.2/dt:end))) 2*max(VTh(ind_plot))]);
    xlim([bins(1) bins(end)]);
    xlab = 'Time [s]';
    ylab = 'mV';
    tt = '';
    fntsz = 15;
    aux.figset(gca,xlab,ylab,tt,fntsz);

end


%-------------------------------------------------------

% 3A. Cluster PSTH

function [rate, x_bins, numfigFlag] = PlotPopRates(bins,firings,N_e,N_i,NcE,N_decision,numfig)
    
    numfigFlag = 0;
    MinRate = 10;
    BinSize = diff(bins(1:2));
    fntsz = 15;
    cols = aux.distinguishable_colors(numel(NcE)+1);
    binwidth = 8;
    cusumNcE = [0 cumsum(NcE)'];
    rate = zeros(numel(NcE),numel(bins)+1);
    for clu = 1:numel(NcE)
        ind_clu = (firings(:,2)>cusumNcE(clu) & firings(:,2)<=cusumNcE(clu+1));
        spikes = firings(ind_clu,1);
        x = histc(spikes,bins);
        [tout, z] = aux.gaussfilt(bins,x,binwidth,0);
        x_bins = tout;
        rate(clu,:) = z/(NcE(clu)*BinSize);
        plot(x_bins-diff(x_bins(1:2))/2,rate(clu,:),'color',cols(clu,:),'linewidth',1.5);
        hold on;
    end
    % unclustered subpop
    if N_e - cusumNcE(end) > 0 % if there's an unclustered population
        ind_unclu = (firings(:,2)>cusumNcE(end) & firings(:,2)<N_e);
        spikes = firings(ind_unclu,1);
        x = histc(spikes,bins);
        [tout, z] = aux.gaussfilt(bins,x,binwidth,0);
        x_bins = tout;
        rate = [rate; z'/((N_e-cusumNcE(end))*BinSize)];
        plot(x_bins-diff(x_bins(1:2))/2,rate(clu+1,:),'color',cols(clu+1,:),'linewidth',1.5,'linestyle','--');
    end
    % inhibitory neurons
    ind_inh = (firings(:,2) > N_e & firings(:,2) <= (N_e + N_i));
    spikes = firings(ind_inh,1);
    x = histc(spikes,bins);
    [tout, z] = aux.gaussfilt(bins,x,binwidth,0);
    x_bins = tout;
    rate = [rate; z'/((N_i)*BinSize)];
    plot(x_bins-diff(x_bins(1:2))/2,rate(clu+2,:),'color','r','linewidth',2,'linestyle',':');

    line([x_bins(1)-diff(x_bins(1:2))/2 x_bins(end)],[MinRate MinRate],'linestyle','--','color','k');
    hold on;
    xlim([bins(1) bins(end)]);
    xlab = 'Time [s]';
    ylab = 'Firing rate [spks/s]';
    tt = '';
    aux.figset(gca,xlab,ylab,tt,fntsz);
    ylim([0 max(max(rate))+10]);
    
    ind_dec = firings(:,2) > (N_e + N_i) & firings(:,2) <= (N_e + N_i + 2*N_decision);
    if sum(ind_dec) ~= 0
        figure(numfig+1); clf; numfigFlag = 1;
        ind_L = firings(:,2) > (N_e + N_i) & firings(:,2) <= (N_e + N_i + N_decision);
        spikes_L = firings(ind_L,1);
        x = histc(spikes_L,bins);
        [tout, z] = aux.gaussfilt(bins,x,binwidth,0);
        x_bins = tout;
        rate = [rate; z'/((N_decision)*BinSize)];
        plot(x_bins-diff(x_bins(1:2))/2,rate(clu+3,:),'color','g','linewidth',2,'linestyle',':');
        hold on;
        ind_R = firings(:,2) > (N_e + N_i + N_decision) & firings(:,2) <= (N_e + N_i + 2*N_decision);
        spikes_R = firings(ind_R,1);
        x = histc(spikes_R,bins);
        [tout, z] = aux.gaussfilt(bins,x,binwidth,0);
        x_bins = tout;
        rate = [rate; z'/((N_decision)*BinSize)];
        plot(x_bins-diff(x_bins(1:2))/2,rate(clu+4,:),'color','c','linewidth',2,'linestyle',':');
        xlim([bins(1) bins(end)]);
        xlab = 'Time [s]';
        ylab = 'Firing rate [Hz]';
        tt = '';
        legend({'Left Decision Cluster','Right Decision Cluster'});
        aux.figset(gca,xlab,ylab,tt,fntsz);
        figure(numfig); hold on;
    end
end


function PlotPopRaster(rate,x_bins,bins)


    Cmin = min(min(rate));
    Cmax = max(max(rate));
    imagesc(x_bins,1:size(rate,1),rate); axis xy;
    caxis([Cmin, Cmax]);
    % colormap gray;
    % colormap(1-colormap);
    xlim([bins(1) bins(end)]);
    t = colorbar; get(t,'ylabel');
    set(get(t,'ylabel'),'String', 'Firing rate [spks/s]');
    hold off
    xlab = 'Time [s]';
    ylab = 'Population index';
    tt = '';
    fntsz = 15;
    aux.figset(gca,xlab,ylab,tt,fntsz);

end

end