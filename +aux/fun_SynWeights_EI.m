% function [J, params]=SynWeights(params)
%
% OUTPUT
%       J                =synaptic matrix
%       OUT.popsize      =array of dim # pops, how many exc neurons in each
%                         pop
%       OUT.clustermatrix=matrix of dimension # pops x # clusters, each row shows
%                   which clusters belong to that pop (1s and 0s)
%
% Luca Mazzucato February 2016


function [J, params, J_decision] = fun_SynWeights_EI(paramsfile,tempfile)

% rng(99995);

% LOAD PARAMETERS

isChangeJDOnly = false;
isUniform = false;
isRegularizeInputNums = false;
isOptimizeSelectivity = true;

if ~isChangeJDOnly
    
params = load(paramsfile);
aux.v2struct(params);
Network = params.Network;
aux.v2struct(Network);
Sim = params.Sim;

% if isUniform
%     rng(99996);
% end

%-----------------------
% PARAMETERS VALUES
%-----------------------
% numfig = 1;
Next = N_e; % external units
% CLUSTERS
Q = p; % number of clusters
%-----------------------
% SYNAPTIC WEIGHTS
%-----------------------
% WEIGHTS
% depression
gam = 1/(2-f*(Q+1));%
Jminus = 1.-gam*f*(Jplus-1.);
params.Jminus = Jminus;

% INHIBITORY CLUSTERS
if any(strcmp(clusters,'EI'))
    %     JminusEI = 1.-gam*fI*(JplusEI-1.);
    JplusEI = 1/(1/p+(1-1/p)/factorEI);
    JminusEI = JplusEI/factorEI;
    params.JminusEI = JminusEI;
    jei_out = -JminusEI*Jei; % inter-cluster
    jei_in = -JplusEI*Jei; % intra-cluster
    fprintf('JplusEI = %0.03g; JminusEI = %0.03g\n',JplusEI,JminusEI);
end
if any(strcmp(clusters,'IE'))
    %     JminusIE = 1.-gam*fI*(JplusIE-1.);
    JplusIE = 1/(1/p+(1-1/p)/factorIE);
    JminusIE = JplusIE/factorIE;
    params.JminusIE = JminusIE;
    jie_out = JminusIE*Jie; % inter-cluster
    jie_in = JplusIE*Jie; % intra-cluster
    fprintf('JplusIE = %0.03g; JminusIE = %0.03g\n',JplusIE,JminusIE);
end
if any(strcmp(clusters,'II'))
    JplusII = factorII;
    JminusII = 1.-gam*fI*(JplusII-1.);
    params.JminusII = JminusII;
    jii_out = -JminusII*Jii; % inter-cluster
    jii_in = -JplusII*Jii; % intra-cluster
    fprintf('JplusII = %0.03g; JminusII = %0.03g\n',JplusII,JminusII);
end

%
jee = Jee;
jee_out = Jminus*Jee; % inter-cluster potentiation (strength of connection to neurons outside cluster)
jee_in = Jplus*Jee; % intra-cluster potentiation (strength of connection to neurons in same cluster)
jei = -Jei;
jie = Jie;
jii = -Jii;
% connection probability
pee = Cee/N_e;
pie = Cie/N_e;
pei = Cei/N_i;
pii = Cii/N_i;
pext = Cext/Next;
peeout = pee;
peein = pee;
peiout = pei;
peiin = pei;
pieout = pie;
piein = pie;
piiout = pii;
piiin = pii;
%
fprintf('  --- Jplus = %0.03g, Jminus = %0.03g\n',Jplus,Jminus);

%----------------------------
% SYNAPTIC MATRIX
%----------------------------
% if delta>0, generate a distribution of synaptic weights with mean J and
% variance delta^2 J^2
% peeout = pee;
% peein = pee;

%-------------------------------------------
% CLUSTERS INVOLVING ONLY EXCITATORY UNITS
%-------------------------------------------
% check #clusters and coding level are consistent
NcUnits = round(f*N_e);    % number of Exc units per cluster
Numbg = round(N_e*(1-f*p)); % number of background (i.e. non-selective) Exc units
switch Network.clust
    case 'hom'
        popsize = repmat(NcUnits,Q,1);
    case 'het'
        Nc = [];
        clust_std = Network.clust_std;
        while (sum(Nc)-(N_e-Numbg))~=0 || any(Nc<0)
            Nc = round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
        end
        popsize = Nc; % array of cluster sizes
        if any(sum(popsize)-(N_e-Numbg))
            fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
        end
end
cusumNcE = [0 cumsum(popsize)'];

%----------------------
% EE CLUSTERS (E --> E)
%----------------------
%----------------------
% EE weights (E --> E)
%----------------------
% background units (if present), if not, override in next line
JEE = (jee*(ones(N_e)+delta*randn(N_e,N_e))).*(rand([N_e,N_e])<peeout);
JEI = (jei*(ones(N_e,N_i)+deltaEI*randn(N_e,N_i))).*(rand([N_e,N_i])<pei);
JIE = (jie*(ones(N_i,N_e)+deltaIE*randn(N_i,N_e))).*(rand([N_i,N_e])<pie);
JII = (jii*(ones(N_i)+delta*randn(N_i,N_i))).*(rand([N_i,N_i])<pii);
%clustermatrix = eye(Q);
if strcmp(Network.clust,'het') || strcmp(Network.clust,'hom')
    if isUniform
        sparseEEinter = randn([popsize(1),popsize(1)])<peeout;
        sparseEEintra = randn([popsize(1),popsize(1)])<peein;
        for clu1 = 1:14
            for clu2 = 1:14
                if clu1 == clu2
                    JEE(cusumNcE(clu1)+1:cusumNcE(clu1+1),cusumNcE(clu2)+1:cusumNcE(clu2+1)) = ...
                        jee_in*ones(popsize(clu1),popsize(clu2)).*sparseEEintra;
                else
                    JEE(cusumNcE(clu1)+1:cusumNcE(clu1+1),cusumNcE(clu2)+1:cusumNcE(clu2+1)) = ...
                        jee_out*ones(popsize(clu1),popsize(clu2)).*sparseEEinter;
                end
            end
        end
    else
        % clustered units: inter-cluster weights
        JEE(1:cusumNcE(Q+1),1:cusumNcE(Q+1)) = ...
            (jee_out*(ones(cusumNcE(Q+1))+delta*randn(cusumNcE(Q+1), ...
            cusumNcE(Q+1)))).*(rand([cusumNcE(Q+1),cusumNcE(Q+1)])<peeout); % inter-cluster weights
        % intra-cluster higher weights
        for clu = 2:Q+1 
            JEE(1+cusumNcE(clu-1):cusumNcE(clu),1+cusumNcE(clu-1):cusumNcE(clu)) = ...
                (jee_in*(ones(popsize(clu-1))+delta*randn(popsize(clu-1),popsize(clu-1)))).* ...
                (rand([popsize(clu-1),popsize(clu-1)])<peein);
        end
    end
end
clustermatrix = eye(Q);
params.clustermatrix = clustermatrix;
params.popsize = popsize;
params.cusumNcE = cusumNcE;

%-------------------------------------
% CLUSTERS INVOLVING INHIBITORY UNITS
%-------------------------------------
% check # clusters and coding level are consistent
if any([strcmp(clusters,'EI'), strcmp(clusters,'EI'), strcmp(clusters,'II')])
    NcUnits = round(fI*N_i);    %  number of inh units per cluster
    fprintf('  --- Synaptic weights: %d units/cluster \n',NcUnits);
    Numbg = round(N_i*(1-fI*p)); % number of background (i.e. non-selective) inh units
    fprintf('  --- fraction of bg Inh units: %0.03g',Numbg/N_i);
    switch Network.clustEI
        case 'hom'
            popsizeI = repmat(NcUnits,Q,1);
        case 'het'
            Nc = [];
            clust_std = Network.clust_std;
            while (sum(Nc)-(N_i-Numbg))~=0 || any(Nc<0)
                Nc = round(NcUnits+(NcUnits*clust_std)*randn(Q,1));
            end
            popsizeI = Nc; % array of cluster sizes
            if any(sum(popsizeI)-(N_i-Numbg))
                fprintf('\n---ERROR: Heterogeneous clusters: Problem with cluster numbers\n');
            end
    end
    cusumNcI = [0 cumsum(popsizeI)'];
    
    %----------------------
    % EI CLUSTERS (I --| E)
    %----------------------
    %----------------------
    % EI weights (I --| E)
    %----------------------
    if any(strcmp(Network.clusters,'EI'))
        % background units (if present), if not, override in next line
        if strcmp(Network.clustEI,'het') || strcmp(Network.clustEI,'hom')
            % clustered units: inter-cluster weights
            if isUniform
                sparseEIinter = randn([popsize(1),popsizeI(1)])<peiout;
                sparseEIintra = randn([popsize(1),popsizeI(1)])<peiin;
                for clu1 = 1:14
                    for clu2 = 1:14
                        if clu1 == clu2
                            JEI(cusumNcE(clu1)+1:cusumNcE(clu1+1),cusumNcI(clu2)+1:cusumNcI(clu2+1)) = ...
                                jei_in*ones(popsize(clu1),popsizeI(clu2)).*sparseEIintra;
                        else
                            JEI(cusumNcE(clu1)+1:cusumNcE(clu1+1),cusumNcI(clu2)+1:cusumNcI(clu2+1)) = ...
                                jei_out*ones(popsize(clu1),popsizeI(clu2)).*sparseEIinter;
                        end
                    end
                end
            else
                JEI(1:cusumNcE(Q+1),1:cusumNcI(Q+1)) = ...
                    (jei_out*(ones(cusumNcE(Q+1),cusumNcI(Q+1))+deltaEI*randn(cusumNcE(Q+1), ...
                    cusumNcI(Q+1)))).*(rand([cusumNcE(Q+1),cusumNcI(Q+1)])<peiout); % inter-cluster weights
                % intra-cluster higher weights
                for clu = 2:Q+1 
                    JEI(1+cusumNcE(clu-1):cusumNcE(clu),1+cusumNcI(clu-1):cusumNcI(clu)) = ...
                        (jei_in*(ones(popsize(clu-1),popsizeI(clu-1))+deltaEI*randn(popsize(clu-1),popsizeI(clu-1)))).*...
                        (rand([popsize(clu-1),popsizeI(clu-1)])<peiin);
                end
            end
        end
    end
    
    %----------------------
    % IE CLUSTERS (E --> I)
    %----------------------
    %----------------------
    % IE weights (E --> I)
    %----------------------
    if any(strcmp(Network.clusters,'IE'))
        % background units (if present), if not, override in next line
        if strcmp(Network.clustIE,'het') || strcmp(Network.clustIE,'hom')
            if isUniform
                sparseIEinter = randn([popsizeI(1),popsize(1)])<pieout;
                sparseIEintra = randn([popsizeI(1),popsize(1)])<piein;
                for clu1 = 1:14
                    for clu2 = 1:14
                        if clu1 == clu2
                            JIE(cusumNcI(clu1)+1:cusumNcI(clu1+1),cusumNcE(clu2)+1:cusumNcE(clu2+1)) = ...
                                jie_in*ones(popsizeI(clu1),popsize(clu2)).*sparseIEintra;
                        else
                            JIE(cusumNcI(clu1)+1:cusumNcI(clu1+1),cusumNcE(clu2)+1:cusumNcE(clu2+1)) = ...
                                jie_out*ones(popsizeI(clu1),popsize(clu2)).*sparseIEinter;
                        end
                    end
                end
            else
                % clustered units: inter-cluster weights
                JIE(1:cusumNcI(Q+1),1:cusumNcE(Q+1)) = ...
                    (jie_out*(ones(cusumNcI(Q+1),cusumNcE(Q+1))+deltaIE*randn(cusumNcI(Q+1), ...
                    cusumNcE(Q+1)))).*(rand([cusumNcI(Q+1),cusumNcE(Q+1)])<pieout); % inter-cluster weights
                % intra-cluster higher weights
                for clu = 2:Q+1 
                    JIE(1+cusumNcI(clu-1):cusumNcI(clu),1+cusumNcE(clu-1):cusumNcE(clu)) = ...
                        (jie_in*(ones(popsizeI(clu-1),popsize(clu-1))+deltaIE*randn(popsizeI(clu-1),popsize(clu-1)))).*...
                        (rand([popsizeI(clu-1),popsize(clu-1)])<piein);
                end
            end
        end
    end
    
    %----------------------
    % II CLUSTERS (I --| I)
    %----------------------
    %----------------------
    % II weights (I --| I)
    %----------------------
    if any(strcmp(Network.clusters,'II'))
        % background units (if present), if not, override in next line
        if strcmp(Network.clustII,'het') || strcmp(Network.clustII,'hom')
            if isUniform
                sparseIIinter = randn([popsizeI(1),popsizeI(1)])<piiout;
                sparseIIintra = randn([popsizeI(1),popsizeI(1)])<piiin;
                for clu1 = 1:14
                    for clu2 = 1:14
                        if clu1 == clu2
                            JII(cusumNcI(clu1)+1:cusumNcI(clu1+1),cusumNcI(clu2)+1:cusumNcI(clu2+1)) = ...
                                jii_in*ones(popsizeI(clu1),popsizeI(clu2)).*sparseIIintra;
                        else
                            JII(cusumNcI(clu1)+1:cusumNcI(clu1+1),cusumNcI(clu2)+1:cusumNcI(clu2+1)) = ...
                                jii_out*ones(popsizeI(clu1),popsizeI(clu2)).*sparseIIinter;
                        end
                    end
                end
            else
                % clustered units: inter-cluster weights
                JII(1:cusumNcI(Q+1),1:cusumNcI(Q+1)) = ...
                    (jii_out*(ones(cusumNcI(Q+1),cusumNcI(Q+1))+deltaII*randn(cusumNcI(Q+1), ...
                    cusumNcI(Q+1)))).*(rand([cusumNcI(Q+1),cusumNcI(Q+1)])<piiout); % inter-cluster weights
                % intra-cluster higher weights
                for clu = 2:Q+1 
                    JII(1+cusumNcI(clu-1):cusumNcI(clu),1+cusumNcI(clu-1):cusumNcI(clu)) = ...
                        (jii_in*(ones(popsizeI(clu-1),popsizeI(clu-1))+deltaII*randn(popsizeI(clu-1),popsizeI(clu-1)))).*...
                        (rand([popsizeI(clu-1),popsizeI(clu-1)])<piiin);
                end
            end
        end
    end
    params.popsizeI = popsizeI;
    params.cusumNcI = cusumNcI;
end

JEI(JEI>0) = 0;
JIE(JIE<0) = 0;
JII(JII>0) = 0;
JEE(JEE<0) = 0;
J = [JEE JEI; JIE JII];
J = J-diag(diag(J)); % eliminate self-couplings
fprintf('  --- New synaptic weights set: done...\n');
fprintf('      Overall: Jee = %g; Jie = %g; Jei = %g; Jii = %g \n',jee,jie,jei,jii);
fprintf('      Var[J] = (Jx%0.03g)^2\n',delta);

if isRegularizeInputNums
    J = aux.regularizeInputNums(J,N_e,N_i,cusumNcE,cusumNcI);
    fprintf('  --- Inputs in J have been regularized.\n'); 
end

if strcmp(Network.clust,'het')
    for i = 1:Q
        a(i) = sum(popsize(clustermatrix(:,i)>0));
    end
    fprintf('  --- Clusters size: mean = %0.03g neurons/cluster, sd/mean = %0.03g\n',mean(a),std(a)/mean(a));
end
fprintf('  --- Fraction of bg Exc units: %0.03g\n',N_e/Numbg);

% spike thresholds for each population
theta = [params.theta_e, params.theta_i];
Theta = zeros(1,numel(params.popsize)+2);
Theta(1:numel(params.popsize)+1) = theta(1);
Theta(end) = theta(2);
params.Theta = Theta;
fprintf('Spike thresholds Theta calculated for each population and stored in params\n');

%%
%-----------------------
% DECISION CLUSTERS
%-----------------------

% New option: run a baseline test and assign stimulus-selective clusters in
% a way that minimizes their variance in firing rate
if isOptimizeSelectivity
    fprintf('Running initial simulation to obtain baseline firing rates...\n');
    ParamsRun = params;
    ParamsRun.TrialStimuli = {};
    ParamsRun.stimuli = {};
    ParamsRun.Ext.Mu = ParamsRun.Mu;
    ParamsRun.J = J;
    ParamsRun.J_decision = zeros(2*ParamsRun.N_decision,ParamsRun.N_e+ParamsRun.N_i+2*ParamsRun.N_decision);
    [firings, ~] = aux.fun_LIF_SIM(ParamsRun);
    params.prelimFirings = firings;
    fprintf('Initial simulation complete. Optimizing cluster selectivity...\n');
    
    % Calculate E cluster firing rates 
    clusterFRs = NaN(Q,1); 
    for clust = 1:Q    
        ids = (cusumNcE(clust)+1):cusumNcE(clust+1); 
        clusterFRs(clust) = sum(ismember(firings(:,2),ids))/length(ids)/(ParamsRun.Sim.t_End - ParamsRun.Sim.t_Start);
    end 
    params.clusterFRs = clusterFRs; 
    
    %
    figure(998); clf;
    plot(firings(firings(:,2)<=ParamsRun.N_e,1),firings(firings(:,2)<=ParamsRun.N_e,2),'k.','markersize',1); hold on;
    plot(firings(firings(:,2)>ParamsRun.N_e & firings(:,2)<=ParamsRun.N_e+ParamsRun.N_i,1),firings(firings(:,2)>ParamsRun.N_e & firings(:,2)<=ParamsRun.N_e+ParamsRun.N_i,2),'r.','markersize',1); hold on;
    for clust = 1:ParamsRun.p
        % Mark E cluster boundaries
        yline(ParamsRun.cusumNcE(clust)+0.5,'k'); hold on;
        yline(ParamsRun.cusumNcE(clust+1)+0.5,'k'); hold on;
        text(3.05,(ParamsRun.cusumNcE(clust)+ParamsRun.cusumNcE(clust+1)+1)/2,['E' num2str(clust)]);
    end
    text(3.05,(ParamsRun.cusumNcE(end)+ParamsRun.N_e+1)/2,'EB');
    for clust = 1:ParamsRun.p
        % Mark I cluster boundaries
        yline(ParamsRun.N_e+ParamsRun.cusumNcI(clust)+0.5,'k'); hold on;
        yline(ParamsRun.N_e+ParamsRun.cusumNcI(clust+1)+0.5,'k'); hold on;
        text(3.05,(2*ParamsRun.N_e+ParamsRun.cusumNcI(clust)+ParamsRun.cusumNcI(clust+1)+1)/2,['I' num2str(clust)]);
    end
    text(3.05,(ParamsRun.N_e+ParamsRun.cusumNcI(end)+ParamsRun.N_e+ParamsRun.N_i+1)/2,'IB');
    yline(0.5,'b'); hold on;
    yline(ParamsRun.cusumNcE(end)+0.5,'b');
    yline(ParamsRun.N_e+0.5,'b');
    yline(ParamsRun.N_e+ParamsRun.cusumNcI(end)+0.5,'b');
    yline(ParamsRun.N_e+ParamsRun.N_i+0.5,'b');
    ylim([0, ParamsRun.N_e+ParamsRun.N_i]); xlim([ParamsRun.Sim.t_Start, ParamsRun.Sim.t_End]);
    xlab = 'Time [s]';
    ylab = 'Neurons';
    tt = 'Baseline Simulation';
    aux.figset(gca,xlab,ylab,tt,15);
    ax = gca; ax.set('TickDir','out');
    text(-1.65,2000,['Jee: ' num2str(ParamsRun.Jee/ParamsRun.Scale)]);
    text(-1.65,1950,['Jii: ' num2str(ParamsRun.Jii/ParamsRun.Scale)]);
    text(-1.65,1900,['Jie: ' num2str(ParamsRun.Jie/ParamsRun.Scale)]);
    text(-1.65,1850,['Jei: ' num2str(ParamsRun.Jei/ParamsRun.Scale)]);
    text(-1.65,1800,['Jplus: ' num2str(ParamsRun.Jplus)]);
    text(-1.65,1750,['factorII: ' num2str(ParamsRun.Network.factorII)]);
    text(-1.65,1700,['factorIE: ' num2str(ParamsRun.Network.factorIE)]);
    text(-1.65,1650,['factorEI: ' num2str(ParamsRun.Network.factorEI)]);
    text(-1.65,1600,['ni ext: ' num2str(ni_ext)]);
    text(-1.65,1550,['regInputs?: ' num2str(isRegularizeInputNums)]);
%     f = gcf; f.WindowState = 'maximized';
%     HOME = pwd; addpath(genpath(HOME)); cd('New Plots'); cd('Network'); cd('sim'); cd('parameterTuning');
%     d = dir; name = 1; 
%     if numel({d.name}) > 2
%         myCell = {d.name};
%         myCell = myCell(cellfun(@(x) contains(x,'.jpg'), myCell));
%         name = max(cellfun(@(x) str2num(x(1:end-4)), myCell)) + 1;
%     end
%     saveas(f,[num2str(name) '.jpg']); cd(HOME);
%     f.WindowState = 'normal';
%     figure(999); histogram(clusterFRs);
%     disp(clusterFRs);
%     fprintf(['Baseline firing rates: ' num2str(mean(clusterFRs)) ' +/- ' num2str(std(clusterFRs)) '\n']);
%     fprintf(['Baseline firing rate range: ' num2str(min(clusterFRs)) ' - ' num2str(max(clusterFRs)) '\n']);
%     pause;
    %
    lClust = 13;
    rClust = 14;
%     [~,lClust] = min(abs(clusterFRs - mean(clusterFRs)));
%     temp = abs(clusterFRs - clusterFRs(lClust)); temp(lClust) = inf;
%     [~,rClust] = min(temp);
%     dClusts = [lClust rClust]; params.dClusts = dClusts;
%     
%     results = aux.optimizeSelectivity(clusterFRs,dClusts); 
%     params.selOptRes = results;
%     %[~,ind] = min([results{9,:}].^2.*[results{10,:}]);
%     [~,ind] = min([results{9,:}]); % assign clusters to minimize L vs. R baseline variance ==> average 80% performance
%     %[~,ind] = min([results{7,:}].*[results{8,:}]); % assign clusters to minimize overall variance of means and of variances
%     params.Stimulus.feat(2).selective = ismember((1:Q)',results{1,ind}); % sucrose clusters
%     params.Stimulus.feat(3).selective = ismember((1:Q)',results{2,ind}); % quinine clusters
%     params.Stimulus.feat(4).selective = ismember((1:Q)',results{3,ind}); % maltose clusters
%     params.Stimulus.feat(5).selective = ismember((1:Q)',results{4,ind}); % octaacetate clusters
%     fprintf('Optimal clusters selected\n');
    
    % Increase strength from S,Q --> L and M,O --> R within-network clusters 
    for ll = [find(params.Stimulus.feat(2).selective) ; find(params.Stimulus.feat(3).selective)]'
        J(cusumNcE(lClust)+1:cusumNcE(lClust+1),cusumNcE(ll)+1:cusumNcE(ll+1)) = ...
            1.025*J(cusumNcE(lClust)+1:cusumNcE(lClust+1),cusumNcE(ll)+1:cusumNcE(ll+1));
    end
    for rr = [find(params.Stimulus.feat(4).selective) ; find(params.Stimulus.feat(5).selective)]'
        J(cusumNcE(rClust)+1:cusumNcE(rClust+1),cusumNcE(rr)+1:cusumNcE(rr+1)) = ...
            1.025*J(cusumNcE(rClust)+1:cusumNcE(rClust+1),cusumNcE(rr)+1:cusumNcE(rr+1));
    end

end
    % row 5: allGroupErr; 6: LVRErr; 7: varOfMeans; 8: varOfVars; 9: varOfMeans_LR; 10: varOfVars_LR 
%     powers = [1 0 0 0;
%               0 1 0 0;
%               0 0 1 0;
%               0 0 0 1;
%               1 1 1 1;
%               1 1 0 0;
%               0 0 1 1;
%               2 1 2 1;
%               1 2 1 2;
%               1 1 2 2];
%     paramCell = cell(size(powers,1),1);
%     J_decisionCell = cell(size(powers,1),1);
%     for row = 1:size(powers,1)
%         [~,ind] = min([results{7,:}].^powers(row,1).*[results{8,:}].^powers(row,2).*[results{9,:}].^powers(row,3).*[results{10,:}].^powers(row,4));
%         params.Stimulus.feat(2).selective = ismember((1:Q)',results{1,ind}); % sucrose clusters
%         params.Stimulus.feat(3).selective = ismember((1:Q)',results{2,ind}); % quinine clusters
%         params.Stimulus.feat(4).selective = ismember((1:Q)',results{3,ind}); % maltose clusters
%         params.Stimulus.feat(5).selective = ismember((1:Q)',results{4,ind}); % octaacetate clusters
%         %fprintf('Optimal clusters selected\n');
%         J_decision = aux.getJ_decision(params,N_e,N_i,N_decision,Q,jBackground_e,jBackground_i,pBackgroundE,pBackgroundI,...
%         cusumNcE,cusumNcI,jDefault_e,jDefault_i,pDefaultE,pDefaultI,jHigh_e,pPrefE,popsize,jdd_in,jdd_out,pSelf,pCross);
%         paramCell{row} = params;
%         J_decisionCell{row} = J_decision;
%     end
    
else
    
fprintf('\nUsing parameters and J from \n    %s \nwith J_decision parameters specified by \n    %s\n\n',tempfile,paramsfile);   
PARAMS = load(paramsfile); % constructed with this simulation
cd('data'); template = load(tempfile); cd('../'); % old results
J = template.J; % use old J
% use old params, but overwrite J_decision parameters
params = template.params; 
params.N_decision = PARAMS.N_decision;
params.jBackground_e = PARAMS.jBackground_e; params.jBackground_i = PARAMS.jBackground_i;
params.pBackgroundE = PARAMS.pBackgroundE; params.pBackgroundI = PARAMS.pBackgroundI;
params.jDefault_e = PARAMS.jDefault_e; params.jDefault_i = PARAMS.jDefault_i;
params.pDefaultE = PARAMS.pDefaultE; params.pDefaultI = PARAMS.pDefaultI;
params.jHigh_e = PARAMS.jHigh_e; params.pPrefE = PARAMS.pPrefE;
params.jHigh_i = PARAMS.jHigh_i; params.pPrefI = PARAMS.pPrefI;
params.jdd_in = PARAMS.jdd_in; params.jdd_out = PARAMS.jdd_out;
params.pSelf = PARAMS.pSelf; params.pCross = PARAMS.pCross;
params.tau_decision = PARAMS.tau_decision; params.tausyn_decision = PARAMS.tausyn_decision;
% unpack all necessary parameters
N_e = params.N_e; N_i = params.N_i; N_decision = params.N_decision;
Q = params.p;
jBackground_e = params.jBackground_e; jBackground_i = params.jBackground_i;
pBackgroundE = params.pBackgroundE; pBackgroundI = params.pBackgroundI;
cusumNcE = params.cusumNcE; cusumNcI = params.cusumNcI;
jDefault_e = params.jDefault_e; jDefault_i = params.jDefault_i;
pDefaultE = params.pDefaultE; pDefaultI = params.pDefaultI;
jHigh_e = params.jHigh_e; pPrefE = params.pPrefE;
popsize = params.popsize;
jdd_in = params.jdd_in; jdd_out = params.jdd_out;
pSelf = params.pSelf; pCross = params.pCross;

end

J_decision = aux.getJ_decision(params,N_e,N_i,N_decision,Q,jBackground_e,jBackground_i,pBackgroundE,pBackgroundI,...
    cusumNcE,cusumNcI,jDefault_e,jDefault_i,pDefaultE,pDefaultI,jHigh_e,pPrefE,popsize,jdd_in,jdd_out,pSelf,pCross);

% % Background input
% J_D = ...
%     [jBackground_e*ones(2*N_decision,N_e).*(rand(2*N_decision,N_e) < pBackgroundE), ...
%      jBackground_i*ones(2*N_decision,N_i).*(rand(2*N_decision,N_i) < pBackgroundI)];
% % Default clustered input
% J_D(1:2*N_decision, 1:cusumNcE(Q+1)) = ...
%     jDefault_e*ones(2*N_decision, cusumNcE(Q+1)).*(rand(2*N_decision, cusumNcE(Q+1)) < pDefaultE);
% J_D(1:2*N_decision, (N_e+1):(N_e+cusumNcI(Q+1))) = ...
%     jDefault_i*ones(2*N_decision, cusumNcI(Q+1)).*(rand(2*N_decision, cusumNcI(Q+1)) < pDefaultI);
% % OPTION 1: Preferential input wired by neuron
% %
% % selectiveSNs = Stimulus.feat(2).ind_inc;
% % selectiveQNs = Stimulus.feat(3).ind_inc;
% % selectiveMNs = Stimulus.feat(4).ind_inc;
% % selectiveONs = Stimulus.feat(5).ind_inc;
% % selectiveLeftNeurons = [];
% % selectiveRightNeurons = [];
% % for i = 1:1400
% %     if (ismember(i,selectiveSNs) || ismember(i,selectiveQNs)) && ~(ismember(i,selectiveMNs) || ismember(i,selectiveONs))
% %         selectiveLeftNeurons = [selectiveLeftNeurons i];
% %     end
% %     if (ismember(i,selectiveMNs) || ismember(i,selectiveONs)) && ~(ismember(i,selectiveSNs) || ismember(i,selectiveQNs))
% %         selectiveRightNeurons = [selectiveRightNeurons i];
% %     end
% % end
% % J_D(1:N_decision,selectiveLeftNeurons) = ...
% %     2.2*jHigh_e*ones(N_decision,length(selectiveLeftNeurons)).*(rand(N_decision,length(selectiveLeftNeurons)) < pPrefE);
% % J_D((N_decision+1):2*N_decision,selectiveRightNeurons) = ...
% %     jHigh_e*ones(N_decision,length(selectiveRightNeurons)).*(rand(N_decision,length(selectiveRightNeurons)) < pPrefE);
% %
% % OPTION 2: Preferential input wired by cluster
% %
% LeftClusters = [find(params.Stimulus.feat(2).selective') find(params.Stimulus.feat(3).selective')];
% for i = LeftClusters
%     J_D(1:N_decision,1+cusumNcE(i):cusumNcE(i+1)) = ...
%         0.20*jHigh_e*ones(N_decision,popsize(i)).*(rand(N_decision,popsize(i)) < pPrefE);
% end
% RightClusters = [find(params.Stimulus.feat(4).selective') find(params.Stimulus.feat(5).selective')];
% for i = RightClusters
%     J_D((N_decision+1):2*N_decision,1+cusumNcE(i):cusumNcE(i+1)) = ...
%         0.20*jHigh_e*ones(N_decision,popsize(i)).*(rand(N_decision,popsize(i)) < pPrefE);
% end
% %
% % OLD OPTIONS
% %
% % % Preferential excitatory clustered input
% % % S --> L
% % selectiveClusters = find(Stimulus.feat(2).selective == 1)';
% % for clu = selectiveClusters
% %     J_D(1:N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % % Q --> L
% % selectiveClusters = find(Stimulus.feat(3).selective == 1)';
% % for clu = selectiveClusters
% %     J_D(1:N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % % M --> R
% % selectiveClusters = find(Stimulus.feat(4).selective == 1)';
% % for clu = selectiveClusters
% %     J_D((N_decision+1):2*N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % % O --> R
% % selectiveClusters = find(Stimulus.feat(5).selective == 1)';
% % for clu = selectiveClusters
% %     J_D((N_decision+1):2*N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % %----------------------------
% % % Suc: E,I --> L; E,I --> R
% % %----------------------------
% % % Background input
% % JED_S = ...
% %     [jBackground_e*ones(2*N_decision,N_e).*(rand(2*N_decision,N_e) < pBackgroundE), ...
% %      jBackground_i*ones(2*N_decision,N_i).*(rand(2*N_decision,N_i) < pBackgroundI)];
% % % Default clustered input
% % JED_S(1:2*N_decision, 1:cusumNcE(Q+1)) = ...
% %     jDefault_e*ones(2*N_decision, cusumNcE(Q+1)).*(rand(2*N_decision, cusumNcE(Q+1)) < pDefaultE);
% % JED_S(1:2*N_decision, (N_e+1):(N_e+cusumNcI(Q+1))) = ...
% %     jDefault_i*ones(2*N_decision, cusumNcI(Q+1)).*(rand(2*N_decision, cusumNcI(Q+1)) < pDefaultI);
% % % Preferential excitatory clustered input
% % selectiveClusters = find(Stimulus.feat(2).selective == 1)';
% % for clu = selectiveClusters
% %     JED_S(1:N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % selectiveClusters = find(Stimulus.feat(2).selective == -1)';
% % for clu = selectiveClusters
% %     JED_S(1:N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_i*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefI);
% % end
% % %-------------------------
% % % Qui: E,I --> L; E,I --> R
% % %-------------------------
% % % Background input
% % JED_Q = ...
% %     [jBackground_e*ones(2*N_decision,N_e).*(rand(2*N_decision,N_e) < pBackgroundE), ...
% %      jBackground_i*ones(2*N_decision,N_i).*(rand(2*N_decision,N_i) < pBackgroundI)];
% % % Default clustered input
% % JED_Q(1:2*N_decision, 1:cusumNcE(Q+1)) = ...
% %     jDefault_e*ones(2*N_decision, cusumNcE(Q+1)).*(rand(2*N_decision, cusumNcE(Q+1)) < pDefaultE);
% % JED_Q(1:2*N_decision, (N_e+1):(N_e+cusumNcI(Q+1))) = ...
% %     jDefault_i*ones(2*N_decision, cusumNcI(Q+1)).*(rand(2*N_decision, cusumNcI(Q+1)) < pDefaultI);
% % % Preferential excitatory clustered input
% % selectiveClusters = find(Stimulus.feat(3).selective == 1)';
% % for clu = selectiveClusters
% %     JED_Q(1:N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % selectiveClusters = find(Stimulus.feat(3).selective == -1)';
% % for clu = selectiveClusters
% %     JED_Q(1:N_decision, 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_i*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefI);
% % end
% % %---------------------------
% % % Mal: E,I --> L; E,I --> R
% % %---------------------------
% % % Background input
% % JED_M = ...
% %     [jBackground_e*ones(2*N_decision,N_e).*(rand(2*N_decision,N_e) < pBackgroundE), ...
% %      jBackground_i*ones(2*N_decision,N_i).*(rand(2*N_decision,N_i) < pBackgroundI)];
% % % Default clustered input
% % JED_M(1:2*N_decision, 1:cusumNcE(Q+1)) = ...
% %     jDefault_e*ones(2*N_decision, cusumNcE(Q+1)).*(rand(2*N_decision, cusumNcE(Q+1)) < pDefaultE);
% % JED_M(1:2*N_decision, (N_e+1):(N_e+cusumNcI(Q+1))) = ...
% %     jDefault_i*ones(2*N_decision, cusumNcI(Q+1)).*(rand(2*N_decision, cusumNcI(Q+1)) < pDefaultI);
% % % Preferential excitatory clustered input
% % selectiveClusters = find(Stimulus.feat(4).selective == 1)';
% % for clu = selectiveClusters
% %     JED_M((N_decision+1):(2*N_decision), 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % selectiveClusters = find(Stimulus.feat(4).selective == -1)';
% % for clu = selectiveClusters
% %     JED_M((N_decision+1):(2*N_decision), 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_i*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefI);
% % end
% % %---------------------------
% % % Oct: E,I --> L; E,I --> R
% % %---------------------------
% % % Background input
% % JED_O = ...
% %     [jBackground_e*ones(2*N_decision,N_e).*(rand(2*N_decision,N_e) < pBackgroundE), ...
% %      jBackground_i*ones(2*N_decision,N_i).*(rand(2*N_decision,N_i) < pBackgroundI)];
% % % Default clustered input
% % JED_O(1:2*N_decision, 1:cusumNcE(Q+1)) = ...
% %     jDefault_e*ones(2*N_decision, cusumNcE(Q+1)).*(rand(2*N_decision, cusumNcE(Q+1)) < pDefaultE);
% % JED_O(1:2*N_decision, (N_e+1):(N_e+cusumNcI(Q+1))) = ...
% %     jDefault_i*ones(2*N_decision, cusumNcI(Q+1)).*(rand(2*N_decision, cusumNcI(Q+1)) < pDefaultI);
% % % Preferential excitatory clustered input
% % selectiveClusters = find(Stimulus.feat(5).selective == 1)';
% % for clu = selectiveClusters
% %     JED_O((N_decision+1):(2*N_decision), 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_e*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefE);
% % end
% % selectiveClusters = find(Stimulus.feat(5).selective == -1)';
% % for clu = selectiveClusters
% %     JED_O((N_decision+1):(2*N_decision), 1+cusumNcE(clu):cusumNcE(clu+1)) = ...
% %         jHigh_i*ones(N_decision,popsize(clu)).*(rand(N_decision,popsize(clu)) < pPrefI);
% % end
% %
% %----------------------------
% % L --> L
% JLL = jdd_in*ones(N_decision,N_decision).*(rand(N_decision,N_decision) < pSelf);
% JLL = JLL - diag(diag(JLL));
% % R --> R
% JRR = jdd_in*ones(N_decision,N_decision).*(rand(N_decision,N_decision) < pSelf);
% JRR = JRR - diag(diag(JRR));
% % R --| L
% JRL = jdd_out*ones(N_decision,N_decision).*(rand(N_decision,N_decision) < pCross);
% % L --| R
% JLR = jdd_out*ones(N_decision,N_decision).*(rand(N_decision,N_decision) < pCross);
% %----------------------------
% % % Construct weight matrices 
% % J_decision_S = [JED_S , [JLL ; JLR], [JRL ; JRR]];
% % J_decision_Q = [JED_Q , [JLL ; JLR], [JRL ; JRR]];
% % J_decision_M = [JED_M , [JLL ; JLR], [JRL ; JRR]];
% % J_decision_O = [JED_O , [JLL ; JLR], [JRL ; JRR]];
% % J_decision = struct();
% % J_decision.S = J_decision_S;
% % J_decision.Q = J_decision_Q;
% % J_decision.M = J_decision_M;
% % J_decision.O = J_decision_O;
% J_decision = [J_D, [JLL ; JLR], [JRL ; JRR]];

%%
% params.popsize = popsize;
% params.clustermatrix = clustermatrix;

%-------------------
% PLOT weight matrix
%-------------------
% figure(1000); clf; 
% subplot(2,1,1);
% colormap(parula); %xy=J; fun_colormapLim;
% imagesc(J);
% aux.figset(gca,'neurons','neurons','weights',10);
% colorbar;
% subplot(2,1,2);
% lambda = eig(J);
% plot(real(lambda),imag(lambda),'.');
% aux.figset(gca,'Re(\lambda)','Im(\lambda)','eig(weights)',10);
% % saveas(gcf,fullfile('data','weights.pdf'),'pdf');
% 
% if any(strcmp(clusters,'decision'))
%     figure(1001); clf; 
%     colormap(parula); 
%     subplot(2,2,1);
%     imagesc(J_decision.S);
%     aux.figset(gca,'neurons','neurons','J_{Decision} (Sucrose)',10);
%     colorbar;
%     subplot(2,2,2);
%     imagesc(J_decision.M);
%     aux.figset(gca,'neurons','neurons','J_{Decision} (Maltose)',10);
%     colorbar;
%     subplot(2,2,3);
%     imagesc(J_decision.Q);
%     aux.figset(gca,'neurons','neurons','J_{Decision} (Quinine)',10);
%     colorbar;
%     subplot(2,2,4);
%     imagesc(J_decision.O);
%     aux.figset(gca,'neurons','neurons','J_{Decision} (Octaacetate)',10);
%     colorbar;
% end

