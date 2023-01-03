function [stimulus_save, params] = fun_stim(params)

% unpack vars
p = params.p;
Stimulus = params.Stimulus;
stimuli = params.stimuli;
Sim = params.Sim;

% for each event, create external current and stim properties in .Ext 
% and add clusters selective to the event .Stimulus.feat(n).StimClust
stimulus_save = struct('Ext',[],'Stimulus',[]);

% select stimuli
temp_Stimulus = struct('input',Stimulus.input);
indfeat = zeros(1,numel(stimuli));
for ev = 1:numel(stimuli)
    fprintf('Stimulus %s\n',stimuli{ev});
    % match current stimuli to features in Stimulus
    indfeat(ev) = find(cell2mat(arrayfun(@(x) strcmp(stimuli{ev},x(:).name), ...
        Stimulus.feat(:),'uniformoutput',false)));
end
fprintf('\n');
if ~isempty(indfeat)
    temp_Stimulus.feat(1:numel(indfeat)) = Stimulus.feat(indfeat);
    for n = 1:numel(indfeat)
        if any(strcmp(temp_Stimulus.feat(n).name,{'Sucrose','Maltose','Quinine','Octaacetate'}))
            sclust_inc = [];
            sclust_dec = [];
            if ~isempty(temp_Stimulus.feat(n).selective)
                sclust_inc = find(temp_Stimulus.feat(n).selective(:,1)==1);
                sclust_dec = find(temp_Stimulus.feat(n).selective(:,1)==-1);
            end
            temp_Stimulus.feat(n).StimClust_Inc = sclust_inc;
            temp_Stimulus.feat(n).StimClust_Dec = sclust_dec;
        else
            sclust = [];
            if ~isempty(temp_Stimulus.feat(n).selective)
                sclust = find(temp_Stimulus.feat(n).selective(1,:));
            end
            temp_Stimulus.feat(n).StimClust = sclust;
        end
    end
end
Stimulus = temp_Stimulus;
Ext = struct('Mu',[]);

% LOAD PARAMETERS
fieldNames = {'Sim','Network','p','popsize','clustermatrix','N_e','N_i','Cext','Jee_ext','Jie_ext','ni_ext','tau_e','tau_i','fieldNames'};
aux.v2struct(params,fieldNames);
cusumNcE = [0 cumsum(popsize)'];
Tseq = Sim.t_Start:Sim.dt_step:Sim.t_End; 

if ~isempty(stimuli)
    feat = Stimulus.feat;
    nstim = numel(feat); % number of stimuli in current trials
    stim = repmat(struct('profile',[],'ind',[],'interval',[]),1,nstim);
    temp_ind = repmat(struct('ind',[]),1,nstim); % stores indices for mixed cue (see below)
    for n = 1:nstim
        % stimulus interval
        interv = feat(n).interval;
        % stimulus profile
        if ~any(strcmp(feat(n).name,{'Sucrose','Quinine','Maltose','Octaacetate'}))
            Gain = feat(n).gain;
            if ~isempty(strfind(feat(n).name,'gauss'))
                Gain = 1; % with gaussian stim set profile to peak at 1, then multiply each profile by gaussian with SD feat(n).gain for each neuron in feat(n).gauss
            end
            Profile = feat(n).profile;
            Profile = @(t)Profile(t-interv(1));
            MaxThInput = max(abs(Profile(Tseq(Tseq>interv(1) & Tseq<interv(2)))));
            Profile = @(t)Gain*Profile(t)/MaxThInput;
            stim(n).profile = @(t)Profile(t); % fraction increase above baseline
        else
            Gain_inc = feat(n).gain_inc;
            Gain_dec = feat(n).gain_dec;
            if ~isempty(strfind(feat(n).name,'gauss'))
                Gain_inc = 1; % with gaussian stim set profile to peak at 1, then multiply each profile by gaussian with SD feat(n).gain for each neuron in feat(n).gauss
                Gain_dec = -1;
            end
            Profile_inc = feat(n).profile_inc;
            Profile_inc = @(t)Profile_inc(t-interv(1));
            MaxThInput_inc = max(abs(Profile_inc(Tseq(Tseq>interv(1) & Tseq<interv(2)))));
            Profile_inc = @(t)Gain_inc*Profile_inc(t)/MaxThInput_inc;
            stim(n).profile_inc = @(t)Profile_inc(t);
            Profile_dec = feat(n).profile_dec;
            Profile_dec = @(t)Profile_dec(t-interv(1));
            MaxThInput_dec = max(abs(Profile_dec(Tseq(Tseq>interv(1) & Tseq<interv(2)))));
            Profile_dec = @(t)Gain_dec*Profile_dec(t)/MaxThInput_dec;
            stim(n).profile_dec = @(t)Profile_dec(t);
        end
        
        % selective neurons
        if any(strcmp(feat(n).name,{'Sucrose','Quinine','Maltose','Octaacetate'}))
            StimClust_Inc = Stimulus.feat(n).StimClust_Inc; % clusters activated by current stimulus
    %             a = randperm(numel(StimClust));
    %             StimClust_Inc = StimClust(a(1:round(feat(n).fracInc*numel(StimClust))));
            StimClust_Dec = Stimulus.feat(n).StimClust_Dec;
            if isfield(Stimulus.feat(n),'ind_inc') && isfield(Stimulus.feat(n),'ind_dec')
                ind_inc = Stimulus.feat(n).ind_inc;
                ind_dec = Stimulus.feat(n).ind_dec;
            else
                % units selective to stimulus
                ind_inc = []; % indices of stim sel units that increase
                ind_dec = [];
                switch feat(n).selectivity
                    case 'mixed'
                        for c = StimClust_Inc'
                            pop_ind = find(clustermatrix(:,c));
                            for k = 1:numel(pop_ind)
                                ind_inc = [ind_inc cusumNcE(pop_ind(k))+1:cusumNcE(pop_ind(k)+1)]; % stim selective units
                            end
                        end
                        for c = StimClust_Dec'
                            pop_ind = find(clustermatrix(:,c));
                            for k = 1:numel(pop_ind)
                                ind_dec = [ind_dec cusumNcE(pop_ind(k))+1:cusumNcE(pop_ind(k)+1)]; % stim selective units
                            end
                        end
                    case 'exc'
                        ind = 1:N_e;
                    otherwise
                        ind = 1:N_e;
                end
                % sparsify
                a_inc = randperm(numel(ind_inc));
                temp_ind(n).ind_inc = ind_inc;
                ind_inc = ind_inc(a_inc(1:round(feat(n).connectivity*numel(ind_inc))));
                a_dec = randperm(numel(ind_dec));
                temp_ind(n).ind_dec = ind_dec;
                ind_dec = ind_dec(a_dec(1:round(feat(n).connectivity*numel(ind_dec))));
                % gaussian stimulus, draw from randn
                if ~isempty(strfind(feat(n).name,'gauss'))
                    stim(n).gauss_inc = feat(n).gain_inc*randn(numel(ind_inc),1);
                    stim(n).gauss_dec = feat(n).gain_dec*randn(numel(ind_dec),1);
                end
            end
            %
            stim(n).ind_inc = ind_inc;
            stim(n).ind_dec = ind_dec;
            stim(n).interval = interv;
            stim(n).name = feat(n).name;
            stim(n).StimClust_Inc = StimClust_Inc;
            stim(n).StimClust_Dec = StimClust_Dec;
            stim(n).selectivity = feat(n).selectivity;
            % Additional gain modifier to ensure targeted neurons respond
            % nonhomogeneously
            DiffGain_Inc = rand(numel(ind_inc),1);
            seg1 = DiffGain_Inc < 19/25;
            seg2 = DiffGain_Inc >= 19/25 & DiffGain_Inc < 21/25;
            seg3 = DiffGain_Inc >= 21/25;
            DiffGain_Inc(seg1) = 1; % 0.25;
            DiffGain_Inc(seg2) = 1; % 0.5;
            DiffGain_Inc(seg3) = 1; % 0.75;
            stim(n).DiffGain_Inc = DiffGain_Inc;
            
            DiffGain_Dec = rand(numel(ind_dec),1);
            seg1 = DiffGain_Dec < 10/28;
            seg2 = DiffGain_Dec >= 10/28 & DiffGain_Dec < 19/28;
            seg3 = DiffGain_Dec >= 19/28 & DiffGain_Dec < 26/28;
            seg4 = DiffGain_Dec >= 26/28;
            DiffGain_Dec(seg1) = 1; % 1;
            DiffGain_Dec(seg2) = 1; % 0.7;
            DiffGain_Dec(seg3) = 1; % 0.5;
            DiffGain_Dec(seg4) = 1; % 0.25;
            stim(n).DiffGain_Dec = DiffGain_Dec;
        else
            StimClust = Stimulus.feat(n).StimClust; % clusters activated by current stimulus
            % units selective to stimulus
            ind = []; % indices of stim sel units
            switch feat(n).selectivity
                case 'mixed'
                    for c = StimClust
                        pop_ind = find(clustermatrix(:,c));
                        for k = 1:numel(pop_ind)
                            ind = [ind cusumNcE(pop_ind(k))+1:cusumNcE(pop_ind(k)+1)]; % stim selective units
                        end
                    end
                case 'exc'
                    ind = 1:N_e;
                otherwise
                    ind = 1:N_e;
            end
            % sparsify
            a = randperm(numel(ind));
            temp_ind(n).ind = ind;
            ind = ind(a(1:round(feat(n).connectivity*numel(ind))));
            % gaussian stimulus, draw from randn
            if ~isempty(strfind(feat(n).name,'gauss'))
                stim(n).gauss = feat(n).gain*randn(numel(ind),1);
            end
            %
            stim(n).ind = ind;
            stim(n).interval = interv;
            stim(n).name = feat(n).name;
            stim(n).StimClust = StimClust;
            stim(n).selectivity = feat(n).selectivity;
        end
    end
    Ext.stim = stim;
end 
Ext.Mu = params.Mu;

stimulus_save.Ext = Ext;
stimulus_save.Stimulus = temp_Stimulus;
