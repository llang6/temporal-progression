function CLASSIFIED_LABEL = classifyState(STATE,MY_CELL,varargin)
%CLASSIFIED_LABEL = classifyState(STATE,MY_CELL)
%
% The input MY_CELL is critical for this function. It is a cell with very
% specific entries. The 1st and 3rd columns are vectors of state IDs over
% trials for states that were decoded, and the 2nd and 4th columns are
% total numbers of trials. The rows are for different trial types --
% Sucrose, Quinine, Maltose, Octaacetate in that order -- and the first two
% columns are for correct trials while the last two are for incorrect
% trials.
%
% Given the MY_CELL input and a particular STATE to classify, classifyState
% will return a label for the state, according to a set of rules defined by
% the CLASSIFICATION_MODE. 
%
% The mode used for classifying states was *mode 3*, which here is the
% default. Look in the switch/case block under case 3 for classification
% details.
%
% Optional additional arguments can be given as 'Name',Value pairs:
%     'classificationMode': a number from 1 through 5 that changes
%                           classification rules (default 3)
%     'stricterSubclassification': true/false flag for imposing an
%                                  additional criterion on
%                                  subclassification of Exclusive
%                                  Decision-coding states (if true, session
%                                  must include at least 10 incorrectly
%                                  performed trials of each cued direction
%                                  for subclassification to proceed)
%                                  (default false)
%     'isVerbose': true/false flag for outputting results to command window
%                  (default false)
%
% -LL
%

% parse input
if isempty(varargin)
    % default parameters
    CLASSIFICATION_MODE = 3;
    IMPOSE_STRICTER_CRITERIA = false;
    IS_VERBOSE = false;
else
    ind = find(cellfun(@(x)strcmpi(x,'classificationMode'),varargin),1);
    if ~isempty(ind), CLASSIFICATION_MODE = varargin{ind+1};
    else, CLASSIFICATION_MODE = 3; end
    ind = find(cellfun(@(x)strcmpi(x,'isVerbose'),varargin),1);
    if ~isempty(ind), IS_VERBOSE = varargin{ind+1};
    else, IS_VERBOSE = false; end
    ind = find(cellfun(@(x)strcmpi(x,'stricterSubclassification'),varargin),1);
    if ~isempty(ind), IMPOSE_STRICTER_CRITERIA = varargin{ind+1};
    else, IMPOSE_STRICTER_CRITERIA = false; end
end

% set of decoded states
decodedStates = unique([MY_CELL{:,[1,3]}]);

% classify
switch CLASSIFICATION_MODE
%-----------------------------------------------------------------------------------------------------------------------    
    case 1 
        % simplest way: split 4 categories by cued direction and by taste
        % quality, run a chi-squared for each

        if ~ismember(STATE,decodedStates)
            CLASSIFIED_LABEL = 'Not decoded';
        else
            contTab_Dec_Corr = [sum([MY_CELL{[1,2],1}]==STATE), sum([MY_CELL{[1,2],2}]);
                                sum([MY_CELL{[3,4],1}]==STATE), sum([MY_CELL{[3,4],2}])];
            contTab_Dec_Inc = [sum([MY_CELL{[1,2],3}]==STATE), sum([MY_CELL{[1,2],4}]);
                               sum([MY_CELL{[3,4],3}]==STATE), sum([MY_CELL{[3,4],4}])];
            contTab_Tas_Corr = [sum([MY_CELL{[1,3],1}]==STATE), sum([MY_CELL{[1,3],2}]);
                                sum([MY_CELL{[2,4],1}]==STATE), sum([MY_CELL{[2,4],2}])];
            prefLR_Corr = diff(contTab_Dec_Corr(:,1)./contTab_Dec_Corr(:,2));
            prefLR_Inc = diff(contTab_Dec_Inc(:,1)./contTab_Dec_Inc(:,2));
            prefSwitchLR = prefLR_Corr*prefLR_Inc < 0;
            prefRetainLR = prefLR_Corr*prefLR_Inc > 0;
            isDecisionCode = fun.chi2test(contTab_Dec_Corr) < 1 - (1 - 0.05)^(1/length(decodedStates));
            isQualityCode = fun.chi2test(contTab_Tas_Corr) < 1 - (1 - 0.05)^(1/length(decodedStates));
            if isDecisionCode
                if isQualityCode
                    CLASSIFIED_LABEL = 'Dual-coding';
                else
                    if IMPOSE_STRICTER_CRITERIA
                        additionalCriteria = all(contTab_Dec_Inc(:,2)>=10);
                    else
                        additionalCriteria = true;
                    end
                    if prefSwitchLR && additionalCriteria
                        CLASSIFIED_LABEL = 'Action-coding';
                    elseif prefRetainLR && additionalCriteria
                        CLASSIFIED_LABEL = 'Cue-coding';
                    else
                        CLASSIFIED_LABEL = 'Exclusive Decision-coding';
                    end
                end
            elseif isQualityCode
                CLASSIFIED_LABEL = 'Exclusive Quality-coding';
            else
                CLASSIFIED_LABEL = 'Non-coding';
            end
        end
%-----------------------------------------------------------------------------------------------------------------------    
    case 2 
        % same as case 1 except we first run a 4-category chi-squared and
        % use all 6 pairwise post-hoc Marascuilo tests to identify the
        % signature of taste ID-coding states, then exclude any of these
        % from the normal analysis

        if ~ismember(STATE,decodedStates)
            CLASSIFIED_LABEL = 'Not decoded';
        else
            contTab_ID_Corr = [sum(MY_CELL{1,1}==STATE), MY_CELL{1,2};
                               sum(MY_CELL{2,1}==STATE), MY_CELL{2,2};
                               sum(MY_CELL{3,1}==STATE), MY_CELL{3,2};
                               sum(MY_CELL{4,1}==STATE), MY_CELL{4,2}];
            contTab_Dec_Corr = [sum([MY_CELL{[1,2],1}]==STATE), sum([MY_CELL{[1,2],2}]);
                                sum([MY_CELL{[3,4],1}]==STATE), sum([MY_CELL{[3,4],2}])];
            contTab_Dec_Inc = [sum([MY_CELL{[1,2],3}]==STATE), sum([MY_CELL{[1,2],4}]);
                               sum([MY_CELL{[3,4],3}]==STATE), sum([MY_CELL{[3,4],4}])];
            contTab_Tas_Corr = [sum([MY_CELL{[1,3],1}]==STATE), sum([MY_CELL{[1,3],2}]);
                                sum([MY_CELL{[2,4],1}]==STATE), sum([MY_CELL{[2,4],2}])];
            prefLR_Corr = diff(contTab_Dec_Corr(:,1)./contTab_Dec_Corr(:,2));
            prefLR_Inc = diff(contTab_Dec_Inc(:,1)./contTab_Dec_Inc(:,2));
            prefSwitchLR = prefLR_Corr*prefLR_Inc < 0;
            prefRetainLR = prefLR_Corr*prefLR_Inc > 0;
            if fun.chi2test(contTab_ID_Corr) < 1-(1-0.05)^(1/length(decodedStates))
                res = fun.Marascuilo(contTab_ID_Corr, 1-(1-0.05)^(1/length(decodedStates)));
                allPostSig = false(1,4);
                for group = 1:4
                    [inds,~] = find(res(:,1:2)==group);
                    allPostSig(group) = all(res(inds,3));
                end
                if sum(allPostSig)==1
                    CLASSIFIED_LABEL = 'Taste ID-coding';
                    return
                end
            end
            isDC = fun.chi2test(contTab_Dec_Corr) < 1 - (1 - 0.05)^(1/length(decodedStates));
            isQC = fun.chi2test(contTab_Tas_Corr) < 1 - (1 - 0.05)^(1/length(decodedStates));
            if isDC && isQC
                CLASSIFIED_LABEL = 'Dual-coding';
            elseif ~isDC && ~isQC
                CLASSIFIED_LABEL = 'Non-coding';
            elseif isQC && ~isDC
                CLASSIFIED_LABEL = 'Exclusive Quality-coding';
            elseif isDC && ~isQC
                if IMPOSE_STRICTER_CRITERIA
                    additionalCriteria = all(contTab_Dec_Inc(:,2)>=10);
                else
                    additionalCriteria = true;
                end
                if prefSwitchLR && additionalCriteria
                    CLASSIFIED_LABEL = 'Action-coding';
                elseif prefRetainLR && additionalCriteria
                    CLASSIFIED_LABEL = 'Cue-coding';
                else
                    CLASSIFIED_LABEL = 'Exclusive Decision-coding';
                end
            end
        end
%-----------------------------------------------------------------------------------------------------------------------        
    case 3
        % same as case 2, but the 4-way chi-squared test used for
        % identifying ID-coding states must actually be passed or the state
        % is 'Non-coding'
        
        if ~ismember(STATE,decodedStates)
            % not decoded --> 'Not decoded'
            CLASSIFIED_LABEL = 'Not decoded';
        else
            % was decoded --> classify the state
            contTab_ID_Corr = [sum(MY_CELL{1,1}==STATE), MY_CELL{1,2};
                               sum(MY_CELL{2,1}==STATE), MY_CELL{2,2};
                               sum(MY_CELL{3,1}==STATE), MY_CELL{3,2};
                               sum(MY_CELL{4,1}==STATE), MY_CELL{4,2}];
            contTab_Dec_Corr = [sum([MY_CELL{[1,2],1}]==STATE), sum([MY_CELL{[1,2],2}]);
                                sum([MY_CELL{[3,4],1}]==STATE), sum([MY_CELL{[3,4],2}])];
            contTab_Dec_Inc = [sum([MY_CELL{[1,2],3}]==STATE), sum([MY_CELL{[1,2],4}]);
                               sum([MY_CELL{[3,4],3}]==STATE), sum([MY_CELL{[3,4],4}])];
            contTab_Tas_Corr = [sum([MY_CELL{[1,3],1}]==STATE), sum([MY_CELL{[1,3],2}]);
                                sum([MY_CELL{[2,4],1}]==STATE), sum([MY_CELL{[2,4],2}])];
            prefLR_Corr = diff(contTab_Dec_Corr(:,1)./contTab_Dec_Corr(:,2));
            prefLR_Inc = diff(contTab_Dec_Inc(:,1)./contTab_Dec_Inc(:,2));
            prefSwitchLR = prefLR_Corr*prefLR_Inc < 0;
            prefRetainLR = prefLR_Corr*prefLR_Inc > 0;
            if fun.chi2test(contTab_ID_Corr) >= 1-(1-0.05)^(1/length(decodedStates))
                % failed prelim chi-squared --> 'Non-coding'
                CLASSIFIED_LABEL = 'Non-coding';
            else
                % passed prelim chi-squared, now check for ID-coding...
                res = fun.Marascuilo(contTab_ID_Corr, 1-(1-0.05)^(1/length(decodedStates)));
                allPostSig = false(1,4);
                for group = 1:4
                    [inds,~] = find(res(:,1:2)==group);
                    allPostSig(group) = all(res(inds,3));
                end
                if sum(allPostSig)==1
                    % passed ID-coding test --> 'ID-coding'
                    CLASSIFIED_LABEL = 'Taste ID-coding';
                else
                    % failed ID-coding test --> check for Decision- and Quality-coding
                    isDC = fun.chi2test(contTab_Dec_Corr) < 1 - (1 - 0.05)^(1/length(decodedStates));
                    isQC = fun.chi2test(contTab_Tas_Corr) < 1 - (1 - 0.05)^(1/length(decodedStates));
                    if isDC && isQC
                        % Decision- and Quality-coding --> 'Dual-coding'
                        CLASSIFIED_LABEL = 'Dual-coding';
                    elseif ~isDC && ~isQC
                        % not Decision- nor Quality-coding --> 'Non-coding'
                        CLASSIFIED_LABEL = 'Non-coding';
                    elseif isQC && ~isDC
                        % only Quality-coding --> 'Exclusive Quality-coding'
                        CLASSIFIED_LABEL = 'Exclusive Quality-coding';
                    elseif isDC && ~isQC
                        % only Decision-coding --> try to subclassify
                        if IMPOSE_STRICTER_CRITERIA
                            % must have at least 10 incorrectly performed trials of
                            % each cued direction
                            additionalCriteria = all(contTab_Dec_Inc(:,2)>=10);
                        else
                            % no additional criteria (true so you automatically pass)
                            additionalCriteria = true;
                        end
                        if prefSwitchLR && additionalCriteria
                            % switches preference --> 'Action-coding'
                            CLASSIFIED_LABEL = 'Action-coding';
                        elseif prefRetainLR && additionalCriteria
                            % retains preference --> 'Cue-coding'
                            CLASSIFIED_LABEL = 'Cue-coding';
                        else
                            % no preference is ambiguous ... can't subclassify,
                            % stick with default --> 'Exclusive Decision-coding'
                            CLASSIFIED_LABEL = 'Exclusive Decision-coding';
                        end
                    end
                end
            end
        end
%-----------------------------------------------------------------------------------------------------------------------        
    case 4 
        % use only the 4-category chi-squared test and the 6 pairwise
        % post-hoc Marascuilo tests to match all signatures with coding
        % classification
        
        if ~ismember(STATE,decodedStates)
            CLASSIFIED_LABEL = 'Not decoded';
        else
            contTab_ID_Corr = [sum(MY_CELL{1,1}==STATE) MY_CELL{1,2};
                               sum(MY_CELL{2,1}==STATE) MY_CELL{2,2};
                               sum(MY_CELL{3,1}==STATE) MY_CELL{3,2};
                               sum(MY_CELL{4,1}==STATE) MY_CELL{4,2}];
            prefSwitchLR = (sum([MY_CELL{[3,4],1}]==STATE)/sum([MY_CELL{[3,4],2}]) - sum([MY_CELL{[1,2],1}]==STATE)/sum([MY_CELL{[1,2],2}]))*(sum([MY_CELL{[3,4],3}]==STATE)/sum([MY_CELL{[3,4],4}]) - sum([MY_CELL{[1,2],3}]==STATE)/sum([MY_CELL{[1,2],4}])) < 0;
            prefRetainLR = (sum([MY_CELL{[3,4],1}]==STATE)/sum([MY_CELL{[3,4],2}]) - sum([MY_CELL{[1,2],1}]==STATE)/sum([MY_CELL{[1,2],2}]))*(sum([MY_CELL{[3,4],3}]==STATE)/sum([MY_CELL{[3,4],4}]) - sum([MY_CELL{[1,2],3}]==STATE)/sum([MY_CELL{[1,2],4}])) > 0;
            if fun.chi2test(contTab_ID_Corr) >= 1-(1-0.05)^(1/length(decodedStates))
                CLASSIFIED_LABEL = 'Non-coding';
                return
            else
                res = fun.Marascuilo(contTab_ID_Corr, 1-(1-0.05)^(1/length(decodedStates)));
                match = double(res(:,3));
                if IS_VERBOSE
                    fprintf('\nPost-hoc result: ');
                    for j = 1:length(match), fprintf('%i',match(j)); end
                    fprintf('\n');
                end
                if all(match==0)
                    CLASSIFIED_LABEL = 'Non-coding';
                    return
                else
                    idTemplate = [1 1 1 1 1 1 1 1 0 1 0 0 0 1 0 0 ; ...
                                  1 1 1 1 0 1 0 0 1 1 1 1 0 0 1 0 ; ...
                                  1 1 1 1 0 0 1 0 0 0 1 0 1 1 1 1 ; ...
                                  0 1 0 0 1 1 1 1 1 1 1 1 0 0 0 1 ; ...
                                  0 0 1 0 1 1 1 1 0 0 0 1 1 1 1 1 ; ...
                                  0 0 0 1 0 0 0 1 1 1 1 1 1 1 1 1];
                    if ismember(match',idTemplate','rows')
                        CLASSIFIED_LABEL = 'Taste ID-coding';
                    elseif all(match==[0 1 1 1 1 0]')
                        CLASSIFIED_LABEL = 'Exclusive Decision-coding';
                        if IMPOSE_STRICTER_CRITERIA
                            additionalCriteria = all(contTab_Dec_Inc(:,2)>=10);
                        else
                            additionalCriteria = true;
                        end
                        if prefSwitchLR && additionalCriteria
                            CLASSIFIED_LABEL = 'Action-coding';
                        elseif prefRetainLR && additionalCriteria
                            CLASSIFIED_LABEL = 'Cue-coding';
                        end
                    elseif all(match==[1 0 1 1 0 1]')
                        CLASSIFIED_LABEL = 'Exclusive Quality-coding';
                    else
                        CLASSIFIED_LABEL = 'Other';
                    end
                end
            end
        end
%-----------------------------------------------------------------------------------------------------------------------        
end

if IS_VERBOSE
    fprintf('\nP(State %i | Sucrose, Correct) = %.2f\n',STATE,sum(MY_CELL{1,1}==STATE)/MY_CELL{1,2});
    fprintf('P(State %i | Quinine, Correct) = %.2f\n',STATE,sum(MY_CELL{2,1}==STATE)/MY_CELL{2,2});
    fprintf('P(State %i | Maltose, Correct) = %.2f\n',STATE,sum(MY_CELL{3,1}==STATE)/MY_CELL{3,2});
    fprintf('P(State %i | Octaacetate, Correct) = %.2f\n',STATE,sum(MY_CELL{4,1}==STATE)/MY_CELL{4,2});
    fprintf('P(State %i | Sucrose, Incorrect) = %.2f\n',STATE,sum(MY_CELL{1,3}==STATE)/MY_CELL{1,4});
    fprintf('P(State %i | Quinine, Incorrect) = %.2f\n',STATE,sum(MY_CELL{2,3}==STATE)/MY_CELL{2,4});
    fprintf('P(State %i | Maltose, Incorrect) = %.2f\n',STATE,sum(MY_CELL{3,3}==STATE)/MY_CELL{3,4});
    fprintf('P(State %i | Octaacetate, Incorrect) = %.2f\n',STATE,sum(MY_CELL{4,3}==STATE)/MY_CELL{4,4});
end
    
end