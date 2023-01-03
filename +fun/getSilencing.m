function SILENCING_CURVE = getSilencing(START,DURATION,STRENGTH,WARMUP,TOTAL_T,DT,varargin)
%SILENCING_CURVE = getSilencing(START,DURATION,STRENGTH,WARMUP,TOTAL_T,DT,...) 
%
%Creates a square pulse function that is (1 + STRENGTH/100) from START [s]
%     to START + DURATION [s], and is 1 otherwise
%WARMUP, TOTAL_T, and DT are the time parameters we need to form the trial
%     time vector. Note that these are in [ms] and the silencing parameters
%     are converted to [ms] inside this code
%The value of this function will be a multiplicative modifier of baseline
%     external current to all inhibitory neurons, so e.g. a STRENGTH value of
%     100 means a 100% increase in baseline
%
%Note: if 'center' is given as the last argument to the function, the first
%     argument (START) will be used as the *center* time of the pulse, rather
%     than as the start time of the pulse (which is the default behavior)
%
% -LL
%
timeV = 0:DT:(WARMUP+TOTAL_T);
if isempty(varargin)
    % default behavior: use first argument as silencing start time
    start = START;
elseif strcmp(varargin{1},'center')
    % optional behavior: use first argument as silencing center time
    start = START - DURATION/2;
else
    error('Improper syntax. The last argument can only be ''center'' or omitted.');
end
SILENCING_CURVE = ones(1,length(timeV));
SILENCING_CURVE((timeV>=WARMUP+1000*start)&(timeV<=WARMUP+1000*start+1000*DURATION)) = 1 + STRENGTH/100;
end
