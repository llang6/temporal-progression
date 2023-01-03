function CURVE = getActionGate(T_START,DURATION,FLOOR,CEILING,WARMUP,TOTAL_T,DT)
%CURVE = getActionGate(T_START,DURATION,FLOOR,CEILING,WARMUP,TOTAL_T,DT) 
%
%Creates a ramp function that is FLOOR before time T_START [s], then ramps up
%     to CEILING over the course of DURATION [s]
%WARMUP, TOTAL_T, and DT are the time parameters we need to form the trial
%     time vector. Note that these are in [ms] and the ramp parameters are
%     converted to [ms] inside this code
%
% -LL
%
timeV = 0:DT:(WARMUP+TOTAL_T);
CURVE = zeros(1,length(timeV));
CURVE(timeV>=WARMUP+1000*T_START) = min(1/(1000*DURATION)*(0:DT:(TOTAL_T-1000*T_START)),1);
CURVE = (CEILING-FLOOR)*CURVE + FLOOR;
end