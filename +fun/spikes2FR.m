function [T_VECTOR,FR_VECTOR] = spikes2FR(SPIKE_DATA,BIN_EDGES,INDS)
%[T_VECTOR,FR_VECTOR] = spikes2FR(SPIKE_DATA,BIN_EDGES,INDS) 
%
%Calculates the average firing rate over time (in bins defined by BIN_EDGES
%    [s]) of neurons (indexed by INDS) whose spike times [ms] and indices
%    are contained in SPIKE_DATA
%Assumes SPIKE_DATA is formatted as in simulation, i.e. a 2-column matrix
%    with the first column indicating spike times (in ms) and the second
%    column indicating which neuron spiked at that time
%INDS is the indices (i.e., values in column 2 of SPIKE_DATA) of all
%    neurons you want to average over (could be the indices of all neurons
%    in a particular cluster, the index of a single neuron, etc.)
%
% -LL
%
binSize = diff(BIN_EDGES(1:2));
spikes = SPIKE_DATA(ismember(SPIKE_DATA(:,2),INDS),:);
FR_VECTOR = zeros(1,length(BIN_EDGES)-1);
for i = 1:length(BIN_EDGES)-1
    FR_VECTOR(i) = sum(spikes(:,1)/1000>=BIN_EDGES(i)&spikes(:,1)/1000<BIN_EDGES(i+1))/length(INDS)/binSize;
end
T_VECTOR = BIN_EDGES(1:end-1)+binSize/2;
end