function CURVE_BINARIZED = threshold(CURVE,THRESHOLD_VALUE,MIN_TIME,BIN_SIZE)
%CURVE_BINARIZED = threshold(CURVE,THRESHOLD_VALUE,MIN_TIME,BIN_SIZE) 
%
%Converts CURVE into a logical vector indicating which entries are above
%    THRESHOLD_VALUE for at least MIN_TIME
%BIN_SIZE is the amount of time between consecutive values in CURVE (must
%    be same units as MIN_TIME)
%
% -LL
%
curve_binarized = CURVE>THRESHOLD_VALUE;
onsets = find(diff([0 curve_binarized 0])==1); 
offsets = find(diff([0 curve_binarized 0])==-1);
inds = find(BIN_SIZE*(offsets-onsets)<MIN_TIME);
if ~isempty(inds)
    for ind = inds
        curve_binarized(onsets(ind):min(offsets(ind)-1,length(curve_binarized))) = 0;
    end
end
CURVE_BINARIZED = curve_binarized;
end