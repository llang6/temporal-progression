function P = chi2test(x,varargin)
% function P = chi2test(x,yates_correction);
%
% proportion test for N proportions.
%
% It performs the chi^2 test for Nx2 contingency table x = matrix(N,2) and
% returns the P-value. This is also called the proportion test for the N
% proportions success(i)/tot.count(i), for i=1,2,...,N.
% It can also perform yates correction for small samples (should be used
% for bins with less than 5 counts, however be aware that it tends to
% over-correct). Requires statistics toolbox (chi2ddf function). 
%
% INPUT:
% - x: Nx2 contingency table where x(:,1)=successes and x(:,2)=tot.counts
% - yates_correction: if ==1, performs Yates correction for small samples
%   (default: 0).
% 
% OUTPUT:
% - P: the P value of the test.
% 
% NOTE: this is better than 'chisquarecont' and 'prop_test', which only
% apply to 2x2 tables (but on those tables all these functions give the
% same results). 
%
% INFO: Based on a lecture by Sanjeev Sanjeev Pillai
% http://jura.wi.mit.edu/bio/education/hot_topics/pdf/matlab.pdf
% Also, check out his other stuff such as visualization of complex data:
% http://jura.wi.mit.edu/bio/education/hot_topics/
% Yates correction taken from Laurie's prop_test (see).
%
% NOTE: in Pillai's version, x(:,2)=failures (and not tot.counts as here).
%
% GLC Nov 2013 - revised Nov 2014 (added Yates correction)

correct=0; % default
if nargin>1 correct=varargin{1}; end

x(:,2)=x(:,2)-x(:,1); % now x(:,2)= # failures
e=sum(x')'*sum(x)/sum(sum(x));
X2=(x-e).^2./e;
if correct X2=(abs(x-e)-0.5).^2./e; end % Yates correction
X2=sum(sum(X2));
df=prod(size(x)-[1,1]);
P=1-chi2cdf(X2,df);