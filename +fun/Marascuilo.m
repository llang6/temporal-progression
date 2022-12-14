function RES = Marascuilo(CONT_TAB,ALPHA)
%RES = Marascuilo(CONT_TAB,ALPHA)
%
% Marascuilo post-hoc test (compares all categories pair-wise after a
% significant chi-squared test)
% 
% Input: 
% CONT_TAB: an nx2 contingency table where each row is a different category 
% and the columns are number of "successes," then total number of trials
% ALPHA: significance level for test (if not provided, default to 0.05)
% 
% Returns results as an (nchoose2)x3 matrix where each row is a pairing
% between categories and columns tell you the categories being compared (in
% terms of their row indices in the original contingency table) and the
% respective result for that pairing (yes/no significance at alpha level)
%
% This code follows the procedure explained in Methods section of La Camera
% and Richmond 2008 ("Modeling the Violation of Reward Maximization...",
% PLOS Comp Bio)
%
% -LL
%
if nargin < 2
    ALPHA = 0.05;
end
num_cats = size(CONT_TAB,1);
pairs = nchoosek(1:num_cats,2);
num_pairs = size(pairs,1);
RES = NaN(num_pairs,3);
for n = 1:num_pairs
    RES(n,1) = pairs(n,1);
    RES(n,2) = pairs(n,2);
    x_i = CONT_TAB(pairs(n,1),1);
    n_i = CONT_TAB(pairs(n,1),2);
    f_i = x_i/n_i;
    x_j = CONT_TAB(pairs(n,2),1);
    n_j = CONT_TAB(pairs(n,2),2);
    f_j = x_j/n_j;
    Mij = abs(f_i-f_j); % we want the C.I. on this distance
    % find the critical value
    step = 0.001;
    X2crit = step;
    while chi2cdf(X2crit,num_cats-1,'upper') > ALPHA
        X2crit = X2crit + step;
    end
    % calculate C.I. and compare to critical value
    Mij_hat = sqrt(X2crit)*sqrt(f_i*(1-f_i)/n_i + f_j*(1-f_j)/n_j);
    if Mij > Mij_hat
        RES(n,3) = true;
    else
        RES(n,3) = false;
    end
end