function VARIABLE = loadVar(FILE)
%VARIABLE = loadVar(FILE)
%
% Takes the name of a '.mat' file in your current directory as input and,
% if there is only one variable saved inside, loads the variable and saves
% it to the workspace as the output.
%
% This function was designed to make loading simple files easier, since you
% don't need to know what the name of the variable inside the file
% originally was -- you only need to know the name of the file it was
% saved to.
%
% -LL
%
S = load(FILE);
fields = fieldnames(S);
if length(fields)~=1, error('Input file has multiple variables saved inside. This function only works for files containing a single variable.'); end
VARIABLE = S.(fields{1});
end