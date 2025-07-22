function dataTable_ = getPigeon_blankDataTable(numRows, varargin)
% function dataTable_ = getPigeon_blankDataTable(options)
%
% Varargin is <name> <type> pairs to add as columns
%
% Returns empty dataTable with the following default columns:
%   1. subjectIndex
%   2. blockIndex
%   3. trialNumber
%   4. bound (pigeon position in OL)
%   5. RT (response time, including ndt, in number of steps)
%   6. DT (response time, not including ndt, in number of steps)
%   7. choice (0/1)
%   8. correct (correct=1,error=0)
%   9. coinCount
%   10. snr
%   11. steps (cell)

% Defaults
variableNames = ...
    {'subjectIndex', 'blockIndex', 'trialNumber', 'bound', 'absBound', ...
    'RT', 'DT', 'changePoint', 'choice', 'correct', 'coinCount', 'snr', 'steps', 'congruence'};
variableTypes = cat(2, repmat("double", 1, length(variableNames)-2), "cell", "cell");

% Add optional columns
for nn = 1:2:nargin-2
    variableNames = cat(2,variableNames,varargin{nn});
    variableTypes = cat(2,variableTypes,varargin{nn+1});
end

% Make the table
dataTable_ = table('Size', [numRows, length(variableNames)], ...
    'VariableTypes', variableTypes, ...
    'VariableNames', variableNames);