function Lg_ = getPigeon_goodTrialArray(dataTable, options)
% function Lg_ = getPigeon_goodTrialArray(dataTable, options)
%
% Standard selection criteria (to ensure consistency across analyses)

arguments
    dataTable = [];
    options.DT = 1;
    options.trialNumber = 10;
    options.blockIndex = [];
end

Lg_ = true(size(dataTable,1),1);

% Min DT
if ~isempty(options.DT)
    Lg_ = Lg_ & (dataTable.DT >= options.DT);
end

% Min trial number
if ~isempty(options.trialNumber)
    Lg_ = Lg_ & (dataTable.trialNumber >= options.trialNumber);
end

% Block index
if ~isempty(options.blockIndex)
    Lg_ = Lg_ & ismember(dataTable.blockIndex, options.blockIndex);
end

