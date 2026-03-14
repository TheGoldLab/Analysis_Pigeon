function performanceSummary_ = getPigeon_performanceSummary(dataTable, options)
% function performanceSummaryTable = getPigeon_performanceSummary(dataTable, options)
%
% Input is dataTable with columns:
%
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
%
% Output is performanceSummary_ matrix with dimensions:
%   1. subject index
%   2. block index
%   3. SNR index
%   4. accuracy / RT mean / RT STD / n

arguments
    dataTable;
    options.blocks = 'all' % all or vector
    options.snrs = 'all' % all or vector
end

% Parse subjects
subjectIndices = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjectIndices);

% Parse blocks
if isnumeric(options.blocks)
    blockIndices = options.blocks;
else
    blockIndices = 1:3;
end
numBlocks = length(blockIndices);

% Parse SNRs
absSNR = abs(dataTable.snr); % get array of abs values
if isnumeric(options.snrs)
    absSNRs = options.snrs;
else
    absSNRs = nonanunique(absSNR);
end
numSNRs = length(absSNRs);

% Standard selection array
Lg = getPigeon_goodTrialArray(dataTable);

% Set up the summary data matrix    
performanceSummary_ = zeros(numSubjects,numBlocks,numSNRs,4); % last is %cor/rt mean/rt std/n

% Loop through the subjects
for ss = 1:length(subjectIndices)
    Ls = Lg & dataTable.subjectIndex==subjectIndices(ss);

    % Loop through the blocks
    for bb = 1:numBlocks
        Lsb = Ls & dataTable.blockIndex==bb; %Lbound(:,bb);

        % Loop through the SNRS
        for nn = 1:numSNRs
            Lsnr = Lsb & absSNR==absSNRs(nn);
            if any(Lsnr)
                performanceSummary_(ss,bb,nn,:) = [ ...
                    sum(dataTable.correct(Lsnr))./sum(Lsnr), ...
                    mean(dataTable.DT(Lsnr)), ...
                    std(dataTable.DT(Lsnr)), ...
                    sum(Lsnr)];
            else
                performanceSummary_(ss,bb,nn,:) = [nan, nan, nan, 0];
            end
        end
    end
end