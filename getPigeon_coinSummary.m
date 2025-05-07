function coinCounts_ = getPigeon_coinSummary(dataTable, options)
% function coinCounts_ = getPigeon_coinSummary(dataTable, options)
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
% Output is matrix of coin counts per subject (dim 1) and block (dim 2)

arguments
    dataTable;
    options.blocks = 'all' % all or vector
    options.splitBySNR = true;
    options.useGoodTrials = false;
end

% Parse subjects
subjectIndices = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjectIndices);

% Parse blocks
if isnumeric(options.blocks)
    blockIndices = options.blocks;
else
    blockIndices = unique(dataTable.blockIndex);
end
numBlocks = length(blockIndices);

% Split by SNR
if options.splitBySNR
    absSNR = abs(dataTable.snr);
    snrs = nonanunique(absSNR);
    numSNRs = length(snrs);
    Lsnrs = false(size(dataTable,1), numSNRs);
    for rr = 1:numSNRs
        Lsnrs(:,rr) = absSNR == snrs(rr);
    end
else
    numSNRs = 1;
    Lsnrs = true(size(dataTable(:,1),1),1);
end

% Standard selection array
if options.useGoodTrials
    Lg = getPigeon_goodTrialArray(dataTable);
else
    Lg = true(size(dataTable,1),1);
end

% Set up return matrix
coinCounts_ = zeros(numSubjects,numBlocks, numSNRs);

% Loop through the subjects
for ss = 1:length(subjectIndices)
    Ls = Lg & dataTable.subjectIndex==subjectIndices(ss);

    % Loop through the blocks
    for bb = 1:length(blockIndices)
        Lsb = Ls & dataTable.blockIndex == blockIndices(bb);

        % Loop through the SNRs
        for rr = 1:numSNRs
            Lsnr = Lsb & Lsnrs(:,rr);

            if any(Lsnr)
                coinCounts_(ss,bb,rr) = dataTable.coinCount(find(Lsnr, 1, 'last'));
            else
                % disp([ss bb rr])
            end
        end
    end
end