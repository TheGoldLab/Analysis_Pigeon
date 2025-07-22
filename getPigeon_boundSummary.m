function boundSummary_ = getPigeon_boundSummary(dataTable, options)
% function boundSummary_ = getPigeon_boundSummary(dataTable, options)
%
%   dataTable is:
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
% Returns:
%   boundSummary_ .. matrix with dims:
%       subjects
%       blocks
%       snrs
%       rt bins
%       bound avg, bound std, n, bound avg zscore

arguments
    dataTable;
    options.correctOnly = false;
    options.blocks = 'all' % all or vector
    options.splitBySNR = true;
    options.maxRT = 10; % scalar, 'all', or 'medsplit'
end

% Collect some useful variables
subjectIndices = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjectIndices);
absBounds = abs(dataTable.bound);

% parse blocks
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

% Check for RT dependence (sort by RT bins/median split RT/none)
if isnumeric(options.maxRT)
    RTbins = 1:options.maxRT;
elseif strcmp(options.maxRT, 'medsplit')
    RTbins = 1:2;
else %if strcmp(options.RTs, 'all')
    RTbins = 1;
end
numRTs = length(RTbins);

% Last is: bound avg, bound std, n
boundSummary_ = nan(numSubjects,numBlocks,numSNRs,numRTs,3);

% Selection array
if options.correctOnly
    Lgood = dataTable.correct==1;
else
    Lgood = dataTable.correct>=0;
end
% Lgood = Lgood & isfinite(dataTable.bound);
Lgood = Lgood & dataTable.bound~=0;

% for z-scored bounds
if isnumeric(options.maxRT)
    zbounds = nan(size(dataTable,1),1);
end

% Loop through the subjects
for ss = 1:numSubjects
    Lsubject = Lgood & dataTable.subjectIndex == subjectIndices(ss);

    % Loop through the blocks
    for bb = 1:numBlocks
        Lsb = Lsubject & dataTable.blockIndex==blockIndices(bb);

        % Loop through the SNRs
        for nn = 1:numSNRs
            Lsnr = Lsb & Lsnrs(:,nn);            

            % Per RT type, as specified
            if isnumeric(options.maxRT)

                % z-score per subject/block/snr
                zbounds(Lsnr) = median(absBounds(Lsnr)); %zscore(absBounds(Lsnr));
               
                % Per RT bin
                for rr = 1:numRTs
                    Lrt = Lsnr & dataTable.DT == RTbins(rr);
                    boundSummary_(ss,bb,nn,rr,:) = [ ...
                        mean(absBounds(Lrt), 'omitnan'), ...
                        std(absBounds(Lrt), 'omitnan'), ...
                        sum(Lrt)];
                end

            elseif strcmp(options.maxRT, 'medsplit')

                % Median split
                medRT = median(dataTable.RT(Lsnr));
                Lrts = [ ...
                    Lsnr & dataTable.RT <= medRT, ...
                    Lsnr & dataTable.RT > medRT];
                for rr = 1:2
                    boundSummary_(ss,bb,nn,rr,:) = [ ...
                        mean(absBounds(Lrts(:,rr)), 'omitnan'), ...
                        std(absBounds(Lrts(:,rr)), 'omitnan'), ...
                        sum(Lrts(:,rr))];
                end

            else %if strcmp(options.RTs, 'all')

                % all at once
                % Lsnr = Lsnr & dataTable.DT>1;% & dataTable.DT<6;
                boundSummary_(ss,bb,nn,1,:) = [ ...
                    mean(absBounds(Lsnr), 'omitnan'), ...
                    std(absBounds(Lsnr), 'omitnan'), ...
                    sum(Lsnr)];
            end
        end
    end
end

