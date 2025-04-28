function accuracy_ = getPigeon_accuracySummary(dataTable)
% function accuracy_ = getPigeon_accuracySummary(dataTable)
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

subjectIndices = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjectIndices);
numBlocks = 3;
numSNRs = 2;
absSNR = abs(dataTable.snr);

% Check for task type (OL or MX)
% if isscalar(nonanunique(absSNR(dataTable.blockIndex == 1)))
%     Lbound = [...
%         dataTable.blockIndex==1 | dataTable.blockIndex==4, ...
%         dataTable.blockIndex==2 | dataTable.blockIndex==5, ...
%         dataTable.blockIndex==3 | dataTable.blockIndex==6, ...
%         ];
% else
%     Lbound = [...
%         dataTable.blockIndex==1, ...
%         dataTable.blockIndex==2, ...
%         dataTable.blockIndex==3, ...
%         ];
% end

% Set up the data matrix
accuracy_ = zeros(numSubjects,numBlocks,numSNRs,2); % last is %cor/n

% Loop through the subjects
for ss = 1:length(subjectIndices)
    Ls = dataTable.subjectIndex==subjectIndices(ss);

    % Loop through the blocks
    for bb = 1:numBlocks

        Lsb = Ls & dataTable.blockIndex==bb; %Lbound(:,bb);
        snrs = nonanunique(absSNR(Lsb));

        for nn = 1:length(snrs)
            Lsnr = Lsb & absSNR==snrs(nn);
            accuracy_(ss,bb,nn,:) = [ ...
                sum(dataTable.correct(Lsnr))./sum(Lsnr), sum(Lsnr)];
        end
    end
end