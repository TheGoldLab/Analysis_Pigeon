function coinCounts_ = getPigeon_coinSummary(dataTable)
% function coinCounts_ = getPigeon_coinSummary(dataTable)
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
blockIndices = nonanunique(dataTable.blockIndex);
numBlocks = length(blockIndices);
coinCounts_ = zeros(numSubjects,numBlocks);

for ss = 1:length(subjectIndices)
    Ls = dataTable.subjectIndex==subjectIndices(ss);
    for bb = 1:length(blockIndices)
        Lsb = Ls & dataTable.blockIndex == blockIndices(bb);
        coinCounts_(ss,bb) = dataTable.coinCount(find(Lsb, 1, 'last'));
    end
end