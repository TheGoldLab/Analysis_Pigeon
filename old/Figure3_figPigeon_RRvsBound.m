function Figure3_figPigeon_RRvsBound(dataTable, num)
% function figPigeon_boundByTime(dataTable, num)
%
% Figure:

arguments
    dataTable
    num = 3
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {3,3};
hts     = 5;
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.8, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);
SNRs = nonanunique(abs(dataTable.snr));
numSNRs = length(SNRs);

bData = nan(numSubjects,numBlocks, numSNRs,3); % save [mean(bound) std(bound) rewardRate] for each subject/block

Lg = dataTable.DT>2 & dataTable.trialNumber > 10;
for bb = 1:numBlocks
    for ns = 1:numSNRs
        for ss = 1:numSubjects
            Lsb = Lg & abs(dataTable.snr) == SNRs(ns) & ...
                dataTable.subjectIndex==subjects(ss) & ...
                dataTable.blockIndex==blocks(bb);
            bounds = abs(dataTable.bound(Lsb));
            coinCount = dataTable.coinCount(Lsb);
            if ~isempty(coinCount)
                bData(ss,bb,ns,:) = [mean(bounds) std(bounds)/sqrt(length(bounds)) coinCount(end)];
            end
        end
    end
end

%% Plotz
% Per block
%gry = 0.9.*ones(1,3);
wt = 0.99.*ones(1,3);

for ns = 1:numSNRs
    for bb = 1:numBlocks
    
        % Set axes
        axes(axs(3*(ns-1)+bb)); cla reset; hold on;
    
        meanBound = bData(:,bb,ns,1);
        semBound = bData(:,bb,ns,2);
        rewardRate = bData(:,bb,ns,3)/600;
    
        scatter(meanBound, rewardRate, 'ko', 'MarkerFaceColor', wt);
        title(sprintf("Block %d", bb));
        ylabel("Reward Rate");
        if bb == 2
            xlabel("Bound")
        end
    end
end
