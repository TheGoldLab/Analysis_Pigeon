function Figure4_figPigeon_boundBlockTransition(dataTable, num)
% function figPigeon_boundByTime(dataTable, num)
%
% Figure:

arguments
    dataTable
    num = 4
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {1,1};
hts     = 9;
[axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
blocks = nonanunique(dataTable.blockIndex);
numBlocks = length(blocks);
% numTransitions = numBlocks-1;
SNRs = nonanunique(abs(dataTable.snr));
numSNRs = length(SNRs);

numSteps = 20;
bData = nan(numSteps*numSNRs*2, numSubjects, numSNRs); % last is mean/sem/n bound

for ss = 1:numSubjects
    for ns = 1:numSNRs
        tData = [];
        for bb = 1:numBlocks-1
            Lsn = dataTable.subjectIndex==subjects(ss) & abs(dataTable.snr) == SNRs(ns);
            Lsb =  Lsn & dataTable.blockIndex==blocks(bb);
            Lsbn = Lsn & dataTable.blockIndex==blocks(bb+1);
            bounds = abs(dataTable.bound(Lsb));
            % bounds = (bounds - mean(bounds))/(std(bounds)/length(bounds));
            nBounds = abs(dataTable.bound(Lsbn));
            % nBounds = (nBounds - mean(nBounds))/(std(nBounds)/length(nBounds));
            nData = [(bounds(end-numSteps:end) - mean(bounds(end-numSteps:end)))/(std(bounds(end-numSteps:end))/sqrt(length(bounds(end-numSteps:end))));
                (nBounds(1:numSteps) - mean(nBounds(1:numSteps)))/(std(nBounds(1:numSteps))/sqrt(length(nBounds(1:numSteps))))];
            nData(round(size(nData,1)/2),:,:) = [];
            tData = [tData; nData];
        end
        bData(:,ss,ns) = tData;
    end
end



%% Plotz
% Per block
%gry = 0.9.*ones(1,3);
wt = 0.99.*ones(1,3);


for ns = 1:numSNRs
    % Set axes
    axes(axs(ns)); cla reset; hold on;

    meanBySubject = mean(bData(:,:,ns), 2, "omitnan");
    semBySubject = std(bData(:,:,ns), 0, 2, 'omitnan')/sqrt(numSubjects);

    % trialNums = 1:80;
    trialNums = (1:81);
    meanBySubject = [meanBySubject(1:round(size(meanBySubject,1)/2)); NaN; meanBySubject(round(size(meanBySubject,1)/2)+1:end)];
    semBySubject  = [semBySubject(1:round(size(semBySubject,1)/2)); NaN; semBySubject(round(size(semBySubject,1)/2)+1:end)];
    upper = meanBySubject + semBySubject;
    lower = meanBySubject - semBySubject;


    segments = {(1:40), (42:81)};
    colorshade = [0.7 0.7 0.9];
    for s = 1:numel(segments)
        seg = segments{s};
        upper = meanBySubject(seg) + semBySubject(seg);
        lower = meanBySubject(seg) - semBySubject(seg);
        fill([trialNums(seg), fliplr(trialNums(seg))], ...
             [upper', fliplr(lower')], ...
             colorshade, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    end

    plot(trialNums, meanBySubject, 'k', "LineWidth", 2);
    % grandMean = mean(meanBySubject, "omitnan");

    xline(20, 'k', "Label", "Block 1 \rightarrow 2"); ylim([-5 5]);
    xline(40, 'k', "Label", "Trial +20, Block 2", "LabelHorizontalAlignment","left"); ylim([-5 5]);
    xline(42, 'k', "Label", "Trial -20, Block 3", "LabelHorizontalAlignment","right"); ylim([-5 5]);
    xline(60, 'k', "Label", "Block 2 \rightarrow 3"); ylim([-5 5]);
    % Labels
    xticks([0 10 20 30 51 61 71]);
    xticklabels({'-20', '-10', '0', '+10', '-10', '0', '+10'});

    % plot([-numSteps numSteps], [grandMean grandMean], 'k:');
    ylabel("Normalized Bound");
    if bb == 2
        xlabel("Trial Number (relative to block transition)")
    end
    xlim([trialNums(1) trialNums(end)]);
end


% function Figure4_figPigeon_boundBlockTransition(dataTable, num)
% % function figPigeon_boundByTime(dataTable, num)
% %
% % Figure:
% 
% arguments
%     dataTable
%     num = 4
% end
% 
% %% Set up figure
% %
% wid     = 17.6; % total width
% cols    = {1,1};
% hts     = 9;
% [axs,~] = getPLOT_axes(num, wid, hts, cols, 1.3, 0.5, [], 'Pigeons', true);
% set(axs,'Units','normalized');
% 
% % Collect data per subject/block
% subjects = nonanunique(dataTable.subjectIndex);
% numSubjects = length(subjects);
% blocks = nonanunique(dataTable.blockIndex);
% numBlocks = length(blocks);
% numTransitions = numBlocks-1;
% SNRs = nonanunique(abs(dataTable.snr));
% numSNRs = length(SNRs);
% 
% numSteps = 20;
% bData = nan(numSteps*numSNRs*2, numSubjects, numSNRs); % last is mean/sem/n bound
% 
% for ss = 1:numSubjects
%     for ns = 1:numSNRs
%         tData = [];
%         for bb = 1:numBlocks-1
%             Lsn = dataTable.subjectIndex==subjects(ss) & abs(dataTable.snr) == SNRs(ns);
%             Lsb =  Lsn & dataTable.blockIndex==blocks(bb);
%             Lsbn = Lsn & dataTable.blockIndex==blocks(bb+1);
%             bounds = abs(dataTable.bound(Lsb));
%             % bounds = (bounds - mean(bounds))/(std(bounds)/length(bounds));
%             nBounds = abs(dataTable.bound(Lsbn));
%             % nBounds = (nBounds - mean(nBounds))/(std(nBounds)/length(nBounds));
%             nData = [(bounds(end-numSteps:end) - mean(bounds(end-numSteps:end)))/(std(bounds(end-numSteps:end))/sqrt(length(bounds(end-numSteps:end))));
%                 (nBounds(1:numSteps) - mean(nBounds(1:numSteps)))/(std(nBounds(1:numSteps))/sqrt(length(nBounds(1:numSteps))))];
%             nData(round(size(nData,1)/2),:,:) = [];
%             tData = [tData; nData];
%         end
%         bData(:,ss,ns) = tData;
%     end
% end
% 
% 
% 
% %% Plotz
% % Per block
% %gry = 0.9.*ones(1,3);
% wt = 0.99.*ones(1,3);
% 
% 
% for ns = 1:numSNRs
%     % Set axes
%     axes(axs(ns)); cla reset; hold on;
% 
%     meanBySubject = mean(bData(:,:,ns), 2, "omitnan");
%     semBySubject = std(bData(:,:,ns), 0, 2, 'omitnan')/sqrt(numSubjects);
% 
%     trialNums = 1:80;
%     upper = meanBySubject + semBySubject;
%     lower = meanBySubject - semBySubject;
% 
% 
% 
%     fill([trialNums, fliplr(trialNums)], [upper', fliplr(lower')], [0.7 0.7 0.9], ...
%         'EdgeColor', 'none', 'FaceAlpha', 0.8);
% 
%     plot(trialNums, meanBySubject, 'k', "LineWidth", 2);
%     grandMean = mean(meanBySubject, "omitnan");
% 
%     xline(20, 'k', "Label", "Block 1 \\rightarrow 2"); ylim([-5 5]);
%     % xline(40, 'k', "Label", "Block Transition"); ylim([-5 5]);
%     xline(60, 'k', "Label", "Block 2 \\rightarrow 3"); ylim([-5 5]);
%     % Labels
%     xticks([10 30 50 70]);
%     xticklabels({'-10', '10', '-10', '10'});
% 
%     plot([-numSteps numSteps], [grandMean grandMean], 'k:');
%     title(sprintf("Block %d \\rightarrow %d", bb, bb+1));
%     ylabel("Normalized Bound");
%     if bb == 2
%         xlabel("Trial Number (relative to block transition)")
%     end
%     xlim([trialNums(1) trialNums(end)]);
% end
