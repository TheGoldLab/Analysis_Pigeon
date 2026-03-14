function Figure03_boundSummary(dataTable, options)
% function Figure03_boundSummary(dataTable, options)
%
% dataTable is dataTableMX, with columns:
%   1. subjectIndex
%   2. blockIndex (1:3 are mixed SNR, 4:6 are avg SNR)
%   3. trialNumber
%   4. bound (pigeon position in OL)
%   5. RT (response time, including ndt, in number of steps)
%   6. DT (response time, not including ndt, in number of steps)
%   7. choice (0/1)
%   8. correct (correct=1,error=0)
%   9. coinCount
%   10. snr
%   11. steps (cell)
arguments
    dataTable
    options.exampleSubject = 16;
    options.blockIndex = 2;
    options.figureNumber = 3;
end

%% Set up figure
%
wid     = 8.5; % total width
cols    = {1,1,2,1};
hts     = [4 4 2 4];
[axs,~] = getPLOT_axes(options.figureNumber, wid, hts, cols, 1.3, 1.5, [], 'Pigeons', true);
%set(axs,'Units','normalized');

% colors
gr = ones(3,1).*0.5;
wh = ones(3,1).*0.99;

%% Collect data per subject/snr
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
numTimeBins = 10;
abounds = abs(dataTable.bound); % absolute values of bounds
zbounds = abounds; % z-scored (abs) bounds, per subject
dbounds = [nan; diff(abounds)]; % delta bounds
boundByDTData = nan(numSubjects,numTimeBins);
boundVsDTData = nan(numSubjects,4);  % rows are mean/sem bound, mean/sem DT
boundBySideData = nan(numSubjects,2); % diff of medians, p for diff
regressionData = nan(numSubjects, 2); % boundVsDT, deltaBoundVsBound

aSNRs = abs(dataTable.snr);
Lg = isfinite(dataTable.bound) & aSNRs==min(aSNRs) & dataTable.blockIndex==options.blockIndex;
Lgr = dataTable.DT >= 2;
Lrl = [Lg&dataTable.choice==0 Lg&dataTable.choice==1];
for ss = 1:numSubjects
    Ls = Lg & dataTable.subjectIndex == subjects(ss);

    % Collect z-scored bound per DT
    zbounds(Ls) = zscore(abounds(Ls));
    for dd = 1:numTimeBins
        Ld = Ls & dataTable.DT==dd;
        if any(Ld) %sum(Ld) > 2
            boundByDTData(ss,dd) = mean(zbounds(Ld));
        end
    end

    % Compare left/right
    leftBounds = -dataTable.bound(Ls&Lrl(:,1));
    rightBounds = dataTable.bound(Ls&Lrl(:,2));
    boundBySideData(ss,:) = [ ...
        median(rightBounds)-median(leftBounds), ...
        ranksum(rightBounds, leftBounds)];

    % Possibly show example
    if ss == options.exampleSubject
        axes(axs(1)); cla reset; hold on;
        xs = dataTable.DT(Ls);
        ys = abounds(Ls);
        plot(xs, ys, 'ko', 'MarkerFaceColor', wh);
        for dd = 1:15
            plot(dd,mean(ys(xs==dd)),'ro', 'MarkerFaceColor', 'r')
        end
        axis([1 10 0 1]);
        ylabel('Bound')
    end

    % Collect mean/sem bound and DT
    boundVsDTData(ss,:) = [ ...
        mean(abounds(Ls), 'omitnan'), ...
        sem(abounds(Ls)), ...
        mean(dataTable.DT(Ls), 'omitnan'), ...
        sem(dataTable.DT(Ls)), ...
        ];

    % regression bound vs RT
    Lreg = Ls & Lgr;
    if sum(Lreg) > 5
        xs = ones(sum(Lreg),2);
        xs(:,2) = dataTable.DT(Lreg);
        b = regress(abounds(Lreg), xs);
        regressionData(ss,1) = b(2);
    end

    % regression deltaBound vs bound
    if sum(Ls) > 5
        xs = ones(sum(Ls),2);
        xs(:,2) = abounds(Ls);
        deltaBounds = dbounds(Ls);
        b = regress(deltaBounds(2:end), xs(2:end,:));
        regressionData(ss,2) = b(2);
    end
end

%% Second row is bound summary
axes(axs(2)); cla reset; hold on;
plot(boundByDTData(:,:)', '-', 'Color', gr);
plot(mean(boundByDTData(:,:), 'omitnan'), 'r-', 'LineWidth', 2)
axis([1 10 -2 2])
xlabel('DT (steps)')
ylabel('Final position (z-score)')

%% Third row is regression histograms
xlabels = {'Slope Bound Vs DT', 'Slope Delta Bound Vs Bound'};
hax = {-0.04:0.008:0.04, 0:0.15:2};
offset = [0 1];
for hh = 1:2 % For each histogram
    axes(axs(2+hh)); cla reset; hold on;
    histogram(regressionData(:,hh), hax{hh});
    fprintf('median [IQR]=%.2f [%.2f %.2f], p=%.3f\n', ...
        prctile(regressionData(:,hh), 50), ...
        prctile(regressionData(:,hh), 25), ...
        prctile(regressionData(:,hh), 75), ...
        signrank(regressionData(:,hh)-offset(hh)))
    plot(prctile(regressionData(:,hh), 50).*[1 1], [0 20], 'r:')
    xlabel(xlabels{hh})
    if hh==1
        ylabel('Count')
    end
end

%% Fourth row is bound vs DT per subject
% plot data
axes(axs(5)); cla reset; hold on;
errorbar(boundVsDTData(:,3), boundVsDTData(:,1), ...
    0.5.*boundVsDTData(:,2), 0.5.*boundVsDTData(:,2), ...
    0.5.*boundVsDTData(:,4), 0.5.*boundVsDTData(:,4), ...
    'ko', 'MarkerFaceColor', wh);
[R,P] = corr(boundVsDTData(:,3), boundVsDTData(:,1),'type', 'Spearman');
fprintf('Spearman=%.3f,p=%.3f\n',R,P)
axis([0 10 0 0.75])
xlabel('DT (steps)')
ylabel('Bound')

% 
% 
% % 
% % 
% % 
% % %% Get summary data per subject
% % %
% % % Set up data tables (one for mixed SNR, one for average SNR)
% % variableNames = ...
% %     {'DTMedian', 'DT25', 'DT75', 'DTMean', 'DTSEM', ...
% %     'boundMedian', 'bound25', 'bound75', 'boundMean', 'boundSEM', ...
% %     'BoundVsDT_B0', 'BoundVsDT_B1', 'BoundVsDT_p', ...
% %     'DeltaBoundVsBound_B0', 'DeltaBoundVsBound_B1', 'DeltaBoundVsBound_p'};
% % 
% % % Use standard selection criteria
% % Lg = getPigeon_goodTrialArray(dataTable, 'DT', 0);
% % 
% % % make two tables, average and mixed
% % tableNames = {'mixedSNR', 'averageSNR'};
% % for tt = 1:2
% % 
% %     % Make the table
% %     summaryTables.(tableNames{tt}) = ...
% %         makeTableOfDoubles(numSubjects, variableNames);
% % 
% %     % select appropriate block
% %     Lb = Lg & dataTable.blockIndex == options.block + (tt-1)*3;
% % 
% %     % Loop through the subjects
% %     for ss = 1:numSubjects
% % 
% %         % Get subject-specific data
% %         Lsb = Lb & dataTable.subjectIndex==subjects(ss);
% % 
% %         % Summarize DT
% %         DT = dataTable.DT(Lsb);
% %         summaryTables.(tableNames{tt})(ss, {'DTMedian', 'DT25', 'DT75', 'DTMean', 'DTSEM'}) = ...
% %             {prctile(DT,50) prctile(DT,25) prctile(DT,75) mean(DT,'omitnan') nanse(DT)};
% % 
% %         % Summarize bound
% %         Lsb = Lb & dataTable.subjectIndex==subjects(ss);
% %         bound = abs(dataTable.bound(Lsb));
% %         summaryTables.(tableNames{tt})(ss, {'boundMedian', 'bound25', 'bound75', 'boundMean', 'boundSEM'}) = ...
% %             {prctile(bound,50) prctile(bound,25) prctile(bound,75) mean(bound,'omitnan') nanse(bound)};
% % 
% %         %% Linear regressions
% % 
% %         % 1. bound vs RT
% %         Lreg = dataTable.DT(Lsb)>=3;
% %         if sum(Lreg) > 5
% %             xs = ones(sum(Lreg),2);
% %             xs(:,2) = DT(Lreg);
% %             [b,~,~,~,stats] = regress(bound(Lreg), xs);
% %             summaryTables.(tableNames{tt})(ss, {'BoundVsDT_B0', 'BoundVsDT_B1', 'BoundVsDT_p'}) = ...
% %                 {b(1) b(2) stats(3)};
% %         end
% % 
% %         % 2. deltaBound vs bound
% %         if length(bound) > 5
% %             deltaBound = diff(bound);
% %             xs = ones(sum(Lsb),2);
% %             xs(:,2) = bound;
% %             [b,~,~,~,stats] = regress(deltaBound, xs(1:end-1,:));
% %             summaryTables.(tableNames{tt})(ss, {'DeltaBoundVsBound_B0', 'DeltaBoundVsBound_B1', 'DeltaBoundVsBound_p'}) = ...
% %                 {b(1) b(2) stats(3)};
% %         end
% % 
% %         %% Possibly show example
% %         if ismember(ss, options.exampleSubject)
% % 
% %             % DT vs Bound
% %             axes(axs((tt-1)*2+1)); cla reset; hold on;
% %             scatter(DT(Lreg), bound(Lreg), 'ko', ...
% %                 'SizeData', 20, 'MarkerFaceColor', wt);
% %             lsline
% %             scatter(DT, bound, 'ko', ...
% %                 'SizeData', 20, 'MarkerFaceColor', wt);
% %             xlabel('RT');
% %             ylabel('Bound');
% % 
% %             % deltaBound vs Bound
% %             axes(axs((tt-1)*2+2)); cla reset; hold on;
% %             scatter(bound(1:end-1), deltaBound, 'ko', ...
% %                 'SizeData', 20, 'MarkerFaceColor', wt);
% %             lsline
% %             plot([0 1.0], [0 0], 'k:');
% %             xlabel('Bound');
% %             ylabel('\Delta bound');
% %         end
% %     end
% % end
% 
% %% Column 1 is bound vs DT
% axes(axs(end-2)); cla reset; hold on;
% 
% % Possibly plot CIs from sims
% if options.showSimBounds
%     simBounds = 0.01:0.01:0.5;
%     numBounds = length(simBounds);
%     simBoundData = nan(numBounds, 3); % 2.5/50/97.5 prctiles
%     specs.blocks = 2;
%     specs.numSubjects = 1;
%     for bb = 1:numBounds
%         simTable = getPigeon_simulatedDataTable(specs, ...
%             'generativeMean',   [0.05 0.15], ...
%             'numTrials',        100000, ...
%             'NDTMin',           0, ...
%             'NDTmax',           0, ...
%             'boundMean',        simBounds(bb), ...
%             'boundSTD',         0.1);
%         simBoundData(bb,:) = prctile(simTable.DT(simTable.blockIndex==2), [25 50 75]);
%     end
%     % plot median/IQR from sims
%     plot(simBounds, simBoundData(:,1), 'r-')
%     % plot(simBounds, simBoundData(:,2), 'r-', 'LineWidth', 2)
%     plot(simBounds, simBoundData(:,3), 'r-')
% end
% 
% % plot data
% errorbar(summaryTables.mixedSNR.boundMean, summaryTables.mixedSNR.DTMean, ...
%     0.5.*summaryTables.mixedSNR.DTSEM, 0.5.*summaryTables.mixedSNR.DTSEM, ...
%     0.5.*summaryTables.mixedSNR.boundSEM, 0.5.*summaryTables.mixedSNR.boundSEM, ...
%     'ko', 'MarkerFaceColor', wt);
% xlabel('Bound')
% ylabel('DT')
% axis([0 0.75 0 10])
% 
% %% Column 2 slope of BoundVsDT average vs mixed
% axes(axs(end-1)); cla reset; hold on;
% plot([-0.2 0.2], [-0.2 0.2], 'k:')
% plot([-0.2 0.2], [0 0], 'k:')
% plot([0 0], [-0.2 0.2], 'k:')
% plot(summaryTables.averageSNR.BoundVsDT_B1, ...
%     summaryTables.mixedSNR.BoundVsDT_B1, ...
%     'ko', 'MarkerFaceColor', wt);
% axis([-0.05 0.15 -0.05 0.15])
% title(sprintf('slope of BoundVsDT,\np_c=%.2f,p_a=%.2f,p_m=%.2f', ...
%     signrank(summaryTables.averageSNR.BoundVsDT_B1, ...
%     summaryTables.mixedSNR.BoundVsDT_B1), ...
%     signrank(summaryTables.averageSNR.BoundVsDT_B1), ...
%     signrank(summaryTables.mixedSNR.BoundVsDT_B1)));
% xlabel('Average SNR')
% ylabel('Mixed SNR')
% 
% %% Column 3 slope of DeltaBoundVsBound average vs mixed
% axes(axs(end)); cla reset; hold on;
% plot([-2.5 0], [-1 -1], 'k:')
% plot([-1 -1], [-2.5 0], 'k:')
% plot(summaryTables.averageSNR.DeltaBoundVsBound_B1, ...
%     summaryTables.mixedSNR.DeltaBoundVsBound_B1, ...
%     'ko', 'MarkerFaceColor', wt);
% axis([-1.5 0 -2.5 0])
% title(sprintf('slope of DeltaBoundVsBound\np=%.2f', ...
%     signrank(summaryTables.averageSNR.DeltaBoundVsBound_B1+1), ...
%     signrank(summaryTables.mixedSNR.DeltaBoundVsBound_B1+1)));
% xlabel('Average SNR')
% ylabel('Mixed SNR')
% 
% 
% 