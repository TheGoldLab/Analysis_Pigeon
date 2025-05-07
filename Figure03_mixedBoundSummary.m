function Figure03_mixedBoundSummary(dataTable, options)
% function Figure03_mixedBoundSummary(dataTable, options)
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
    options.exampleSubject = 2
    options.block = 2
    options.figureNumber = 3
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {4,3};
hts     = [3 5];
[axs,~] = getPLOT_axes(options.figureNumber, wid, hts, cols, 1.3, 1.5, [], 'Pigeons', true);
set(axs,'Units','normalized');
wt = ones(3,1).*0.99;

%% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);

%% Get summary data per subject
%
% Set up data tables (one for mixed SNR, one for average SNR)
variableNames = ...
    {'DTMedian', 'DT25', 'DT75', 'DTMean', 'DTSEM', ...
    'boundMedian', 'bound25', 'bound75', 'boundMean', 'boundSEM', ...
    'BoundVsDT_B0', 'BoundVsDT_B1', 'BoundVsDT_p', ...
    'DeltaBoundVsBound_B0', 'DeltaBoundVsBound_B1', 'DeltaBoundVsBound_p'}; 

% Use standard selection criteria
Lg = getPigeon_goodTrialArray(dataTable);

% make two tables
tableNames = {'mixedSNR', 'averageSNR'};
for tt = 1:2

    % Make the table
    summaryTables.(tableNames{tt}) = ...
        makeTableOfDoubles(numSubjects, variableNames);
    
    % select appropriate block
    Lb = Lg & dataTable.blockIndex == options.block + (tt-1)*3;

    % Loop through the subjects
    for ss = 1:numSubjects

        % Get subject-specific data
        Lsb = Lb & dataTable.subjectIndex==subjects(ss);

        % Summarize DT
        DT = dataTable.DT(Lsb);
        summaryTables.(tableNames{tt})(ss, {'DTMedian', 'DT25', 'DT75', 'DTMean', 'DTSEM'}) = ...
            {prctile(DT,50) prctile(DT,25) prctile(DT,75) mean(DT,'omitnan') nanse(DT)};

        % Summarize bound
        bound = abs(dataTable.bound(Lsb));
        summaryTables.(tableNames{tt})(ss, {'boundMedian', 'bound25', 'bound75', 'boundMean', 'boundSEM'}) = ...
            {prctile(bound,50) prctile(bound,25) prctile(bound,75) mean(bound,'omitnan') nanse(bound)};

        % Linear regressions
        xs = ones(sum(Lsb),2);

        % 1. bound vs RT
        xs(:,2) = DT;
        [b,~,~,~,stats] = regress(bound, xs);
        summaryTables.(tableNames{tt})(ss, {'BoundVsDT_B0', 'BoundVsDT_B1', 'BoundVsDT_p'}) = ...
            {b(1) b(2) stats(3)};

        % 2. deltaBound vs bound
        deltaBound = diff(bound);
        xs(:,2) = bound;
        [b,~,~,~,stats] = regress(deltaBound, xs(1:end-1,:));
        summaryTables.(tableNames{tt})(ss, {'DeltaBoundVsBound_B0', 'DeltaBoundVsBound_B1', 'DeltaBoundVsBound_p'}) = ...
            {b(1) b(2) stats(3)};

        % Possibly show example
        if ismember(ss, options.exampleSubject)

            % DT vs Bound
            axes(axs((tt-1)*2+1)); cla reset; hold on;
            scatter(DT, bound, 'ko', ...
                'SizeData', 20, 'MarkerFaceColor', wt);
            lsline
            xlabel('RT');
            ylabel('Bound');
            
            % deltaBound vs Bound
            axes(axs((tt-1)*2+2)); cla reset; hold on;
            scatter(bound(1:end-1), deltaBound, 'ko', ...
                'SizeData', 20, 'MarkerFaceColor', wt);
            lsline
            plot([0 1.0], [0 0], 'k:');
            xlabel('Bound');
            ylabel('\Delta bound');
        end
    end
end

%% Column 1 is bound vs DT
axes(axs(end-2)); cla reset; hold on;
errorbar(summaryTables.mixedSNR.boundMean, summaryTables.mixedSNR.DTMean, ...
    0.5.*summaryTables.mixedSNR.DTSEM, 0.5.*summaryTables.mixedSNR.DTSEM, ...
    0.5.*summaryTables.mixedSNR.boundSEM, 0.5.*summaryTables.mixedSNR.boundSEM, ...
    'ko', 'MarkerFaceColor', wt);
xlabel('Bound')
ylabel('DT')

%% Column 2 slope of BoundVsDT average vs mixed
axes(axs(end-1)); cla reset; hold on;
plot([-0.05 0.2], [-0.05 0.2], 'k:')
plot(summaryTables.averageSNR.BoundVsDT_B1, ...
    summaryTables.mixedSNR.BoundVsDT_B1, ...
    'ko', 'MarkerFaceColor', wt);
axis([-0.05 0.2 -0.05 0.2])
title(sprintf('slope of BoundVsDT,p=%.2f', ...
    signrank(summaryTables.averageSNR.BoundVsDT_B1, ...
    summaryTables.mixedSNR.BoundVsDT_B1)))
xlabel('Average SNR')
ylabel('Mixed SNR')

%% Column 3 slope of DeltaBoundVsBound average vs mixed
axes(axs(end)); cla reset; hold on;
plot([-2.5 0], [-1 -1], 'k:')
plot([-1 -1], [-2.5 0], 'k:')
plot(summaryTables.averageSNR.DeltaBoundVsBound_B1, ...
    summaryTables.mixedSNR.DeltaBoundVsBound_B1, ...
    'ko', 'MarkerFaceColor', wt);
axis([-1.5 0 -2.5 0])
title('slope of DeltaBoundVsBound')
xlabel('Average SNR')
ylabel('Mixed SNR')
