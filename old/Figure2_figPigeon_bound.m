function Figure2_figPigeon_bound(dataTable, options)
% function figPigeon_bound(dataTable, options)
%
% dataTable columns:
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
arguments
    dataTable
    options.exampleSubject = [1, 2, 3, 4]
    options.block = 2
    options.figureNumber = 2
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {4,4,2};
ht      = 5;
[axs,~] = getPLOT_axes(options.figureNumber, wid, ht, cols, 1.3, 1.5, [], 'Pigeons', true);
set(axs,'Units','normalized');
wt = ones(3,1).*0.99;

%% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);

%% Get summary data per subject
%
% Set up data table
variableNames = ...
    {'DTMedian', 'DT25', 'DT75', 'DTMean', 'DTSEM', ...
    'boundMedian', 'bound25', 'bound75', 'boundMean', 'boundSEM', ...
    'BoundVsDT_B0', 'BoundVsDT_B1', 'BoundVsDT_p', ...
    'DeltaBoundVsBound_B0', 'DeltaBoundVsBound_B1', 'DeltaBoundVsBound_p'}; 
variableTypes = cell(size(variableNames));
[variableTypes{:}] = deal('double');
summaryTable = table('Size', [numSubjects length(variableNames)], ...
    'VariableTypes', variableTypes, ...
    'VariableNames', variableNames);

% Use standard selection criteria
Lg = getPigeon_goodTrialArray(dataTable,'blockIndex',2);

% Loop through the subjects
exampleIndex = 1;
for ss = 1:numSubjects

    % Get subject-specific data
    Lsb = Lg & dataTable.subjectIndex==subjects(ss);

    % Summarize DT
    DT = dataTable.DT(Lsb);
    summaryTable(ss, {'DTMedian', 'DT25', 'DT75', 'DTMean', 'DTSEM'}) = ...
        {prctile(DT,50) prctile(DT,25) prctile(DT,75) mean(DT,'omitnan') nanse(DT)};

    % Summarize bound
    bound = abs(dataTable.bound(Lsb));
    summaryTable(ss, {'boundMedian', 'bound25', 'bound75', 'boundMean', 'boundSEM'}) = ...
        {prctile(bound,50) prctile(bound,25) prctile(bound,75) mean(bound,'omitnan') nanse(bound)};

    % Linear regressions
    xs = ones(sum(Lsb),2);

    % 1. bound vs RT
    xs(:,2) = DT;
    [b,~,~,~,stats] = regress(bound, xs);
    summaryTable(ss, {'BoundVsDT_B0', 'BoundVsDT_B1', 'BoundVsDT_p'}) = ...
        {b(1) b(2) stats(3)};
    
    % 2. deltaBound vs bound
    deltaBound = diff(bound);
    xs(:,2) = bound;
    [b,~,~,~,stats] = regress(deltaBound, xs(1:end-1,:));
    summaryTable(ss, {'DeltaBoundVsBound_B0', 'DeltaBoundVsBound_B1', 'DeltaBoundVsBound_p'}) = ...
        {b(1) b(2) stats(3)};

    % Possibly show example
    if ismember(ss, options.exampleSubject)

        % DT vs Bound
        axes(axs(exampleIndex)); cla reset; hold on;
        scatter(DT, bound, 'ko', ...
            'SizeData', 20, 'MarkerFaceColor', wt);
        lsline
        xlabel('RT');
        if exampleIndex == 1
            ylabel("Bound");
        end

        % deltaBound vs Bound
        axes(axs(exampleIndex+4)); cla reset; hold on;
        scatter(bound(1:end-1), deltaBound, 'ko', ...
            'SizeData', 20, 'MarkerFaceColor', wt);
        lsline
        plot([0 1.0], [0 0], 'k:');
        xlabel("Bound");
        if exampleIndex == 1
            ylabel('\Delta bound');
        end

        exampleIndex = exampleIndex + 1;
    end
end

%% Column 1 is bound vs DT
axes(axs(9)); cla reset; hold on;
errorbar(summaryTable.DTMean, summaryTable.boundMean, ...
    0.5.*summaryTable.boundSEM, 0.5.*summaryTable.boundSEM, ...
    0.5.*summaryTable.DTSEM, 0.5.*summaryTable.DTSEM, ...
    'ko', 'MarkerFaceColor', wt);
xlabel('DT')
ylabel('Bound')

%% Column 2 is slope of BoundVsDT vs slope of DeltaBoundVsBound
axes(axs(10)); cla reset; hold on;
plot(summaryTable.BoundVsDT_B1, summaryTable.DeltaBoundVsBound_B1, ...
    'ko', 'MarkerFaceColor', wt);
xlabel('slope of BoundVsDT')
ylabel('slope of DeltaBoundVsBound')
