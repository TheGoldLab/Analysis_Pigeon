function Figure04_RRvsBound(dataTable, options)
% function Figure04_RRvsBound(dataTable, options)
%
% Figure

arguments
    dataTable
    options.blocks = 1:3
    options.showSimulations = true
    options.simBounds = 0:0.01:0.75;
    options.numSimSubjects  = 100;
    options.figureNumber = 5
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {3,3};
hts     = 5;
[axs,~] = getPLOT_axes(options.num, wid, hts, cols, 1.8, 0.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Collect bound data
%   boundSummary_ .. matrix with dims:
%       subjects
%       blocks
%       snrs
%       rt bins
%       bound avg, bound std, n
boundSummary = getPigeon_boundSummary(dataTable, ...
    'blocks',       options.blocks, ...
    'maxRT',        'all');

% Collect coin data
coinSummary = getPigeon_coinSummary(dataTable, ...
    'blocks',       options.blocks);

%% Plotz
% Colors
wt = 0.99.*ones(1,3);
gr = 0.75.*ones(1,3);
rp = [1 0.9 0.9];

%% Top row is RR vs bound Per block
%
numBlocks = size(boundSummary,2);
for bb = 1:numBlocks

    % Set axes
    axes(axs(bb)); cla reset; hold on;

    % simulate using different bounds 
    if options.showSimulations
        specs.numSubjects = options.numSimSubjects;
        specs.blocks = options.blocks(bb);
        numSimBounds = length(options.simBounds);
        simBoundData = nan(numSimBounds,3);
        for ss = 1:numSimBounds
            if mod(ss,10) == 0
                fprintf('simulating bounds %d of %d\n', ss, numSimBounds)
            end

            simTable = getPigeon_simulatedDataTable(specs, ...
                'boundMean', options.simBounds(ss));
            simCoinSummary = getPigeon_coinSummary(simTable, ...
                'blocks',       options.blocks(bb));
            simBoundData(ss,:) = prctile(simCoinSummary, [5 50 95])./600;
        end

        % Plot as ribbon
        h=patch([options.simBounds flip(options.simBounds)], ...
            [simBoundData(:,1); flip(simBoundData(:,3))]', rp);
        set(h, 'LineStyle', 'none')
        plot(options.simBounds, simBoundData(:,2), '-', ...
            'Color', 'r', 'LineWidth', 2);
    end

    % Summarize bound & RR data
    meanBound = boundSummary(:,bb,1,1,1);
    stdBound = boundSummary(:,bb,1,1,2);
    rewardRate = coinSummary(:,bb,1)./600;

    % plot with horizontal error bars
    plot([meanBound-stdBound meanBound+stdBound]', [rewardRate rewardRate]', ...
         '-', 'Color', gr, 'LineWidth', 0.5);
    plot(meanBound, rewardRate, 'ko', 'MarkerFaceColor', wt);
    title(sprintf('Block %d', bb));
    axis([0 0.7 -0.8 0.4])
    if bb == 1
        xlabel('Bound')
        ylabel('Reward Rate');
    end
end

%% Bottom row is all three bound comparisons
% 
comparisons = [1 2; 2 3; 1 3];
numComparisons = size(comparisons, 1);
for cc = 1:numComparisons

    % Set axes
    axes(axs(3+cc)); cla reset; hold on;

    % Plot it
    plot([0 0.75], [0 0.75], 'k:');
    xs = boundSummary(:,comparisons(cc,1),1,1,1);
    ys = boundSummary(:,comparisons(cc,2),1,1,1);
    plot(xs, ys, 'ko', 'MarkerFaceColor', wt);
    title(sprintf('p=%.2f', ranksum(xs,ys)))
    xlabel(sprintf('block %d', comparisons(cc,1)));
    ylabel(sprintf('block %d', comparisons(cc,2)));    
end
