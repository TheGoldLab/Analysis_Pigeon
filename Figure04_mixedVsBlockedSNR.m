function Figure03_mixedVsBlockedSNR(dataTableMX, dataTableOL, options)
% function Figure03_mixedVsBlockedSNR(dataTableMX, dataTableBL, options)
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
    dataTableMX
    dataTableOL
    options.showRR = true;
    options.maxRT = 10;
    options.minN = 4;
    options.numRRReps=100;
    options.scaleCoins = true;
    options.block = 2
    options.figureNumber = 4
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {3,3};
hts     = [4 4];
[axs,~] = getPLOT_axes(options.figureNumber, wid, hts, cols, 1.3, 1.5, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Colors
wt = ones(3,1).*0.99;
gr = ones(3,1).*0.5;

%  Set up rt axis for column 3 plotz
rtAxis = 2:options.maxRT;
numRT = length(rtAxis);

% Collect summary data & plot
types = {'MX', 'OL'};
titles = {'mixed', 'blocked'};
for tt = 1:length(types)

    %% Collect data
    thisTable = eval(['dataTable' types{tt}]);

    % Bound summary, collapsed across all RTs
    boundSummary = getPigeon_boundSummary(thisTable, ...
        'blocks',       options.block, ...
        'maxRT',        'all');

    % Coin count summary
    coinSummary = getPigeon_coinSummary(thisTable, ...
        'blocks',       options.block, ...
        'splitBySNR',   false);

    % RR Matrix
    rrBounds = 0.01:0.05:0.8;
    RRMatrix = getPigeon_RRMatrix( ...
        'blockType',    types{tt}, ...
        'blocks',       options.block,  ...
        'bounds',       rrBounds, ...
        'numReps',      options.numRRReps);

    %% Plotz
    % --- first column ---
    %
    % hiSNR vs loSNR bound
    %
    axes(axs((tt-1)*3+1)); cla reset;  hold on;
    imagesc(rrBounds, flip(rrBounds), flipud(RRMatrix(:,:,1)));
    set(gca,'YDir', 'normal')
    axis(rrBounds([1 end 1 end]))
    plot([0 1], [0 1], 'k:')

    if options.scaleCoins

        % scale marker size by coin count
        maxCoins = max(coinSummary);
        for ss = 1:size(boundSummary,1)
            h = plot(boundSummary(ss,1,1,1,1), boundSummary(ss,1,2,1,1), 'ko', ...
                'MarkerFaceColor', wt, 'MarkerSize', abs(coinSummary(ss)/maxCoins*12+0.1));
            if coinSummary(ss) <= 0
                set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'r')
            end
        end
    else

        % Show using same marker
        plot(boundSummary(:,1,1,1,1), boundSummary(:,1,2,1,1), 'ko', 'MarkerFaceColor', wt);
    end

    %plot(medianBounds{xx}(:,bb,1), medianBounds{xx}(:,bb,2), 'ko', ...
    %    'MarkerFaceColor', wtc);
    if tt == 2
        xlabel('loSNR bound');
        ylabel('hiSNR bound');
    end

    title(sprintf('%s: p=%.3f', titles{tt}, signrank(boundSummary(:,1,1,1,1), boundSummary(:,1,2,1,1))))
    axis([0 0.75 0 0.75])

    % --- second column ---
    %
    % hiSNR vs loSNR RT (with std error bars)

    % performance summary is matrix with dimensions:
    %   1. subject index
    %   2. block index
    %   3. SNR index
    %   4. accuracy / RT mean / RT STD / n
    performanceSummary = getPigeon_performanceSummary(thisTable, ...
        'blocks',       options.block);

    % plotz
    axes(axs((tt-1)*3+2)); cla reset;  hold on;
    plot([0 10], [0 10], 'k:')
    % plot std errobars
    meanLo = performanceSummary(:,1,1,2);
    stdLo  = performanceSummary(:,1,1,3);
    meanHi = performanceSummary(:,1,2,2);
    stdHi  = performanceSummary(:,1,2,3);
    plot([meanLo meanLo]', [meanHi-stdHi meanHi+stdHi]', '-', 'Color', gr)
    plot([meanLo-stdLo meanLo+stdLo]', [meanHi meanHi]', '-', 'Color', gr)
    plot(meanLo, meanHi, 'ko', 'MarkerFaceColor', wt);
    axis([0 10 0 10])

    % --- third column ---
    %
    % hiSNR minus loSNR bound per RT

    % Bound summary, per RT (up to max)
    boundSummaryPerRT = getPigeon_boundSummary(thisTable, ...
        'blocks',       options.block, ...
        'maxRT',        options.maxRT);

    % plotz
    axes(axs((tt-1)*3+3)); cla reset;  hold on;
    plot(rtAxis([1 end]), [0 0], 'k:');

    % Compute as difference in bound for two SNRs (per RT)
    ys = squeeze(diff(boundSummaryPerRT(:,1,1:2,rtAxis,1),[],3));
    ns = squeeze(boundSummaryPerRT(:,1,1:2,rtAxis,3));
    Ln = squeeze(ns(:,1,:))>=options.minN & squeeze(ns(:,2,:))>=options.minN;
    xax = repmat(rtAxis,size(ys,1),1) + normrnd(0,0.1,size(ys));
    plot(xax(Ln), ys(Ln), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', wt);
    for rr = 1:length(rtAxis)
        diffs = ys(Ln(:,rr),rr);
        diffs = diffs(isfinite(diffs));
        if any(diffs)
            h=plot(rtAxis(rr), median(diffs), 'ro', 'MarkerSize', 7);
            if signrank(diffs) < 0.05
                set(h, 'MarkerFaceColor', 'r');
            end
        end
    end
    axis([rtAxis([1 end]) -0.5 0.5])
end
