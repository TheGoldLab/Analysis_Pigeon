function Figure05_RRvsBound(dataTable, options)
% function Figure05_RRvsBound(dataTable, options)
%
% Figure

arguments
    dataTable
    options.blocks = 1:3;
    options.showSimulations = true;
    options.simBounds = 0:0.01:0.75;
    options.generativeMeans = [0.05 0.15];
    options.numSimSubjects  = 100;
    options.figureNumber = 5
end

%% Set up figure
%
wid     = 14; % total width
cols    = {4,4,4,4};
hts     = 3;
[axs,~] = getPLOT_axes(options.figureNumber, wid, hts, cols, 1.8, 0.5, [], 'Pigeons', true);
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

% Collect coin data (per subject, block, snr)
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

    % For each SNR
    for rr = 1:2

        % Summarize bound & RR data
        meanBound = boundSummary(:,bb,rr,1,1);
        stdBound = boundSummary(:,bb,rr,1,2);
        rewardRate = coinSummary(:,bb,rr)./600;

        % Plot RR vs bound per subject
        axes(axs((bb-1)*4+(rr-1)*2+1)); cla reset; hold on;

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
                    'generativeMean', options.generativeMeans(rr), ...
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

        % plot with horizontal error bars
        plot([meanBound-stdBound meanBound+stdBound]', [rewardRate rewardRate]', ...
            '-', 'Color', gr, 'LineWidth', 0.5);
        plot(meanBound, rewardRate, 'ko', 'MarkerFaceColor', wt);
        title(sprintf('Block %d', bb));
        axis([0 0.7 -0.8 0.4])
        if bb == 3
            xlabel('Bound')
            ylabel('Reward Rate');
        end
    end
end


%% Add bound comparisons 1/2 and 2/3
%
comparisons = [1 2; 2 3];
for cc = 1:size(comparisons, 1)
    for rr = 1:2

        % Set axes
        axes(axs((cc-1)*4+(rr-1)*2+2)); cla reset; hold on;

        % Plot it
        plot([0 0.75], [0 0.75], 'k:');
        xs = boundSummary(:,comparisons(cc,1),rr,1,1);
        ys = boundSummary(:,comparisons(cc,2),rr,1,1);
        plot(xs, ys, 'ko', 'MarkerFaceColor', wt);
        %title(sprintf('p=%.2f', ranksum(xs,ys)))
        if cc==2
            xlabel(sprintf('block %d', comparisons(cc,1)));
            ylabel(sprintf('block %d', comparisons(cc,2)));
        end
        [R,P] = corr(xs, ys, 'Type', 'Spearman');
        fprintf('block %d vs %d: median diff p=%.2f, corr %.2f (p=%.2f)\n', ...
            comparisons(cc,1), comparisons(cc,2), ranksum(xs,ys), R, P)
    end
end

%% Exponential fits to transition from block 1 to 2
% Set up data matrix, fit parameters
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);
aSNRs = abs(dataTable.snr);
SNRs = nonanunique(aSNRs);
numSNRs = length(SNRs);
bounds = abs(dataTable.bound);
Lg = isfinite(dataTable.bound);
Lblock = [dataTable.blockIndex==1 dataTable.blockIndex==2];
tauData = nan(numSubjects,2); % Save taus for each snr, transition 1->2
fitFcn = @(b,x) b(1)-b(2).*exp(-x./b(3));
fitOpts = optimset('MaxFunEvals', 10000);

% Loop through the subjects,snrs
for ss = 1:numSubjects
    Lsub = Lg & dataTable.subjectIndex==subjects(ss);
    for rr = 1:numSNRs
        Lsnr = Lsub & aSNRs==SNRs(rr);
        Lb1 = Lblock(:,1) & Lsnr;
        Lb2 = Lblock(:,2) & Lsnr;
        if sum(Lb1) > 10 && sum(Lb2) > 10 && ranksum(bounds(Lb1), bounds(Lb2)) < 0.05

            % Fit tau
            tax = (0:sum(Lb2)-1)';
            theseBounds = bounds(Lb2);
            b0 = [mean(theseBounds(1:3))-1 mean(theseBounds(end-10:end))./mean(theseBounds(1:10)) 15];

            fits = fmincon(@(b) sum((theseBounds-fitFcn(b,tax)).^2), ...
                b0,[],[],[],[],[-3 -10 1],[3 10 100],[], fitOpts);
            % cla reset; hold on;
            % plot(tax, theseBounds, 'ko');
            % plot(tax, fitFcn(fits, tax), 'r-')
            % title(sprintf('%.2f', fits(3)))
            % r = input('next')
            tauData(ss,rr) = fits(3);
        end
    end
end
axes(axs(13)); cla reset; hold on;
plot([0 100], [0 100], 'k:')
plot(tauData(:,1), tauData(:,2), 'ko')
Lg = isfinite(tauData(:,1)) & isfinite(tauData(:,2));
[R,P] = corr(tauData(Lg,1), tauData(Lg,2), 'Type', 'Spearman');
title(sprintf('r=%.2f, p=%.2f', R, P))
xlabel('Tau, low SNR')
ylabel('Tau, high SNR')

set(axs([10 12 14:16]), 'Visible', 'off')