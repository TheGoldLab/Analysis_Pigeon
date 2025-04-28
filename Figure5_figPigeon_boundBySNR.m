function Figure5_figPigeon_boundBySNR(dataTableMX, dataTableOL, options)
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
    dataTableMX
    dataTableOL
    options.simDataTables={};
    options.scaleCoins=true;
    options.showRR=true;
    options.numRRReps=100;
end

%% Set up figure
wid     = 17.6; % total width
cols    = {3,3};
hts     = [3.5 3.5];
[axs,~] = getPLOT_axes(1, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
set(axs,'Units','normalized');

% matrices with dims subjects,blocks,snrs,rts,data (bound avg/std/n):
bounds = { ...
    getPigeon_boundSummary(dataTableMX, 'maxRT', 'all'), ...
    getPigeon_boundSummary(dataTableOL, 'maxRT', 'all')};    
coinCounts = { ...
    getPigeon_coinSummary(dataTableMX), ...
    getPigeon_coinSummary(dataTableOL)};

% Possibly standardize by simulated data
if ~isempty(options.simDataTables)
    for xx = 1:2
        simBounds = getPigeon_boundSummary(options.simDataTables{xx}, 'maxRT', 'all');
        bounds{xx}(:,:,2,:,1) = bounds{xx}(:,:,2,:,1) - diff(simBounds(:,:,:,:,1),[],3);
    end
end

%% Plotz
%
% Possibly show RR function
if options.showRR
    getPigeon_RRMatrix('blockType', 'MX', 'axes', axs(1:3), 'numReps', options.numRRReps);
    getPigeon_RRMatrix('blockType', 'OL', 'axes', axs(4:6), 'numReps', options.numRRReps);
end

titles = {'mixed', 'blocked'};
wtc = ones(3,1).*0.99;
for bb = 1:3
    for xx = 1:2
        axes(axs((xx-1)*3+bb)); hold on;
        plot([0 1], [0 1], 'k:', 'LineWidth', 3)

        if options.scaleCoins

            % scale marker size by coin count
            maxCoins = max(coinCounts{xx}(:,bb));
            for ss = 1:size(bounds{xx},1)
                h = plot(bounds{xx}(ss,bb,1,1,1), bounds{xx}(ss,bb,2,1,1), 'ko', ...
                    'MarkerFaceColor', wtc, 'MarkerSize', abs(coinCounts{xx}(ss,bb)/maxCoins*12+0.1));
                if coinCounts{xx}(ss,bb) <= 0
                    set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'r')
                end
            end
        else

            % Show using same marker
            plot(bounds{xx}(:,bb,1,1,1), bounds{xx}(:,bb,2,1,1), 'ko', 'MarkerFaceColor', wtc);
        end

        %plot(medianBounds{xx}(:,bb,1), medianBounds{xx}(:,bb,2), 'ko', ...
        %    'MarkerFaceColor', wtc);
        if bb == 1 
            ylabel('hiSNR bound');
        end
        xlabel('loSNR bound');

        title(sprintf('%s: p=%.3f', titles{xx}, signrank(bounds{xx}(:,bb,1,1,1), bounds{xx}(:,bb,2,1,1))))
        axis([0 0.75 0 0.75])
    end
end
