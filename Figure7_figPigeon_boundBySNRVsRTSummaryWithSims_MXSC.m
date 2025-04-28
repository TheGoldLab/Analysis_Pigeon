function Figure7_figPigeon_boundBySNRVsRTSummaryWithSims_MXSC(dataTableMXSC, simTableMXSC, options)
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
    dataTableMXSC
    simTableMXSC
    options.minN=2;
    options.maxRT=30;
    options.numSims = 1;
end

% rt axis
rtAxis = 2:options.maxRT;
numRTs = length(rtAxis);
lagAxis = -4:4;
numLags = length(lagAxis);
numBlocks = length(unique(dataTableMXSC.blockIndex));
numSubjects = length(unique(dataTableMXSC.subjectIndex));

%% Set up figure
wid     = 17.6; % total width
cols    = {2,2};
hts     = [4.5 4.5];
[axs,~] = getPLOT_axes(1, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
set(axs,'Units','normalized');


% matrices with dims subjects,blocks,snrs,rts,data (bound avg/std/n):
data = getPigeon_boundSummary(dataTableMXSC, 'maxRT', options.maxRT);
simData = getPigeon_boundSummary(simTableMXSC, 'maxRT', options.maxRT);

%% Plotz
titles = {'fixed SNRdiff', 'mixed SNRdiff', 'loHi change', 'hiLo change'};
wtc = 0.99.*ones(3,1);
for bb = 2:numBlocks
    axes(axs(bb-1)); hold on;
    title(titles{bb-1});
    if bb == 2
        plot(rtAxis([1 end]), [0 0], 'k:');

        % compute simulated differences
        ys = squeeze(simData(:,2,2,rtAxis,1)) - squeeze(simData(:,1,1,rtAxis,1));
        ns = squeeze(simData(:,bb,mod(bb,3),rtAxis,3));
        Ln = squeeze(ns(:,:))>=options.minN;
        ys(~Ln) = nan;
        simBounds(:,:,:,1) = ys;
        % Compute as difference in bound for two SNRs (per RT)
        ys = squeeze(data(:,2,2,rtAxis,1)) - squeeze(data(:,1,1,rtAxis,1));
        ns = squeeze(data(:,bb,mod(bb,3),rtAxis,3));
        Ln = squeeze(ns(:,:))>=options.minN;
        if options.numSims > 0
            ys = ys - squeeze(mean(simBounds,1,'omitnan'));
            % ys = squeeze(mean(simBounds(:,bb,:,:,xx),4,'omitnan'));
            % Ln = isfinite(ys);
        end
        xax = repmat(rtAxis,size(ys,1),1) + normrnd(0,0.1,size(ys));
        plot(xax(Ln), ys(Ln), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', wtc);
        for rr = 1:numRTs
            diffs = ys(Ln(:,rr), rr);
            diffs = diffs(isfinite(diffs));
            if any(diffs)
                h=plot(rtAxis(rr), median(diffs), 'ro', 'MarkerSize', 7);
                if signrank(diffs) < 0.05
                    set(h, 'MarkerFaceColor', 'r');
                end
            end
        end
        axis([rtAxis([1 options.maxRT/2]) -0.5 0.5])
    elseif bb == 3 % mixed snr block
        plot(rtAxis([1 end]), [0 0], 'k:', 'LineWidth', 3);

        % compute simulated differences
        ys = squeeze(diff(simData(:,bb,1:2,rtAxis,1),[],3));
        ns = squeeze(simData(:,bb,1:2,rtAxis,3));
        Ln = squeeze(ns(:,:,1,:))>=options.minN & squeeze(ns(:,:,2,:))>=options.minN;
        ys(~Ln) = nan;
        simBounds(:,:,:,1) = ys;
        % Compute as difference in bound for two SNRs (per RT)
        ys = squeeze(diff(data(:,bb,1:2,rtAxis,1),[],3));
        ns = squeeze(data(:,bb,1:2,rtAxis,3));
        Ln = squeeze(ns(:,1,:))>=options.minN & squeeze(ns(:,2,:))>=options.minN;
        if options.numSims > 0
            ys = ys - squeeze(mean(simBounds,1,'omitnan'));
            % ys = squeeze(mean(simBounds(:,bb,:,:,xx),4,'omitnan'));
            % Ln = isfinite(ys);
        end
        xax = repmat(rtAxis,size(ys,1),1) + normrnd(0,0.1,size(ys));
        plot(xax, ys, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', wtc);
        for rr = 1:numRTs
            diffs = ys(Ln(:,rr), rr);
            diffs = diffs(isfinite(diffs));
            if any(diffs)
                h=plot(rtAxis(rr), median(diffs), 'ro', 'MarkerSize', 7);
                if signrank(diffs) < 0.05
                    set(h, 'MarkerFaceColor', 'r');
                end
            end
        end
        axis([rtAxis([1 options.maxRT/2]) -0.5 0.5])
    else % changepoint blocks
        plot(lagAxis([1 end]), [0 0], 'k:', 'LineWidth',3);
        diffs = zeros(numSubjects,length(lagAxis));
        for ss = 1:numSubjects
            Lsb = (dataTableMXSC.subjectIndex == ss) & (dataTableMXSC.blockIndex == bb);
            changePoint = unique(dataTableMXSC.changePoint(Lsb));
            if changePoint > abs(lagAxis(1)) && changePoint < options.maxRT - lagAxis(end)
                ys = squeeze(simData(ss,bb,mod(bb,3),changePoint + lagAxis,1))';
                ns = squeeze(simData(ss,bb,mod(bb,3),changePoint + lagAxis,3))';
                % Ln = squeeze(ns(:,:,mod(bb,3),:))>=options.minN;
                % ys(~Ln) = nan;
                simBounds = ys;
    
                ys = squeeze(data(ss,bb,mod(bb,3),changePoint + lagAxis,1))';
                ns = squeeze(data(ss,bb,mod(bb,3),changePoint + lagAxis,3))';
                % Ln = squeeze(ns(:,1,:))>=options.minN & squeeze(ns(:,2,:))>=options.minN;
                if options.numSims > 0
                    ys = ys - simBounds;
                    % ys = squeeze(mean(simBounds(:,bb,:,:,xx),4,'omitnan'));
                    % Ln = isfinite(ys);
                end
                xax = lagAxis + normrnd(0,0.1,size(lagAxis));
                plot(xax, ys, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', wtc);
                diffs(ss, :) = ys;
            end
        end
        if changePoint > abs(lagAxis(1)) && changePoint < options.maxRT - lagAxis(end)
            for ll = 1:numLags
                d = diffs(:,ll);
                d = d(isfinite(d));
                if any(d)
                    h=plot(lagAxis(ll), median(d), 'ro', 'MarkerSize', 7);
                    if signrank(d) < 0.05
                        set(h, 'MarkerFaceColor', 'r');
                    end
                end
            end
        end

        axis([lagAxis([1 end]) -0.5 0.5])
    end
    xlabel('RT')
    if mod(bb,2) == 0
        ylabel('Bound diff re: Sims')
    end
end
