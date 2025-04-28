function Figure6_figPigeon_boundBySNRVsRTSummaryWithSims(dataTableMX, dataTableOL, options)
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
    options.minN=3;
    options.maxRT=10;
    options.numSims = 20;
end

% rt axis
rtAxis = 2:options.maxRT;
numRT = length(rtAxis);
numBlocks = 3;

%% Set up figure
wid     = 17.6; % total width
cols    = {3,3};
hts     = [3.5 3.5];
[axs,~] = getPLOT_axes(1, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
set(axs,'Units','normalized');

% matrices with dims subjects,blocks,snrs,rts,data (bound avg/std/n):
data = { ...
    getPigeon_boundSummary(dataTableMX, 'maxRT', options.maxRT), ...
    getPigeon_boundSummary(dataTableOL, 'maxRT', options.maxRT)};

% Get sim data
if options.numSims > 0
    numSubjects = size(data{1},1);
    simBounds = nan(numSubjects, numBlocks, numRT, options.numSims, 2);
    tables = {dataTableMX, dataTableOL};

    for mm = 1:options.numSims
        fprintf('Sim %d/%d\n', mm, options.numSims)

        for xx = 1:2
            % returns subj,block,snr,rtbin,mean/std/n
            boundSummary = getPigeon_boundSummary( ...
                getPigeon_simulatedDataTable(tables{xx}, 'boundType', 'true'), ...
                'maxRT', options.maxRT);

            ys = squeeze(diff(boundSummary(:,:,1:2,rtAxis,1),[],3));
            ns = squeeze(boundSummary(:,:,1:2,rtAxis,3));
            Ln = squeeze(ns(:,:,1,:))>=options.minN & squeeze(ns(:,:,2,:))>=options.minN;
            ys(~Ln) = nan;
            simBounds(:,:,:,mm,xx) = ys;
        end
    end
end


%% Plotz
%
titles = {'mixed', 'blocked'};
wtc = 0.99.*ones(3,1);
for bb = 1:numBlocks
    for xx = 1:2
        axes(axs((xx-1)*numBlocks+bb)); cla reset; hold on;
        plot(rtAxis([1 end]), [0 0], 'k:');

        % Compute as difference in bound for two SNRs (per RT)
        ys = squeeze(diff(data{xx}(:,bb,1:2,rtAxis,1),[],3));
        ns = squeeze(data{xx}(:,bb,1:2,rtAxis,3));
        Ln = squeeze(ns(:,1,:))>=options.minN & squeeze(ns(:,2,:))>=options.minN;
        if options.numSims > 0
            ys = ys - squeeze(mean(simBounds(:,bb,:,:,xx),4,'omitnan'));
            % ys = squeeze(mean(simBounds(:,bb,:,:,xx),4,'omitnan'));
            % Ln = isfinite(ys);
        end
        xax = repmat(rtAxis,size(ys,1),1) + normrnd(0,0.1,size(ys));
        plot(xax(Ln), ys(Ln), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', wtc);
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
        if bb==1
            % title(titles{xx})
            ylabel('Bound diff re: SNR')
        end
        xlabel('RT')
    end
end