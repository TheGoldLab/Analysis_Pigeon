function Figure06_boundBlockTransition(dataTableOL, dataTableMX, options)
% function Figure06_boundBlockTransition(dataTableOL, dataTableMX, options)
%
% Figure:

arguments
    dataTableOL
    dataTableMX
    options.stepsToShow = 40;
    options.blocks = 1:2;
    options.figureNumber = 6
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {2,2,2};
hts     = 5;
[axs,~] = getPLOT_axes(options.num, wid, hts, cols, 1.6, 1.6, [], 'Pigeons', true);
set(axs,'Units','normalized');

% For fits
fitFcn = @(b,x) b(1)-b(2).*exp(-x./b(3));
fitOpts = optimset('MaxFunEvals', 10000);

% For plotz
wt = 0.99.*ones(1,3);
gr = 0.8.*ones(1,3);
xax = -options.stepsToShow:options.stepsToShow-1;

% Loop through the two dataTables
titles = {'Mixed' 'Blocked'};
for dd = 1:2
    if dd == 1
        dataTable = dataTableMX;
    else
        dataTable = dataTableOL;
    end

    % Collect data per subject/block/snr
    subjects = nonanunique(dataTable.subjectIndex);
    numSubjects = length(subjects);
    aSNRs = abs(dataTable.snr);
    SNRs = nonanunique(aSNRs);
    numSNRs = length(SNRs);

    % Set up data matrix
    summaryData = nan(numSubjects, options.stepsToShow*2, numSNRs, 2);
    fitData = nan(numSubjects,2,3); % Save taus for each snr

    % Loop through the subjects,snrs
    for ss = 1:numSubjects
        Lsub = dataTable.subjectIndex==subjects(ss);
        for rr = 1:numSNRs
            Lsnr = Lsub & aSNRs==SNRs(rr);

            % get z-scored abs(bound) from sequential blocks
            Lb = Lsnr & ismember(dataTable.blockIndex, options.blocks);
            blockIndex = dataTable.blockIndex(Lb);
            bound = zscore(abs(dataTable.bound(Lb)));
            % avoid last trial of first block, which seems a little wonky
            block1Bound = bound(blockIndex==1);
            block1Bound = block1Bound(max(1,end-options.stepsToShow):end-1);
            block2Bound = bound(blockIndex==2);
            block2Bound = block2Bound(1:min(options.stepsToShow-1,length(block2Bound)));

            % Save the data
            summaryData(ss,options.stepsToShow-length(block1Bound)+1:options.stepsToShow,rr) = block1Bound;
            summaryData(ss,options.stepsToShow+1:options.stepsToShow+length(block2Bound),rr) = block2Bound;

            % Fit tau to second block if there is learning
            tax = (0:length(block2Bound)-1)';
            fitData(ss,rr,:) = fmincon(@(b) sum((block2Bound-fitFcn(b,tax)).^2), ...
                [mean(block2Bound(1:3))-1 mean(block2Bound(end-10:end))./mean(block2Bound(1:10)) 15],[],[],[],[],[-3 -10 1],[3 10 100],[], fitOpts);

            % subplot(2,1,rr); cla reset; hold on;
            % plot(tax,block2Bound,'ko');
            % plot(tax, fitFcn(fitData(ss,rr,:), tax), 'r-')
            % title(sprintf('b1=%.2f, b2=%.2f, tau=%.2f', fitData(ss,rr,1), fitData(ss,rr,2), fitData(ss,rr,3)))
        end
        %r=input('next')
    end

    %% Plotz
    % Per block
    for rr = 1:numSNRs

        % get median/iqr
        prcts = prctile(summaryData(:,:,rr), [25 50 75]);

        % means = mean(summaryData(:,:,rr), 'omitnan');
        % sems = std(summaryData(:,:,rr), 'omitnan')/sqrt(numSubjects);
        % ys = [means-sems flip(means+sems)];

        % Set axes
        axes(axs((rr-1)*2+dd)); cla reset; hold on;

        % Plot it
        plot(xax([1 end]), [0 0], 'k:');
        plot([0 0], [-2 2], 'k:');
        xs = [xax flip(xax)];
        ys = [prcts(1,:) flip(prcts(3,:))];
        Lg = isfinite(ys);
        h=patch(xs(Lg), ys(Lg), gr);
        set(h, 'LineStyle', 'none')
        plot(xax, prcts(2,:), 'k-', 'LineWidth', 2)
        axis([xax(1) xax(end) -1.5 1.5])
        if rr==1
            title(titles{dd})
        end
    end
    xlabel('trial number re: block transition')
    ylabel('z-scored abs(bound)')

    % Show taus from fits
    axes(axs(4+dd)); cla reset; hold on;
    plot([0 100], [0 100], 'k:')
    Lboth = fitData(:,1,3)<90 & fitData(:,2,3)<90;
    Leither = fitData(:,1,3)<90 | fitData(:,2,3)<90;
    plot(fitData(:,1,3), fitData(:,2,3), 'ko', 'MarkerFaceColor', wt);
    plot(fitData(Leither,1,3), fitData(Leither,2,3), 'ko', 'MarkerFaceColor', gr);
    plot(fitData(Lboth,1,3), fitData(Lboth,2,3), 'ko', 'MarkerFaceColor', 'k');
    title(sprintf('p=%.2f', ranksum(fitData(Lboth,1,3), fitData(Lboth,2,3))))
    if dd == 1
        xlabel('Block 1 tau')
        ylabel('Block 2 tau')
    end
end    
