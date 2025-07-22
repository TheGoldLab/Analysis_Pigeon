function boundBiasCorrection = makePigeon_biasCorrections(options)
% function boundBiasCorrection = makePigeon_biasCorrections(options)
%
% Corrects biases in estimating bounds by fitting linear functions to
%   true vs extracted bounds per RT, from simulations

arguments
    options.blockIndex = 2;
    options.numTrials = 1000000;
    options.simBounds = (0.01:0.05:0.76)';
    options.generativeMeans = [0.05 0.1 0.15];
    options.maxRT = 12;
    options.NDTMin = 1;
    options.NDTMax = 3;
    options.showOutput = false;
    options.filename = 'boundBiasCorrection';
end

% Set up simulations
specs.blocks = options.blockIndex;
specs.numSubjects = 1;
numMeans = length(options.generativeMeans);
numSimBounds = length(options.simBounds);
numRTs = options.maxRT;
boundSummary = nan(numSimBounds, numRTs, numMeans);
correctedBoundSummary = nan(numSimBounds, numRTs, numMeans);

% Loop through the simulated obounds
for ss = 1:numSimBounds
    disp(options.simBounds(ss))
    for mm = 1:numMeans

        % Get simulated bounds without any correction
        simTable = getPigeon_simulatedDataTable(specs, ...
            'generativeMean',   options.generativeMeans(mm), ...
            'numTrials',        options.numTrials, ...
            'stepsPerBlock',    options.numTrials, ...
            'NDTMin',           options.NDTMin, ...
            'NDTmax',           options.NDTMax, ...
            'boundMean',        options.simBounds(ss), ...
            'boundSTD',         0, ...
            'correctBias',      []);

        %   boundSummary_ .. matrix with dims:
        %       subjects
        %       blocks
        %       snrs
        %       rt bins
        %       bound avg, bound std, n
        boundSummaryByRT = getPigeon_boundSummary(simTable, 'maxRT', options.maxRT);
        boundSummary(ss,:,mm) = squeeze(boundSummaryByRT(1,1,1,:,1));

        if options.showOutput
            % If we are showing the output, we will compare to the
            % corrected bounds
            correctedDataTable = getPigeon_biasCorrections(simTable, 'boundBiasCorrection.mat');
            correctedBoundSummaryByRT = getPigeon_boundSummary(correctedDataTable, 'maxRT', options.maxRT);
            correctedBoundSummary(ss,:,mm) = squeeze(correctedBoundSummaryByRT(1,1,1,:,1));
        end
    end
end

% Possibly show figure
if options.showOutput
    wid     = 14; % total width
    cols    = repmat({3},1,options.maxRT);
    hts     = 0.8;
    [axs,~] = getPLOT_axes('BiasCorrection', wid, hts, cols, 0.4, 0.4, [], 'Pigeons', true);
    set(axs,'Units','normalized');

    for xx = 1:length(axs)
        axes(axs(xx)); cla reset; hold on;
        plot([-0.5 1], [-0.5 1], 'k:')
        axis([0 1 0 1])        
        if xx < length(axs)-2
            set(gca,'XTickLabel', '');
        end
        if mod(xx-1,3)
            set(gca,'YTickLabel', '');
        else
            title(sprintf('DT=%d',(xx-1)/3+1 ));
        end        
    end
    axes(axs(end-2));
    xlabel('Measured bound')
    ylabel('True bound')
end

% lookup table is really matrix of linear regression coefficients
boundBiasCorrection = nan(numRTs,2,3); % intercept/slope, 2 means
Xax = cat(2, ones(numSimBounds,1), options.simBounds);
for mm = 1:3
    for rr = 1:numRTs
        Lg = isfinite(boundSummary(:,rr,mm));
        X = cat(2, ones(sum(Lg),1), boundSummary(Lg,rr,mm));
        boundBiasCorrection(rr,:,mm) = regress(options.simBounds(Lg), X);
        if options.showOutput
            axes(axs((rr-1)*3+mm));
            plot(boundSummary(Lg,rr,mm), options.simBounds(Lg), 'k:', 'LineWidth', 2)
            plot(options.simBounds, Xax*boundBiasCorrection(rr,:,mm)', 'r-')
            plot(correctedBoundSummary(Lg,rr,mm), options.simBounds(Lg), 'k-', 'LineWidth', 2)
        end
    end
end

if ~isempty(options.filename)
    save([options.filename '.mat'], options.filename);
end

% % plot actual bound (y) vs measured bound (x) for each RT
% for rr = 1:size(boundSummary,2)
%     for mm=1:2
%         subplot(2,1,mm); cla reset; hold on;
%         plot(boundSummary(:,rr,mm), options.simBounds, 'k-')
%         % plot(boundSummary(:,rr,mm), options.simBounds, 'k-')
%         plot([-1 2], [-1 2], 'k:')
%         %axis([-0.2 1.8 0 0.8])
%     end
%     r = input('next')
% end
%
