function boundBiasCorrection = makePigeon_biasCorrections(options)
% function boundBiasCorrection = makePigeon_biasCorrections(options)
%
% Corrects biases in estimating bounds by fitting linear functions to
%   true vs extracted bounds per RT, from simulations
% 
% Remember to run it twice if you change params and want to plot it -- 

arguments
    options.blockIndex = 2;
    options.numTrials = 1000000;
    options.simBounds = (0.01:0.05:0.71)';
    options.generativeMeans = [0.05 0.1 0.15];
    options.generativeSTD = 0.15;
    options.maxStepsPerTrial = 50;
    options.boundMax = 0.75;
    options.maxRT = 12;
    options.NDTMin = 1;
    options.NDTMax = 3;
    options.showOutput = true;
    options.maxOutputRows = 4;
    options.maxOutputCols = 3;
    options.filename = 'boundBiasCorrection';
    options.save = true;
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
            'generativeSTD',    options.generativeSTD, ...
            'maxStepsPerTrial', options.maxStepsPerTrial, ...
            'numTrials',        options.numTrials, ...
            'stepsPerBlock',    options.numTrials, ...
            'NDTMin',           options.NDTMin, ...
            'NDTmax',           options.NDTMax, ...
            'boundMean',        options.simBounds(ss), ...
            'boundSTD',         0, ...
            'boundMax',         options.boundMax, ...
            'correctBias',      []);

        %   boundSummary_ .. matrix with dims:
        %       subjects
        %       blocks
        %       snrs
        %       rt bins
        %       bound avg, bound std, n
        boundSummaryByRT = getPigeon_boundSummary(simTable, 'maxRT', options.maxRT);
        boundSummary(ss,:,mm) = squeeze(boundSummaryByRT(1,1,1,:,1));

        if options.showOutput && exist([options.filename '.mat'], 'file')
            % If we are showing the output, we will compare to the
            % corrected bounds
            correctedDataTable = getPigeon_biasCorrections(simTable, [options.filename '.mat']);
            correctedBoundSummaryByRT = getPigeon_boundSummary(correctedDataTable, 'maxRT', options.maxRT);
            correctedBoundSummary(ss,:,mm) = squeeze(correctedBoundSummaryByRT(1,1,1,:,1));
        end
    end
end

% Possibly show figure
% Columns are: 
if options.showOutput
    options.maxOutputCols = min(options.maxOutputCols, numMeans);
    options.maxOutputRows = min(options.maxOutputRows, numRTs);

    wid     = 18; % total width
    cols    = repmat({options.maxOutputCols},1,options.maxOutputRows);
    hts     = 5; %0.8;
    [axs,~] = getPLOT_axes('BiasCorrection', wid, hts, cols, 0.4, 0.4, [], 'Pigeons', true);
    % set(axs,'Units','normalized');

    for xx = 1:length(axs)
        axes(axs(xx)); cla reset; hold on;
        plot([0 options.boundMax], [0 options.boundMax], 'k:')
        axis([0 options.boundMax 0 options.boundMax])        
        if xx < length(axs)-options.maxOutputCols-1
            set(gca,'XTickLabel', '');
        end
        if mod(xx-1,options.maxOutputCols)
            set(gca,'YTickLabel', '');
        else
            title(sprintf('DT=%d',(xx-1)/options.maxOutputCols+1 ));
        end
        if xx == 1
            title(sprintf('Mean=%.2f\nDT=1', options.generativeMeans(xx)))
        elseif xx <= options.maxOutputCols
            title(sprintf('Mean=%.2f', options.generativeMeans(xx)))
        end
    end
    axes(axs(end-2));
    xlabel('Measured bound')
    ylabel('True bound')
    axIndex = 1; % counter
end

% lookup table is really matrix of linear regression coefficients
boundBiasCorrection = nan(numRTs,2,numMeans); % intercept/slope, 2 means
Xax = cat(2, ones(numSimBounds,1), options.simBounds);
for rr = 1:numRTs
    for mm = 1:numMeans
        Lg = isfinite(boundSummary(:,rr,mm));
        X = cat(2, ones(sum(Lg),1), boundSummary(Lg,rr,mm));
        boundBiasCorrection(rr,:,mm) = regress(options.simBounds(Lg), X);
        if options.showOutput && mm<=options.maxOutputCols && rr<=options.maxOutputRows
            axes(axs(axIndex));
            axIndex = axIndex + 1;
            % Thick red is uncorrected bound
            plot(boundSummary(Lg,rr,mm), options.simBounds(Lg), 'r-', 'LineWidth', 2)
            % Dashed blue is fit to bound
            plot(options.simBounds, Xax*boundBiasCorrection(rr,:,mm)', 'b--')
            % Thick green is corrected
            plot(correctedBoundSummary(Lg,rr,mm), options.simBounds(Lg), 'g-', 'LineWidth', 2)
        end
    end
end

if ~isempty(options.filename) && options.save
    snrs = options.generativeMeans./options.generativeSTD;
    save([options.filename '.mat'], 'boundBiasCorrection', 'snrs');
end


% for figure
wid     = 6.5; % total width
cols    = {1};
hts     = 5;
[axs,~] = getPLOT_axes('BiasCorrection', wid, hts, cols, 1.3, 1.5, [], 'Pigeons', true);
axes(axs); cla reset; hold on;
Xax = cat(2, ones(numSimBounds,1), options.simBounds);
plot([0 options.boundMax], [0 options.boundMax], 'k:')
axis([0 options.boundMax 0 options.boundMax])
xlabel('Simulated (true) bound')
ylabel('Measured bound')
mm = 1;
rr = 2;
Lg = isfinite(boundSummary(:,rr,mm));
X = cat(2, ones(sum(Lg),1), options.simBounds(Lg));
bbc = regress(boundSummary(Lg,rr,mm), X);
% Thick gray is uncorrected bound
plot(options.simBounds(Lg), boundSummary(Lg,rr,mm), '-', 'LineWidth', 2, 'Color', 0.8.*ones(3,1));
% Blue is fit to bound
plot(options.simBounds, Xax*bbc, 'b-')
% Thick black is corrected
plot(options.simBounds(Lg), correctedBoundSummary(Lg,rr,mm), 'k-', 'LineWidth', 2)

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
