function dataTable_ = getPigeon_simulatedDataTable(dataTable, options)
%
% Returns:
%   dataTable_ as in getPigeon_dataTable
%
% Arguments:
%   property/value pairs, see below

% Arguments
arguments
    dataTable
    options.generativeMean = 0.05;
    options.generativeSTD = 0.15;
    options.numTrials = 600;
    options.maxStepsPerTrial = 50;
    options.boundMean = 0;
    options.boundSTD = 0;
    options.boundSlope = 0;
    options.boundMax = 0.75;
    options.NDTMin = 2;
    options.NDTmax = 3;
    options.lapseRate = 0;
    options.boundType = 'given'; %'fixed', 'true', 'var', 'varBySNR'
    options.stepsPerBlock = 600;
    options.P2Psteps = 1;
    options.P2Pcoins = 0;
    options.coinsGainedPerCorrect = 1;
    options.coinsLostPerError = 0;
    options.stepsLostPerError = 0;
    options.bounds = [];
    options.correctBias = 'boundBiasCorrection.mat';
end

% hard-coding in definitions of blocks 1-6 in Alice's data
blockDefaults = struct( ...
    'stepsPerBlock',            repmat({options.stepsPerBlock}, 1, 6), ...
    'P2Psteps',                 {1 1 1 1 1 1}, ...
    'P2Pcoins',                 {0 0 0 0 0 0}, ...
    'coinsGainedPerCorrect',    {1 1 1 1 1 1}, ...
    'coinsLostPerError',        {0 4 0 0 4 0}, ...
    'stepsLostPerError',        {0 0 30 0 0 30});

% Make blockSpecs -- matrix of option structs to run simulations
%   rows are subjects
%   columns are blocks
% Check to build from existing data Table
if nargin >= 1 && istable(dataTable)
    subjects = nonanunique(dataTable.subjectIndex);
    numSubjects = length(subjects);
    blocks = nonanunique(dataTable.blockIndex);
    numBlocks = length(blocks);
    Lgood = dataTable.trialNumber & dataTable.RT>=0; % non-zero RT
    blockSpecs = repmat(options, numSubjects, numBlocks);    

    % Loop through the blocks and subjects
    for bb = 1:numBlocks

        % Add block defaults
        for ff = fieldnames(blockDefaults)'
            [blockSpecs(:,bb).(ff{:})] = deal(blockDefaults(blocks(bb)).(ff{:}));
        end

        % Loop through subject data
        Lb = Lgood & dataTable.blockIndex==blocks(bb);
        aSNRs = abs(dataTable.snr);
        for ss = 1:numSubjects
            Lbs = Lb & dataTable.subjectIndex==subjects(ss);
            SNRs = nonanunique(aSNRs(Lbs));
            blockSpecs(ss,bb).generativeMean = SNRs .* ...
                blockSpecs(ss,bb).generativeSTD;
            nSNRs = length(SNRs);
            % parse bound type
            switch blockSpecs(ss,bb).boundType
                case 'given'
                    blockSpecs(ss,bb).boundSTD = 0;
                case 'fixed'
                    absBounds = abs(dataTable.bound(Lbs));
                    blockSpecs(ss,bb).boundMean = median(absBounds,'omitnan');
                    blockSpecs(ss,bb).boundSTD = 0;
                case 'true'
                    abounds = abs(dataTable.bound(Lbs));
                    abounds(~isfinite(abounds)) = 0;
                    blockSpecs(ss,bb).bounds = abounds;
                case 'var'
                    absBounds = abs(dataTable.bound(Lbs));
                    blockSpecs(ss,bb).boundMean = median(absBounds,'omitnan');
                    blockSpecs(ss,bb).boundSTD = std(absBounds,'omitnan');
                case 'changepoint'
                    Lpc = Lbs & (dataTable.DT < dataTable.changePoint);
                    blockSpecs(ss,bb).changePoint = nonanunique(dataTable.changePoint(Lbs));
                    blockSpecs(ss,bb).boundMean = NaN(options.maxStepsPerTrial,1);
                    blockSpecs(ss,bb).boundMean(1:blockSpecs(ss,bb).changePoint-1) = median(dataTable.absBound(Lpc),'omitnan');
                    blockSpecs(ss,bb).boundMean(blockSpecs(ss,bb).changePoint:end) = median(dataTable.absBound(~Lpc),'omitnan');
                    blockSpecs(ss,bb).boundSTD = 0;
                    % save both snrs
                    blockSpecs(ss,bb).snrs = unique(aSNRs);
                case 'counterfactual'
                    Lpc = Lbs & (dataTable.RT > dataTable.changePoint);
                    blockSpecs(ss,bb).boundMean = median(dataTable.absBound(Lpc),'omitnan');
                    blockSpecs(ss,bb).boundSTD = 0;
                otherwise % 'varBySNR'
                    % Bound computed separately per snr
                    blockSpecs(ss,bb).boundMean = nan(nSNRs,1);
                    blockSpecs(ss,bb).boundSTD = nan(nSNRs,1);
                    for mm = 1:nSNRs
                        Lm = Lbs & aSNRs == SNRs(mm);
                        absBounds = abs(dataTable.bound(Lm));
                        blockSpecs(ss,bb).boundMean(mm) = median(absBounds,'omitnan');
                        blockSpecs(ss,bb).boundSTD(mm) = std(absBounds,'omitnan');
                    end
            end
        end
    end
elseif nargin >= 1 && isstruct(dataTable)
    % expects struct with fields:
    %   numSubjects
    %   blocks
    blocks = nonanunique(dataTable.blocks);
    numBlocks = length(blocks);
    blockSpecs = repmat(options, dataTable.numSubjects, numBlocks);
    for bb = 1:numBlocks
        % Add block defaults
        for ff = fieldnames(blockDefaults)'
            [blockSpecs(:,bb).(ff{:})] = deal(blockDefaults(blocks(bb)).(ff{:}));
        end
    end
end

% Set up the data table
%   Note that we're using the same fields as from getDataTable, with the
%   addition of "trueBound" -- bound is computed empirically, trueBound is,
%   duh, the real one used in the simulations
dataTable_ = getPigeon_blankDataTable(sum([blockSpecs.numTrials]), 'trueBound', "double");
tableIndex = 1;
argNames = {'generativeMean', 'generativeSTD', 'numTrials', ...
    'maxStepsPerTrial', 'boundMean', 'boundSTD', ...
    'boundSlope', 'boundMax', 'NDTMin', ...
    'NDTmax', 'lapseRate', 'bounds'};
args = cell(1,length(argNames)*2);
args(1:2:end) = argNames;

% Now loop through the specs and do the simulations
for ss = 1:size(blockSpecs,1) % Per subject
    for bb = 1:size(blockSpecs,2) % Per block

        % Make the arg list for the simulations
        for aa = 1:length(argNames)
            args{aa*2} = blockSpecs(ss,bb).(argNames{aa});
        end

        % Do the sim
        [choices, rts, ~, ~, ~, snrs, steps] = getPigeon_simulatedData(args{:});

        % Update the data Table
        % jig changed to match real data, where steps always include zero
        % start
        %         stepCounts(:) = cumsum(repmat(blockArgs(bb).P2Psteps,simArgs(bb).numTrials,1) + ...
        %             steps + blockArgs(bb).stepsLostPerError .* (choices==0));
        stepCounts = cumsum(rts + blockSpecs(ss,bb).stepsLostPerError .* (choices==0));
        coinCounts = cumsum(repmat(-blockSpecs(ss,bb).P2Pcoins,blockSpecs(ss,bb).numTrials,1) - ...
            blockSpecs(ss,bb).coinsLostPerError .* (choices==0) + ...
            blockSpecs(ss,bb).coinsGainedPerCorrect .* (choices==1));
        trialCount = find(stepCounts>=blockSpecs(ss,bb).stepsPerBlock,1);

        % Add for final aborted trial
        if stepCounts > blockSpecs(ss,bb).stepsPerBlock
            choices(trialCount) = nan;
            coinCounts(trialCount) = coinCounts(trialCount-1);
        end

        % Update the table
        inds = tableIndex:tableIndex+trialCount-1;
        if isempty(inds)
            disp(length(inds))
        end
        if inds(end) > height(dataTable_)
            dataTable_ = cat(1, dataTable_, ...
                getPigeon_blankDataTable(10000, 'trueBound', "double"));
        end

        % Save stuff
        dataTable_.subjectIndex(inds) = ss;
        dataTable_.blockIndex(inds) = blocks(bb);
        dataTable_.trialNumber(inds) = 1:trialCount;
        dataTable_.RT(inds) = rts(1:trialCount);
        dataTable_.choice(inds) = choices(1:trialCount);
        dataTable_.correct(inds) = choices(1:trialCount);
        dataTable_.coinCount(inds) = coinCounts(1:trialCount);
        dataTable_.snr(inds) = snrs(1:trialCount);

        % Get bound, etc in the standard way
        [dataTable_.bound(inds), ...
            dataTable_.DT(inds), ...
            dataTable_.RT(inds), ...
            dataTable_.congruence{inds(1)}] = ...
            getPigeon_bounds(steps(1:trialCount), choices(1:trialCount), ...
            'SNR', snrs(1:trialCount));
        dataTable_.steps(inds) = steps(1:trialCount);

        % update index
        tableIndex = inds(end)+1;
    end
end

% Remove fluff
dataTable_ = dataTable_(1:tableIndex-1,:);

% Possibly correct biases
if exist(options.correctBias, 'file')
    dataTable_ = getPigeon_biasCorrections(dataTable_, options.correctBias);
end
