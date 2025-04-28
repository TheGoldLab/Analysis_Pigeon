function RRMatrix_ = getPigeon_RRMatrix(options)
% function RRMatrix_ = getPigeon_RRMatrix(options)
%
% simulated RR for block vs mixed snr, different bounds
%
% dataTable columns:
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
    options.bounds = 0.01:0.05:0.8;
    options.boundSTD = 0.1;
    options.numReps = 100;
    options.gMeans = [0.05 0.15];
    options.gSTD = 0.15;
    options.blocks = 1:3;
    options.blockType = 'OL'; % 'OL' or 'MX'
    options.axes = [];
end

numBounds = length(options.bounds);
numBlocks = length(options.blocks);
RRMatrix_ = zeros(numBounds, numBounds, numBlocks);
simParams = struct('numSubjects', options.numReps, 'blocks', options.blocks);

if strcmpi(options.blockType, 'OL')

    % Blockwise SNR
    for b1 = 1:numBounds

        disp(b1)

        % low SNR block
        loSNRTable = getPigeon_simulatedDataTable(simParams, ...
            'generativeMean',   options.gMeans(1), ...
            'boundMean',        options.bounds(b1), ...
            'boundSTD',         options.boundSTD);
        loCoinCounts = getPigeon_coinSummary(loSNRTable);

        % high SNR block
        hiSNRTable = getPigeon_simulatedDataTable(simParams, ...
            'generativeMean',   options.gMeans(2), ...
            'boundMean',        options.bounds(b1), ...
            'boundSTD',         options.boundSTD);
        hiCoinCounts = getPigeon_coinSummary(hiSNRTable);

        % Save reward rates
        for bb = 1:3 % for each block
            RRMatrix_(:,b1,bb) = RRMatrix_(:,b1,bb,1) + mean(loCoinCounts(:,bb));
            RRMatrix_(b1,:,bb) = RRMatrix_(b1,:,bb,1) + mean(hiCoinCounts(:,bb));
        end
    end
else

    % Mixed block SNR
    for b1 = 1:numBounds
        disp([b1 numBounds])        
        for b2 = 1:numBounds

            % mixed SNR block
            mxSNRTable = getPigeon_simulatedDataTable(simParams, ...
                'generativeMean',   options.gMeans, ...
                'boundType',        'varBySNR', ...
                'boundMean',        options.bounds([b1 b2]), ...
                'boundSTD',         [options.boundSTD options.boundSTD]);

            % Save reward rates
            RRMatrix_(b1,b2,:) = mean(getPigeon_coinSummary(mxSNRTable));
        end
    end
end

if options.axes
    for bb = 1:3
        axes(options.axes(bb)); cla reset; hold on;
        imagesc(options.bounds, flip(options.bounds), flipud(RRMatrix_(:,:,bb)));
        set(gca,'YDir', 'normal')
        axis(options.bounds([1 end 1 end]))
        % xlabel('low SNR')
        % ylabel('high SNR')
    end
end
