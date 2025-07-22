% playpen for testing bound estimates

specs.blocks = 2;
specs.numSubjects = 1;
maxRT = 15;
generativeMeans = [0.05 0.15];
numMeans = length(generativeMeans);
simBounds = 0.0:0.05:0.6;
numSimBounds = length(simBounds);
simBoundSummary = nan(numSimBounds, numMeans, 2); % last is mean/std
simBoundSummaryByRT = nan(numSimBounds, numMeans, maxRT+1); % last is mean by RT

ndts = nan(numSimBounds,2);
for ss = 1:numSimBounds
    disp(simBounds(ss))
    for mm = 1:numMeans

        simTable = getPigeon_simulatedDataTable(specs, ...
            'generativeMean',   generativeMeans(mm), ...
            'numTrials',        100000, ...
            'stepsPerBlock',    100000, ...
            'NDTMin',           1, ...
            'NDTmax',           2, ...
            'boundMean',        simBounds(ss), ...
            'boundSTD',         0);
        ndts(ss,mm) = mean(simTable.RT-simTable.DT);

        %   boundSummary_ .. matrix with dims:
        %       subjects
        %       blocks
        %       snrs
        %       rt bins
        %       bound avg, bound std, n
        % summary across all RTs
        boundSummary = getPigeon_boundSummary(simTable, 'maxRT', 'all');
        simBoundSummary(ss,mm,:) = squeeze(boundSummary(1,1,1,1,1:2));

        % summary by RT
        boundSummaryByRT = getPigeon_boundSummary(simTable, 'maxRT', maxRT);
        simBoundSummaryByRT(ss,mm,:) = squeeze(boundSummaryByRT(1,1,1,:,1));
    end
end

figure
co = {'r' 'b'};
% all RTs
subplot(2,1,1); cla reset; hold on;
plot([0 0.8], [0 0.8], 'k:');
for mm = 1:2
    plot(simBounds, simBoundSummary(:,mm,1), '-', 'Color', co{mm});
end

% diff in estimated bound (hi-lo snr) as a function of RT, per bound
subplot(2,1,2); cla reset; hold on;
plot([0 maxRT], [0 0], 'k:');    
tax = 0:maxRT;
diffs = squeeze(diff(simBoundSummaryByRT,[],2));
for bb = 1:numSimBounds
    h=plot(tax, diffs(bb,:), '.-', 'Color', ones(3,1).*bb/numSimBounds);
end
set(h, 'Color', 'r')

% % As a function of bound, per RT
% subplot(2,1,2); cla reset; hold on;
% plot([0 0.8], [0 0.8], 'k:');
% for mm = 1:2
%     for rr = 1:maxRT
%         h=plot(simBounds, simBoundSummaryByRT(:,mm,rr), '-', 'Color', co{mm});
%         set(h, 'LineWidth', 0.25)
%         if rr==maxRT
%             set(h, 'LineWidth', 3);
%         end
%     end
% end
