function Figure02_performanceSummary(dataTable, options)
% function Figure02_performanceSummary(dataTable, options)
%
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
    dataTable
    options.block=2    
    options.figureNumber=2
end

%% Set up figure
wid     = 11.7; % total width
cols    = {2,2};
hts     = [4,4];
[axs,~] = getPLOT_axes(options.figureNumber, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Useful variables
numSubjects = length(unique(dataTable.subjectIndex));
aSNRs = abs(dataTable.snr);
uniqueSNRs = unique(aSNRs);
numSNRs = length(uniqueSNRs);
congruenceSummary = nan(numSubjects, 4, numSNRs); % numNDTs
performanceSummary = nan(numSubjects, numSNRs, 2);  % accuracy/DT

% Loop through the subjects, snrs to collect data
for ss = 1:numSubjects
    Ls = dataTable.subjectIndex == ss; 
    congruences = dataTable.congruence{find(Ls,1)};
    congruenceSummary(ss,:,1) = congruences(:,1);
    congruenceSummary(ss,:,2) = congruences(:,end);
    
    for aa = 1:numSNRs
        Lsb = Ls & dataTable.blockIndex == options.block & aSNRs == uniqueSNRs(aa);

        % save p(Correct), DT
        performanceSummary(ss,aa,:) = [sum(dataTable.correct(Lsb)==1)/sum(dataTable.correct(Lsb)>=0) ...
            mean(dataTable.DT(Lsb), 'omitnan')];
    end
end

%% Top row is congruence between pigeon positions and choices, per NDT
% Separately per SNR
for xx = 1:2
    axes(axs(xx)); cla reset; hold on;
    plot(congruenceSummary(:,:,xx)', '-', 'Color', 0.5.*ones(3,1));
    plot(median(congruenceSummary(:,:,xx), 'omitnan'), 'k-', 'LineWidth', 2);
    if xx == 1
        ylabel('Congruence')
        xlabel('NDT (steps)')
    end
end

%% Bottom row is accuracy vs NDT per subject
% Separately per SNR
for xx = 1:2
    axes(axs(2+xx)); cla reset; hold on;
    plot(performanceSummary(:,xx,2), performanceSummary(:,xx,1), 'ko', ...
        'MarkerFaceColor', 0.99.*ones(3,1));
    plot(0, median(performanceSummary(:,xx,1)), 'rd', 'MarkerFaceColor', 'r')
    plot(median(performanceSummary(:,xx,2)), 0.4, 'rd', 'MarkerFaceColor', 'r')
    fprintf('SNR %d: Accuracy = %.2f [%.2f %.2f], mean NDT = %.2f [%.2f %.2f]\n', ...
        xx, ...
        prctile(performanceSummary(:,xx,1), 50), ...
        prctile(performanceSummary(:,xx,1), 25), ...
        prctile(performanceSummary(:,xx,1), 75), ...
        prctile(performanceSummary(:,xx,2), 50), ...
        prctile(performanceSummary(:,xx,2), 25), ...
        prctile(performanceSummary(:,xx,2), 75) ...
        );
    axis([0 10 0.4 1]);
    if xx == 1
        ylabel('Accuracy');
        xlabel('Mean DT (steps)');
    end
end
fprintf('Accuracy diff hi vs lo SNR, p=%.2f\n', signrank( ...
    performanceSummary(:,1,1), performanceSummary(:,1,1)))
fprintf('DT diff hi vs lo SNR, p=%.2f\n', signrank( ...
    performanceSummary(:,1,2), performanceSummary(:,1,2)))
% 
% 
% SNRs = nonanunique(dataTable.snr);
% aSNRs = 
% numSNRs = 
% summaryData = nan(numSubjects, length(SNRs), 2);
% ndtSummary
% performanceSummary = nan(numSubjects, )
% % Color
% gr = 0.5.*ones(3,1);
% 
% % Loop through the subjects, snrs to collect data
% for ss = 1:numSubjects
%     for aa = 1:length(SNRs)
%         Lsb = dataTable.subjectIndex == ss & dataTable.blockIndex == options.block & dataTable.snr == SNRs(aa);
%         % save p(choose Right)
%         summaryData(ss,aa,1) = sum(dataTable.choice(Lsb))/length(dataTable.choice(Lsb));
%         % save mean RT
%         summaryData(ss,aa,2) = mean(dataTable.RT(Lsb));
%     end
%     axes(axs(1)); hold on;
%     plot(SNRs, summaryData(ss,:,1), 'Color', gr);
%     axes(axs(2)); hold on;
%     plot(SNRs, summaryData(ss,:,2), 'Color', gr);
% end
% axes(axs(1));
% plot(SNRs, mean(summaryData(:,:,1),1), 'k-', 'LineWidth', 3);
% ylabel('p(choose Right)');
% axes(axs(2));
% plot(SNRs, mean(summaryData(:,:,2),1), 'k-', 'LineWidth', 3);
% ylabel('RT');
% xlabel('Signed SNR');
% 
% 