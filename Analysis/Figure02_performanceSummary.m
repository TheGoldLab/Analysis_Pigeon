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
    options.exampleSession=10;
    options.figureNumber=2
end

%% Set up figure
wid     = 8.5; % total width
cols    = {1,1,1};
hts     = [4,4,4];
[axs,~] = getPLOT_axes(options.figureNumber, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
%set(axs,'Units','normalized');

% Useful variables
numSubjects = length(unique(dataTable.subjectIndex));
aSNRs = abs(dataTable.snr);
congruenceSummary = nan(numSubjects, 5); % numNDTs
performanceSummary = nan(numSubjects, 2);  % accuracy/DT
Lgood = dataTable.DT > 0 & dataTable.blockIndex == options.block & aSNRs == min(aSNRs);

% Colors
gr = ones(3,1).*0.8;
wh = ones(3,1).*0.99;

% Loop through the subjects to collect data
for ss = 1:numSubjects

    % Save congruences
    Ls = dataTable.subjectIndex == ss; 
    congruenceSummary(ss,:) = dataTable.congruence{find(Ls,1)}(:,1);
    
    % save p(Correct), DT
    Lsg = Lgood & Ls;
    performanceSummary(ss,:) = [sum(dataTable.correct(Lsg)==1)/sum(dataTable.correct(Lsg)>=0) ...
        median(dataTable.DT(Lsg), 'omitnan')];
end

%% Top row is example session trajectories
%
Ls = Lgood & dataTable.subjectIndex==options.exampleSession;
exampleTable = dataTable(Ls,{'choice', 'correct', 'steps'});
axes(axs(1)); cla reset; hold on;
% draw trajectories
for tt = 1:size(exampleTable,1)
    steps = exampleTable.steps{tt};
    plot(0:length(steps)-1, steps, '-', 'Color', gr);
end
% draw points
for tt = 1:size(exampleTable,1)
    steps = exampleTable.steps{tt};
    h=plot(length(steps)-1,steps(end),'ko','MarkerSize',5);
    if exampleTable.choice(tt) == 0
        set(h, 'MarkerFaceColor', 'r')
    else
        set(h, 'MarkerFaceColor', 'g')
    end
end
plot([0 12], [0 0], 'k:');

%% Middle row is congruence between pigeon positions and choices, per NDT
%
axes(axs(2)); cla reset; hold on;
xax = 0:4;
plot(xax,congruenceSummary(:,:)', '-', 'Color', gr);
plot(xax, median(congruenceSummary, 'omitnan'), 'k-', 'LineWidth', 2);
ylabel('Congruence')
xlabel('NDT (steps)')

%% Bottom row is accuracy vs NDT per subject
%
axes(axs(3)); cla reset; hold on;
plot([0 10], [0.5 0.5], 'k:')
plot(performanceSummary(:,2), performanceSummary(:,1), 'ko', ...
    'MarkerFaceColor', wh);
plot(0, median(performanceSummary(:,1)), 'rd', 'MarkerFaceColor', 'r')
plot(median(performanceSummary(:,2)), 0.4, 'rd', 'MarkerFaceColor', 'r')
fprintf('Accuracy = %.2f [%.2f %.2f], mean NDT = %.2f [%.2f %.2f]\n', ...
    prctile(performanceSummary(:,1), 50), ...
    prctile(performanceSummary(:,1), 25), ...
    prctile(performanceSummary(:,1), 75), ...
    prctile(performanceSummary(:,2), 50), ...
    prctile(performanceSummary(:,2), 25), ...
    prctile(performanceSummary(:,2), 75) ...
    );
axis([0 10 0.4 1]);
ylabel('Accuracy');
xlabel('Median DT (steps)');
[Rlo,Plo] = corr(performanceSummary(:,1), performanceSummary(:,2), 'type', 'Spearman');
fprintf('Spearman acc vs DT=%.2f,p=%.4f\n', Rlo, Plo)

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