function Figure02_psychoChronometricCurves(dataTable, options)
% function Figure02_psychoChronometricCurves(dataTable, options)
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
wid     = 6.5; % total width
cols    = {1,1};
hts     = [5 5];
[axs,~] = getPLOT_axes(options.num, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
set(axs,'Units','normalized');

% Useful variables
numSubjects = length(unique(dataTable.subjectIndex));
SNRs = nonanunique(dataTable.snr);
summaryData = nan(numSubjects, length(SNRs), 2);

% Color
gr = 0.5.*ones(3,1);

% Loop through the subjects, snrs to collect data
for ss = 1:numSubjects
    for aa = 1:length(SNRs)
        Lsb = dataTable.subjectIndex == ss & dataTable.blockIndex == options.block & dataTable.snr == SNRs(aa);
        % save p(choose Right)
        summaryData(ss,aa,1) = sum(dataTable.choice(Lsb))/length(dataTable.choice(Lsb));
        % save mean RT
        summaryData(ss,aa,2) = mean(dataTable.RT(Lsb));
    end
    axes(axs(1)); hold on;
    plot(SNRs, summaryData(ss,:,1), 'Color', gr);
    axes(axs(2)); hold on;
    plot(SNRs, summaryData(ss,:,2), 'Color', gr);
end
axes(axs(1));
plot(SNRs, mean(summaryData(:,:,1),1), 'k-', 'LineWidth', 3);
ylabel('p(choose Right)');
axes(axs(2));
plot(SNRs, mean(summaryData(:,:,2),1), 'k-', 'LineWidth', 3);
ylabel('RT');
xlabel('Signed SNR');

    