function Figure1_figPigeon_PsychoChronometricCurves(dataTableOL)
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
    dataTableOL
end

%% Set up figure
wid     = 13; % total width
cols    = {1,1};
hts     = [6.5 6.5];
[axs,~] = getPLOT_axes(1, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
set(axs,'Units','normalized');

blockNum = 2;
numSubjects = length(unique(dataTableOL.subjectIndex));
aSNRs = nonanunique(dataTableOL.snr);
sData = nan(numSubjects, length(aSNRs), 2);

for ss = 1:numSubjects
    for aa = 1:length(aSNRs)
        Lsb = dataTableOL.subjectIndex == ss & dataTableOL.blockIndex == blockNum & dataTableOL.snr == aSNRs(aa);
        % save p(choose Right)
        sData(ss,aa,1) = sum(dataTableOL.choice(Lsb))/length(dataTableOL.choice(Lsb));
        % save mean RT
        sData(ss,aa,2) = mean(dataTableOL.RT(Lsb));
    end
    axes(axs(1)); hold on;
    plot(aSNRs, sData(ss,:,1), "Color", [0.5 0.5 0.5]);
    axes(axs(2)); hold on;
    plot(aSNRs, sData(ss,:,2), "Color", [0.5 0.5 0.5]);
end
axes(axs(1));
plot(aSNRs, mean(sData(:,:,1),1), "Color", [0 0 0], "LineWidth", 3);
ylabel("p(choose Right)");
axes(axs(2));
plot(aSNRs, mean(sData(:,:,2),1), "Color", [0 0 0], "LineWidth", 3);
ylabel("RT");
xlabel("Signed SNR");

    