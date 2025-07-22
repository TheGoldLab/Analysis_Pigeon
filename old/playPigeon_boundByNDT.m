% playpen for testing bound estimates per NDT

% Real data
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
%   12. congruence (per snr)
% Get real/simulated data
dataTables = cell(2,1);
dataTables{1} = getPigeon_dataTable('taskType', 'OL'); % 'OL' or 'MX'
dataTables{2} = getPigeon_simulatedDataTable(dataTables{1}, 'boundType', 'true');

% collect and show bound data
subjectIndices = unique(dataTables{1}.subjectIndex);
numSubjects = length(subjectIndices);
uniqueSNRs = unique(abs(dataTables{1}.snr));
numSNRs = length(uniqueSNRs);
numRTs = 10;
numNDTs = 3;
boundData = nan(numSubjects,numNDTs,numSNRs,2); % last is real/sim
blockIndex = 2;

% Loop throught the subjects
for dd = 1:2 % real/sim
    for ss = 1:numSubjects
        Ls = dataTables{dd}.subjectIndex == subjectIndices(ss) & dataTables{dd}.blockIndex==blockIndex;
        for mm = 1:2 % snrs
            % get the data
            Lm = Ls & abs(dataTables{dd}.snr)==uniqueSNRs(mm);
            bounds = abs(dataTables{dd}.bound(Lm));
            steps = dataTables{dd}.steps(Lm);
            numTrials = sum(Lm);
            boundPerTrial = nan(numTrials,numNDTs);
            for tt = 1:numTrials
                trialSteps = flip(steps{tt});
                for nn = 1:numNDTs
                    if length(trialSteps)>=nn+1
                        boundPerTrial(tt,nn) = mean(trialSteps(nn:nn+1).*sign(trialSteps(nn)));
                    end
                end
            end
            boundData(ss,:,mm,dd) = mean(boundPerTrial,'omitnan');
        end
    end
end

wh = ones(3,1).*0.99;
for dd = 1:2
    for nn = 1:numNDTs
        subplot(2,numNDTs,(dd-1)*numNDTs+nn); cla reset; hold on;
        plot([0 0.8], [0 0.8], 'k:');
        plot(boundData(:,nn,1,dd), boundData(:,nn,2,dd), 'ko', 'MarkerFaceColor', wh);
    end
end
