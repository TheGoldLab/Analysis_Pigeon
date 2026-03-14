%% Script to generate figures
%

%% Possibly make bias correction
%
makePigeon_biasCorrections('showOutput', true);%, 'filename', ''); % just make the figure

%% Real data
%
dataTableOL = getPigeon_dataTable('taskType', 'OL'); % 'OL' or 'MX'
dataTableMX = getPigeon_dataTable('taskType', 'MX'); 
dataTableMXall = getPigeon_dataTable('taskType', 'MX', 'combineSNR', false); 

%% Simulated data
% 
% Uses median/STD bounds from each subject for simulations
simDataTableOL = getPigeon_simulatedDataTable(dataTableOL, 'boundType', 'true');
simDataTableMX = getPigeon_simulatedDataTable(dataTableMX, 'boundType', 'true');
simDataTableMXall = getPigeon_simulatedDataTable(dataTableMXall, 'boundType', 'true');

%% Figures
%

% Figure 2: Summary of accuracy and DT
Figure02_performanceSummary(dataTableOL)
% Figure02_performanceSummary(simDataTableOL,'figureNumber','2S')

% Figure 3: summary of performance using block 2
Figure03_boundSummary(dataTableOL)
% Figure03_mixedBoundSummary(simDataTableOL, 'blockIndex', 2,'figureNumber','2S')

% Figure 4: summary of bounds vs RR in different blocks
Figure04_RRvsBound(dataTableOL);

% Figure 5: summary of snr-dependent bounds in fixed/mixed blocks
Figure05_mixedVsBlockedSNR(dataTableMX, dataTableOL); %, 'showRR', false);
% Figure05_mixedVsBlockedSNR(simDataTableMX, simDataTableOL, 'showRR', false);

% Figure 6: change-point summary
Figure06_cpSummary();


%% compare sampling noise, real vs simulated bounds
%
subjects = nonanunique(dataTableOL.subjectIndex);
numSubjects = length(subjects);
stds = nan(numSubjects, 2);

% Real
aboundsReal = abs(dataTableOL.bound); % absolute values of bounds
LgReal = isfinite(aboundsReal) & dataTableOL.blockIndex==2 & dataTableOL.DT > 2 & dataTableOL.DT < 8;

% Sim
aboundsSim = abs(simDataTableOL.bound); % absolute values of bounds
LgSim = isfinite(aboundsSim) & simDataTableOL.blockIndex==2 & simDataTableOL.DT > 2 & simDataTableOL.DT < 8;

for ss = 1:numSubjects
    LsReal = LgReal &    dataTableOL.subjectIndex == subjects(ss);
    LsSim  = LgSim  & simDataTableOL.subjectIndex == subjects(ss);
    stds(ss,:) = [std(aboundsReal(LsReal)), std(aboundsSim(LsSim))];
end

Lg = isfinite(stds(:,1)) & isfinite(stds(:,2));
signrank(stds(Lg,1), stds(Lg, 2))