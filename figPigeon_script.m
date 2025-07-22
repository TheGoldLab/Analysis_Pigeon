%% Script to generate figures
%

%% Possibly make bias correction
%
makePigeon_biasCorrections('showOutput', true, 'filename', ''); % just make the figure

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

% Figure 3: summary of performance using mixed block 2
Figure03_mixedBoundSummary(dataTableOL, 'blockIndex', 2)

% Figure 4: summary of snr-dependent bounds in fixed/mixed blocks
Figure04_mixedVsBlockedSNR(dataTableMX, dataTableOL); %, 'showRR', false);

% Figure 5: summary of bounds vs RR in different blocks
Figure05_RRvsBound(dataTableOL);

