%% Script to generate figures
%

% Real data
%
dataTableOL = getPigeon_dataTable('taskType', 'OL'); % 'OL' or 'MX'
dataTableMX = getPigeon_dataTable('taskType', 'MX'); 
dataTableMXall = getPigeon_dataTable('taskType', 'MX', 'combineSNR', false); 

% Simulated data
% 
% Includes random non-decision time of 1-3 steps
% Uses median/STD bounds from each subject for simulations
simDataTableOL = getPigeon_simulatedDataTable(dataTableOL, 'boundType', 'true');
simDataTableMX = getPigeon_simulatedDataTable(dataTableMX, 'boundType', 'true');

% Figure 2: psychometric and chronometric curves
Figure02_psychoChronometricCurves(dataTableMX)

% Figure 3: summary of performance using mixed block 2
Figure03_mixedBoundSummary(dataTableMXall, 'block', 2)

% Figure 4: summary of snr-dependent bounds in fixed/mixed blocks
Figure04_mixedVsBlockedSNR(dataTableMX, dataTableOL);

% Figure 4.1: summary of snr-dependent bounds in fixed/mixed blocks with
% SIM DATA
Figure04_mixedVsBlockedSNR(simDataTableMX, simDataTableOL, 'figureNumber', 4.1);

% Figure 5: summary of bounds vs RR in different blocks
Figure05_RRvsBound(dataTableOL);

% Figure 6: time course of transition block 1->2 for OL and MX blocks
Figure06_boundBlockTransition(dataTableOL, dataTableMX);
