function [dataTable_] = getPigeon_dataTable_MXSC(options)
% function [dataTable_] = getPigeon_dataTable(options)
%
% Returns dataTable with the following columns:
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
    options.taskType = 'MXSC'    % 'OL', 'MX', 'SC', 'MXSC', or 'PD'
    options.dataDir = ...
        fullfile('/Users', 'ishankalburge', 'Documents', 'penn', 'Pigeon', 'Data');
    options.blocks = 'all';
    options.combineSNR = true; % default to 3 blocks MX/OL snrs
end

% Set up the data table
dataTable_ = getPigeon_blankDataTable(0);

% Loop through the file list
files = dir(fullfile(options.dataDir, ['Pigeon_' options.taskType], 'prolificcsvs', '*.csv'));
saveParticipantID = cell(length(files),1);
saveCoinCount = zeros(length(files),1);
for ff = 1:length(files)

    fprintf('% 1d: %s\n\n', ff, fullfile(files(ff).folder, files(ff).name))

    % Read data into Table
    tmpTable = readtable(fullfile(files(ff).folder, files(ff).name));

    % Save participant ID for printing later
    saveParticipantID{ff} = string(tmpTable.PROLIFIC_PID{1});

    % Use only good entries from selected blocks
    Lgood = isfinite(tmpTable.trial_number);
    if isnumeric(options.blocks) && sum(options.blocks) ~= sum(1:6)
        Lgood = Lgood & ismember(tmpTable.block_number, options.blocks);
    end

    % Get the columns we want
    if strcmp(options.taskType, 'OL') || contains(options.taskType, 'MX')
        tmpTable = tmpTable(Lgood, ...
            {'block_number', 'trial_number', 'pigeon_steps', 'choice', ...
            'correct' 'steps_count', 'step_mean_val', 'step_std_val', 'coins_count'});
    else
        indices = find(Lgood);
        tmpTable = cat(2, ...
            tmpTable(indices, {'block_number', 'trial_number', 'predefined_bound', 'correct'}), ...
            tmpTable(indices+1, {'pigeon_steps', 'choice', 'coins_count'}));
    end

    % Add the data
    numTrials = size(tmpTable, 1);
    tableToAppend = getPigeon_blankDataTable(numTrials);
    tableToAppend.subjectIndex(1:end) = ff;
    tableToAppend.blockIndex(1:end) = tmpTable.block_number(1:end);
    tableToAppend.trialNumber(1:end) = tmpTable.trial_number(1:end);
    tableToAppend.snr(1:end) = tmpTable.step_mean_val(1:end)./tmpTable.step_std_val(1:end);
    
    % default changepoint for SC tasks
    if contains(options.taskType, 'SC')
        tableToAppend.changePoint(1:end) = 0;
    end

    % Loop through the trials to fill in the blanks
    for tt = 1:numTrials

        % Get steps
        tableToAppend.steps{tt} = str2num(tmpTable.pigeon_steps{tt});

        % Get bound, rt for PD task
        if strcmp(options.taskType, 'PD')
            tableToAppend.bound(tt) = tmpTable.predefined_bound(tt);
            tableToAppend.RT(tt) = length(tableToAppend.steps(tt));
        end

        % Choice/outcome
        if ismember(tmpTable.choice{tt}, {'right', 'left'})
            tableToAppend.choice(tt) = double(strcmp(tmpTable.choice{tt}, 'right'));
            tableToAppend.correct(tt) = double(strcmp(tmpTable.choice{tt}, ...
                tmpTable.correct{tt}));
        end

        % Coin count
        tableToAppend.coinCount(tt) = tmpTable.coins_count(tt);
    end

    

    % Get NDT, bounds for on-line tasks
    if strcmp(options.taskType, 'OL') || contains(options.taskType, 'MX')
        [tableToAppend.bound(1:end), ...
            tableToAppend.DT(1:end), ...
            tableToAppend.RT(1:end)] = ...
            getPigeon_bounds(tableToAppend.steps, tableToAppend.choice);
    end
    tableToAppend.absBound(1:end) = abs(tableToAppend.bound(1:end));


    % save change-point for SC task
    if contains(options.taskType, 'SC')
        Lblock  = tableToAppend.blockIndex == 1;
        Hblock  = tableToAppend.blockIndex == 2;
        LhBlock = tableToAppend.blockIndex == 4;
        HlBlock = tableToAppend.blockIndex == 5;
        tableToAppend.changePoint(LhBlock) = round(mean(tableToAppend.RT(Lblock)))*ones(size(tableToAppend.changePoint(LhBlock)));
        tableToAppend.changePoint(HlBlock) = round(mean(tableToAppend.RT(Hblock)))*ones(size(tableToAppend.changePoint(HlBlock)));
    end

    % add the data
    dataTable_ = cat(1, dataTable_, tableToAppend);


    % save final scores / bonuses for printing later
    saveCoinCount(ff) = max(tableToAppend.coinCount);
end


% Possibly set to 3 blocks of mixed or blocked snr
if length(unique(dataTable_.blockIndex))== 6 && all(unique(dataTable_.blockIndex)==(1:6)') && options.combineSNR

    % check for block type 
    if isscalar(nonanunique(abs(dataTable_.snr(dataTable_.blockIndex == 1))))
        % blocked, change block number 4:6 to 1:3
        for bb = 1:3
            dataTable_.blockIndex(dataTable_.blockIndex==bb+3)=bb;
        end
    else
        % mixed, remove blocks 4:6
        dataTable_ = dataTable_(dataTable_.blockIndex<=3,:);
    end
end

% for bulk approval on Prolific
for ff = 1:length(files)
    fprintf("%s,\n", saveParticipantID{ff});
end

% for bulk bonus on Prolific
for ff = 1:length(files)
    % give bonus of $2.00 if score equals or exceeds 100, $1.00 if score equals or exceeds 80
    % (for this specific set of task conditions, i.e. 1000 total steps, 50
    % steps per trial, 1 coin gained for correct and 1 lost for incorrect)
    bonus = 1.00 * (saveCoinCount(ff) >= 100) + 1.00 * (saveCoinCount(ff) >= 80);
    if bonus > 0
        fprintf("%s,%.1f\n", saveParticipantID{ff}, bonus);
    end
end
