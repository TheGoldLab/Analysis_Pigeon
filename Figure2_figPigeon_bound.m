function Figure2_figPigeon_bound(dataTable, exampleSubject, block, num)
% function figPigeon_bound(dataTable, exampleSubject, block, num)
%
% dataTable columns:
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
    exampleSubject = [1, 2, 3, 4]
    block = 3
    num = 2
end

%% Set up figure
%
wid     = 17.6; % total width
cols    = {4,4,2};
ht      = 5;
[axs,~] = getPLOT_axes(num, wid, ht, cols, 1.3, 1.5, [], 'Pigeons', true);
set(axs,'Units','normalized');
wt = ones(3,1).*0.99;

%% Collect data per subject/block
subjects = nonanunique(dataTable.subjectIndex);
numSubjects = length(subjects);

% Need to find the appropriate selection criteria -- right now these
% are pretty arbitrary
Lg = dataTable.DT>1 & dataTable.trialNumber>30 & dataTable.blockIndex==block;

%% Column 1 is bound vs DT
% Loop through each subject
bData = nan(numSubjects,3); % 2: median bound/slope beta/slope p
sData = nan(numSubjects,4); % mean-bound std-bound mean-DT std-DT
for ss = 1:numSubjects
    % Get subject-specific data
    Lsb = Lg & dataTable.subjectIndex==subjects(ss);

    if ismember(ss, exampleSubject)
        axes(axs(ss)); cla reset; hold on;
        scatter(dataTable.DT(Lsb), abs(dataTable.bound(Lsb)), 'ko', ...
            'SizeData', 20, 'MarkerFaceColor', wt);
        lsline
        xlabel("RT")
        if ss == 1
            ylabel("Bound");
        end
        %disp(ss)
        %r = input('next')
    end

    if sum(Lsb) > 4

        [B,~,~,~,STATS] = regress(abs(dataTable.bound(Lsb)), ...
            cat(2, ones(sum(Lsb),1), dataTable.DT(Lsb)));

        bData(ss,:) = [median(abs(dataTable.bound(Lsb))), B(2), STATS(3)];
        sData(ss,:) = [mean(abs(dataTable.bound(Lsb))), std(abs(dataTable.bound(Lsb))/length(dataTable.bound(Lsb))), ... 
                       mean(abs(dataTable.DT(Lsb))), std(abs(dataTable.DT(Lsb)))/length(dataTable.DT(Lsb))];
    end
end
axes(axs(9)); cla reset; hold on;
% plot([0 0.8], [0 0], 'k:');
errorbar(sData(:,3), sData(:,1), sData(:,2), sData(:,2), sData(:,4), sData(:,4), 'ko', 'MarkerFaceColor', wt);
xlabel("RT"); ylabel("RT")

% axes(axs(10)); cla reset; hold on;
% plot([0 0.8], [0 0], 'k:');
% plot(bData(:,1), bData(:,2), 'ko', 'MarkerFaceColor', wt);

%% Column 2 is Pigeon step size vs steps before DT
% Loop through each subject
nSteps = 10;
sData = nan(numSubjects,nSteps); % mean step size
vData = nan(numSubjects,2); % var bound, var final step size
for ss = 1:numSubjects
    % Get subject-specific data
    Lsb = Lg & dataTable.subjectIndex==subjects(ss);

    steps = dataTable.steps(Lsb);
    ndts = dataTable.RT(Lsb) - dataTable.DT(Lsb);
    dirs = dataTable.choice(Lsb) == dataTable.correct(Lsb);
    nTrials = length(ndts);
    tmpData = nan(nTrials,nSteps);
    for ii = 1:nTrials
        trialSteps = steps{ii}(1:end-ndts(ii));
        if dirs(ii) == 0
            trialSteps = -trialSteps;
        end        
        tmpData(ii,1:length(trialSteps)-1) = fliplr(diff(trialSteps));
    end
    for tt = 1:nSteps
        sData(ss,tt) = mean(tmpData(:,tt), 'omitnan');
    end

    % Var bound, final step
    vData(ss,:) = [ ...
        var(dataTable.bound(Lsb), 'omitnan'), ...
        var(diff(tmpData(:,1:2),[],2),'omitnan')];

    % if ss == exampleSubject
    %     axes(axs(2)); cla reset; hold on;
    %     plot([0 11], [0 0], 'k:');
    %     plot(sData(ss,:,1), 'ko', 'MarkerSize', 5)
    % end
end
% axes(axs(11)); cla reset; hold on;
% plot([0 11], [0 0], 'k:');
% boxplot(sData(:,:,1))
% ylim([-0.2 0.2])

Lgv = isfinite(vData(:,1)) & isfinite(vData(:,2));
vb = prctile(vData(Lgv,1),[25 50 75]);
vs = prctile(vData(Lgv,2),[25 50 75]);
fprintf('Median [IQR] Variance of bound = %.2f [%.2f %.2f], of step = %.2f [%.2f %.2f]\n', ...
    vb(2), vb(1), vb(3), vs(2), vs(1), vs(3))

%% Column 3 is regression to mean
% Loop through each subject
rData = nan(numSubjects, 3); % last is median bound, slope/intercept
for ss = 1:numSubjects
    % Get subject-specific data
    Lsb = Lg & dataTable.subjectIndex==subjects(ss);

    bound = abs(dataTable.bound(Lsb));
    if length(bound) > 30
        bound = bound(20:end);
        rData(ss,1) = median(bound);
        deltaBound = diff(bound);
        bound = bound(1:end-1);
        rData(ss,2:3) = regress(deltaBound, ...
            cat(2, ones(size(bound)), bound));

        if ismember(ss, exampleSubject)
            axes(axs(ss+length(exampleSubject))); cla reset; hold on;
            plot([0 0.8], [0 0], 'k:');
            plot(bound, deltaBound, 'ko', ...
                'MarkerSize', 5, 'MarkerFaceColor', wt);
            lsline
            xlabel("Bound");
            if ss == 1
                ylabel("\Delta Bound");
            end
            %disp(ss)
            %r = input('next')
        end
    end
end

axes(axs(10)); cla reset; hold on;
plot([0 0.8], [-1 -1], 'k:')
plot(rData(:,1), rData(:,3), 'ko', 'MarkerFaceColor', wt);
axis([0 0.8 -2 0])
xlabel("Bound");
ylabel("\beta_{\Delta bound}");
%ylim([-0.2 0.2])

