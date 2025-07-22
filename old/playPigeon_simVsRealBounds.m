dataTableOL = getPigeon_dataTable('taskType', 'OL'); % 'OL' or 'MX'
dataTableMX = getPigeon_dataTable('taskType', 'MX');
dataTableMXall = getPigeon_dataTable('taskType', 'MX', 'combineSNR', false);
dataTableOLall = getPigeon_dataTable('taskType', 'OL', 'combineSNR', false); % 'OL' or 'MX'

% Simulated data
%
% Includes random non-decision time of 1-3 steps
% Uses median/STD bounds from each subject for simulations
simDataTableOL = getPigeon_simulatedDataTable(dataTableOL, 'boundType', 'true');
simDataTableMX = getPigeon_simulatedDataTable(dataTableMX, 'boundType', 'true');
simDataTableMXall = getPigeon_simulatedDataTable(dataTableMXall, 'boundType', 'true');
simDataTableOLall = getPigeon_simulatedDataTable(dataTableOLall, 'boundType', 'true');

%% Compare per RT
%   boundSummary_ .. matrix with dims:
%       subjects
%       blocks
%       snrs
%       rt bins
%       bound avg, bound std, n, mean zbound
dataOLSummaryByRT = getPigeon_boundSummary(dataTableOLall);
simOLSummaryByRT = getPigeon_boundSummary(simDataTableOLall);

figure
wh = 0.99.*ones(1,3);
blocks = [2 5];
snrs = [1 2];
for bb = 1:length(blocks)
    for rr = 1:10
        subplot(2,10,(bb-1)*10+rr); cla reset; hold on;
        plot([0 1], [0 1], 'k:')
        xs = dataOLSummaryByRT(:,blocks(bb),snrs(bb),rr,1);
        ys = simOLSummaryByRT(:,blocks(bb),snrs(bb),rr,1);
        plot(xs, ys, 'ko', 'MarkerFaceColor', wh);
        Lg = isfinite(xs) & isfinite(ys);
        if sum(Lg) > 3
            title(sprintf('mr=%.2f,ms=%.2f,p=%.3f', ...
                mean(xs(Lg)), mean(ys(Lg)), signrank(xs(Lg), ys(Lg))))
        end
    end
end

%% Compare overall
dataOLSummary = getPigeon_boundSummary(dataTableOLall, 'maxRT', 'all');
simOLSummary = getPigeon_boundSummary(simDataTableOLall, 'maxRT', 'all');

%figure
wh = 0.99.*ones(1,3);
blocks = [2 5];
for bb = 1:length(blocks)
    subplot(2,1,bb); cla reset; hold on;
    plot([0 0.8], [0 0.8], 'k:')
    xs = dataOLSummary(:,blocks(bb),bb,1,1);
    ys = simOLSummary(:,blocks(bb),bb,1,1);
    plot(xs, ys, 'ko', 'MarkerFaceColor', wh);
    Lc = dataOLSummary(:,blocks(bb),bb,1,4) > 0.8;
    plot(xs(Lc), ys(Lc), 'ko', 'MarkerFaceColor', 'k');
    Lg = isfinite(xs) & isfinite(ys);
    title(sprintf('p=%.3f', signrank(xs(Lc&Lg), ys(Lc&Lg))))
end



bb=2;
for ss = 1:size(dataOLSummary,1)
    plot(dataOLSummary(ss,blocks(bb),bb,1,1), simOLSummary(ss,blocks(bb),bb,1,1), 'go', 'MarkerFaceColor', 'g');
    disp(ss)
    r = input('next')
end

bi = 5; % block index
for xx = 1:2
    if xx == 1
        dt = dataTableOLall;
    else
        dt = simDataTableOLall;
    end
    figure
    subjects = unique(dt.subjectIndex);
    numSubjects = length(subjects);
    for ss = 1:numSubjects
        subplot(6,10,ss)
        Lsub = dt.subjectIndex==subjects(ss) & dt.blockIndex==bi;
        h=histogram(abs(dt.bound(Lsub)), 0:0.05:0.8);
        if median(dt.congruence(Lsub))<0.8
            set(h, 'FaceColor', 'r')
        end
        title(sprintf('%d', ss))
    end
end
