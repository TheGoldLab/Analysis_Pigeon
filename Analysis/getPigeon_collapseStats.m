%% Compare collapse, blocked vs mixed
dataTableMXall = getPigeon_dataTable('taskType', 'MX', 'combineSNR', false); 

%% Collect data per subject/block
blocks = [2 5];
numBlocks = length(blocks);
subjects = nonanunique(dataTableMXall.subjectIndex);
numSubjects = length(subjects);
abounds = abs(dataTableMXall.bound); % absolute values of bounds
regressionData = nan(numSubjects, numBlocks);
blockData = nan(numSubjects,2);
Lg = isfinite(dataTableMXall.bound) & dataTableMXall.DT >= 2;

for ss = 1:numSubjects
    Ls = Lg & dataTableMXall.subjectIndex == subjects(ss);
    for bb = 1:numBlocks
        Lb = Ls & dataTableMXall.blockIndex==blocks(bb);

         % regression bound vs RT
         if sum(Lb) > 5
             xs = ones(sum(Lb),2);
             xs(:,2) = dataTableMXall.DT(Lb);
             b = regress(abounds(Lb), xs);
             regressionData(ss,bb) = b(2);
         end

         % Save median bound
         blockData(ss,bb) = median(abounds(Lb));
    end
end

% Compare bounds
Lg = isfinite(blockData(:,1)) & isfinite(blockData(:,2));
signrank(blockData(Lg,1), blockData(Lg,2))

% Compare slopes
Lg = isfinite(regressionData(:,1)) & isfinite(regressionData(:,2));
signrank(regressionData(Lg,1), regressionData(Lg,2))
disp(prctile(regressionData(Lg,1), [50 25 75]))
disp(prctile(regressionData(Lg,2), [50 25 75]))

% Compare slopes to zero
signrank(regressionData(Lg,1))
signrank(regressionData(Lg,2))
