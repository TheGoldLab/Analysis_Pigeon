function Figure06_cpSummary(options)
% function Figure02_performanceSummary(options)
%
% Notes:
% The generative mean was 0.01 and 0.02 for low and high SNR conditions 
% respectively, and the generative std was always 0.05. So, for instance, 
% in LH change trials, we would go from 0.01 to 0.02 at the CP.
arguments
    options.makeRR=false;
    options.showRR=true;
    options.figureNumber=6
end

wid     = 8.5; % total width
cols    = {2,2,2,2};
hts     = [4,2,4,4];
[axs,~] = getPLOT_axes(options.figureNumber, wid, hts, cols, 2.0, 1.0, [], 'Pigeons', true);
% set(axs,'Units','normalized');

% Get the data
load('boundSummaries_wCP.mat');

% Get the RR matrices
if options.makeRR
    RRMatrixLowHigh = getPigeon_RRMatrix( ...
        'numReps',      500, ...
        'gMeans',       [0.01 0.02], ...
        'gSTD',         0.05, ...
        'blockType',    'CP', ...
        'blocks',       unique(boundSummaries.CPLowHigh));

    RRMatrixHighLow = getPigeon_RRMatrix( ...
        'numReps',      500, ...
        'gMeans',       [0.02 0.01], ...
        'gSTD',         0.05, ...
        'blockType',    'CP', ...
        'blocks',       unique(boundSummaries.CPHighLow));

    pmats = {RRMatrixLowHigh, RRMatrixHighLow};
    save('RRMats.mat', 'pmats')
else
    load('RRMats.mat');
end

bounds = 0.01:0.05:0.8;
cpTimes = {unique(boundSummaries.CPLowHigh), unique(boundSummaries.CPHighLow)};

fprintf('LowHigh, median [IQR] cp times = %.1f [%.1f %.1f]\n', ...
    prctile(boundSummaries.CPLowHigh, 50), ...
    prctile(boundSummaries.CPLowHigh, 25), ...
    prctile(boundSummaries.CPLowHigh, 75));
fprintf('HighLow, median [IQR] cp times = %.1f [%.1f %.1f]\n', ...
    prctile(boundSummaries.CPHighLow, 50), ...
    prctile(boundSummaries.CPHighLow, 25), ...
    prctile(boundSummaries.CPHighLow, 75));

%% Row 1: RR functions
inds = [3 5];
co = {0.7.*ones(3,1) zeros(3,1)};
titles = {'LowHigh', 'HighLow'};
xx = 7; % cp time index
for bb = 1:2
    axes(axs(bb)); cla reset; hold on;
    title(titles{bb})
    for ii = 1:length(inds)
        ys = nanrunmean(pmats{bb}(inds(ii),:,xx),1);
        mxi = find(ys==max(ys),1);
        plot(bounds, ys, '-', 'Color', co{ii}, 'LineWidth', 2);
        plot(bounds([mxi mxi]), [min(ys) max(ys)], '-', 'Color', co{ii});
        plot(bounds(inds(ii)), ys(inds(ii)), 'o', ...
            'Color', co{ii}, 'MarkerFaceColor', co{ii}, 'MarkerSize', 15);
    end
    if bb == 1
        xlabel('Bound')
        ylabel('Reward Rate (coins/trial)')
    end
end

%% Row 2: example RR gradients
%% Row 3: post vs pre bound cp trials, color coded by optimal direction
titles = {'LowHigh', 'HighLow'};
for tt = 1:length(titles)

    % Get data
    xs = boundSummaries.(titles{tt})(:,1);
    ys = boundSummaries.(titles{tt})(:,2);
    Lg = isfinite(xs) & isfinite(ys);

    % For CP trials, find optimal direction and color code
    cpt = boundSummaries.(['CP' titles{tt}]);
    % LoDirs = false(length(cpt),1);
    grads = nan(length(cpt),1);
    for pp = 1:length(cpt)
        if(Lg(pp))
            % Get RR matrix for the given cp time
            rmat = pmats{tt}(:,:,cpt(pp)==cpTimes{tt});
            % Get RR vs bound for the given pre-cp bound
            prei = find(bounds<=xs(pp),1,'last');
            rrs = rmat(prei,:);           
            maxi = find(rrs == max(rrs), 1);

            % Find gradient to RR max
            if maxi == prei
                grads(pp) = 0;
            % elseif maxi < prei
            %     grads(pp) = (rrs(prei-1)-rrs(prei))./(bounds(prei-1)-bounds(prei));
            % else % maxi > prei
            %     grads(pp) = (rrs(prei+1)-rrs(prei))./(bounds(prei+1)-bounds(prei));
            else
                grads(pp) = (rrs(maxi)-rrs(prei))./(bounds(maxi)-bounds(prei));
            end

            if (tt==1 && pp==53) || (tt==2 && pp==52)
                axes(axs(tt)+2); cla reset; hold on;
                rrs = nanrunmean(rrs,1);
                plot(bounds, rrs, 'k-');
                plot(bounds(prei), rrs(prei), 'bx')
                plot(bounds(maxi), rrs(maxi), 'rx')
                posti = find(bounds<=ys(pp),1,'last');
                plot(bounds(posti), rrs(posti), 'gx')
            end

            % find bound associated with max RR, take diff from pre-cp bound
            %rdiff = bounds(find(rrs == max(rrs), 1))-xs(pp);
            % mark if change post-cp is in the direction of max RR
            %LoDirs(pp) = sign(ys(pp)-xs(pp)) == sign(rdiff);% && abs(rdiff) > 0.005;
        end
    end

    % Plot data
    axes(axs(tt)+4); cla reset; hold on;
    % plot([0 0.5], [0 0], 'k:')
    % yv = abs(ys - xs);
    % yv(~LoDirs) = -yv(~LoDirs);
    % plot(xs(Lg), yv(Lg), 'ko', 'MarkerFaceColor', 0.99.*ones(3,1))
    plot([-50 300], [0 0], 'k:')
    plot(grads(Lg), ys(Lg)-xs(Lg), 'ko', 'MarkerFaceColor', 0.99.*ones(3,1))
    [R,P] = corr(grads(Lg), ys(Lg)-xs(Lg), 'type', 'Spearman');
    title(sprintf('%s: r=%.3f, p=%.4f', titles{tt}, R, P));
    axis([-50 300 -0.3 0.3])
    if tt == 1
        xlabel('RR gradient')
        ylabel('Bound change')
    end

    % plot(xs(Lg), ys(Lg)-xs(Lg), 'ko', 'MarkerFaceColor', 0.99.*ones(3,1))
    % plot(xs(Lg&LoDirs), ys(Lg&LoDirs)-xs(Lg&LoDirs), 'ko', 'MarkerFaceColor', 'k');
    
    % Binomial test
    % p = binopdf(sum(LoDirs), sum(Lg), 0.5);
    % title(sprintf('%s: p=%.4f', titles{tt}, p));
    % title(sprintf('%s: p=%.4f', titles{tt}, signrank(yv(Lg))));
    % axis([0 0.48 -0.3 0.3])
    % if tt == 1
    %     xlabel('Bound pre CP')
    %     ylabel('Bound post CP')
    % end
end

%% Plot post-cp diff, Low->High vs High->Low, Low->Low vs High->High
tags = {'LowHigh', 'HighLow'; 'Low', 'High'};
for tt = 1:size(tags,2)
    xs = boundSummaries.(tags{tt,1})(:,1);
    ys = boundSummaries.(tags{tt,2})(:,1);
    Lg = isfinite(xs) & isfinite(ys);
    axes(axs(tt+6)); cla reset; hold on;

    % check to show RR
    if tt==2 && options.showRR
        RRMatrix = getPigeon_RRMatrix( ...
            'numReps',      500, ...
            'gMeans',       [0.01 0.02], ...
            'gSTD',         0.05, ...
            'blocks',       7);

        cla reset; hold on;
        imagesc(bounds, flip(bounds), flipud(RRMatrix), ...
            [0 max(RRMatrix(:))]);
        set(gca,'YDir', 'normal')
        rr=RRMatrix(:,:,1);
        %[i,j] = find(rr>=0.95.*max(rr(:)));
        %plot(bounds(j), bounds(i), 'kx')
        [i,j] = find(rr==max(rr(:)),1);
        plot(bounds(j), bounds(i), 'rx')
        axis(bounds([1 end 1 end]))
    end
    plot([0 0.8], [0 0.8], 'k:')
    plot(xs(Lg), ys(Lg), 'ko', 'MarkerFaceColor', 0.99.*ones(3,1))
    title(sprintf('p=%.3f', signrank(xs(Lg), ys(Lg))))
    xlabel(tags{tt,1})
    ylabel(tags{tt,2})
    axis([0 0.6 0 0.7])
end

% %% To show RR matrices
% bounds = 0.01:0.05:0.8;
% for pp = 1:15
%     for xx = 1:2
%         rmat = pmats{xx}(:,:,pp);
%         rmat = imgaussfilt(rmat,1);
%         % rmat(rmat<0)=0;
%         subplot(1,2,xx); cla reset; hold on;
%         plot(1,1,'k.');
%         imagesc(bounds,bounds,rmat')
%         %[i,j] = find(rmat>=0.9.*max(rmat(:)));
%         %plot(bounds(i),bounds(j),'ro','MarkerFaceColor','r','MarkerSize',14)
%         [i,j] = find(rmat==max(rmat,[],2));
%         plot(bounds(i),bounds(j),'ro','MarkerFaceColor','r','MarkerSize',14)
%         plot([0 0.8], [0 0.8], 'k--', 'LineWidth', 2)
%         plot([0 0.4], [0.4 0.4], 'k--', 'LineWidth', 2)
%         plot([0.4 0.4], [0 0.4], 'k--', 'LineWidth', 2)
%         axis([0 0.8 0 0.8])
%         xlabel('Pre bound')
%         ylabel('Post bound')
%         title(sprintf('%d', pp))
%     end
%     r = input('next')
% end

% co = {'r' 'g' 'b' 'k' 'c' 'm' 'y'};
% ps = 5:9;
% nps = length(ps);
% for xx = 1:2
%     subplot(1,2,xx); cla reset; hold on;
%     plot([0 0.8], [0 0], 'k:');
%     for pp = 1:nps
%         rmat = pmats{xx}(:,:,ps(pp));
%         rmat = imgaussfilt(rmat,1);
%         [i,j] = find(rmat==max(rmat,[],2));
%         plot(bounds, bounds(j)-bounds(i), '-', 'Color', co{pp});
%     end
%     axis([0 0.8 -0.8 0.8])
% end


    

% 
% 
% 
% %% Plot post vs pre bound cp trials, color coded by optimal direction
% titles = {'LowHigh', 'HighLow', 'Low', 'High'};
% for tt = 1:length(titles)
% 
%     % Get data
%     xs = boundSummaries.(titles{tt})(:,1);
%     ys = boundSummaries.(titles{tt})(:,2);
%     Lg = isfinite(xs) & isfinite(ys);
% 
%     % Plot data
%     subplot(2,2,tt); cla reset; hold on;
%     plot([0 0.5], [0 0.5], 'k:')
%     plot(xs(Lg), ys(Lg), 'ko', 'MarkerFaceColor', 0.99.*ones(3,1))
% 
%     % For CP trials, find optimal direction and color code
%     if tt <= 2
%         cpt = boundSummaries.(['CP' titles{tt}]);
%         LoDirs = false(length(cpt),1);
%         for pp = 1:length(cpt)
%             if(Lg(pp))
%                 rrs = pmats{tt}(find(bounds(bounds<=0.5)<=xs(pp),1,'last'),:,find(cpt(pp)==cpTimes{tt},1));
%                 LoDirs(pp) = sign(ys(pp)-xs(pp)) == sign(bounds(find(rrs == max(rrs), 1))-xs(pp));
%             end
%         end
%         plot(xs(Lg&LoDirs), ys(Lg&LoDirs), 'ko', 'MarkerFaceColor', 'r');
%         p = binopdf(sum(LoDirs), sum(Lg), 0.5);
%     else
%         p = binopdf(sum(xs>ys), sum(Lg), 0.5);
%     end
% 
%     % Binomial test
%     title(sprintf('%s: p=%.3f', titles{tt}, p)); ...
% end
% subplot(2,2,3);
% xlabel('Bound pre CP')
% ylabel('Bound post CP')
% 
% 
% 
% 
% 
% 
% %% Junk below
% vals = { ...
%     'Low', 1, 'High', 1; ...
%     'Low', 2, 'High', 2; ...
%     'LowHigh', 1, 'HighLow', 1; ...
%     'LowHigh', 2, 'HighLow', 2; ...
%     };
% for xx = 1:4
%     subplot(3,2,xx+2); cla reset; hold on;
%     xs = boundSummaries.(vals{xx,1})(:,vals{xx,2});
%     ys = boundSummaries.(vals{xx,3})(:,vals{xx,4});  
%     Lg = isfinite(xs) & isfinite(ys);
%     plot([0 0.5], [0 0.5], 'k:')
%     plot(xs, ys, 'ko');
%     title(sprintf('%s (%d) vs %s (%d), p=%.3f', ...
%         vals{xx,1}, vals{xx,2}, vals{xx,3}, vals{xx,4}, ...
%         signrank(xs(Lg), ys(Lg))))
% end
% 
% 
% 
% 
% 
% 
% xs = boundSummaries.LowHigh(:,2);
% ys = boundSummaries.HighLow(:,2);
% cla reset; hold on;
% plot(xs, ys, 'ko')
% 
% %% Good plot of early vs late bound
% tags = {'LowHigh', 'Low', 'HighLow', 'High'};
% X = ones(60,2);
% for tt = 1:length(tags)
%     xs = boundSummaries.(tags{tt})(:,1);
%     ys = boundSummaries.(tags{tt})(:,2);
%     Lg = isfinite(xs) & isfinite(ys);
%     X(:,2) = xs;
%     [b,bint,r,rint,stats] = regress(ys(Lg), X(Lg,:));
%     disp([tt signrank(xs(Lg), ys(Lg))])
% 
%     subplot(4,1,tt); cla reset; hold on;
%     title(sprintf('%s: b=%.2f [%.2f %.2f]', tags{tt}, b(2), bint(2,1), bint(2,2)))
%     plot([0 0.5], [0 0.5], 'k:')
%     plot(X(Lg,2), ys(Lg), 'ko', 'MarkerFaceColor', 0.99.*ones(3,1))
%     X(:,2) = linspace(0,0.5,60);
%     plot(X(:,2), X*b, 'r-')
% end
% xlabel('Early bound')
% ylabel('Late bound')
% 
% 
% 
% % Compute peaks
% peaks = cell(2,1);
% for mm = 1:2
%     peaks{mm} = nan(size(pmats{mm},3),2);
%     for pp = 1:size(pmats{mm},3)
%         rmat = pmats{mm}(:,:,pp);
%         maxv = max(rmat(:));
%         [i,j] = find(rmat>=0.9.*maxv);
%         peaks{mm}(pp,:) = [mean(i), mean(j)];
%     end
% end
% 
% titles = {'LowHigh', 'HighLow'};
% 
% for tt = 1:length(titles)
%     xs = boundSummaries.(titles{tt})(:,1);
%     ys = boundSummaries.(titles{tt})(:,2);
%     Lg = isfinite(xs) & isfinite(ys);
%     X(:,2) = xs;
% 
%     subplot(2,1,tt); cla reset; hold on;
%     title(sprintf('%s: b=%.2f [%.2f %.2f]', titles{tt}, b(2), bint(2,1), bint(2,2)))
%     plot([0 0.5], [0 0.5], 'k:')
%     plot(X(Lg,2), ys(Lg), 'ko', 'MarkerFaceColor', 0.99.*ones(3,1))
% 
%     cpt = boundSummaries.(['CP' titles{tt}]);
%     for pp = 1:length(cpt)
%         rrs = pmats{tt}(find(bounds(bounds<=0.5)<=xs(pp),1,'last'),:,find(cpt(pp)==cpTimes{tt},1));
%         if max(rrs) > 0
%             br = bounds(rrs>=0.9.*max(rrs));
%         else
%             br = bounds(rrs>=1.1.*max(rrs));
%         end            
%         plot(repmat(xs(pp),1,length(br)), br, 'r.')        
%     end
% end
% 
% for pp = 1:16
%     for mm = 1:2
%         subplot(2,1,mm); cla reset; hold on;
%         rmat = pmats{mm}(:,:,pp);
%         rmat(rmat<0)=0;
%         imagesc(rmat);
%         maxv = max(rmat(:));
%         [i,j] = find(rmat>=0.9.*maxv);
%         plot(j,i,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 15);
%         [i,j] = find(rmat==maxv);
%         plot(j,i,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15)
%         axis([1 16 1 16])
%         title(sprintf('%s: %d', titles{mm}, cpTimes{mm}(pp)))
%         xlabel('Bound after CP')
%         ylabel('Bound before CP')
%     end
%     r = input('next')
% end
% 
% 
% xlabel('Early bound')
% ylabel('Late bound')
% 
% 
% 
% 
% 
% %%% junk below
% for mm = 1:2
%     rmat = pmats{mm};
%     rmat(rmat<0)=0;
%     for pp = 1:size(rmat,3)
%         cla reset; hold on;
%         imagesc(rmat(:,:,pp));
%         maxv = max(reshape(rmat(:,:,pp),[],1));       
%         [i,j] = find(rmat(:,:,pp)>=0.9.*maxv);
%         plot(j,i,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 15);
%         [i,j] = find(rmat(:,:,pp)==maxv);
%         plot(j,i,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15)
%         axis([1 16 1 16])
%         title(sprintf('%d', cpTimes{mm}(pp)))
%         xlabel('Bound after CP')
%         ylabel('Bound before CP')
%         r = input('next')
%     end
% end
% 
% 
% for pp = 1:size(RRMatrixLowHigh,3)
%     cla reset; hold on;
%     imagesc(RRMatrixLowHigh(:,:,pp))
%     [i,j] = find(RRMatrixLowHigh(:,:,pp)==max(reshape(RRMatrixLowHigh(:,:,pp),[],1)),1)
%     plot(j,i,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15)
%     axis([1 16 1 16])
%     title(sprintf('%d', cpTimes(pp)))
%     xlabel('Bound after CP')
%     ylabel('Bound before CP')
%     r = input('next')
% end
% 
% cla reset; hold on;
% for bb = 1:16
%     plot(RRMatrixLowHigh(bb,:,pp), 'k-')
%     r = input('next')
% end
% 
% 
% cpTimes = unique(boundSummaries.CPHighLow);
% for pp = 1:size(RRMatrixHighLow,3)
%     cla reset; hold on;
%     imagesc(RRMatrixHighLow(:,:,pp))
%     [i,j] = find(RRMatrixHighLow(:,:,pp)==max(reshape(RRMatrixHighLow(:,:,pp),[],1)),1)
%     plot(j,i,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15)
%     axis([1 16 1 16])
%     title(sprintf('%d', cpTimes(pp)))
%     xlabel('Bound after CP')
%     ylabel('Bound before CP')
%     r = input('next')
% end
