function [bounds_, dts_, rts_, congruences_] = getPigeon_bounds(steps, choices, options)%, trueBound)
% function [bounds_, dts_, rts_, congruences_, boundMinMax_] = getPigeon_bounds(steps, choices, maxNDT, trueBound)
%
% Arguments:
%   steps ... cell array of steps
%   choices ... array of choices
%   maxNDT ... max non-decision time to consider
%
% Returns:
%   bounds_ ... array of bounds
%   dts_ ... array of decision times
%   rts_ ... array of rts
%   congruences_ ... array of fraction congruent choices (wrt pigeon position) for
%                       each NDT considered

%% Check args
arguments
    steps
    choices
    options.maxNDT = 4;
    options.minRT = 2; % min steps needed to compute bound
    options.SNR = [];
end

%% Make matrix of bounds, defined as the midpoint between sequential steps,
%   for different NDTs
numTrials = size(steps,1);
boundMatrix = zeros(numTrials,options.maxNDT); % assume default zero bound (guess)
rts_ = nan(numTrials,1);
for ii = 1:numTrials

    % Parse bound based on maxNDT last steps
    rts_(ii) = length(steps{ii});

    % Get final steps to parse bound
    finalSteps = fliplr(steps{ii}(max(rts_(ii)-options.maxNDT,1):end));

    % bounds are midpoints
    if length(finalSteps) > 1
        boundMatrix(ii,1:length(finalSteps)-1) = finalSteps(2:end)-diff(finalSteps)./2;
    end
end

% Check for congrudences by SNR and NDT
if isempty(options.SNR)
    options.SNR = zeros(numTrials,1); % dummy array
end
uniqueSNRs = unique(options.SNR);
numSNRs = length(uniqueSNRs);
choices = choices*2-1;
NDTs = nan(numTrials,1);
bounds_ = nan(numTrials,1);
congruences_ = nan(options.maxNDT,numSNRs);

%% Compute congruences
for ss = 1:numSNRs
    Ls = options.SNR == uniqueSNRs(ss);
    % Loop through NDTs
    for ii = 1:options.maxNDT
        Lc = Ls & boundMatrix(:,ii)~=0;
        if sum(Lc)>5
            % Congruence between choice and pigeon position at given NDT
            congruences_(ii,ss) = sum( ...
                sign(boundMatrix(Lc,ii)) == choices(Lc))./sum(Lc);
        end
    end
    % max congruence determines NDT
    % NDT of zero means that the congruence was highest between step end-1
    %   and end
    % NDT of one means that the congruence was highest between step end-2
    %   and end-1
    % etc
    ci = find(congruences_(:,ss)==max(congruences_(:,ss)),1);
    if ~isempty(ci)
        NDTs(Ls) = max(0, min(rts_(Ls)-2, ci-1));
    end
end

% Get bounds from NDTs
% Loop through the trials
for bi = find(isfinite(NDTs))'
    bounds_(bi) = boundMatrix(bi,NDTs(bi)+1);
end

% remove bad bounds
bounds_(rts_<options.minRT) = nan;

% Set decision time per trial
dts_ = max(0, rts_-NDTs-1);

%% Old version
% %% Check args
% arguments
%     steps
%     choices
%     options.maxNDT = 4;
%     options.minSteps = 2; % min steps needed to compute bound
%     options.SNR = [];
%     options.NDTBySNR = true;
%     options.boundBySNR = false;
%     options.generativeSTD = 0.1;
% end
% 
% %% Make matrix of bounds, defined as the midpoint between sequential steps,
% %   for different NDTs
% numTrials = size(steps,1);
% boundMatrix = zeros(numTrials,options.maxNDT+1,2); % assume default zero bound (guess)
% rts_ = nan(numTrials,1);
% for ii = 1:numTrials
% 
%     % Parse bound based on maxNDT last steps
%     rts_(ii) = length(steps{ii});
% 
%     % Get final steps to parse bound
%     finalSteps = fliplr(steps{ii}(max(rts_(ii)-options.maxNDT-1,1):end));
% 
%     % Save final step sizes/positions
%     if length(finalSteps) > 1
%         boundMatrix(ii,1:length(finalSteps)-1,1) = -diff(finalSteps);
%         boundMatrix(ii,1:length(finalSteps)-1,2) = finalSteps(1:end-1);
%         % boundMatrix(ii,1:length(finalSteps)) = finalSteps(1:end-1)+diff(finalSteps)./2;
%     end
% end
% 
% % Check for SNR
% if ~isempty(options.SNR) && options.NDTBySNR
%     SNR = abs(options.SNR);
%     uniqueSNRs = unique(options.SNR);
%     numSNRs = length(uniqueSNRs);
% else
%     SNR = zeros(numTrials,1);
%     uniqueSNRs = 0;
%     numSNRs = 1;
% end
% 
% % Determine non-decision time (NDT) from choice congruence
% % return arrays of bounds, congruences per trial
% choices = choices*2-1;
% NDTs = nan(numTrials,1);
% bounds_ = nan(numTrials,1);
% congruences_ = nan(options.maxNDT,numSNRs);
% 
% %% Compute congruences
% % Start by possibly looping through SNRs
% for ss = 1:numSNRs
%     Ls = SNR == uniqueSNRs(ss);
%     % Loop through NDTs
%     for ii = 1:options.maxNDT
%         Lc = Ls & boundMatrix(:,ii,2)~=0;
%         if any(Lc)
%             % Congruence between choice and pigeon position at given NDT
%             congruences_(ii,ss) = sum( ...
%                 sign(boundMatrix(Lc,ii,2)) == choices(Lc))./sum(Lc);
%             % congruences_(ii,ss) = sum( ...
%             %     sign(boundMatrix(Lc,ii,1)) == choices(Lc) & ...
%             %     sign(boundMatrix(Lc,ii,2)) == choices(Lc))./sum(Lc);
%         end
%     end
%     % max congruence determines NDT
%     ci = find(congruences_(:,ss)==max(congruences_(:,ss)),1);
%     NDTs(Ls) = ci;
% 
%     % Check if pigeon moved backwards; if so, add 1 step NDT
%     % Lskipped = sign(boundMatrix(:,ci,1)) ~= sign(boundMatrix(:,ci,2));
%     % NDTs(Ls & Lskipped) = NDTs(Ls & Lskipped) + 1;
% end
% 
% % check that ndts are not > rts
% NDTs(NDTs)
% 
% %% Compute bound
% if ~isempty(options.SNR) && options.boundBySNR
% 
%     % Compute separeately per SNR by taking expectation of jump
%     % distribution
%     mus = abs(options.SNR).*0.1;
% 
%     % Loop through the trials
%     for ii = 1:numTrials
% 
%         % Check if NDT is too large, if so, decrement
%         if length(steps{ii})<(NDTs(ii)+1) && length(steps{ii}) > 1
%             % disp(sprintf('NDT is %d', length(steps{ii})-1))
%             NDTs(ii) = length(steps{ii})-1;
%         end
% 
%         % first position + expected value of step between first and second
%         % positions (after accounting for ndt)
%         % first position = boundMatrix(ii,ndts(ii)+1,2)
%         % step size = boundMatrix(ii,ndts(ii),1)
%         if boundMatrix(ii,NDTs(ii),1) > 0 && boundMatrix(ii,NDTs(ii),2) > 0
%             bounds_(ii) = boundMatrix(ii,NDTs(ii)+1,2) + ...
%                 quadgk(@(x) x.*normpdf(x,mus(ii),0.1)./...
%                 ((normcdf(boundMatrix(ii,NDTs(ii),1),mus(ii),0.1) - ...
%                 normcdf(0,mus(ii),0.1))),0,boundMatrix(ii,NDTs(ii),1));
%         elseif boundMatrix(ii,NDTs(ii),1) < 0 && boundMatrix(ii,NDTs(ii),2) < 0
%             bounds_(ii) = boundMatrix(ii,NDTs(ii)+1,2) + ...
%                 quadgk(@(x) x.*normpdf(x,mus(ii),0.1)./...
%                 ((normcdf(0,mus(ii),0.1) - ...
%                 normcdf(boundMatrix(ii,NDTs(ii),1),mus(ii),0.1))), boundMatrix(ii,NDTs(ii),1),0);
%         else
%             bounds_(ii) = mean(boundMatrix(ii,NDTs(ii):NDTs(ii)+1,2));
%         end
%     end
% 
% else % not by SNR
% 
%     % Just take midpoints
%     % Loop through the trials
%     for ii = 1:numTrials
%         bounds_(ii) = mean(boundMatrix(ii,NDTs(ii):NDTs(ii)+1,2));
%     end
% end
% 
% % Set decision time per trial
% dts_ = max(0, rts_-NDTs-1);
% 
% % if nargin >= 4 && ~isempty(trueBound)
% %     cla reset; hold on;
% %     fracs = (trueBound-boundMatrix(:,2,2))./boundMatrix(:,1,1);
% %     %     histogram(fracs(choices==1),0:0.05:1)
% %     Lc = choices ==1;
% %     boxplot(fracs(Lc), rts_(Lc))
% %     plot([1 50], [0.5 0.5], 'k-')
% % end
% % probably shouldn't happen
% % bounds_(~isfinite(choices)) = nan;
% 
