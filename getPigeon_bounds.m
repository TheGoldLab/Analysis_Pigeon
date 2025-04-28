function [bounds_, dts_, rts_, congruences] = getPigeon_bounds(steps, choices, maxNDTs)%, trueBound)
% function [bounds_, dts_, rts_, congruences_, boundMinMax_] = getPigeon_bounds(steps, choices, maxNDTs, trueBound)
%
% Arguments:
%   steps ... cell array of steps
%   choices ... array of choices
%   maxNDTs ... how many non-decision times to consider
%
% Returns:
%   bounds_ ... array of bounds
%   dts_ ... array of decision times
%   rts_ ... array of rts
%   congruences_ ... fraction congruent choices (wrt pigeon position) for
%                       each NDT considered

% Check arg
if nargin < 3 || isempty(maxNDTs)
    maxNDTs = 5;
end

% Make matrix of bounds, defined as the midpoint between sequential steps,
%   for different NDTs
numTrials = size(steps,1);
boundMatrix = zeros(numTrials,maxNDTs,2); % assume default zero bound (guess)
rts_ = nan(numTrials,1);
for ii = 1:numTrials

    % Parse bound based on maxNDTs last steps
    rts_(ii) = length(steps{ii});

    % Get final steps to parse bound
    finalSteps = fliplr(steps{ii}(max(rts_(ii)-maxNDTs-1,1):end));

    % Save final step sizes/positions
    if length(finalSteps) > 1
        boundMatrix(ii,1:length(finalSteps)-1,1) = -diff(finalSteps);
        boundMatrix(ii,1:length(finalSteps)-1,2) = finalSteps(1:end-1);
        % boundMatrix(ii,1:length(finalSteps)) = finalSteps(1:end-1)+diff(finalSteps)./2;
    end
end

% Measure choice congruence
congruences = nan(maxNDTs,1);
choices = choices*2-1;
for ii = 1:maxNDTs-1
    Lc = boundMatrix(:,ii,2)~=0;
    if any(Lc)
        congruences(ii) = sum( ...
            sign(boundMatrix(Lc,ii,1)) == choices(Lc) & ...
            sign(boundMatrix(Lc,ii,2)) == choices(Lc))./sum(Lc);
    end
end
ndt = find(congruences==max(congruences),1)-1;
ndt = ndt*(ndt>0) + 1 * (ndt<=0);

% Set bound per trial
ndt_idx = ndt;
bounds_ = boundMatrix(:,ndt_idx,2) - 0.5.*(boundMatrix(:,ndt_idx,1));
Lbump = sign(boundMatrix(:,ndt_idx,1)) ~= sign(boundMatrix(:,ndt_idx,2)); 
bounds_(Lbump) = boundMatrix(Lbump,ndt_idx+1,2) - 0.5.*(boundMatrix(Lbump,ndt_idx+1,1));
bounds_(~isfinite(choices)) = nan;

% if nargin >= 4 && ~isempty(trueBound)
%     cla reset; hold on;
%     fracs = (trueBound-boundMatrix(:,2,2))./boundMatrix(:,1,1);
%     %     histogram(fracs(choices==1),0:0.05:1)
%     Lc = choices ==1;
%     boxplot(fracs(Lc), rts_(Lc))
%     plot([1 50], [0.5 0.5], 'k-')
% end

% Set decision time per trial
dts_ = max(0, rts_-ndt);
