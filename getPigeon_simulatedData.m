function [choices, rts, DV, bounds, ndts, snrs, steps] = getPigeon_simulatedData(options)
%
% Returns:
%   choices ... 0/1
%   rts     ... values (really indices into DV)
%   DV      ... decision variable (step matrix, rows are trials, cols are
%                   steps)
%   bounds  ... array of bounds per trial
%   ndts    ... array of ndts per trial
%   snrs    ... matrix of mean/std per trial
%   steps   ... cell array of step arrays
%
% Check arguments
arguments
    options.generativeMean      double = 0.05   % generative mean
    options.generativeSTD       double = 0.15   % generative std
    options.numTrials           double = 1000000 % Simulated trials per SNR
    options.maxStepsPerTrial    double = 50     % Total possible samples per "trial"
    options.boundMean           double = 0      % bound mean/start value
    options.boundSTD            double = 0      % bound variability (trial-to-trial)
    options.boundSlope          double = 0      % Time-dependent slope
    options.boundMax            double = 0.75;  % max, fixed by task geometry
    options.NDTMin              double = 1;     % non-decision time
    options.NDTmax              double = 3;     % non-decision time
    options.lapseRate           double = 0;     % fraction "guess" trials
    options.bounds              double = [];
end

% Make the DV per trial
% Might be given multiple generativeMeans to interleave randomly
nMeans = length(options.generativeMean);
iMeans = randi(nMeans,options.numTrials,1);
options.generativeMean = options.generativeMean(:);
generativeMean = repmat(options.generativeMean(iMeans),1,options.maxStepsPerTrial);

% DV starts at zero
DV = cat(2, zeros(options.numTrials, 1), cumsum(normrnd(...
    generativeMean, options.generativeSTD),2));

% make array of bound means -- might be one per snr
if ~isempty(options.bounds)
    % if given an array of bounds, choose randomly
    boundMean = options.bounds(randi(length(options.bounds), options.numTrials, 1));
else
    % if given parameters, parse mean/std
    if length(options.boundMean) > 1 && length(options.boundMean) == nMeans
        boundMean = options.boundMean(iMeans);
        boundSTD = options.boundSTD(iMeans);
    else
        boundMean = repmat(options.boundMean(1), options.numTrials, 1);
        boundSTD = repmat(options.boundSTD(1), options.numTrials, 1);
    end

    % Check for variable bound
    LvariableBound = boundSTD > 0;
    if any(LvariableBound)
        % variable bound
        boundMean(LvariableBound) = max(0.03, normrnd( ...
            boundMean(LvariableBound), ...
            boundSTD(LvariableBound)));
    end
end
boundMatrix = repmat(boundMean(:), 1, options.maxStepsPerTrial+1);
if options.boundSlope~=0
    boundMatrix = boundMatrix .* repmat(cat(2,linspace(1,options.boundSlope,10), ...
        options.boundSlope.*ones(1,size(boundMatrix,2)-10)), options.numTrials, 1);
end

% Make array of non-decision times
if options.NDTmax > options.NDTMin
    % uniform random between min and max
    ndts = randi(options.NDTmax-options.NDTMin+1,options.numTrials,1)-options.NDTMin+1;
else
    % fixed value
    ndts = options.NDTMin*ones(options.numTrials,1);
end

% To save
choices = nans(options.numTrials, 1);
rts = repmat(options.maxStepsPerTrial, options.numTrials, 1);
bounds = nans(options.numTrials, 1);
snrs = generativeMean(:,1)./options.generativeSTD;
steps = cell(options.numTrials,1);

% disp(sprintf('looping through %d trials', options.numTrials))
% Loop through the trials to find first crossings

% Look for bound crossings
aDV = abs(DV);
Lcrossed = false(options.numTrials,1);
for ss = 1:options.maxStepsPerTrial
    % Check at each step
    Lrt = ~Lcrossed & aDV(:,ss)>=boundMatrix(:,ss);
    if any(Lrt)
        bounds(Lrt) = boundMatrix(Lrt,ss);
        Lcrossed = Lcrossed | Lrt;
        if ss == 1
            % Should only happen if bound=0; give NDT of 0 or 1
            rts(Lrt) = ss+randi(2,sum(Lrt),1)-1;
            choices(Lrt) = randi(2,sum(Lrt),1)-1;
        else
            % Otherwise add random ndt
            rts(Lrt) = ss+ndts(Lrt);
            choices(Lrt) = double(DV(Lrt,ss)>0);
        end
    end
end

% Save steps
rts = min(rts, options.maxStepsPerTrial);
for tt = 1:options.numTrials
    steps{tt} = DV(tt,1:rts(tt));
end

% Add lapses
if options.lapseRate > 0
    numLapseTrials = ceil(options.lapseRate*options.numTrials);
    choices(randi(options.numTrials, numLapseTrials)) = ...
        randi(2, numLapseTrials) - 1;
end
