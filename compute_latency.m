function latencyTime = compute_latency(spTimes, evTimes, method, vizBool, figStr)
%   latencyTime = compute_latency(spTimes, evTimes, method)
%
% Computes the latency at which activity in time-varying bins differs from
% a baseline bin. Written to operate on one neuron at a time. 
%
% INPUT:
%   spTimes - vector of spike times for a given neuron 
%   * These do not need to be segmented into trials. They are simply all 
%   spike times of a neuron within a session)
%
%   evTimes - vector of event times for which you want to compute latency
%   * These also do not need to be segmented into trials. It might be the 
%   case that each event happens once in a trial (e.g. target time in a
%   guided saccade task) but it might also be the case that you get
%   mulitple events within a trial (e.g. stimulus onset time in a serial
%   presentation of stimuli).
%
%   method - string. Needs to correspond to one of the cases of the switch
%       loop below (eg. 'white2017' / 'yuKatz2023' / ...)
%
%   vizBool - (optional) boolean to visualize (ie plot a figure)
%
%   figStr - (optional) string that appears on the figure
%
% OUTPUT:
%   latencyTime - in s.

% 20231124 lnk wrote it. Contact leor.katz@nih.gov for praise/repraoch


switch method
    case 'white2017'
        % parameters by: White, Kan, Levy, Itti, and Munoz, 2017, PNAS.
        p.binSize       = 0.005; % in s
        p.binIncrement  = 0.001; % in s
        p.winLatency    = [0 0.1]; % temporal window within which to detect latency relative to event time, in s
        p.winBaseline   = [-0.1 0]; % temporal window within which to compute baseline activity relative to event time, in s
        p.alphaVal      = 0.05; % for statistical significance testing
        p.kBins         = 1; % number of consecutive bins that need to cross alphaVal in order to decalre "onset"

    case 'yuKatz2023'
        % parameters by: Yu, Katz, Quaia, Messinger, Krauzlis, 2023, bioRxiv.
        p.binSize       = 0.02; % in s
        p.binIncrement  = 0.001; % in s
        p.winLatency    = [0 0.1]; % temporal window within which to detect latency relative to event time, in s
        p.winBaseline   = [-0.1 0]; % temporal window within which to compute baseline activity relative to event time, in s
        p.alphaVal      = 0.05; % for statistical significance testing
        p.kBins         = 20; % number of consecutive bins that need to cross alphaVal in order to decalre "onset"
end


%% computation:
tic;
% computation method may differ across methods. If you add a method that
% requires a different computation, please add it as its onw case
switch method
    case {'white2017', 'yuKatz2023'}
        % get spike rate for baseline window:
        spRateBl = binSpTimes(spTimes, evTimes, p.winBaseline, range(p.winBaseline)) ./ range(p.winBaseline); % convert to rate
        % define bins for within which to obtain spike rate:
        binEdges    = 0:p.binIncrement:(range(p.winLatency) - p.binSize);
        binCenters  = binEdges(1:end) + p.binSize ./ 2;
        nBins       = numel(binCenters);
        % init:
        spRateLat   = cell(nBins,1);
        pVal        = nan(nBins,1);
        isGreater   = nan(nBins,1);
        % get spike count for latency windows (in every bin, for ever increment):
        for iB = 1:nBins
            tStart  = p.winLatency(1) + binEdges(iB);
            tEnd    = tStart + p.binSize;
            spRateLat{iB} = binSpTimes(spTimes, evTimes, [tStart tEnd], range([tStart tEnd])) ./ range([tStart tEnd]);
            % compute pVal vs. baseline:
            [~, pVal(iB)] = ttest(spRateBl, spRateLat{iB}); 
            %      pVal(ik) = ranksum(spRateBl, spRateLat{ik}); % rank sum is giving me trouble... lots of false 0's... 
            % is binned activity greater than baseline?
            isGreater(iB) = mean(spRateLat{iB}) > mean(spRateBl);
        end         

        % only use pVals for which activity is *greater* than baseline:
        pVal(isGreater==0) = nan;

        % threshold pValues to detect which are smaller than alphaVal:
        thresh = pVal < p.alphaVal;
        % find first element for which thresholded pVal crosses below 
        % alphaVal for k bins:
        a = thresh;
        b = ones(p.kBins,1);
        latencyPtr  = find(conv(a,b) == sum(b),1) - (sum(b)-1);
        if ~isempty(latencyPtr)
            latencyTime = binCenters(latencyPtr);
        else
            latencyTime  = nan;
        end
end
disp(['Latency detected in ' num2str(round(toc,1)) ' seconds']);


%% viz?
% set vizBool to true if you wish to visualize the latency detection with
% respect to mean firing of the neuron. Helpful for sanity-checking and
% refining a detection method.
if ~exist('vizBool', 'var')
    vizBool = false;
end
if ~exist('figStr', 'var')
    figStr = [];
end
if vizBool
    figure;
    annotation('textbox', [.02 .98 .98 .02], 'string', figStr, 'linestyle', 'none', 'fontsize', 8, 'Interpreter', 'none');
    %% plot the mean firing rate relative to baseline:
    subplot(211); hold on; title('Mean firing rate');
    xlim(p.winLatency)
    xlabel('Time')
    ylabel('Mean response')
    % mean firing rate:
    m = cellfun(@(x) mean(x), spRateLat);
    plot(binCenters, m, 'lineWidth', 2);
    % baseline mean +/- 1 SEM:
    [mBl, ciBl] = bootMeanAndCI(spRateBl, 1e3, 2.5, 97.5);
    line(xlim, [mBl mBl], 'Color', 'k', 'LineStyle', '-');
    line(xlim, [ciBl(1) ciBl(1)], 'Color', 'k', 'LineStyle', '--');
    line(xlim, [ciBl(2) ciBl(2)], 'Color', 'k', 'LineStyle', '--');
    % onset latency:
    if ~isempty(latencyPtr)
        plot(binCenters(latencyPtr), m(latencyPtr), 'o', 'Color', 'r', 'LineWidth', 2);
        line([binCenters(latencyPtr) binCenters(latencyPtr)], ylim, 'Color', 'r');
    end
    
    %% plot pValue over time:
    subplot(212); hold on; title('pValue relative to criterion');
    xlabel('Time')
    ylabel('pVal')
    % plot pValue:
    plot(binCenters, pVal, 'lineWidth', 2);
    xlim(p.winLatency)
    % criterion:
    line(xlim, [p.alphaVal, p.alphaVal], 'Color', 'k', 'LineStyle', '--');
    % onset latency:
    if ~isempty(latencyPtr)
        plot(binCenters(latencyPtr), pVal(latencyPtr), 'o', 'Color', 'r', 'LineWidth', 2);
        line([binCenters(latencyPtr) binCenters(latencyPtr)], ylim, 'Color', 'r');
    end
end

end

%%
function [spcnt, bcenters,validEvents] = binSpTimes(sptimes, ev, win, binSize)
% bin spike times around event
% [spcnt, bcenters,validEvents]=binSpTimes(sptimes, ev, win, binSize)


if nargin < 4
    dsptimes=diff(sptimes);
    if nargin < 3
        win = [-mean(dsptimes) 10*mean(dsptimes)];
        if nargin < 2
            help binSpCounts
            return
        end
    end
    binSize = min(dsptimes);
end
sptimes=sptimes(:);

be = win(1):binSize:win(2);
bcenters = be(1:end-1)+binSize/2;

ev=ev(:);
binfun = @(t) (t == 0) + ceil(t/binSize);

validEvents = ~isnan(ev);
nTrials=numel(ev);

nEvents = sum(validEvents);

sbn = [];
str = [];
% assert(nEvents<2e3, 'too many events for this to run fast!')
for kEvent=1:nTrials
    if ~validEvents(kEvent)
        continue
    end
    spo = sptimes(sptimes > ev(kEvent) + be(1) & sptimes < ev(kEvent) + be(end))- ev(kEvent);
    if ~isempty(spo)
    sbn = [sbn; binfun(spo- be(1))]; %#ok<AGROW>
    str = [str; ones(numel(spo),1)*kEvent]; %#ok<AGROW>
    end
end

spcnt = full(sparse(str, sbn, 1, nTrials, binfun(be(end)-be(1))));
spcnt(~validEvents,:)=nan;

end

%%
function [bootMean, bootCI, bootDstMeans] = bootMeanAndCI(dst, nBoots, lPrc, uPrc)
%   [bootMean, bootCI] = bootMeanAndCI(dst, nBoots, lPrc, uPrc)
%
% bootstrasp a given dstribution to compute a mean and confidence interval
% of given prcentil (default is inner 68.2%, 1 SEM)
%
% Inputs:
%   dst     - the dstribution of numbers.
%   [nBoots] - number of boots (default is 1000)
%   [lprc]   - lower percentile (default is 15.9)
%   [uprc]   - upper percentile (default is 84.1)
% Outputs:
%   bootMean        - mean of means of bootstrapped dstributions
%   bootCI          - [2,1] lower and upper values for confidence interval
%   bootDstMeans   - [nBoots,1] the means for all bootstrapped dsts.
%
% see also bootMeanAndErr

% 201404 lnk wrote it

%%
if nargin < 3
    lPrc = 15.9;
    uPrc = 84.1;
    if nargin < 2
        nBoots = 1000;
    end
end

if ~isempty(dst)
    nValues = numel(dst);
    
    
    bootDst         = RandSample(dst, [nValues, nBoots]);
    bootDstMeans    = nanmean(bootDst,1);
    
    bootMean        = nanmean(bootDstMeans);
    bootCI          = prctile(bootDstMeans, [lPrc, uPrc]);
else
    bootMean        = nan;
    bootCI          = [nan,nan];
end

end


