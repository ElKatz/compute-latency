%% main latency computation script 
% Simple example script for how to use compute_latency.m to compute visual
% onset latencies of neurons and plot them. 

%% please input details of dataset:
dataCollectedBy = 'Katz';
brainAreaName   = 'SC';
taskName        = 'FNFLGN';
datasetString   = [dataCollectedBy, ' ', brainAreaName, ' during ' taskName];

%% please enter the lag between your event time and actual monitor refresh:
% in Krauzlis section, using Viewpixx monitors and Plexon omniplex, we have
% a 10m lag between when we strobe a flip and when the monitor actually
% presentes it. 
lag = 0.01; % in secdonds. 


%% compute latencies for each neuron in a dataset 
% latency can be computed in various methods. In this script we compute it
% using several. 
methodList  = {'white2017', 'yuKatz2023'}; % see compute_latency.m for more details about method
nMethods    = numel(methodList);
vizBool     = true; % set this to true if you wish to visuzlize each neuron 

% - I organize my data in a struct array called 'R' (of length nNeurons). 
% In each element of R I spike times (in 'R.times'). 
% - I only use neurons that have a clear visual RF, and that have a visual
% response to stimuli in their RF.
nNeurons    = numel(R);

% loop over latency computation methods to compute latency:
latencyTime = cell(nMethods,1);
for iM = 1:nMethods
    latencyTime{iM} = nan(nNeurons,1);
    hW          = waitbar(0/nNeurons, 'Computing latency for each and every neuron..');
    disp('Latency computation has commenced')
    for ii = 1:nNeurons
        % for each neuron, get all event times (eg times at which stimulus
        % appeared):
        evTimes = XXXX;
        % compute latency:
        latencyTime{iM}(ii) = compute_latency(R(ii).times, evTimes, methodList{iM}, vizBool, [datasetString, ', Neuron #' num2str(ii) ' - ' methodList{iM}]);
        % update waitbar:
        waitbar(ii/nNeurons, hW);
    end
    close(hW)
    disp('Done!')
    
    % correct for lag:
    latencyTime{iM} = latencyTime{iM} - lag;
end

%% plot:

hFig = figure; hold on;
for iM = 1:nMethods
    subplot(nMethods, 1, iM); hold on
    h = cdfplot(latencyTime{iM});
    title(['method - ' methodList{iM}])
    set(h, 'Color','k', 'LineWidth', 2);
    xlabel('Time from stimulus onset')
    ylabel('Cumulative proportion of neurons')
    xlim([0 0.1])
    ylim([0 1])
    set(gca, 'XTick', 0:.025:.1, 'YTick', 0:.25:1);
    % median:
    md = nanmedian(latencyTime{iM});
    line([0 md], [.5 .5], 'Color', 'r', 'LineStyle', '--')
    line([md md], [0 .5], 'Color', 'r', 'LineStyle', '-')
    text(md+0.001, 0.4, num2str(md,2), 'fontSize', 10, 'Color', 'r')
    % n:
    nDetected = sum(~isnan(latencyTime{iM}));
    text(md+0.001, 0.3, ['Latency detected for ' num2str(nDetected) ' of ' num2str(nNeurons) ' neurons'], 'fontSize', 10, 'Color', 'k');
end
annotation('textbox', [.02 .98 .98 .02], 'string', datasetString, 'linestyle', 'none', 'fontsize', 10, 'Interpreter', 'none');

% save:
dirFigs = '~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Projects/2022 - scFaces/Figures';
figSz = [8, 4*nMethods];
set(hFig, 'PaperSize', figSz, 'PaperPosition', [0 0 figSz]);
saveas(hFig, fullfile(dirFigs, ['visualOnsetLatency_' datasetString '.pdf']))



