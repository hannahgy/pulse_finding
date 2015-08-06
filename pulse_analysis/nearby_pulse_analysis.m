%% Nearby pulse analysis

pulseOI = pulse(1:5);
name = 'wt';

clear neighbor_definition
time_windows = 60; % seconds
neighbor_definition.temporal.def = ...
    @(time_diff,tau) abs(time_diff) < tau & time_diff > 0;

neighbor_definition.temporal.windows = time_windows;
neighbor_definition.spatial.def = 'identity';
% neighbor_definition.spatial.threshold = 8;

pulseOI.find_near_fits(neighbor_definition);
nearIDs = cat(1,pulseOI.getFits.nearIDs);

% Convert to number of pulses
num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);

%%
% MC stackID

entries = {'Ratcheted (stereotyped)','Ratcheted (weak)','Ratcheted (delayed)','Un-ratcheted','Stretched'};

clear o
o.Nboot = 50;
o.timewindows = time_windows;
o.neighbor_def = neighbor_definition;
o.monte_carlo = 'simulation';
o.filter = 'on';
% o.savepath = ['~/Desktop/simulated pulses/mc_stackID_' name, '_', ...
%     'iter_', num2str(n), '_' ...
%     neighb_str, '_', neighbor_defition.spatial.def, ...
%     '_Nboot', num2str(o.Nboot), '_', o.monte_carlo, '_neighborfilt_', o.filter ...
%     , '_k' num2str(num_clusters)];

MC = monte_carlo_pulse_location(pulse(1:5), o);

% end

%% Select correct timing

% select dataset
% MC = MC_wt;

% MC = filter_mc(MC_control, [11 12 13 14 15] );

tau = 1;    % neighborhood time window
clear opt temporal_bins
temporal_bins(1,:) = [-Inf];
temporal_bins(2,:) = [Inf];

opt.normalize = 'off';
opt.breakdown = 'off';
opt.xlim = [1 4];
% opt.normalize = [5.06 5.00 5.29 5.01];

zscores_wt = plot_mc_results(MC,tau,temporal_bins,opt);

%% Visualize raw distributions

bins = 0:20;
idx = MC.random_cell(1).origin_labels';
RC_num_near = cat(3,MC.random_cell.num_near);

figure,

for i = 1:num_clusters
    
    subplot(1,num_clusters,i);
    
    Nrc = hist(flat(RC_num_near(idx == i,window,:)),bins);
    Nemp = hist(MC.empirical.num_near(idx == i,window),bins);
    plot(bins,cat(1,(Nrc)/sum(Nrc),(Nemp)/sum(Nemp))')
    xlim([-1 10])
    title(behaviors{i})
    
end

xlabel('Number of neighbors')
ylabel('Probability')

Nemp = hist(MC.empirical.num_near(:,window),bins);
Nrc = hist(RC_num_near(:,window),bins);
figure,
plot(bins,cat(1,Nrc/sum(Nrc),Nemp/sum(Nemp))');
xlim([-1 15])
