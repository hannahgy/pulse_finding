%% Stability analysis of FCM clustering

Ks2try = 2:10;

[RI,randRI] = pulse.fcm_stabiliy(Ks2try);

%%

num_clusters = 3;

pulse.fcm_cluster(num_clusters,'corrected_area_norm',3);

%%

clear cluster*

% for i = 1:num_clusters
%     
%     eval(['cluster' num2str(i) ' = fits([fits.cluster_label] == ' num2str(i) ');']);
%     
% end
% 
% behaviors = {'Ratcheted',...
%     'Un-ratcheted', ...
%     'Unconstricting'};
% 
% colors = {'b','m','r'};
