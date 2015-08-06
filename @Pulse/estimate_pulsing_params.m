function [freq,neighbor_count] = estimate_pulsing_params(pulse)
%ESTIMATE_PULSING_PARAMETERS Returns frequency and neighbor cell count
% distributions for a given set of Pulses.
%
% USAGE: [freq,nCount] = pulses.estimate_pulsing_params;
% 
% INPUT: pulse - array of pulses (of the same genotype!)
% 
% OUTPUT: freq - functional handle of gamma distribution with parameters
%            estimated from all Pulses
%         neighbor_count - histogram of neighbor counts
%
% xies@mit.edu

cells = [pulse.cells]; fits = [pulse.fits];

% Check inputs
if numel([fits.cluster_label]) ~= numel(fits)
    error('All pulses must have a cluster label');
end

% Estimate frequency via gamma distribution fit
f = pulse.get_frequency;
[phat,pci] = gamfit([f{:}]);
freq.fun = @(x) gamcdf(x,phat(1),phat(2));

% Estimate # of neighboring cells
num_neighbors = zeros(1,numel(fits));
index = 0;

for i = 1:numel(pulse)
    
    f = pulse(i).fits;
    N = pulse(i).cells.get_adjacency_matrix;
    
    for j = 1:numel(f)
        
        index = index + 1;
        this_fit = f(j);
        this_conn = N(this_fit.cellID,:,this_fit.center_frame);
        
        num_neighbors(index) = numel(this_conn(this_conn > 0));
        
    end
    
end

neighbor_count = cell(1,max([fits.cluster_label]));
for i = 1:max([fits.cluster_label])
    neighbor_count{i} = ...
        hist(num_neighbors([fits.cluster_label] == i),0:10) ...
        /numel(find([fits.cluster_label]==i));
end

end