function X = fcm_cluster(pulse,k,datafield,max_nan,varargin)
%FCM_CLUSTER Uses fuzzy c-means to cluster a given datafield in
% the fit_array. In order to standardize the cluster naming
% schematic, user needs to input an order vector.
%
% USAGE: A = pulse.fcm_cluster(5,'corrected_area_norm')
%        A = pulse.fcm_cluster(5,'corrected_area_norm',3)
%
% INPUT: pulse - array of Pulse to be clustered (passed by reference)
%        k - the number of seeding clusters
%        datafield - the data you want to cluster. Procedure
%           will normalize by taking the Z score within each
%           pulse
%        max_nan - maximum number of tolerated NaN
%
% OUTPUT: A - area matrix used for clustering
% xies@mit.edu

fits = [pulse.fits];

% clear previous labels & weights
[fits.cluster_label] = deal([]);
[fits.cluster_weight] = deal([]);

filtered = fits(...
    cellfun(@(x) numel(x(isnan(x))),{fits.(datafield)}) < max_nan );

X = cat(1,filtered.(datafield));
if nargin > 4
    X = cat(X,varargin{1});
end
X = bsxfun(@minus,X,nanmean(X));
X = bsxfun(@rdivide,X,nanstd(X,[],2));

X(isnan(X)) = 0;

[cluster_centroid,U] = fcm(X,k,[2 1e3 1e-5 1]);
[max_prob, labels] = max(U);

for i = 1:k
    subplot(k,1,i)
    plot (nanmean( X(labels==i,:)) );
end
display('Enter the label order: 1-Stereotyped, 2-Early, 3-Delayed, 4-Unratcheted, 5-Stretched')
order = input(':');
revorder = reverse_index(order);

labels = revorder(labels);
U = U(revorder,:);

% store labels

% fits = set_field(fits,[filtered.fitID], 'cluster_label', labels);
% fits = set_field(fits,[filtered.fitID], 'cluster_weight', U);
for i = 1:numel(labels)
    filtered(i).cluster_label = labels(i);
    filtered(i).cluster_weight = max_prob(i);
%     [fits([fits.fitID] == filtered(i).fitID).cluster_label] = ...
%         deal( labels(i) );
%     [fits([fits.fitID] == filtered(i).fitID).cluster_weight] = ...
%         deal( max_prob(i) );
end

% deal with non-clustered fits (label = 6, weight = NaN)
[fits(cellfun(@isempty, {fits.cluster_label} )).cluster_label] = ...
    deal(k+1);
[fits(cellfun(@isempty, {fits.cluster_weight} )).cluster_weight] = ...
    deal( NaN );

%             for i = 1:k
%                 [fits([fits.cluster_label] == i).cluster_label] = deal(revorder(i)*10);
%             end
%             for i = 10:10:k*10
%                 [fits([fits.cluster_label] == i).cluster_label] = deal(i/10);
%             end

end % cluster