function zscores = plot_mc_results(MC,window,temporal_bin,opt)
%PLOT_MC_RESULTS Generates visualization of how resampled data compares
% with empirical counts
% xies@mit Oct 2013

empirical = MC.empirical;
random_cell = MC.random_cell;

% temporal bin bounds
nbins = size(temporal_bin,2);
left = temporal_bin(1,:);
right = temporal_bin(2,:);

num_clusters = max(empirical.origin_labels) - 1;

zscores_cell = zeros(nbins,num_clusters);


for K = 1:nbins
    
    if strcmpi(opt.normalize,'on')
        num_emp = bsxfun(@rdivide,empirical.num_near,empirical.near_angle);
    else
        num_emp = cat(1,empirical.num_near);
    end
    num_cell = {random_cell.num_near};
    near_cell = {random_cell.near_angle};
    
    labels_emp = empirical.origin_labels;
    labels_cell = {random_cell.origin_labels};
    
    target_emp = empirical.target_labels;
    target_cell = {random_cell.target_labels};
    
    %filter by temporal bin
%     filter = @(x) (x.centers > left(K) & x.centers <= right(K));
    
%     num_emp = num_emp( filter(empirical),: );
%     labels_emp = labels_emp( filter(empirical),: );
%     target_emp = target_emp( filter(empirical),:,: );
%     
%     num_cell = num_cell( filter(random_cell(1)),:,: );
%     labels_cell = labels_cell( filter(random_cell(1)) );
%     target_cell = target_cell( filter(random_cell(1)),:,: );
    
    figure(2)
    
    for i = 1:num_clusters
        
        % Distribution of means within a behavior
%         this_count_cell = squeeze( num_cell( labels_cell == i,window,:) );
        if strcmpi(opt.normalize,'on')
            this_count_cell = ...
                cellfun(@(x,y,z) (x(y==i,window)./z(y==i)), ...
                num_cell,labels_cell,near_cell,'UniformOutput',0);
        else
            this_count_cell = ...
                cellfun(@(x,y) (x(y==i,window)), ...
                num_cell,labels_cell,'UniformOutput',0);
        end
        this_count_emp = num_emp( labels_emp == i,window );
        
%         this_count_emp = this_count_emp ./ num_emp;
%         this_count_cell = this_count_cell ./ near_cell( labels_cell== i,: );
        
        mean_of_cell = cellfun(@nanmean, this_count_cell);
        mean_of_emp = nanmean( this_count_emp );
        
        [Nmean_cell,bins] = hist(mean_of_cell,25);
        
        H(i) = subplot(num_clusters,1,i);
        
        h = bar(bins, ...
            Nmean_cell/sum(Nmean_cell), ...
            'LineStyle','None');
        
        set(h(1),'FaceColor','red');
        vline(mean_of_emp,'b');
        xlim(opt.xlim);
        
        if i == 1
            xlabel('Average number of neighbors');
            ylabel('Frequency');
%             legend('Random-cell','Empirical');
        end
        
        zscores_cell(K,i) = ...
            ( mean_of_emp - mean(mean_of_cell) ) / nanstd(mean_of_cell);
        
        if strcmpi(opt.breakdown,'on')
            
            for j = 1:6
                
                num_target_emp = cellfun(@(x) numel(x(x == j)), ...
                    squeeze( target_emp( labels_emp == i,window,:) ) );
                num_target_cell = cellfun(@(x) numel(x(x==j)), ...
                    squeeze( target_cell( labels_cell == i,window,:) ) );
                
                mean_of_emp = mean(num_target_emp);
                mean_of_cell = mean(num_target_cell);
                
                z_target_cell(j) = ...
                    (mean_of_emp - mean(mean_of_cell)) / std(mean_of_cell);
                
                figure(3)
                subplot(1,num_clusters,(K-1)*num_clusters + i);
                h = bar(1:6, z_target_pulse );
                set(h(1),'FaceColor','red')
                xlim([0 7])
                ylim([-3 4])
                
            end
            
        end
        
    end
    
    figure(1)
    subplot(1,5,5);
    h = bar(1:num_clusters, zscores_cell(K,:) ,'LineStyle','None');
    set(h(1),'FaceColor','red');
    title([ num2str(left(K)) ' < center <= ' num2str(right(K)) ]);
    ylim([-1 5]);
    
%     linkaxes(H,'x');
    
    
end

zscores = zscores_cell;

end
