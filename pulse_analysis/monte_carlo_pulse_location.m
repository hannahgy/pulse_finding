function MC = monte_carlo_pulse_location(fits,cells,opt)
%MONTE_CARLO_PULSE_LOCATION Generate spatially random pulse patterns by
% permutation analysis
% 
% USAGE: MC = monte_carlo_pulse_location(fits,cells,cx,cy,neighborID,opt)
%
% INPUT: fits - pulses
% 		 cells - cells
%		 neighborID - neighborhood of each cell
%		 opt - options for MC analysis
%		 	.Nboot - number of iterations
%			.time_windows - windows in time to check
%			.neighbor_def - neighborhood defintion
%			.savepath
%
% OUTPUT: MC
%
% xies@mit.edu Oct 2013

% parse input options structure
Nboot = opt.Nboot;
time_windows = opt.timewindows;
neighbor_def = opt.neighbor_def;

% get cell info
% cx = cat(2,cells.centroid_x); cy = cat(2,cells.centroid_y);

% preallocate variables
random_cell(Nboot).num_near = [];
random_cell(Nboot).origin_labels = [];
random_cell(Nboot).target_labels = [];
random_cell(Nboot).correlation = [];
random_cell(Nboot).centers = [];
random_cell(Nboot).fitID = [];

for j = 1:Nboot
    
    tic
    % make permutations of cell location
    [fits_bs_cell,cells_bs_cell] = cells.monte_carlo_stackID(fits,opt.monte_carlo);
    % get nearby pulses
	fits_bs_cell = fits_bs_cell.find_near_fits(cells_bs_cell,neighbor_def);
    
    % get current cluster for empirical (only once)
    if j == 1
        
        if strcmpi(opt.filter,'on')
            fitsOI = fits.find_non_edge(cells);
        end
        
        this_nearIDs = cat(1,fitsOI.nearIDs);
        % calculate num-neighbors for each pulse
        num_near = cellfun(@(x) numel(x(~isnan(x))), this_nearIDs);
        angles = cat(1,fitsOI.near_angles);
        % calculate number of neighboring cells to each pulse
%         num_cells = fits_bs_cell.
        % origin labels
        labels = [fitsOI.cluster_label]';
        % get all target labels
        target_labels = ...
            cellfun( @(x) [fitsOI.get_fitID(x).cluster_label], ...
            this_nearIDs,'UniformOutput',0);
        
        % get centers
        centers = [fitsOI.center];
        CT = [ fitsOI(~isnan([fitsOI.nearest_neighbor])).center ];
        nn = fitsOI.get_fitID([fitsOI.nearest_neighbor]);
        dt = [nn.center] - CT;
        
        empirical.num_near = num_near;
        empirical.origin_labels = labels;
        empirical.target_labels = target_labels;
        empirical.centers = centers;
		empirical.neighbor_windows = dt;
        empirical.near_angle = angles;
%         empirical.correlation = spatial_correlation(cx,cy,fits,30);
        
    end
    
    % random-cell
    if ~isempty(fits_bs_cell)
        
        if strcmpi(opt.filter,'on')
            fits_bs_cell = fits_bs_cell.find_non_edge(cells_bs_cell);
        end
        
        nearIDs_cell = cat(1,fits_bs_cell.nearIDs);
        angles = cat(1,fits_bs_cell.near_angles);
        % tabulate num-neighbors for each pulse
        num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs_cell);
        % get origin labels
        labels = [fits_bs_cell.cluster_label]';
        % break down all target labels
        target_labels = ...
            cellfun(@(x) [fits_bs_cell.get_fitID(x).cluster_label], ...
            nearIDs_cell,'UniformOutput',0);
        centers = [fits_bs_cell.center];
       
		CT = [fits_bs_cell(~isnan([fits_bs_cell.nearest_neighbor])).center];
        nn = fits_bs_cell.get_fitID([fits_bs_cell.nearest_neighbor]);
        
        dt = [nn.center] - CT;
        
        random_cell(j).num_near = num_near;
        random_cell(j).origin_labels = labels;
        random_cell(j).target_labels = target_labels;
        random_cell(j).centers = centers;
        random_cell(j).fitID = [fits_bs_cell.fitID];
		random_cell(j).neighbor_windows = dt;
        random_cell(j).near_angle = angles;
        
    end
    % Compute correlation function
%     random_cell(j).correlation = spatial_correlation(cx,cy,fits_bs_cell,30);
    
    T = toc;
    display(['Done with ' num2str(j) ' in ' num2str(T) ' seconds.']);
    
end

MC.empirical = empirical;
MC.random_cell = random_cell;

MC.neighbor_def = neighbor_def;
MC.time_windows = time_windows;

if ~isempty(opt.savepath)
	save(opt.savepath,'MC');
    display(['Saved to: ' opt.savepath])
end

end
