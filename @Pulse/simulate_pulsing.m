function pulse_sim = simulate_pulsing(pulse,freqHat)
% Simulates spatially random pulses onto the empirical cell lattice
% using existing FITS as seeds and freqHat to estimate the
% frequency between consecutive pulses within a cell and pcHat
% to estimate the number of pulses within a cell

% set up local stream with /dev/random (UNIX only)
% set up local stream with /dev/random to get random seed
sfd = fopen('/dev/urandom');
seed = fread(sfd, 1, 'uint32');
fclose(sfd);

% check MATLAB version
str = version;
if strcmpi(str(1),'8')
    rng(seed)
else
    % MATLAB's start-up settings for random seeds
    stream = RandStream('mt19937ar','Seed',seed);
    RandStream.setDefaultStream(stream);
end

pulse_sim = pulse.clear;

% Repeat for each embryo
for e = 1:numel(pulse)
    
    % Sort pulses by their center of timing
    fits2sim = pulse_sim(e).fits;
    [~,I] = sort([fits2sim.center],'ascend');
    fits2sim = fits2sim(I);
    % Get all cells in this embryo
    cellsOI = pulse_sim(e).cells;
    Ncells = numel(cellsOI);
    
    % Generate the adjacency matrix as a fcn of time
    A = pulse(e).cells.get_adjacency_matrix;
    
    % keep track with Nframe x Ncell matrix of which cell already pulsed
    % when
    already_pulsed = zeros( max([fits2sim.center_frame]),Ncells );
    
    % Loop over all fits
    for i = 1:numel(fits2sim)
        
        accept = 0;
        this_fit = fits2sim(i); % Should be modifying by reference
        this_fit_original = pulse(e).get_fitID(this_fit.fitID);
        % Find the empirical/original cell this_fit comes from
        cellOI = pulse_sim(e).get_cellID(this_fit_original.cellID);
        frame = this_fit.center_frame;
        % Figure out how many adjacent cells current pulsing cell has
        this_fit.neighbor_cells = sum(A(this_fit_original.cellID,:,frame));
        
%         display(['Simulating fitID: ' num2str(this_fit.fitID)])
        
        % TODO: Corner case NaN is center_frame - need to deal with
        % case ... right now just spits out same cell
        if this_fit.neighbor_cells == 0
%             display('NaN')
            accept_move(this_fit,cellOI);
            already_pulsed(i,cellOI.cellID) = 1;
            continue
        end
        
        % Find the number of neighboing cells to the pulse
        N = cat(2,cellsOI.identity_of_neighbors_all);
        num_neighbors = N(frame,:);
        num_neighbors = cellfun(@(x) numel(x(x > 0)), num_neighbors);
        
        % make sure candidate cells have the same number of
        % neighboring cells (index is the index of cells_in_embryo)
        candidate_range = find(num_neighbors == this_fit.neighbor_cells);
        
        % If there is only a single candidate, automatically
        % accept
        if numel(candidate_range) == 1
%             display('Only 1')
            cellOI = cellsOI(candidate_range);
            accept_move(this_fit,cellOI);
            already_pulsed(frame,candidate_range) = 1;
        else
            
            while ~accept
                
                % Find candidate
                randomID = candidate_range(randi(numel(candidate_range)));
                cellOI = cellsOI(randomID);
                
                % Check that the current cell doesn't already have
                % a pulse at this time
                if already_pulsed(frame,randomID) == 1,
                    accept = 0;
                else
                    % additional check for neighbor equality
                    num_neighbors = cellOI.identity_of_neighbors_all{ frame };
                    num_neighbors = numel( num_neighbors( num_neighbors > 0 ) );
                    
                    if num_neighbors ~= this_fit.neighbor_cells
                        accept = 0;
                    else
                        
                        if cellOI.num_fits == 0
                            % TODO: should we automatically accept if this
                            % is the candidate's first pulse?
                            
                            frame = this_fit.center_frame;
                            
                            % Accept this move
                            % TODO: modify acceptance
%                             display('First')
                            accept = 1;
                            accept_move(this_fit,cellOI);
                            already_pulsed(frame,randomID) = 1;
                            
                        else
                            
                            % If there is already a pulse in cell, then
                            % check for interval between pulses
                            interval = this_fit.center - ...
                                max( [pulse_sim.find_fits_from_cell(cellOI).center] );

                            if interval < 0, keyboard; end
                            
                            % Figure out if input frequency is a histogram or not
                            if isfield(freqHat,'bin') && ~isfield(freqHat,'fun')
                                idx = findnearest(freqHat.bin,interval);
                                p = freqHat.prob(idx);
                            elseif ~isfield(freqHat,'bin') && isfield(freqHat,'fun')
                                p = feval(freqHat.fun,interval);
                            end
                            % How to accept?
                            if rand >= p %rand generates a random uniform number [0,1]
                                % Decline move
                                accept = 0;
                            else
                                % Accept this move
%                                 display('MC move')
                                accept_move(this_fit,cellOI);
                                already_pulsed(frame,randomID) = 1;
                                accept = 1;
                                
                            end % whether to accept based on random number
                            
                        end % Condition on distribution of intervals
                        
                    end % Condition on distribution of number of neighboring cells
                    
                end % make sure cell doesn't already have a pulse
                
            end % while loop for accepting move
            
        end % accept if only 1 candidate
        
    end % Loop over all fits within embryo
    
%     pulse_bs(e).fits = fitsOI;
%     pulse_bs(e).cells = cellsOI;
    pulse_sim(e).bootstrapped = 1;
    
end % Loop over all embryos

    function accept_move(f,c)
        % Output unnecessary since passed by reference
        c.flag_tracked = 1;
        c.flag_fitted = 1;
%         f.stackID = c.stackID;
        f.cellID = c.cellID;
        f.bootstrapped = 1;
        
        if isnan(c.num_fits), c.num_fits = 0; end
        
        c.num_fits = c.num_fits + 1;
        c.fitID = [c.fitID f.fitID];
    end

end % simulate_pulsing