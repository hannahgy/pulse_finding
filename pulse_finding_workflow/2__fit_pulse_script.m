
%%

clear fit_opts
[fit_opts(1:num_embryos).to_fit] = deal('myosin_intensity');
[fit_opts(1:num_embryos).bg] = deal('on');

[fit_opts(1:num_embryos).left_margin] = deal(10);
[fit_opts(1:num_embryos).right_margin] = deal(10);
[fit_opts(1:num_embryos).nan_thresh] = deal(30);
[fit_opts(1:num_embryos).nan_consec_thresh] = deal(4);
[fit_opts(1:num_embryos).end_tol] = deal(30);

[fit_opts(1:num_embryos).alpha] = deal(0.01);
[fit_opts(1:num_embryos).sigma_lb] = deal(8); % in seconds
[fit_opts(1:num_embryos).sigma_ub] = deal(35); % in seconds

%% Perform fitting

fits_raw = fit_gaussians(cells_raw,fit_opts);

%% Instantiate temporary Pulse objects

for i = 1:num_embryos
    
    cellsTmp = cells_raw([cells_raw.embryoID] == i);
    fitsTmp = fits_raw([fits_raw.embryoID] == i);
    
    pulse(i) = Pulse(fitsTmp,cellsTmp,fit_opts(i),input(i));
    
end

%% sub-set of pulses

pulse.retrace( fit_opts );
pulse.measure_fits();
