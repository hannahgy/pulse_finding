function fig = plot_single_pulse(pulse,fitID)
% PLOT_SINGLE_PULSE Plot the myosin and area time-series for a
% single pulse.
%
% USAGE:
%   fit.plot_single_pulse;
%   fits.plot_single_pulse(fitID)

% fit = [pulse.fits];
fit = pulse.get_fitID(fitID);
assert( numel(fit) == 1, 'Fit not found!')
% if numel(fit) == 1
%     fitID = fit.fitID;
% else
%     fit = fit.get_fitID(fitID);
% end

%             fig = figure;
fig = gcf;
x = fit.corrected_time;
[ax,h1,h2] = plotyy( ...
    x, fit.corrected_area_norm, x, fit.corrected_myosin );
xlabel('Pulse time (sec)')

set(ax(1),'YTick',-8:4:8)
ylim(ax(1),[-10 10])

set(ax(1),'YColor',[1 0 1]); set(h1,'Color',[1 0 1]);
ylabel(ax(1),'Apical area (\mum^2)')
set(ax(2),'YColor',[0 0.5 0]); set(h2,'Color',[0 0.5 0]);
ylabel(ax(2),'Myosin intensity (a.u.)')

end % plot_single_pulse