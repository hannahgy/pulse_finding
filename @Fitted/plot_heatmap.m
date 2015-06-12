function plot_heatmap(fits,sortname)
% Uses IMAGESC instead of PCOLOR (PCOLOR is werid)

if nargin < 2, sortname = 'amplitude'; end

fits = sort(fits,sortname);

figure

subplot(1,4,1:2)
[X,Y] = meshgrid( fits(1).corrected_time, 1:numel(fits));
pcolor( X,Y, cat(1,fits.corrected_myosin) );
%             imagesc( ...
%                 fits(1).corrected_time, ...
%                 1:numel(fits), ...
%                 cat(1,fits.corrected_myosin) );
shading flat; axis tight; colorbar;
title('Myosin intensity')
xlabel('Pulse time (sec)');
axis xy
%             colormap(pmkmp(255))

subplot(1,4,3:4)
M = cat(1,fits.corrected_area_norm);
pcolor( X,Y, M);
%             imagesc( ...
%                 fits(1).corrected_time, ...
%                 1:numel(fits), ...
%                 cat(1,fits.corrected_area_norm) );
shading flat; axis tight; colorbar;
caxis( [-8 8] );
title('Local area change');
xlabel('Pulse time (sec)');
axis xy
%             colormap(pmkmp(255))

end % plot_heatmap