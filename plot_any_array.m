Mx = img2mx('Data/gggrx_0660pm_anom_l320.img', 4); %array to plot
% Plot
figure;
imagesc([-180 180], [-90 90], Mx);
set(gca, 'YDir', 'normal');  % so latitude increases upward
%demcmap(Mx, 128); % Using a better colormap for topography visualization
colorbar();

%title('Crustal Thickness');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
