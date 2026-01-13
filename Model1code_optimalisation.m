%% Setup

clear;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

fprintf('Loading models\n')

% Optimization output directories
dataDir   = fullfile(HOME, 'Data', 'Model1optdensity');
resultDir = fullfile(HOME, 'Results', 'Model1optdensity');

% Create folders if they do not exist
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

%% Load  Moon data

%LOLA topography radius 1737400 m
LOLAtop = 'LDEM_4.IMG';
resolution_lola = 4;
f = fopen(LOLAtop,"r","ieee-le");
r_moon = 0.5 * fread(f,[360*resolution_lola Inf],'int16')' + 1737400;
disp(size(r_moon))
fclose(f);

%Grail geoid radius 1378000 m gggrx_1200a_geoid_l660.lbl
GRAILgeoid = 'GRAILgeoid.img';
resolution_GRAIL = 16;
f = fopen(GRAILgeoid,"r","ieee-le");
geoid = fread(f, [360*resolution_lola Inf],'float')' + 1738000;
fclose(f);

%Grail SHA data
%% Gravity data Grail gggrx_0900c_sha.tab
data = readmatrix('GRAILdata2.txt');

[rows, cols] = size(r_moon);

target_res = 1;

%Check if dimensions of the moon are okay with the target resoluation
if mod(rows, target_res) ~= 0 || mod(cols, target_res) ~= 0
    error('Input dimensions are not compatible with the target resolution')
end

%Make the smaller matrix

rows_small = rows / target_res;
cols_small = cols / target_res;
geoid_small = zeros(rows_small, cols_small);

for i = 1:rows_small
    for j = 1:cols_small
        row_index = (i-1)*target_res + 1 : i*target_res;
        col_index = (j-1)*target_res + 1 : j*target_res;
        part = geoid(row_index, col_index);
        geoid_small(i, j) = mean(part(:));
    end
end

fprintf('Starting calculations\n')

topography = r_moon - geoid_small;
nmax = 180;

%remove the header of SHA data
data(1,:) = [];

V = data(:,1:4);
V = [0 0 1 0; V];
Vmax = (nmax + 1)*(nmax + 2)/2;
V = V(1:Vmax, :);
V = sortrows(V,2);

%% Convert to spherical harmonic coefficient matrix
lmax = max(V(:,1));
field = zeros(lmax+1, 2*lmax+1);

for k = 1:size(V,1)
    n = V(k,1);
    m = V(k,2);
    C = V(k,3);
    S = V(k,4);

    field(n+1, lmax+1 + m) = C;
    if m ~= 0
        field(n+1, lmax+1 - m) = S;
    end
end

field(1,:) = 0; %C0,0 is zero

fprintf('Completed calculation of SC\n')

fprintf('Model is being made now\n')

%Crust
D = 38000;
D_topo = -D * ones(size(topography));
D_moon = -D * ones(size(r_moon));


rho_c = 2800;
rho_m = 3300;

Model1 = struct();
Model1.number_of_layers = 2;
Model1.name = 'Two_layered_moon_crust';
Model1.GM = 4.90280011526323e12;
Model1.Re = 1738000;
Model1.geoid = 'none';
Model1.nmax = nmax;
Model1.l1.bound = topography;
Model1.l1.dens = rho_c;
Model1.l2.bound = D_topo;
Model1.l2.dens = rho_m;
Model1.l3.bound = -50000;

fprintf('Model has been made\n')

resolution = 4;
latlon = 1/(2*resolution);

latlim = [(-90+latlon) (90-latlon) (1/resolution)];
lonlim = [(0+latlon) (360-latlon) (1/resolution)];
height = 0;
SHbounds = [2 90];

% Cache file for observed gravity synthesis
cacheFile = fullfile(HOME, 'Data', 'model_data_obs.mat');

if isfile(cacheFile)
    fprintf('Loading cached observed gravity model...\n');
    S = load(cacheFile, 'model_data');
    model_data = S.model_data;
else
    fprintf('Computing observed gravity model...\n');
    model_data = model_SH_synthesis(lonlim, latlim, height, SHbounds, V, Model1);
    save(cacheFile, 'model_data', '-v7.3');
    fprintf('Observed gravity model saved.\n');
end

gravity_data = flipud(model_data.vec.R) * 1e5; %unit is mGal

fprintf('Starting iteration\n')

iteration = 0;
max_iterations = 10;
convergence = 100; % mGal
gravity_dif = 1000000;


while max(abs(gravity_dif(:))) > convergence && iteration < max_iterations
    % Save model boundaries (every iteration)
    %modelFile = fullfile(dataDir, ...
    %    sprintf('Model1_boundaries_iter%02d.mat', iteration));
    save(fullfile(dataDir, sprintf('Model1_bounds_iter_%02d.mat', iteration)), 'Model1')
    fprintf('Computing modeled gravity data...\n');
    % Build SH coefficients for the layered model
    V_model = segment_2layer_model(Model1.l1.bound, Model1.l2.bound, Model1.l3.bound, Model1.l1.dens, Model1.l2.dens, D, Model1);
    % Synthesize gravity field
    V_model_data = model_SH_synthesis(lonlim, latlim, height, SHbounds, V_model, Model1);

    gravity_data_model = flipud(V_model_data.vec.R) * 1e5; % mGal
    lon = V_model_data.grd.lon(1,:);
    lat = V_model_data.grd.lat(:,1);

    fprintf('calculating the difference\n')
    gravity_dif = gravity_data - gravity_data_model;
    disp([max(abs(gravity_dif(:))), min(gravity_dif(:))]);

    fprintf('plotting now\n')

    fig2D = figure('Color','w');
    
    ax1 = subplot(2,2,1);
    imagesc(lon,-lat,gravity_data_model);
    colorbar; colormap(ax1,'parula');
    title('Modeled gravity'); ylabel('Latitude'); xlabel('Longitude');
    set(gca,'YDir','normal')
    
    ax2 = subplot(2,2,2);
    imagesc(lon,-lat,gravity_data);
    colorbar; colormap(ax2,'parula');
    title('Observed gravity'); ylabel('Latitude'); xlabel('Longitude');
    set(gca,'YDir','normal')
    
    ax3 = subplot(2,2,3);
    imagesc(lon,-lat,gravity_dif);
    colorbar; colormap(ax3,'parula');
    title('Residual gravity'); ylabel('Latitude'); xlabel('Longitude');
    set(gca,'YDir','normal')
    
    thickness = Model1.l1.bound - Model1.l2.bound;  % meters

    ax4 = subplot(2,2,4);
    imagesc(lon,-lat,thickness / 1e3);   % meters → km
    colorbar; colormap(ax4,'turbo');
    title('Crust thickness'); ylabel('Latitude'); xlabel('Longitude');
    set(gca,'YDir','normal')
    
    t = sgtitle(sprintf('Model 1 – Iteration %d', iteration));
    t.Color = 'k';
    
    % Force black text for saved figures
    set(findall(fig2D,'Type','text'),'Color','k');
    set(findall(fig2D,'Type','axes'),'XColor','k','YColor','k');
    set(findall(fig2D,'Type','colorbar'),'Color','k');
    
    % White background
    set(fig2D, 'Color', 'w');

    fig2DFile = fullfile(resultDir, sprintf('Model1_iter_%02d.png', iteration));
    print(fig2D, fig2DFile, '-dpng', '-r300');

    Model1.l1.dens = Model1.l1.dens + gravity_dif;
    iteration = iteration + 1;
end



save(fullfile(dataDir,sprintf('Model1_FINAL.mat')), 'Model1')
