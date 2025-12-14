function lowRes = compress_to_1deg(highRes, px)
% COMPRESS_TO_1DEG  Resample high-res equirectangular grid to 1 deg spacing.
%   lowRes = compress_to_1deg(highRes, px)
%   highRes : matrix of size (180*px + 1) x (360*px)  (e.g. 721 x 1440 for px=4)
%   px      : pixels per degree (e.g. 4)
%   lowRes  : output size 180 x 361 (lat x lon) with 1 deg resolution
%
%  Assumes data rows correspond to latitudes from +90 down to -90,
%  and columns to longitudes from -180 increasing eastward.

% Validate sizes
[nRows, nCols] = size(highRes);
expectedRows = 180*px + 1;
expectedCols = 360*px;
if nRows ~= expectedRows || nCols ~= expectedCols
    error('highRes must be %d x %d (rows x cols) for px=%d. Got %d x %d.', ...
          expectedRows, expectedCols, px, nRows, nCols);
end

% Original grid coordinates (deg)
% columns: -180 .. 180 - (1/px)  (360*px samples)
dLon = 1/px;
lon_orig = -180 : dLon : (180 - dLon);          % length = 360*px

% rows: +90 down to -90, total 180*px + 1 samples
dLat = 1/px;
lat_orig = 90 : -dLat : -90;                    % length = 180*px + 1

% Target 1-degree grid:
% longitudes: -180 .. 180 in 1° steps -> 361 samples
lon_tgt = -180 : 1 : 180;                       % length = 361
% latitudes: use 1°-centered values between +90 and -90 to avoid endpoint duplications.
% We choose centers at +89.5, +88.5, ..., -89.5 -> 180 samples
lat_tgt = 89.5 : -1 : -89.5;                    % length = 180

% Create meshgrids for interpolation (interp2 expects X=cols, Y=rows)
[LonOrig, LatOrig] = meshgrid(lon_orig, lat_orig);    % size matches highRes
[LonTgt,  LatTgt]  = meshgrid(lon_tgt,  lat_tgt);     % desired output positions

% Interpolate (bilinear). Convert highRes to double for accuracy.
lowRes = interp2(LonOrig, LatOrig, double(highRes), LonTgt, LatTgt, 'linear');

% If you prefer area-averaging rather than interpolation, replace the interp2 call
% with a block-mean approach (requires exact integer downsampling and careful edge handling).
end