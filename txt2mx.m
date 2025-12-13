function Mx = txt2mx(filename)
% TXT2MX  Read float grid and return Matrix 
%   TopoMx = txt2mx(filename) expects the file to contain numeric
%   values arranged column-wise and automatically determines the dimensions.

Mx = readmatrix(filename);
Mx = Mx / 1000;          % convert to kilometers
end
