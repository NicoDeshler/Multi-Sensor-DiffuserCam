% Make save folders
%dtstamp = datestr(datetime('now'),'mm-dd-yyyy HH_MM_SS');
dest = 'Prototype 2D Recon\2DScenes\';
mkdir(dest);

% Lens File info
lensFile = 'C:\Users\Deshler\Documents\Zemax\Samples\DiffuserCam Setups\RaspberryPi MSDiffuserCams\4 Sensor Setup.zmx';
source_surfNum = 1; 
aper_surfNum = 2;
image_surfNum = 6;
focal_len = 0.8;

% Geometric Image Analysis settings -- (ENSURE THESE AGREE WITH YOUR LENSFILE SETUP)
imaSize = [1,1] * 12.2; % Change coeff to adjust size (imageSize is always square)
pixels = [500,500];  % [pix_per_row, pix_per_col]
u2pix = @(x)round(pixels ./ imaSize .* x);  % Convert lens unit dimensions to pixels

% DiffuserCam reconstruction volume parameters
FOV = pi/2;     % Angular FOV in radians (assumed to be radially symmetric)
DOFmax = 10;     % ILU
DOFmin = 2;     % ILU
DOFrange = DOFmax - DOFmin;

% Multi-Sensor Configuration Parameters (In lens units)
sensor_dim = [2.74, 3.76];       % Dimensions of a single sensor in the array [y-dim, x-dim]
spacing_dim = [5.48, 4.68];
sensor_array = [2, 2];            % Number of rows and columns of small sensors [rows, columns]
percent_aper = 0;                 % Percent of sensor array obscured by the aperture
array_cfg = 'full';               % Array configuration %OPTIONS: 'full', 'checkered', 'striped', 'random'

% Multi-Sensor Configuration Parameters (In pixels)
sensor_pix = u2pix(sensor_dim);

if FOV >= pi || FOV < 0
    error('Field of View exceeds realistic limits.')
end

if percent_aper >= 1 || percent_aper < 0
    error('Percent Aperture is greater than or equal to one. Aperture obscures entire sensing plane.')
end
