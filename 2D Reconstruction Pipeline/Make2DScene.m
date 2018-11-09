function Make2DScene(z_plane, num_sources)
%Simulates a sparsely populated scene with n_srcs point
%sources. Saves an image of the ground truth, the PSF, and ms_measurement/std_measurement.


% Initialize the OpticStudio connection
TheApplication = InitConnection();
if isempty(TheApplication)
    % failed to initialize a connection
else
    try
        BeginApplication(TheApplication,z_plane, num_sources);
        CleanupConnection(TheApplication);
    catch err
        CleanupConnection(TheApplication);
        rethrow(err);
    end
end
end

function app = InitConnection()

import System.Reflection.*;

% Find the installed version of OpticStudio.

% This method assumes the helper dll is in the .m file directory.
% p = mfilename('fullpath');
% [path] = fileparts(p);
% p = strcat(path, '\', 'ZOSAPI_NetHelper.dll' );
% NET.addAssembly(p);

% This uses a hard-coded path to OpticStudio
NET.addAssembly('C:\Users\Deshler\Documents\Zemax\ZOS-API\Libraries\ZOSAPI_NetHelper.dll'); 
success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize();
% Note -- uncomment the following line to use a custom initialization path
% success = ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize('C:\Program Files\OpticStudio\');
if success == 1
    LogMessage(strcat('Found OpticStudio at: ', char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory())));
else
    app = [];
    return;
end

% Now load the ZOS-API assemblies
NET.addAssembly(AssemblyName('ZOSAPI_Interfaces'));
NET.addAssembly(AssemblyName('ZOSAPI'));

% Create the initial connection class
TheConnection = ZOSAPI.ZOSAPI_Connection();

% Attempt to create a Standalone connection

% NOTE - if this fails with a message like 'Unable to load one or more of
% the requested types', it is usually caused by try to connect to a 32-bit
% version of OpticStudio from a 64-bit version of MATLAB (or vice-versa).
% This is an issue with how MATLAB interfaces with .NET, and the only
% current workaround is to use 32- or 64-bit versions of both applications.
app = TheConnection.CreateNewApplication();
if isempty(app)
   HandleError('An unknown connection error occurred!');
end
if ~app.IsValidLicenseForAPI
    HandleError('License check failed!');
    app = [];
end
end


%% Boilerplate Code
% NOTE: The abreviation ILU stands for "In Lens Units" and refers to the
% units of measurement used in the .zmx lens file.
function BeginApplication(TheApplication, z_plane, num_sources)
import ZOSAPI.*;

% OpticStudio session variables
TheSystem = TheApplication.PrimarySystem;
TheAnalyses = TheApplication.PrimarySystem.Analyses;
TheLDE = TheSystem.LDE;
SysData = TheSystem.SystemData.Fields;

% Configure lens setup variables and output directories
run('Make2DScene_settings.m')
outdir = [dest, 'num_srcs_',int2str(num_sources),'\'];
mkdir(outdir);

% Load in lens file
success = TheSystem.LoadFile(lensFile, false);
if ~success
    HandleError('Failed to load lens file. Check provided filepath.')
end


% Surfaces
source_surf = TheLDE.GetSurfaceAt(source_surfNum);    % Point Source to Diffuser spacing
aper_surf = TheLDE.GetSurfaceAt(aper_surfNum);      % System Aperture
image_surf = TheLDE.GetSurfaceAt(image_surfNum);    % Diffuser surface

% Fields
f1 = SysData.GetField(1);


% Set Image Plane Aperture
Rect_Aper = image_surf.ApertureData.CreateApertureTypeSettings(ZOSAPI.Editors.LDE.SurfaceApertureTypes.RectangularAperture);
Rect_Aper.S_RectangularAperture_.XHalfWidth = imaSize(2)/2;
Rect_Aper.S_RectangularAperture_.YHalfWidth = imaSize(1)/2;
image_surf.ApertureData.ChangeApertureTypeSettings(Rect_Aper);


% Dimension assertion
imaAper = 2*[image_surf.ApertureData.CurrentTypeSettings.XHalfWidth,...
            image_surf.ApertureData.CurrentTypeSettings.YHalfWidth];
if imaSize ~= imaAper
    HandleError('Analysis image size and image surface aperture are not equal. Aperture assignment failed.')
end


% Set System Aperture
effective_dim = (sensor_dim + spacing_dim) .* sensor_array - spacing_dim; % Dimensions of the complete array.
effective_pix = u2pix(effective_dim);
Rect_Aper = aper_surf.ApertureData.CreateApertureTypeSettings(ZOSAPI.Editors.LDE.SurfaceApertureTypes.RectangularAperture);
Rect_Aper.S_RectangularAperture_.XHalfWidth = (1-percent_aper)*effective_dim(1,1)/2;
Rect_Aper.S_RectangularAperture_.YHalfWidth = (1-percent_aper)*effective_dim(1,2)/2;
aper_surf.ApertureData.ChangeApertureTypeSettings(Rect_Aper);

% Instantiate Zemax analysis object
geom = TheAnalyses.New_Analysis_SettingsFirst(ZOSAPI.Analysis.AnalysisIDM.GeometricImageAnalysis);

% Set field at z-plane
source_surf.Thickness = z_plane;

% Make/save the PSF

% Center the field
f1.Y = 0;
f1.X = 0;

% Run the Geometric Image Analysis
geom.ApplyAndWaitForCompletion();
geom_results = geom.GetResults();
geom_grid = geom_results.GetDataGrid(0);
geom_data = geom_grid.Values.double;

PSF = geom_data;



%% MEASUREMENTS
% Generate random xy coordinate pairs for the sources on z_plane
sources_yx = rand([num_sources, 2]);
sources_yx(:,1:2) = (z_plane* tan(FOV/2))*(2*(sources_yx(:,1:2)-(1/2)));   % y and x coords

% Make the ground truth image
% Matrix indices corresponding to xy cartesian space positions - Const is a
% lateral position scaling factor for a given z_plane.
const = pixels/imaSize * focal_len/z_plane;
sources_ji = round(sources_yx*const);
ground_truth = zeros(pixels);
lin_idx = sub2ind(size(ground_truth),sources_ji(:,1) + floor(pixels(1)/2), sources_ji(:,2) + floor(pixels(2)/2));
ground_truth(lin_idx) = 1;


% Make the standard and ms measurements
% NOTE: If using x,y coordinate pairs for defining field position, ensure
% that your lens file has the field type set to OBJECT, not RADIAL
measurement = zeros(pixels);
for ps = transpose(sources_yx) 
        % Set lateral and axial position for point source
        f1.Y = ps(1);
        f1.X = ps(2);
        
        % Run the Geometric Image Analysis
        geom.ApplyAndWaitForCompletion();
        geom_results = geom.GetResults();
        geom_grid = geom_results.GetDataGrid(0);
        geom_data = geom_grid.Values.double;

        measurement = measurement + geom_data;
end

% Resize standard measurement to only consider the effective sensor size
padSize = ceil((pixels - effective_pix)/2);
std_measurement = measurement .* PadCropResize(padarray(ones(effective_pix), padSize, 0), pixels);

% Create and apply the multi-sensor configuration mask
spacing_pix = u2pix(spacing_dim);
ms_mask = MSConfig_GetMask(array_cfg, sensor_array, pixels, sensor_pix, spacing_pix);
ms_measurement = measurement .* ms_mask;


% Write ground truth, PSF, and measurements to .PNG files. 
imwrite(ground_truth, [outdir, 'ground_truth.png'])
imwrite(PSF, [outdir, 'PSF.png'])
imwrite(std_measurement, [outdir, 'std_measurement.png'])
imwrite(ms_measurement, [outdir, 'ms_measurement.png'])

end


function LogMessage(msg)
disp(msg);
end

function HandleError(error)
ME = MXException(error);
throw(ME);
end

function  CleanupConnection(TheApplication)
% Note - this will close down the connection.
% If you want to keep the application open, you should skip this step
% and store the instance somewhere instead.
TheApplication.CloseApplication();
end
