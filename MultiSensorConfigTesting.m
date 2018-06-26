function [r] = MultiSensorConfigTesting()

% Initialize the OpticStudio connection
TheApplication = InitConnection();
if isempty(TheApplication)
    % failed to initialize a connection
    r = [];
else
    try
        r = BeginApplication(TheApplication);
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

function [r] = BeginApplication(TheApplication)
r = 0; import ZOSAPI.*;

% OpticStudio session variables
TheSystem = TheApplication.PrimarySystem;
TheAnalyses = TheApplication.PrimarySystem.Analyses;
TheLDE = TheSystem.LDE;
SysData = TheSystem.SystemData.Fields;

% Load in .zmx lens file from provided filepath
lensFile = 'C:\Users\Deshler\Documents\Zemax\Samples\DiffuserCam Setups\Reconstruction Testing Setups\5 x 5 cm sensor - 300 pixels.zmx';
TheSystem.LoadFile(lensFile, false)

% Save Folder
timestamp = datestr(now,'mm-dd-yy');
folder = strcat('C:\Users\Deshler\Documents\Zemax\ZOS-API Projects\MATLABStandaloneApplication\DiffuserCam Analyses', timestamp, '\');   
if ~exist(folder,'dir')
    mkdir(folder,timestamp);
end

% File names for saving workspace variables
saveStack = '5x5 cm Diffuser_300 pix - PSF zStack.mat';
saveScene = '5x5 cm Diffuser_300 pix - Scene.mat';

% Function booleans
makeStack = true;       % Create axial PSF zstack
makeScene = true;        % Create and measure a scene 

% Surfaces
s1 = TheLDE.GetSurfaceAt(1);    % Point Source to Diffuser spacing

% Fields
f1 = SysData.GetField(1);

% Analyses
geom = TheAnalyses.New_Analysis_SettingsFirst(ZOSAPI.Analysis.AnalysisIDM.GeometricImageAnalysis);

% Sensor array specifications
% S_DIM       -   Size of each sensor (measured in .zmx file's lens units)
% S_SPACING   -   Spacing between sensors (measured in .zmx file's lens units)
% PIXELS      -   Output matrix size of Zemax analysis results
% MS_DIM      -   Dimensions of the sensor grid (measured in .zmx file's lens units) 
% MS_SLOTS    -   Integer number of sensor positions per row in the square sensor array
% MS_CFG      -   Binary matrix indicating sensor array configuration

% s_dim = 2 * s7.ApertureData.CurrentTypeSettings.XHalfWidth;     
s_dim = 1;
s_spacing = 0;
pixels = 300;
ms_slots = 5;
ms_dim = (s_dim + s_spacing) * ms_slots;     
ms_cfg = GetSensorConfig('checkerboard', ms_slots);

%{
% Analysis Settings 
% NOTE more settings are adjustable through the OpticStudio UI. 
% Open the lens file through Zemax and modify/save your settings there
% before running the analysis in the standlone application.
cfgFile = 'C:\Users\Deshler\Documents\Zemax\Samples\DiffuserCam Setups\DOF v Sensor Size Experiment Setups\5 x 5 cm sensor.CFG';
geomSettings = geom.GetSettings();


%IMA_FIELD      -   The field size
%IMA_IMAGESIZE  -   The image size
%IMA_IMANAME    -   The image file name
%IMA_NA         -   The numerical aperture
%IMA_KRAYS      -   The number of rays x 1000
%IMA_OUTNAME    -   The output file name
%IMA_SURFACE    -   The surface number

geomSettings.ModifySettings(cfgFile, 'IMA_FIELD', '0');
geomSettings.ModifySettings(cfgFile, 'IMA_SIZE', int2str(ms_dim));
geomSettings.ModifySettings(cfgFile, 'IMA_NAME', 'circle');
geomSettings.ModifySettings(cfgFile, 'IMA_KRAYS', '5000');
%}
        
% Convert lens unit dimensions to pixels
u2pix = @(x)round(pixels/ms_dim * x);

% CFG_MASK filters which pixels from a fully tiled sensor array will
% be recorded given the multi-sensor configuration.
cfg_mask = zeros(pixels);
s_pix = u2pix(s_dim);
pad = round(u2pix(s_spacing)/2);
sensor = padarray(ones(s_pix), [pad, pad], 0);
for i = 1:ms_slots
    for j = 1:ms_slots
       if ms_cfg(i,j) == 1
           d = size(sensor,1);
           r = (i-1)*d + 1; rr = i*d;
           c = (j-1)*d + 1; cc = j*d;
           cfg_mask(r:rr, c:cc) = sensor;
       end  
    end
end

if makeStack
    start_xyz = [0,0,1];
    end_xyz = [0,0,4.2];
    n = 32;
    MS_PSFstack = GeneratePSFstack(start_xyz,end_xyz, n);
    
    path = strcat(folder,'\', timestamp,'\', saveStack);
    save(path, 'MS_PSFstack', 'start_xyz', 'end_xyz', 'n');
end

if makeScene
    FOV = pi/6;
    DOF = 20;
    DOF_min = 2;
    m = 5;
    [PS_xyz, measurement] = GenerateScene(FOV, DOF, DOF_min, m);
    
    path = strcat(folder,'\', timestamp,'\', saveScene);
    save(path, 'PS_xyz', 'measurement');
    imwrite(measurement, strcat(folder,'\', timestamp,'\', 'measurement.png'));
end


    function stack = GeneratePSFstack(start_xyz, end_xyz, n)
    % Point source translates from START_XYZ coordinate to END_XYZ coordinate
    % in N steps while CURR_XYZ tracks the immediate position of the source.
    % Hence, only linear trajectories for the point source are possible.
    stack = zeros(pixels, pixels, n);
    curr_xyz = start_xyz;
    step_vec = (end_xyz - start_xyz)/n; 

    for pos = 1:n
        % Set lateral and axial position of point source
        f1.X = curr_xyz(1);
        f1.Y = curr_xyz(2);
        s1.Thickness = curr_xyz(3);
        
        % Run the Geometric Image Analysis
        geom.ApplyAndWaitForCompletion();
        geom_results = geom.GetResults();
        geom_grid = geom_results.GetDataGrid(0);
        geom_data = geom_grid.Values.double;
        
        %Add masked PSF to stack
        stack(:,:,pos) = cfg_mask.*geom_data;

        % Increment Point Source position
        curr_xyz = curr_xyz + step_vec;
    end
    end

    function [PS_xyz, measurement] = GenerateScene(FOV, DOF, DOF_min,  m)
    % Simulate a measurement generated by an assortment of M randomly
    % positioned point sources.
    measurement = zeros(pixels);
         
    % FOV defines the maximum half angle (in radians) of the scene
    % DOF defines the range of z values visible in the scene
    % DOF_min defines the starting axial distance for the DOF
    PS_xyz = (rand(3,m) * DOF)+ DOF_min;
    PS_xyz(1,:) = (PS_xyz(3,:)* tan(FOV)).*PS_xyz(1,:);
    PS_xyz(2,:) = (PS_xyz(3,:)* tan(FOV)).*PS_xyz(1,:);
    
    % Display PS locations
    figure
    scatter3(PS_xyz(1,:), PS_xyz(2,:), PS_xyz(3,:));

    for ps = PS_xyz

        % Set lateral and axial position for point source
        f1.X = ps(1);
        f1.Y = ps(2);
        s1.Thickness = ps(3);

        % Run the Geometric Image Analysis
        geom.ApplyAndWaitForCompletion();
        geom_results = geom.GetResults();
        geom_grid = geom_results.GetDataGrid(0);
        geom_data = geom_grid.Values.double;

        measurement = measurement + cfg_mask.*geom_data;
    end
    end

end

function tiling = GetSensorConfig(BoardPatn, ms_slots)
tiling = zeros(ms_slots);

switch BoardPatn
    case 'full'
        tiling = ones(ms_slots);
    case 'checkerboard'
        tiling = invhilb(ms_slots) < 0;
    case 'stripes'
        for j = 1:ms_slots
            if mod(j,2) == 0
                tiling(j,:) = 1;
            end
        end
    %{    
    case 'squares'
        tiling = ones(2);
        for j = 1:ceil(ms_slots/2)
            if mod(j,2)
                tiling = padarray(tiling, [1,1], 0, 'both');
            else
                tiling = padarray(tiling, [1,1], 1, 'both');
            end
        end
        %}
    case 'random'
        tiling = randi([0,1],ms_slots);
    
    otherwise
        HandleError('Provided sensor tiling does not exist.')
end

if size(tiling,1) ~= ms_slots 
    HandleError('Invalid sensor tiling. Must be a binary square matrix with size = array_dim.')
end
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
