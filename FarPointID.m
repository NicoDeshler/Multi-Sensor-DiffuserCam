function [r] = FarPointID()

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
import ZOSAPI.*;
r = [];
% OpticStudio session variables
TheSystem = TheApplication.PrimarySystem;
TheAnalyses = TheApplication.PrimarySystem.Analyses;
TheLDE = TheSystem.LDE;
SysData = TheSystem.SystemData.Fields;


DOF_farPTS = zeros(1,5);
for j = 1:5

% Load in .zmx lens file from provided filepath
lens = sprintf('%i x %i cm sensor.zmx',j,j);
lensFile = strcat('C:\Users\Deshler\Documents\Zemax\Samples\DiffuserCam Setups\DOF v Sensor Size Experiment Setups\', lens);
TheSystem.LoadFile(lensFile, false)

% Path for saving analysis results
timestamp = datestr(now,'mm-dd-yy');
folder = strcat('C:\Users\Deshler\Documents\Zemax\ZOS-API Projects\MATLABStandaloneApplication\DiffuserCam Analyses', timestamp, '\');
if ~exist(folder,'dir')
    mkdir(folder,timestamp);
end
saveAs = '1x1cm sensor.mat';    % The filename your results will be saved as
path = strcat(folder,'\', timestamp,'\', saveAs);

% Surfaces
PS = TheLDE.GetSurfaceAt(1);    % Point Source to diffuser spacer surface

% Geom Image analysis
geom = TheAnalyses.New_Analysis_SettingsFirst(ZOSAPI.Analysis.AnalysisIDM.GeometricImageAnalysis);


%% COMPARING SIMILARITY TO PSF AT INFINITY

%{
% Set point source distance to infinity
PS.Thickness = Inf;

% Run the Geometric Image Analysis
geom.ApplyAndWaitForCompletion();
geom_results = geom.GetResults();
geom_grid = geom_results.GetDataGrid(0);
geom_data = geom_grid.Values.double;

% Point source at infinity
PSF_inf = geom_data;
inf_autocorr = sum(PSF_inf(:).*PSF_inf(:));


start_z = 0;
end_z = 10;
n = 10;

threshold = 0.2;
far_point = inf;

while end_z - start_z > 0.1
    step = (end_z - start_z)/n;
    zrange = linspace(start_z,end_z,n);
    for i = 1:n
        %Set point source distance
        PS.Thickness = zrange(i);

        geom.ApplyAndWaitForCompletion();
        geom_results = geom.GetResults();
        geom_grid = geom_results.GetDataGrid(0);
        geom_data = geom_grid.Values.double;

        % PSF at z
        PSF_z = geom_data;
        
        z_xcorr = sum(PSF_inf(:).*PSF_z(:));
        if z_xcorr/inf_autocorr >= threshold
            far_point = zrange(i);
            disp(far_point)
            break
        end
    end

    % update start_z and end_z
    if far_point == inf
        end_z = 10 * end_z;
    else
        start_z = far_point - step;
        end_z = far_point;
    end
end

%}

%% COMPARING ADJACENT PSFs -- 'SIMILARITY' DETERMINED WITH THRESHOLD VALUE
z = 2.5; % in cm -- %starting depth for point source
step = 0.1; % in cm
threshold = 0.7;

% Run geometric analysis
PS.Thickness = z; 
geom.ApplyAndWaitForCompletion();
geom_results = geom.GetResults();
geom_grid = geom_results.GetDataGrid(0);
geom_data = geom_grid.Values.double;

% Set PSF starting conditions
PSF_curr = zeros(1000);
PSF_next = geom_data;

% AutoCorrelation and X-Correlation
PSF_Acorr = dot(PSF_next, PSF_next);
PSF_Xcorr = dot(PSF_curr, PSF_next);
metric = PSF_Xcorr/PSF_Acorr;


while metric < threshold
    
    PSF_curr = PSF_next;
    z = z + step;
    PS.Thickness = z;
    
    geom.ApplyAndWaitForCompletion();
    geom_results = geom.GetResults();
    geom_grid = geom_results.GetDataGrid(0);
    geom_data = geom_grid.Values.double;
    
    % Set next PSF
    PSF_next = geom_data;
    
    % Consider only peripherals of PSF <<<< REMOVE TO CONSIDER FULL PSF
    L = j; % sensor dimension
    PSF_curr = Peripheral(PSF_curr, L);
    PSF_next = Peripheral(PSF_next, L);
    
    % Compute correlations
    PSF_Acorr = sum(PSF_next(:).*PSF_next(:));
    PSF_Xcorr = sum(PSF_next(:).*PSF_curr(:));
    metric = PSF_Xcorr/PSF_Acorr;
    
end
disp(z)
DOF_farPTS(j) = z;
end
figure
plot(DOF_farPTS)
t = sprintf('_Threshold = %d',threshold);
info = strcat('DOF Far Points vs. Sensor Size', t);
title(info)

end


function M = Peripheral(PSF, L)
    pix_cm_ratio = 1000/L;
    aper = round((L*3/4) * pix_cm_ratio);
    if mod(aper,2) == 1
        aper = aper + 1;
    end
    p = (1000 - aper)/2;
    z = zeros(aper);
    peripheral = padarray(z,[p p],1,'both');
    M = peripheral.* PSF; 
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
