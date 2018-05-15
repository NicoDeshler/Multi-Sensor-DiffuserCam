function [ r ] = MATLABStandaloneApplication( args )

if ~exist('args', 'var')
    args = [];
end

% Initialize the OpticStudio connection
TheApplication = InitConnection();
if isempty(TheApplication)
    % failed to initialize a connection
    r = [];
else
    try
        r = BeginApplication(TheApplication, args);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   USER INPUT BEGINS BELOW    %%%%%%%%%%%%%%%%%%%%%%%%%



function [r] = BeginApplication(TheApplication, args)
% BEGIN APPLICATION Run a session of Zemax Optic Studio. 
% User inputs include 
% The lens file directory and save folder directory.

import ZOSAPI.*;

% Initiate Primary System
TheSystem = TheApplication.PrimarySystem;

% File/Folder locations
LensFile = 'C:\Users\Deshler\Documents\Zemax\Samples\Diffuser Setups\MatLab Diffuser\Basic DiffuserCam.zmx';
SaveFolder = 'PUT YOUR SAVE FOLDER HERE';

% Load in the lens file 
TheSystem.LoadFile(LensFile, false);

% Run specified analyses
Analyze(TheApplication);

disp("Analysis complete. Data files saved.")
r = [];

end


function Analyze(TheApplication)
% ANALYZE Runs analyses specified by the user in the function body.

% Session Variables
TheAnalyses = TheApplication.PrimarySystem.Analyses;
TheLDE = TheApplication.PrimarySystem.LDE; 

% Surfaces
Source = TheLDE.GetSurfaceAt(0);
ImageSpacing = TheLDE.GetSurfaceAt(2);
i = 1;

% Image Plane Depth: linspace(max, min, n)
for z = linspace(1,0.1,20) 
    
    % Set image plane depth
    ImageSpacing.Thickness = z;
    disp(z)
    
    % Physical Optics Propogation Analysis
    pop = TheAnalyses.New_Analysis_SettingsFirst(ZOSAPI.Analysis.AnalysisIDM.PhysicalOpticsPropagation);

    % Configure POP Settings | NOTE: Make sure you adjust your POP settings
    % in your lens file using the OpticStudio suite and 'save' the settings.
    % Few of the POP settings are easily adjustable in the standalone connection. 
    pop_settings = pop.GetSettings();
    
    % Run the POP Analysis
    pop.ApplyAndWaitForCompletion();
    pop_results = pop.GetResults();
    pop_grid = pop_results.GetDataGrid(0);
    pop_data = pop_grid.Values.double;
                
    
    % Spot Diagram Analysis
    spot = TheAnalyses.New_StandardSpot();

    % Configure Spot Settings
    spot_settings = spot.GetSettings();
    spot_settings.Pattern = ZOSAPI.Analysis.Settings.Spot.Patterns.Square;

    % Run the Spot Analysis
    spot.ApplyAndWaitForCompletion();
    spot_results = spot.GetResults();
    spot_data = spot_results.SpotData;
    
    
    WriteData(pop_data, z);
    
    imagesc(pop_data)
    Mov(i) = getframe(gcf);
    i = i + 1;
end

WriteMov(Mov);

end

function WriteMov(figures)
v = VideoWriter('CloseCaustics.avi', 'Uncompressed AVI');
v.FrameRate = 15;
open(v)
writeVideo(v, figures)
close(v)
end


% Write data to SaveDir as .mat file. Filename contains source depth info.
function WriteData(pop_data, z)

timestamp = datestr(now,'mm-dd-yy');
% AnalysisDir = strcat('YOUR SAVE FOLDER PATH', timestamp, '\');
AnalysisDir = strcat('C:\Users\Deshler\Documents\MATLAB\DiffuserCam Simulations\Zemax Analysis_', timestamp, '\');
if ~exist(AnalysisDir)
    mkdir(AnalysisDir)
end

% Save raw data
filename = strcat(AnalysisDir, sprintf('pop_data @ %d microns', int64(z*1000)));
% save(filename, 'pop_data')

% Save color-scale image
%colorIm = figure;
%imagesc(pop_data);
%saveas(colorIm, filename, 'png');
%close all

disp("POP data written to folder.")

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


