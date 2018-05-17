function [r] = MATLABStandaloneApplication()

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
% BEGINAPPLICATION Contains all user input variables and handles
% input errors. Executes specified analyses and handles data saving.

import ZOSAPI.*;
TheSystem = TheApplication.PrimarySystem;
TheAnalyses = TheApplication.PrimarySystem.Analyses;
TheLDE = TheApplication.PrimarySystem.LDE; 

% Load in lensfile files from provided locations
lensFile = 'C:\Users\Deshler\Documents\Zemax\Samples\DiffuserCam Setups\Setups with MatLab Diffuser\Basic DiffuserCam.zmx'; %where the lens file is located
saveFolder = 'C:\Users\Deshler\Documents\Zemax\ZOS-API Projects\MATLABStandaloneApplication\DiffuserCam Analyses'; %where analysis results will be stored
if ~exist(saveFolder, 'dir')
    HandleError('Provided save folder does not exist.');
end
TheSystem.LoadFile(lensFile, false)

% The list of user-defined analyses functions
a1 = @PhysicalOpticsProp;
a2 = @SpotDiagram;
a3 = @SampleAnalysis3;
a4 = @SampleAnalysis4;

%Specify which analyses to run here e.g. {a1, a2, ...}
toRun = {a1}; 

% Run the specified analyses and save stack video if 
% visualize is true
visualize = false;
for aN = 1 : size(toRun)
    [data, filename] = toRun{aN}(TheAnalyses, TheLDE);
    saveResults(data, filename, saveFolder);    
    if visualize
        avi = makeVideo(data, filename);
        saveResults(avi, avi.Filename, saveFolder);
    end
end

r = [];
end


function [zStack, filename] = PhysicalOpticsProp(TheAnalyses, TheLDE)
% Surfaces
s0 = TheLDE.GetSurfaceAt(0);
s2 = TheLDE.GetSurfaceAt(2);

% Intantiate the POP Analysis object
pop = TheAnalyses.New_Analysis_SettingsFirst(ZOSAPI.Analysis.AnalysisIDM.PhysicalOpticsPropagation);

% Configure POP Settings | NOTE: Make sure you adjust your POP settings
% in your lens file using the OpticStudio suite and 'save' the settings.
% Few of the POP settings are easily adjustable in the standalone connection. 
pop_settings = pop.GetSettings();


maxDist = 20; %millimeters 
minDist = 0; %millimeters
n = 40; 
i = 1;
zStack = zeros(512,512,n); %Specify POP dimensions

% Image Plane Depth: linspace(max, min, n)
for z = linspace(maxDist, minDist, n) 
    
    % Set image plane depth
    s2.Thickness = z;
    
    % Run the POP Analysis
    pop.ApplyAndWaitForCompletion();
    pop_results = pop.GetResults();
    pop_grid = pop_results.GetDataGrid(0);
    pop_data = pop_grid.Values.double;
        
    % Push POP data onto Z stack
    zStack(:,:,i) = pop_data;
    i = i + 1;
end
filename = sprintf('POPData (Max_%dmm Min_%dmm N_%d)', maxDist, minDist, n);
end


function [data, filename] = SpotDiagram(TheAnalyses, TheLDE)
% Surfaces

% Spot Diagram Analysis
spot = TheAnalyses.New_StandardSpot();

% Configure Spot Settings
spot_settings = spot.GetSettings();
spot_settings.Pattern = ZOSAPI.Analysis.Settings.Spot.Patterns.Square;

% Run the Spot Analysis
spot.ApplyAndWaitForCompletion();
spot_results = spot.GetResults();
spot_data = spot_results.SpotData;

data = spot_data;
filename = 'SpotData';    
end


function saveResults(data, filename, saveFolder)
% SAVERESULTS saves data collected in an analysis to the  
timestamp = datestr(now,'mm-dd-yy');
subdir = strcat(saveFolder, timestamp,'\'); %ADD ANOTHER '\' BEFORE THE TIMESTAMP
if ~exist(subdir, 'dir')
    mkdir(subdir)
end
save(strcat(subdir, filename), 'data');
LogMessage(strcat(filename, 'has been saved.'))
end


function v = makeVideo(data, filename)
frameNum = size(data, 3);
video = zeros(1,frameNum); 
for i = 1 : frameNum
    frame = imagesc(data(:,:,i));
    video(i) = frame;
end    
    
vidName = strcat(filename,'Video');
v = VideoWriter(vidName);
v.FrameRate = 10;
open(v)
try writeVideo(v, video)
catch writeErr
    HandleError(writeErr)
end
close(v)
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
