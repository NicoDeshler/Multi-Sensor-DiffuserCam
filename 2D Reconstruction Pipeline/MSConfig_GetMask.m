function mask = MSConfig_GetMask(tiling_patn, sensor_array, pixels, sensor_pix, spacing_pix)
%MSCONFIG returns a multi-sensor configuration from a database of grid-wise
%sensor placements.
% t_pix = total pixels in zemax analysis
% s_pix = sensor size dimensions in pixels
% spc_pix = sensor spacing dimensions in pixels

% Get binary tile configuration
T = getTiling(tiling_patn, sensor_array);

% Create binary pixel mask indicating the active pixels in the sensor array
mask = [];
temp_row = [];
pad = round(spacing_pix/2);
active = padarray(ones(sensor_pix), pad, 0);
empty = padarray(zeros(sensor_pix), pad, 0);
for i = 1:sensor_array(1)
    for j = 1:sensor_array(2)
       if T(i,j) == 1
           temp_row = horzcat(temp_row, active);
       else 
           temp_row = horzcat(temp_row, empty);
       end
    end
   mask = vertcat(mask,temp_row);
   temp_row = [];
end
mask = PadCropResize(mask, pixels);
%{
half_dif1 = ceil((pixels - size(mask))/2);
fill = max(half_dif1, 0);
mask = padarray(mask, max(fill,0), 0);

half_dif2 = floor(abs(pixels - size(mask))/2);
kill = max(half_dif2, 1); % Flip due to imcrop 'rect' argument: rect= [x_min, y_min, width, height]
mask = imcrop(mask, [fliplr(kill), fliplr(pixels-1)]);
%}
end


function T = getTiling(patn, slots)

T = zeros(slots);

switch patn
    case 'full'
        T = ones(slots);
    case 'checkered'
        T = invhilb(slots) < 0;
    case 'vertStriped'
        oneCols = mod(1:slots(2), 2) == 0;
        T = repmat(oneCols, [slots(1), 1]);
    case 'horzStriped'
        oneRows = mod((1:slots(2))', 2) == 0;
        T = repmat(oneRows, [1, slots(1)]);
    case 'random'
        T = randi([0,1],slots);
    
    otherwise
        HandleError('Provided sensor tiling does not exist.')
end

if size(T) ~= slots 
    HandleError('Invalid sensor tiling. Must be a binary square matrix with size = array_dim.')
end
end
