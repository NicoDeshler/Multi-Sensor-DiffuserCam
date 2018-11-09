function resized = PadCropResize(img, resizeTo)
% Resizes an image by zero-padding in the deficient dimensions and cropping
% in the excess dimensions. The center of the image remains fixed.  
half_dif1 = ceil((resizeTo - size(img))/2);
fill = max(half_dif1, 0);
img = padarray(img, max(fill,0), 0);

half_dif2 = floor(abs(resizeTo - size(img))/2);
kill = max(half_dif2, 1); 
resized = imcrop(img, [fliplr(kill), fliplr(resizeTo-1)]); % Flip due to imcrop 'rect' argument: rect= [x_min, y_min, width, height]

end