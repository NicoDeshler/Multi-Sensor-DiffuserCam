%Load and prepare impulse stack
load('./MScfg PSFzstack and Measurement/5x5 cm Diffuser_300 pix - Scene.mat');

% modify DiffuserCam_settings.m to use measurement and PSF zstack
% Point Source coordinate format:
% x1 | x2 | x3 | x4 | x5 | x6
% y1 | y2 | y3 | y4 | y5 | y6
% z1 | z2 | z3 | z4 | z5 | z6
[xhat, f] = DiffuserCam_main('DiffuserCam_settings.m');
temp = xhat;
num_pts = size(PS_xyz,2);
reconst = zeros(num_pts, 3);

% Brightest voxels are assumed to represent the reconstructions estimation
% of the point source positions in the scene.
for i = 1: num_pts
   
    [y,x,z] = ind2sub(size(temp),find(temp == max(max(max(temp))))); %y,x,z ordered this way since the ind2sub() function returns format (row, col, z)  
    
    %volumetric removal of box = r 
     r = 3;
     %min([11,y,x,z])+1;
    
    temp(y-r:y+r,x-r:x+r,z-r:z+r) = 0;
    reconst(i,1) = x;
    reconst(i,2) = y;
    reconst(i,3) = z;

end
diff_xyz = transpose(PS_xyz) - reconst;
L2norm_vec = sqrt(sum(diff_xyz.^2,2));
meanErr = mean(L2norm_vec);