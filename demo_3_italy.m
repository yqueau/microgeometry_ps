clear
close all

addpath('Toolbox/');

dataset = 'Coin_Italy'; % Dataset that is used
images = 1:15;	% Images that are used
do_write_obj = 1;	% Set to 0 to de-active OBJ writing, to 1 to activate it

%%% Set dataset
% Folder containing the photometric stereo images
images_folder = sprintf('Datasets/%s',dataset);
% File containing the initial lighting and intrinsics
calib_folder = sprintf('Datasets/Calib');


%%% Read calibration files
disp('Reading calibration files');
% Load intrinsics
load(sprintf('%s/K.mat',calib_folder));
if(K(1,3)==0) K = K';end % Matlab uses a transposed definition of intrinsics matrix
% Load light field
i = 1;
ii = images(i);
load(sprintf('%s/light_%02d.mat',calib_folder,ii));
S = zeros(size(S_i,1),size(S_i,2),size(S_i,3),length(images));
S(:,:,:,1) = double(S_i); 
for i = 2:length(images)
	ii = images(i);
	load(sprintf('%s/light_%02d.mat',calib_folder,ii));
	S(:,:,:,i) = double(S_i);
end
clear S_i
% Store in a compact structure
calib.S = S; clear S; % Lighting (nrows x ncols x 3 x nimgs)
calib.K = K; clear K; % Intrinsics (3 x 3)

%%% Read dataset
disp('Reading photometric stereo images');
% Read ambient light image
[nrows,ncols,trois,nimgs] = size(calib.S);
% Read each image and substract ambient light
I = zeros(nrows,ncols,nimgs);
for i = 1:nimgs
	ii = images(i);
	I(:,:,i) = double(imread(sprintf('%s/I_%02d.png',images_folder,ii)));
end
% Store in a compact structure
data.I = I; clear I; % Images (nrows x ncols x nchannels x nimgs)
data.mask = ones(nrows,ncols);

%%% Show the input images
disp('Displaying data');
figure(1001)
subplot(3,1,1)
imshow(uint8(255*data.I(:,:,1)./max(data.I(:))))
title('$$I^1$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,2)
imshow(uint8(255*data.I(:,:,2)./max(data.I(:))))
title('$$I^2$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,3)
imshow(uint8(255*data.mask))
title('$$\Omega$$','Interpreter','Latex','Fontsize',18);
drawnow

%%% Set parameters
disp('Setting parameters');
params.precond = 'cmg'; % Use multigrid preconditioner
params.z0 = 70.6383*ones(size(data.I,1),size(data.I,2)); % Initial depth map: a plane at 70mm from camera

%%% Solve photometric stereo
disp('Solving photometric stereo');
[XYZ,N,rho] = microgeometry(data,calib,params);
% Scale albedo for visualization
rho = rho./max(rho(:));
rho = uint8(255*rho);

%%% Show the result
disp('Displaying results');
figure(1002);
surfl(-XYZ(:,:,1),-XYZ(:,:,2),-XYZ(:,:,3),[0 90]);
shading flat;
colormap gray
view(-220,30);
axis ij
axis image
title('Shape')

figure(1003)
imagesc(rho)
if(size(rho,3)==1)
	colormap gray
end
axis image
axis off
title('Albedo')

figure(1004)
Ndisp = N;
Ndisp(:,:,3) = -Ndisp(:,:,3);
Ndisp = 0.5*(Ndisp+1);
imagesc(Ndisp)
axis image
axis off
title('Normals')
drawnow

%%% Save to an obj file
if(do_write_obj)
	disp('Saving results');
	export_obj(XYZ,N,rho,data.mask,dataset);
end
