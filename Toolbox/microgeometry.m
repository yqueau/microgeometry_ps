function [XYZ,N,rho] = microgeometry(data,calib,params)
%NEAR_PS solves photometric stereo under nearby point light sources
%
%	=== USAGE ===	
%
%	XYZ = EXPORT_OBJ(DATA,CALIB) estimates an
%	NROWS x NCOLS x 3 gridded point cloud, using the images in DATA and
%	the calibration parameters in CALIB
%
%	XYZ = EXPORT_OBJ(DATA,CALIB,PARAMS) uses the algorithm parameters in
%	set in PARAMS
%
%	[XYZ,N] = EXPORT_OBJ(...) also provides an NROWS x NCOLS x 3 normal
%	map
%
%	[XYZ,N,rho] = EXPORT_OBJ(...) also provides an 
%	NROWS x NCOLS x 3 RGB albedo map
%
%	=== DATA ===
%	
%	- Required fields:
%		- DATA.I:	NROWS x NCOLS x NIMGS (graylevel images) mosaiced
%		RGB images, in the form RGGB
%	- Optional fields:
%		- DATA.mask:NROWS x NCOLS binary mask (default: ones(NROWS,NCOLS))
%
%	=== CALIB ===
%	- Required fields:
%		- CALIB.S:	NROWS x NCOLS x 3 x NIMGS light sources
%
%	- Optional fields:
%		- CALIB.K:	3 x 3 intrinsics matrix
%				(default: [], i.e. orthographic)
%
%	=== PARAMETERS ===
%	- Optional fields:
%		- PARAMS.z0: NROWS x NCOLS initial depth map
%			(default: 70.6383*ones(NROWS,NCOLS))
%		- PARAMS.ratio: integer downsampling factor 
%			(default: 1 -- should be 1 or a factor of 4)
%
%	=== CREDITS ===
%
%	This is inspired by the microgrometry PS method in :
%	[1]	"Microgeometry capture and RGB albedo estimation by photometric
%	stereo without demosaicing"
%		Quéau et al.
%		Proc. QCAV 2017, Tokyo, Japan
%
%	Data acquisition was made possible by a technological transfer
%	between the IRIT lab and the Pixience company, supervised by
%	Toulouse Tech Transfer

	disp(' ');
	
	if(nargin<2)
		disp('ERROR: must provide data and calibration');
		return;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Check DATA
	% Check images
	if(~isfield(data,'I'))
		disp('ERROR: images not provided in DATA.I')
		return;
	end
	I = double(data.I); clear data.I;
	[nrows,ncols,nimgs] = size(I);
	% Check mask
	if(~isfield(data,'mask'))
		disp('WARNING: mask not provided in DATA.mask, using full 2D grid')
		data.mask = ones(nrows,ncols);
	end
	mask = double(data.mask)>0; clear data.mask;
	if(size(mask,1)~=nrows | size(mask,2)~=ncols | ndims(mask)>2)
		disp(sprintf('ERROR: mask should be %d x %d',nrows,ncols));
		return;
	end
	clear data

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Check CALIB	
	% Check locations
	if(~isfield(calib,'S'))
		disp('ERROR: sources locations not provided in DATA.S')
		return;
	end
	S = double(calib.S); clear calib.S;
	if(size(S,1)~=nrows | size(S,2)~= ncols | size(S,3)~= 3 | size(S,4) ~= nimgs) 
		disp(sprintf('ERROR: S should be %d x %d x 3 x %d',nrows,ncols,nimgs));
		return;
	end
	% Check intrinsics
	if(~isfield(calib,'K'))
		disp('Warning: intrinsics not provided in CALIB.K - using orthographic projection')
		calib.K = [1,0,0.5*size(I,2);0,1,0.5*size(I,1);0,0,1];
		orthographic = 1;
	else
		orthographic = 0;
	end
	K = double(calib.K); clear calib.K;
	if(size(K,1)~=3 | size(K,2)~= 3 | ndims(K)>2) 
		disp('ERROR: K should be 3 x 3');
		return;
	end
	if(K(1,3)==0)
		K = K';
	end
	if(K(1,1)==0 | K(2,2)==0 | K(1,3)==0 | K(2,3)==0 | K(3,3)~=1 )
		disp('ERROR: intrinsics matrix not in a standard format');
		return;
	end
	clear calib

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Check PARAMS
	if(~exist('params','var')|isempty(params)) params=[]; end;
	% Check initial depth
	if(~isfield(params,'z0'))
		disp('WARNING: initial depth not provided in PARAMS.z0, using default values')
		params.z0 = 70.6383*ones(nrows,ncols);
	end
	z0 = params.z0; clear params.z0;
	if(size(z0,1)~=nrows | size(z0,2) ~= ncols | ndims(z0) > 2)
		disp(sprintf('ERROR: z0 should be %d x %d',nrows,ncols));
	end
	% Check ratio
	if(~isfield(params,'ratio'))
		disp('WARNING: downsampling factor not set in PARAMS.ratio, using default values')
		params.ratio = 1;
	end
	ratio = params.ratio; clear params.ratio;
	if(ratio < 0)
		disp('ERROR: ratio should be a positive integer');
		return;
	end
	if(ratio>1)
		K(1:2,:) = K(1:2,:)./ratio;
		mask = mask(1:ratio:end,1:ratio:end);
		z0 = z0(1:ratio:end,1:ratio:end);
		I = I(1:ratio:end,1:ratio:end,:);
		if(orthographic)
			K(1,1) = 1;
			K(2,2) = 1;
		end
		[nrows,ncols] = size(mask);
	end
	clear ratio
	
	clear params
	disp(' ');

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Prepare data
	% Intrinsics
	fx = K(1,1);
	fy = K(2,2);
	x0 = K(1,3);
	y0 = K(2,3);
	clear K

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Estimate normals and albedo pixelwise
	disp('Albedo and normal estimation')	
	indices_mask = find(mask>0);
	npix_masked = length(indices_mask);		
	if(nimgs>5)
		max_px = 2e6; % hack to avoid out of memory crash
	else
		max_px = 1e7;
	end
	N_ = zeros(npix_masked,3);
	rho = zeros(nrows,ncols);		
			
	for pxmin = 1:max_px:npix_masked
		i_submask = indices_mask(pxmin:min(npix_masked,pxmin+max_px-1));			
		II = [];
		JJ = [];
		KK = [];	
		Big_I = zeros(nimgs*length(i_submask),1);
		for i = 1:nimgs
			I_i = I(:,:,i);
			Big_I(i:nimgs:end-nimgs+i) = I_i(i_submask);

			II = [II;transpose(i:nimgs:nimgs*length(i_submask)-nimgs+i);transpose(i:nimgs:nimgs*length(i_submask)-nimgs+i);...
					 transpose(i:nimgs:nimgs*length(i_submask)-nimgs+i)];
			JJ = [JJ;transpose(1:3:3*length(i_submask)-2);...
					 transpose(2:3:3*length(i_submask)-1);...
					 transpose(3:3:3*length(i_submask))];
			S_i_1 = S(:,:,1,i);S_i_1 = S_i_1(i_submask);
			S_i_2 = S(:,:,2,i);S_i_2 = S_i_2(i_submask);
			S_i_3 = S(:,:,3,i);S_i_3 = S_i_3(i_submask);
			KK = [KK;S_i_1;S_i_2;S_i_3];
		end			
		clear I_i S_i_1 S_i_2 S_i_3 
		Big_S = sparse(II,JJ,KK);
		clear II JJ KK
		M_vect = (Big_S'*Big_S)\(Big_S'*Big_I);
		clear Big_S Big_I	
		rho(i_submask) = sqrt(M_vect(1:3:end-2).^2+M_vect(2:3:end-1).^2+M_vect(3:3:end).^2);	
		N_(pxmin:pxmin+length(i_submask)-1,1) = M_vect(1:3:end-2)./rho(i_submask);
		N_(pxmin:pxmin+length(i_submask)-1,2) = M_vect(2:3:end-1)./rho(i_submask);
		N_(pxmin:pxmin+length(i_submask)-1,3) = M_vect(3:3:end)./rho(i_submask);
	end	
	
	N = zeros(nrows,ncols, 3);
	N(:,:,3) = -1;
	for i=1:3
		Ni = N(:,:,i);
		Ni(indices_mask) = N_(:,i);
		N(:,:,i) = Ni;
	end
	indices_NaN = find(isnan(N(:,:,3)));
	for i=1:2
		Ni = N(:,:,i);
		Ni(indices_NaN) = 0;
		N(:,:,i) = Ni;
	end
	Ni = N(:,:,3);
	Ni(indices_NaN) = -1;
	N(:,:,3) = Ni;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Integrate normals into depth
	disp('Integration des normales')
	[u,v] = meshgrid(1:ncols,1:nrows);
	u = u-x0;
	v = v-y0;
	if(orthographic)
		r = -N(:,:,1)./N(:,:,3);
		s = -N(:,:,1)./N(:,:,3);
		r(mask==0) = 0;
		s(mask==0) = 0;
		r(N(:,:,3)==0) = 0;
		s(N(:,:,3)==0) = 0;
		z = simchony(s,r);

		lambda = sum(z(:).*z0(:))/sum(z(:).^2); % Choisi de facon a minimiser le residu quadratique \sum_\Omega (lambda z - Z_0)^2
		x = lambda*u;
		y = lambda*v;
		z = lambda*z;
		indices_NaN = find(mask==0);
		x(indices_NaN) = NaN;
		y(indices_NaN) = NaN;
		z(indices_NaN) = NaN;
		XYZ = zeros(size(x,1),size(x,2),3);
		XYZ(:,:,1)=x;
		XYZ(:,:,2)=y;
		XYZ(:,:,3)=z;
	else		
		focale = 0.5*(fx+fy);
		denom = (N(:,:,1).*u+N(:,:,2).*v+N(:,:,3)*focale);
		r=-N(:,:,1)./denom;
		s=-N(:,:,2)./denom;
		r((mask==0))=0;
		s((mask==0))=0;
		r((denom==0))=0;
		s((denom==0))=0;
		z = exp(simchony(s,r));

		lambda = sum(z(:).*z0(:))/sum(z(:).^2); % Choisi de facon a minimiser le residu quadratique \sum_\Omega (lambda z - Z_0)^2
		x = lambda*z.*u/focale;
		y = lambda*z.*v./focale;
		z = lambda*z;
		indices_NaN = find(mask==0);
		x(indices_NaN) = NaN;
		y(indices_NaN) = NaN;
		z(indices_NaN) = NaN;
		XYZ = zeros(size(x,1),size(x,2),3);
		XYZ(:,:,1)=x;
		XYZ(:,:,2)=y;
		XYZ(:,:,3)=z;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Albedo demosaicing
	disp('Interpolation des albedos')	
	rho_ = zeros(size(rho,1),size(rho,2),3);	
	rho_(:,:,1) = imresize(rho(1:2:end,1:2:end),size(rho));
	rho_(:,:,2) = 0.5*(imresize(rho(1:2:end,2:2:end),size(rho))+imresize(rho(2:2:end,1:2:end),size(rho)));
	rho_(:,:,3) = imresize(rho(2:2:end,2:2:end),size(rho));
	rho = rho_;
	for ch = 1:3
		rhoch = rho(:,:,ch);
		max_val = median(rhoch(indices_mask))+8*std(rhoch(indices_mask));
		rhoch(rhoch>max_val) = max_val;
		rhoch(rhoch<0) = 0;
		rho(:,:,ch) = rhoch;
	end

end

function [U,rmse] = simchony(p,q,gt)
% An implementation of the method from Simchony et Al for integration 
% of a normal field with natural boundary condition
% Ref : Direct Analytical Methods for Solving Poisson Equation in 
% Computer Vision problems - PAMI 1990
%
% example :
% p = ones(100,100); 
% q = ones(100,100); 
% u = simchony(p,q); 
%
% This performs the least square solution to \nabla u = [p,q], i.e. :
% min \int_\Omega \| \nablua U - [p,q] \|^2
% where \Omega is square and the natural Neumann boundary condition 
% \mu \cdotp (\nabla u -[p,q]) = 0 is used (2nd order derivatives 
% are neglected), where \mu is the outer 
% normal to \Omega. 
% This boundary condition is the one given by the 
% calculus of variations.
%
% % Axis : O->y
%        |
%        x
%
% The problem comes to solve a linear system Ax=b, where A is a bloc 
% Toplitz matrix. Fast solution is provided by Discrete Cosine Transform
%
% Implementation : Yvain Queau
% Universite de Toulouse, IRIT, UMR CNRS 5505
% yvain.queau@enseeiht.fr
% See other codes at http://ubee.enseeiht.fr/photometricstereo/

%~ p = 0.5*(p([2:end end],:)+p);
%~ q = 0.5*(q(:,[2:end end])+q);

% Compute div(p,q)
px = 0.5*(p([2:end end],:)-p([1 1:end-1],:));
qy = 0.5*(q(:,[2:end end])-q(:,[1 1:end-1]));

% Div(p,q) + Boundary Condition
f = px+qy;
f(1,2:end-1) = 0.5*(p(1,2:end-1)+p(2,2:end-1));
f(end,2:end-1) = 0.5*(-p(end,2:end-1)-p(end-1,2:end-1));
f(2:end-1,1) = 0.5*(q(2:end-1,1)+q(2:end-1,2));
f(2:end-1,end) = 0.5*(-q(2:end-1,end)-q(2:end-1,end-1));

f(1,1)=0.5*(p(1,1)+p(2,1)+q(1,1)+q(1,2));
f(end,1)=0.5*(-p(end,1)-p(end-1,1)+q(end,1)+q(end,2));
f(1,end)=0.5*(p(1,end)+p(2,end)-q(1,end)-q(1,end-1));
f(end,end)=0.5*(-p(end,end)-p(end-1,end)-q(end,end)-q(end,end-1));

% Sine transform of f
fsin=dct2(f);

% Denominator
[x,y] = meshgrid(0:size(p,2)-1,0:size(p,1)-1);
denom = (2*cos(pi*x/(size(p,2)))-2) + (2*cos(pi*y/(size(p,1))) - 2);
Z = fsin./(denom);
Z(1,1)=0.5*Z(1,2)+0.5*Z(2,1); %Or whatever...

% Inverse Sine transform :
U=idct2(Z);


	if(nargin>2)% Ground truth available
		moyenne_ecarts=mean(U(:)-gt(:));
		U=U-moyenne_ecarts;
		npix=size(p,1)*size(p,2);
		rmse=sqrt((sum((U(:)-gt(:)).^2))/npix);
	else
		U=U-min(U(:));
	end

end

