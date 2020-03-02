function [geo_im_FOV_mask] = get_geo_im_FOV(basename_dem,dpath_dem,rover_nav_coord,cmmdl,im_size,varargin)
% [geo_im_FOV_mask] = get_geo_im_FOV(basename_dem,dpath_dem,rover_nav,cmmdl,im_size)
%   evaluate FOV of an image on an ortho-georeferenced image using a
%   georeferenced DEM image.
%  INPUTS:
%    basename_dem: basename of DEM image
%    dpaht_dem   : directory path of DEM image
%    rover_nav_coord   
%                : struct containing information regarding the coordinate
%                  of ROVER_NAV frame, having fields below
%                    northing  : north positive
%                    easting   : east positive
%                    elevation : zenith positive
%                    roll
%                    pitch
%                    yaw
%    cmmdl       : struct containing camera model parameters for CAHV
%                  model in the ROVER_NAV frame.
%                  Fields are 'C','A','H','V'
%    im_size     : size of the image [L_im, S_im] for which FOV is
%                  evaluated.
%  OUTPUTS
%    geo_im_FOV_mask: boolean image, true if in the FOV, false otherwise.

is_gpu = false;
precision = 'double';
save_imxy = false;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'GPU'
                is_gpu = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
            case 'SAVE_IMXY'
                save_imxy = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if is_gpu
    gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
end

cmmdl_C = cmmdl.C';
cmmdl_A = cmmdl.A;
cmmdl_H = cmmdl.H;
cmmdl_V = cmmdl.V;
L_im = im_size(1);
S_im = im_size(2);

if is_gpu
    cmmdl_C = gpuArray(cmmdl_C);
    cmmdl_A = gpuArray(cmmdl_A);
    cmmdl_H = gpuArray(cmmdl_H);
    cmmdl_V = gpuArray(cmmdl_V);
end

cmmdl_HV_mat = cat(1,cmmdl_H,cmmdl_V);

% Get rover coordinate
rov_northing  = rover_nav_coord.NORTHING;
rov_easting   = rover_nav_coord.EASTING;
rov_elevation = rover_nav_coord.ELEVATION;
rov_roll      = rover_nav_coord.ROLL;
rov_pitch     = rover_nav_coord.PITCH;
rov_yaw       = rover_nav_coord.YAW;
rov_rot_mat = get_rot_mat(rov_roll,rov_pitch, rov_yaw);

hdr_dem = envihdrreadx(joinPath(dpath_dem,[basename_dem '.hdr']));
demimg_path = joinPath(dpath_dem,[basename_dem '.img']);
fid_demimg = fopen(demimg_path,'r');
L_geo = hdr_dem.lines; S_geo = hdr_dem.samples;
geo_northing   = reshape(hdr_dem.y,L_geo,1);
geo_easting    = reshape(hdr_dem.x,1,S_geo);



if is_gpu
    rov_rot_mat = gpuArray(rov_rot_mat);
    geo_northing = gpuArray(geo_northing);
    geo_easting = gpuArray(geo_easting);
end

switch lower(precision)
    case 'double'
        geo_northing = double(geo_northing);
        geo_easting = double(geo_easting);
    case 'single'
        geo_northing = single(geo_northing);
        geo_easting = single(geo_easting);
end
geo_rov0_northing = geo_northing - rov_northing;
geo_rov0_easting = geo_easting - rov_easting;

% Perform line by line operation
geo_im_FOV_mask = false(L_geo,S_geo,gpu_varargin{:});

% Open image to store coordinate in the image space.
if save_imxy
    hdr_out = hdr_dem;
    hdr_out.interleave = 'bil';
    hdr_out.bands = 2;
    hdr_out.band_names = {'imx','imy'};
    out_file_path = joinPath(dpath_dem, [basename_dem '_imxy.img']);
    envihdrwritex(hdr_out,joinPath(dpath_dem, [basename_dem '_imxy.hdr']));
end

geol_rov0 = zeros(3,S_geo);
geol_rov0(2,:) = geo_rov0_easting;
for l = 1:L_geo
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    % get
    % geol_northing   = repmat(geo_northing(l),[1,S_geo]);
    % geol_easting    = hdr_dem.x;
    geol_elevation = lazyEnviReadl(fid_demimg,hdr_dem,l);
    if is_gpu
        geol_elevation = gpuArray(geol_elevation);
    end
    
    % ---------------------------------------------------------------------
    % northing easting to rover coordinate with no rotation
    geol_rov0(1,:) = geo_rov0_northing(l);
    geol_rov0(3,:) = -geol_elevation - (-rov_elevation);
    % geol_rov0 = cat(1, repmat(geo_rov0_northing(l),[1,S_geo]),...
    %                    geo_rov0_easting,...
    %                   -geol_elevation - (-rov_elevation) ); % [3 x S_geo]
    
    % ---------------------------------------------------------------------
    % rover coordinate (performing rotation)
    geol_rov = rov_rot_mat \ geol_rov0;  % [3 x S_geo]
    
    %==========================================================================
    % Projection of ROVER_NAV to the image xy coordinate
    %==========================================================================
    geol_rov_m_cmmdl_C = geol_rov-cmmdl_C;
    right_dir = (cmmdl_A * geol_rov_m_cmmdl_C) > -1; % safeguard

    geol_im = (cmmdl_HV_mat * geol_rov_m_cmmdl_C) ./ (cmmdl_A * geol_rov_m_cmmdl_C); % 2 x S_geo
    geol_im_FOV = and(and(all(geol_im>-200,1),geol_im(1,:)<S_im+200),geol_im(2,:)<L_im+200);

    geol_im_FOV_mask = and(right_dir,geol_im_FOV);
    
    geo_im_FOV_mask(l,:) = geol_im_FOV_mask;
    
    if save_imxy
        if is_gpu
            geol_im = gather(geol_im);
        end
        lazyEnviWritel( out_file_path,single(geol_im),hdr_out,l,'a');
    end
    
end

fclose(fid_demimg);


end