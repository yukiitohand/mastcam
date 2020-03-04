function [dem_imFOV_mask] = get_dem_imFOV(basename_dem,dpath_dem,rover_nav_coord,cmmdl,im_size,varargin)
% [dem_im_FOV_mask] = get_dem_imFOV(basename_dem,dpath_dem,rover_nav_coord,cmmdl,im_size,varargin)
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
%    dem_imFOV_mask: boolean image, true if in the FOV, false otherwise.

is_gpu = false;
precision = 'double';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'GPU'
                is_gpu = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
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
L_dem = hdr_dem.lines; S_dem = hdr_dem.samples;
dem_geo_northing   = reshape(hdr_dem.y,L_dem,1);
dem_geo_easting    = reshape(hdr_dem.x,1,S_dem);



if is_gpu
    rov_rot_mat = gpuArray(rov_rot_mat);
    dem_geo_northing = gpuArray(dem_geo_northing);
    dem_geo_easting = gpuArray(dem_geo_easting);
end

switch lower(precision)
    case 'double'
        dem_geo_northing = double(dem_geo_northing);
        dem_geo_easting = double(dem_geo_easting);
    case 'single'
        dem_geo_northing = single(dem_geo_northing);
        dem_geo_easting = single(dem_geo_easting);
end
dem_rov0_northing = dem_geo_northing - rov_northing;
dem_rov0_easting = dem_geo_easting - rov_easting;

% Perform line by line operation
dem_imFOV_mask = false(L_dem,S_dem,gpu_varargin{:});

deml_rov0 = zeros(3,S_dem,precision,gpu_varargin{:});
deml_rov0(2,:) = dem_rov0_easting;
for l = 1:L_dem
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    deml_elevation = lazyEnviReadl(fid_demimg,hdr_dem,l);
    if is_gpu
        deml_elevation = gpuArray(deml_elevation);
    end
    
    % ---------------------------------------------------------------------
    % northing easting to rover coordinate with no rotation
    deml_rov0(1,:) = dem_rov0_northing(l);
    deml_rov0(3,:) = -deml_elevation - (-rov_elevation);
    
    % ---------------------------------------------------------------------
    % rover coordinate (performing rotation)
    deml_rov = rov_rot_mat \ deml_rov0;  % [3 x S_geo]
    
    %==========================================================================
    % Projection of ROVER_NAV to the image xy coordinate
    %==========================================================================
    deml_rov_m_cmmdl_C = deml_rov-cmmdl_C;
    right_dir = (cmmdl_A * deml_rov_m_cmmdl_C) > -1; % safeguard

    deml_im = (cmmdl_HV_mat * deml_rov_m_cmmdl_C) ./ (cmmdl_A * deml_rov_m_cmmdl_C); % 2 x S_geo
    deml_imFOV = and(and(all(deml_im>-200,1),deml_im(1,:)<S_im+200),deml_im(2,:)<L_im+200);

    deml_imFOV_mask = and(right_dir,deml_imFOV);
    
    dem_imFOV_mask(l,:) = deml_imFOV_mask;
    
end

fclose(fid_demimg);


end