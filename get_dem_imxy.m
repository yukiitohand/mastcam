function [dem_imxy,hdr_dem_imxy] = get_dem_imxy(basename_dem,dpath_dem,rover_nav_coord,cmmdl,dem_im_FOV_mask,varargin)
% [dem_imxy,hdr_dem_imxy] = get_dem_imxy(basename_dem,dpath_dem,rover_nav_coord,cmmdl,geo_im_FOV_mask,varargin)
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
%    dem_im_FOV_mask: boolean matrix, true if a pixel is in the FOV, zero
%                     other wise.
%   
%  OUTPUTS
%    dem_imxy: [len_vl x lenvs x 2] image storing the image (x,y)
%              coordinate for each pixel of dem image. The image is
%              cropped from the original dem image. The information for
%              offset samples and lines are stored in hdr_dem_geo_imxy
%              information.
%    hdr_dem_imxy: header information, storing line_offset and
%    sample_offset

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
% L_im = im_size(1);
% S_im = im_size(2);

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

valid_lines = find(any(dem_im_FOV_mask',1)); 
len_vl = length(valid_lines);

valid_samples = find(any(dem_im_FOV_mask,1));
len_vs = length(valid_samples);

deml_rov0 = zeros(3,len_vs,precision,gpu_varargin{:});
deml_rov0(2,:) = geo_rov0_easting(valid_samples);

dem_imxy = nan(len_vl,len_vs,2);

hdr_dem_imxy = hdr_dem;
hdr_dem_imxy.interleave = 'bil';
hdr_dem_imxy.bands = 2;
hdr_dem_imxy.lines = len_vl;
hdr_dem_imxy.samples = len_vs;
hdr_dem_imxy.band_names = {'imx','imy'};
hdr_dem_imxy.line_offset = valid_lines(1)-1;
hdr_dem_imxy.sample_offset = valid_samples(1)-1;

dem_im_FOV_mask_vls = dem_im_FOV_mask(valid_lines,valid_samples);
    
if save_imxy
    out_file_path = joinPath(dpath_dem, [basename_dem '_imxy2.img']);
    envihdrwritex(hdr_dem_imxy,joinPath(dpath_dem, [basename_dem '_imxy2.hdr']));
end

for li = 1:len_vl
    l = valid_lines(li);
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    geol_elevation = lazyEnviReadl(fid_demimg,hdr_dem,l);
    if is_gpu
        geol_elevation = gpuArray(geol_elevation);
    end
    
    % ---------------------------------------------------------------------
    % northing easting to rover coordinate with no rotation
    deml_rov0(1,:) = geo_rov0_northing(l);
    deml_rov0(3,:) = -geol_elevation(valid_samples) - (-rov_elevation);
    
    % ---------------------------------------------------------------------
    % rover coordinate (performing rotation)
    deml_rov = rov_rot_mat \ deml_rov0;  % [3 x S_geo]
    
    %==========================================================================
    % Projection of ROVER_NAV to the image xy coordinate
    %==========================================================================
    deml_rov_m_cmmdl_C = deml_rov-cmmdl_C;
    deml_im = (cmmdl_HV_mat * deml_rov_m_cmmdl_C) ./ (cmmdl_A * deml_rov_m_cmmdl_C); % 2 x S_geo
    
    deml_im(1,dem_im_FOV_mask_vls(li,:)==0) = nan;
    deml_im(2,dem_im_FOV_mask_vls(li,:)==0) = nan;
    dem_imxy(li,:,:) = permute(deml_im, [3,2,1]);
    
    if save_imxy
        if is_gpu
            deml_im = gather(deml_im);
        end
        lazyEnviWritel( out_file_path,single(deml_im),hdr_dem_imxy,l,'a');
    end
    
end

fclose(fid_demimg);


end