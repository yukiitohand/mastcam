function [geo_im_FOV_mask] = get_geo_im_FOV(basename_dem,dpath_dem,rover_nav_coord,cmmdl,im_size)
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


cmmdl_C = cmmdl.C;
cmmdl_A = cmmdl.A;
cmmdl_H = cmmdl.H;
cmmdl_V = cmmdl.V;
L_im = im_size(1);
S_im = im_size(2);

% Get rover coordinate
rov_northing  = rover_nav_coord.NORTHING;
rov_easting   = rover_nav_coord.EASTING;
rov_elevation = rover_nav_coord.ELEVATION;
rov_roll      = rover_nav_coord.ROLL;
rov_pitch     = rover_nav_coord.PITCH;
rov_yaw       = rover_nav_coord.YAW;
rov_rot_mat = get_rot_mat(rov_roll,rov_pitch, rov_yaw);

hdr_dem = envihdrreadx(joinPath(dpath_dem,[basename_dem '.hdr']));
L_geo = hdr_dem.lines; S_geo = hdr_dem.samples;

% Perform line by line operation
geo_im_FOV_mask = false(L_geo,S_geo);

for l = 1:L_geo
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    % get
    geol_northing   = repmat(hdr_dem.y(l),[1,S_geo]);
    geol_easting    = hdr_dem.x;
    geol_elevation = lazyEnviReadl(joinPath(dpath_dem,[basename_dem '.img']),hdr_dem,l);
    
    % ---------------------------------------------------------------------
    % northing easting to rover coordinate with no rotation
    geol_rov0 = cat(1, geol_northing - rov_northing,...
                       geol_easting - rov_easting,...
                      -geol_elevation - (-rov_elevation) ); % [3 x S_geo]
    
    % ---------------------------------------------------------------------
    % rover coordinate (performing rotation)
    geol_rov = rov_rot_mat \ geol_rov0;  % [3 x S_geo]
    
    %==========================================================================
    % Projection of ROVER_NAV to the image xy coordinate
    %==========================================================================
    right_dir = (cmmdl_A * (geol_rov-cmmdl_C')) > -1; % safeguard

    geol_im = [cmmdl_H;cmmdl_V] * (geol_rov-cmmdl_C') ./ (cmmdl_A * (geol_rov-cmmdl_C')); % 2 x S_geo
    geol_im_FOV = and(and(all(geol_im>-200,1),geol_im(1,:)<S_im+200),geol_im(2,:)<L_im+200);

    geol_im_FOV_mask = and(right_dir,geol_im_FOV);
    
    geo_im_FOV_mask(l,:) = geol_im_FOV_mask;
end

end