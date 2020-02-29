function [rover_nav_coord] = get_rover_nav_coord(site_id,drive_id,pose_id)
% [rover_nav] = get_rover_nav_coord(site_id,drive_id,pose_id)
% Get rover coordinate from telemetry.csv for give site_id, drive_id,
% pose_id
% INPUTS:
%  site_id: integer, site id
%  drive_id: integer, drive id
%  pose_id: integer, pose id
% OUTPUTS
%  rover_nav_coord: struct
%    having fields: something like
%                       FRAME: 'ROVER'
%                        SITE: 3
%                       DRIVE: 372
%                        POSE: 8
%                   LANDING_X: -0.4610
%                   LANDING_Y: 31.8850
%                   LANDING_Z: 1.2700
%                    NORTHING: -2.7204e+05
%                     EASTING: 8.1468e+06
%     PLANETOCENTRIC_LATITUDE: -4.5895
%       PLANETODETIC_LATITUDE: -4.6438
%                   LONGITUDE: 137.4422
%                   ELEVATION: -4.5020e+03
%              MAP_PIXEL_LINE: 2.1108e+03
%            MAP_PIXEL_SAMPLE: 2.3237e+04
%              DEM_PIXEL_LINE: 528.6200
%            DEM_PIXEL_SAMPLE: 5.8098e+03
%                        ROLL: -1.3700
%                       PITCH: -1.6500
%                         YAW: 116.5000
%                        SCLK: 399622879
%                         SOL: 24
%
% rover frame is referenced from a Site frame:
%  +X: pointing north
%  +Y: east
%  +Z: nadir, pointing down 
% Note that elevation in the telemetry.csv is up positive

global msl_env_vars
localrootDir = msl_env_vars.local_pds_msl_imaging_rootDir;
pds_msl_imaging_URL = msl_env_vars.pds_msl_imaging_URL;

mastcam_rootpath = joinPath(localrootDir,pds_msl_imaging_URL);

% get telemetry data
lbl_telemetry = pds3lblread(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS','telemetry.lbl'));
telemetry = msl_telemetryCSVread(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS','telemetry.csv'),lbl_telemetry);
site_telemetry = cat(1,telemetry.SITE);
drive_telemetry = cat(1,telemetry.DRIVE);
pose_telemetry = cat(1,telemetry.POSE);
loc_telemetry = [site_telemetry drive_telemetry pose_telemetry];

% search matching (site,drive,pose) from the telemetry
i_telemetry = find(all(loc_telemetry == [site_id, drive_id, pose_id],2));
rover_nav_coord = telemetry(i_telemetry);

end

% get ROVER_NAV coordinate
% rov_plc_latitude = rover_nav_coord.PLANETOCENTRIC_LATITUDE;
% rov_longitude = rover_nav_coord.LONGITUDE;
% rov_northing = rover_nav_coord.NORTHING;
% rov_easting = rover_nav_coord.EASTING;
% rov_elevation = rover_nav_coord.ELEVATION;
% rov_roll = rover_nav_coord.ROLL;
% rov_pitch = rover_nav_coord.PITCH;
% rov_yaw = rover_nav_coord.YAW;

% compute rotation matrix
% rov_rot_mat = get_rot_mat(rov_roll,rov_pitch, rov_yaw);
% rov_rot_mat_inv = inv(rov_rot_mat);

% Get rover coordinate
% rover_nav = [];
% rover_nav.northing  = rov_northing; 
% rover_nav.easting   = rov_easting;
% rover_nav.elevation = rov_elevation;
% rover_nav.roll      = rov_roll;
% rover_nav.pitch     = rov_pitch;
% rover_nav.yaw       = rov_yaw;
% rover_nav.plc_latitude = rov_plc_latitude;
% rover_nav.longitude = rov_longitude;