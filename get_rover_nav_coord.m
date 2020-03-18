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

dpath_tempfile = '/Users/yukiitoh/src/matlab/toolbox/msl_toolbox/';

% get telemetry data
msl_telemtery_matfname = 'MSLPLC_telemetry.mat';
msl_telemtery_matfpath = joinPath(dpath_tempfile,msl_telemtery_matfname);
if ~exist(msl_telemtery_matfname)
    lbl_telemetry = pds3lblread(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS','telemetry.lbl'));
    telemetry = msl_telemetryCSVread(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS','telemetry.csv'),lbl_telemetry);
    MSL_telemetry_FRAME = cat(1,{telemetry.FRAME})';
    MSL_telemetry_FRAME = cellfun(@(x) sprintf('% 5s',x),MSL_telemetry_FRAME);
    MSL_telemetry_FRAME = cat(1,MSL_telemetry_FRAME{:});
    MSL_telemetry_SITE = cat(1,telemetry.SITE);
    MSL_telemetry_DRIVE = cat(1,telemetry.DRIVE);
    MSL_telemetry_POSE = cat(1,telemetry.POSE);
    MSL_telemetry_LANDING_X = cat(1,telemetry.LANDING_X);
    MSL_telemetry_LANDING_Y = cat(1,telemetry.LANDING_Y);
    MSL_telemetry_LANDING_Z = cat(1,telemetry.LANDING_Z);
    MSL_telemetry_NORTHING = cat(1,telemetry.NORTHING);
    MSL_telemetry_EASTING = cat(1,telemetry.EASTING);
    MSL_telemetry_PLANETOCENTRIC_LATITUDE = cat(1,telemetry.PLANETOCENTRIC_LATITUDE);
    MSL_telemetry_PLANETODETIC_LATITUDE = cat(1,telemetry.PLANETODETIC_LATITUDE);
    MSL_telemetry_LONGITUDE = cat(1,telemetry.LONGITUDE);
    MSL_telemetry_ELEVATION = cat(1,telemetry.ELEVATION);
    MSL_telemetry_MAP_PIXEL_LINE = cat(1,telemetry.MAP_PIXEL_LINE);
    MSL_telemetry_MAP_PIXEL_SAMPLE = cat(1,telemetry.MAP_PIXEL_SAMPLE);
    MSL_telemetry_DEM_PIXEL_LINE = cat(1,telemetry.DEM_PIXEL_LINE);
    MSL_telemetry_DEM_PIXEL_SAMPLE = cat(1,telemetry.MAP_PIXEL_SAMPLE);
    MSL_telemetry_ROLL = cat(1,telemetry.ROLL);
    MSL_telemetry_PITCH = cat(1,telemetry.PITCH);
    MSL_telemetry_YAW = cat(1,telemetry.YAW);
    MSL_telemetry_SCLK = cat(1,telemetry.SCLK);
    MSL_telemetry_SOL = cat(1,telemetry.SOL);
    
    save(msl_telemtery_matfpath,'MSL_telemetry_FRAME','MSL_telemetry_SITE',...
        'MSL_telemetry_DRIVE','MSL_telemetry_POSE','MSL_telemetry_LANDING_X',...
        'MSL_telemetry_LANDING_Y','MSL_telemetry_LANDING_Z',...
        'MSL_telemetry_NORTHING','MSL_telemetry_EASTING','MSL_telemetry_PLANETOCENTRIC_LATITUDE',...
        'MSL_telemetry_PLANETODETIC_LATITUDE','MSL_telemetry_LONGITUDE',...
        'MSL_telemetry_ELEVATION','MSL_telemetry_MAP_PIXEL_LINE','MSL_telemetry_MAP_PIXEL_SAMPLE',...
        'MSL_telemetry_DEM_PIXEL_LINE','MSL_telemetry_DEM_PIXEL_SAMPLE',...
        'MSL_telemetry_ROLL','MSL_telemetry_PITCH','MSL_telemetry_YAW',...
        'MSL_telemetry_SCLK','MSL_telemetry_SOL');
    
else
    load(msl_telemtery_matfpath,'MSL_telemetry_FRAME','MSL_telemetry_SITE',...
        'MSL_telemetry_DRIVE','MSL_telemetry_POSE','MSL_telemetry_LANDING_X',...
        'MSL_telemetry_LANDING_Y','MSL_telemetry_LANDING_Z',...
        'MSL_telemetry_NORTHING','MSL_telemetry_EASTING','MSL_telemetry_PLANETOCENTRIC_LATITUDE',...
        'MSL_telemetry_PLANETODETIC_LATITUDE','MSL_telemetry_LONGITUDE',...
        'MSL_telemetry_ELEVATION','MSL_telemetry_MAP_PIXEL_LINE','MSL_telemetry_MAP_PIXEL_SAMPLE',...
        'MSL_telemetry_DEM_PIXEL_LINE','MSL_telemetry_DEM_PIXEL_SAMPLE',...
        'MSL_telemetry_ROLL','MSL_telemetry_PITCH','MSL_telemetry_YAW',...
        'MSL_telemetry_SCLK','MSL_telemetry_SOL');
    
end


loc_telemetry = [MSL_telemetry_SITE MSL_telemetry_DRIVE MSL_telemetry_POSE];

% search matching (site,drive,pose) from the telemetry
it = find(all(loc_telemetry == [site_id, drive_id, pose_id],2));


MSL_telemetry_FRAME_sel = cellstr(MSL_telemetry_FRAME(it,:));
MSL_telemetry_FRAME_sel = strip(MSL_telemetry_FRAME_sel);

rover_nav_coord = struct('FRAME',MSL_telemetry_FRAME_sel,...
    'SITE',num2cell(MSL_telemetry_SITE(it)),...
    'DRIVE',num2cell(MSL_telemetry_DRIVE(it)),...
    'POSE',num2cell(MSL_telemetry_POSE(it)),...
    'LANDING_X',num2cell(MSL_telemetry_LANDING_X(it)),...
    'LANDING_Y',num2cell(MSL_telemetry_LANDING_Y(it)),...
    'LANDING_Z',num2cell(MSL_telemetry_LANDING_Z(it)),...
    'NORTHING',num2cell(MSL_telemetry_NORTHING(it)),...
    'EASTING',num2cell(MSL_telemetry_EASTING(it)),...
    'PLANETOCENTRIC_LATITUDE',num2cell(MSL_telemetry_PLANETOCENTRIC_LATITUDE(it)),...
    'PLANETODETIC_LATITUDE',num2cell(MSL_telemetry_PLANETODETIC_LATITUDE(it)),...
    'LONGITUDE',num2cell(MSL_telemetry_LONGITUDE(it)),...
    'ELEVATION',num2cell(MSL_telemetry_ELEVATION(it)),...
    'MAP_PIXEL_LINE',num2cell(MSL_telemetry_MAP_PIXEL_LINE(it)),...
    'MAP_PIXEL_SAMPLE',num2cell(MSL_telemetry_MAP_PIXEL_SAMPLE(it)),...
    'DEM_PIXEL_LINE',num2cell(MSL_telemetry_DEM_PIXEL_LINE(it)),...
    'DEM_PIXEL_SAMPLE',num2cell(MSL_telemetry_DEM_PIXEL_SAMPLE(it)),...
    'ROLL',num2cell(MSL_telemetry_ROLL(it)),...
    'PITCH',num2cell(MSL_telemetry_PITCH(it)),...
    'YAW',num2cell(MSL_telemetry_YAW(it)),...
    'SCLK',num2cell(MSL_telemetry_SCLK(it)),...
    'SOL',num2cell(MSL_telemetry_SOL(it)));

end