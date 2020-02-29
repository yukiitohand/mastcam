%% get geographical information from image label.
if ismac
    mastcam_rootpath = '/Volumes/LaCie/data/pds-imaging.jpl.nasa.gov/data/msl/';
elseif isunix
    mastcam_rootpath = '/Volume2/yuki/mastcam/data/pds-imaging.jpl.nasa.gov/data/msl/';
end
propMASTCAM = create_propMASTCAMbasename('SOL',25,'CAM_CODE','M[R]{1}',...
    'SEQ_ID',123,'PRODUCT_TYPE','[DE]{1}','DATA_PROC_CODE','DRLX',...
    'unique_CDPID',698);

[match_list] = mastcam_product_search(propMASTCAM);
site_id_list = cat(1,match_list.ROVER_MOTION_COUNTER_SITE);
drive_id_list = cat(1,match_list.ROVER_MOTION_COUNTER_DRIVE);
pose_id_list = cat(1,match_list.ROVER_MOTION_COUNTER_POSE);
rsm_mc_list = cat(1,match_list.ROVER_MOTION_COUNTER_RSM);

counter_list = [site_id_list drive_id_list pose_id_list rsm_mc_list];

counter_unique = unique(counter_list,'rows');

if size(counter_unique,1)>1
    fprintf('Several different perspectives are included\n');
    fprintf('(site, drive, pose, mc)\n');
    for i=1:size(counter_unique,1)
        match_idx = find(all(counter_list==counter_unique(i,:),2));
        fprintf('(% 4d, % 5d, % 4d, %2d): ', counter_unique(i,1),...
            counter_unique(i,2),counter_unique(i,3),counter_unique(i,4));
        for im = match_idx
            fprintf('%s,',match_list(im).PRODUCT_ID);
        end
        fprintf('\n');
            
    end
        
else
    fprintf('(site, drive, pose, mc)\n');
    fprintf('(% 4d, % 5d, % 4d, %2d):\n', counter_unique(1),...
            counter_unique(2),counter_unique(3),counter_unique(4));
    fprintf('Selected products:\n');
    for im = 1:length(match_list)
        fprintf('%s\n',match_list(im).PRODUCT_ID);
    end

    basename = match_list(1).PRODUCT_ID;
    subdir = joinPath(match_list(1).VOLUME_ID,match_list(1).PATH_NAME);
    dpath = joinPath(mastcam_rootpath,match_list(1).VOLUME_ID,match_list(1).PATH_NAME);
    
    fnamelist = dir(dpath);
    if isempty(extractMatchedBasename_v2(basename,[{fnamelist.name}],'exact',0))
        pds_msl_imaging_downloader(subdir,'basenameptrn',basename,'dwld',2);
    end
    
    lbl = pds3lblread(joinPath(dpath,[basename '.LBL']));
    hdr = mastcam_extract_imghdr_from_lbl(lbl);
    img_mastcam = envidataread_v2(joinPath(dpath,[basename '.IMG']),hdr);

    cmmdl_A_index = find(strcmpi('A',lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.MODEL_COMPONENT_ID));
    cmmdl_A = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.(sprintf('MODEL_COMPONENT_%1d',cmmdl_A_index));
    cmmdl_C_index = find(strcmpi('C', lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.MODEL_COMPONENT_ID));
    cmmdl_C = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.(sprintf('MODEL_COMPONENT_%1d',cmmdl_C_index));
    cmmdl_H_index = find(strcmpi('H', lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.MODEL_COMPONENT_ID));
    cmmdl_H = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.(sprintf('MODEL_COMPONENT_%1d',cmmdl_H_index));
    cmmdl_V_index = find(strcmpi('V', lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.MODEL_COMPONENT_ID));
    cmmdl_V = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.(sprintf('MODEL_COMPONENT_%1d',cmmdl_V_index));

    site_index = find(strcmpi('SITE',lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.COORDINATE_SYSTEM_INDEX_NAME));
    drive_index = find(strcmpi('DRIVE',lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.COORDINATE_SYSTEM_INDEX_NAME));
    pose_index = find(strcmpi('POSE',lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.COORDINATE_SYSTEM_INDEX_NAME));
    site_id = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.REFERENCE_COORD_SYSTEM_INDEX(site_index);
    drive_id = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.REFERENCE_COORD_SYSTEM_INDEX(drive_index);
    pose_id = lbl.GROUP_GEOMETRIC_CAMERA_MODEL_PARMS.REFERENCE_COORD_SYSTEM_INDEX(pose_index);
    % fprintf('(site,drive,pose) = (%d, %d, %d)\n', site_id, drive_id, pose_id);
end

im_size = [L_im, S_im];
cmmdl = [];
cmmdl.C = cmmdl_C;
cmmdl.A = cmmdl_A;
cmmdl.H = cmmdl_H;
cmmdl.V = cmmdl_V;

%%
%==========================================================================
% Get camera pointing vector for each camera image pixel (x_im, y_im)
% using camera CAHV model
%==========================================================================
% get a [x_im,y_im] in the image coordinate
S_im = lbl.OBJECT_IMAGE.LINE_SAMPLES; L_im = lbl.OBJECT_IMAGE.LINES;
[imxy_rov] = get_3d_pointing_from_CAHV([L_im,S_im],cmmdl);


% imx_im_1d = 1:S_im;
% imy_im_1d = reshape(1:L_im,[],1);
% [imx_im,imy_im] = meshgrid(imx_im_1d,imy_im_1d);
% 
% imxy_im_2d = permute(reshape(cat(2,imx_im,imy_im),[L_im*S_im,2]),[2,3,1]);
% 
% % compute pointing of each pixel in rover coordinate.
% % this is done by matrix inversion.
% HmXAt = cmmdl_H - cmmdl_A .* imxy_im_2d(1,:,:);
% VmYAt = cmmdl_V - cmmdl_V .* imxy_im_2d(2,:,:);
% M = cat(1,HmXAt,VmYAt,repmat([0 0 1],1,1,size(imxy_im_2d,3)));
% h = cat(1,mmx('mult',HmXAt,cmmdl_C'),mmx('mult',VmYAt,cmmdl_C'),repmat([1],1,1,size(imxy_im_2d,3)));
% % p_xy = pagefun(@ldivide,M,h);
% imxy_rov = permute(mmx('backslash',M,h),[1,3,2]);
% 
% % normalization
% imxy_rov = normalizevec(imxy_rov,1,'normtype',2);
% % check the direction by taking innerdots with camera axis.
% is_lookback = (cmmdl_A * imxy_rov)<0;
% imxy_rov(:,is_lookback) = -imxy_rov(:,is_lookback);
    
%%
%==========================================================================
% GET ROVER_NAV coordinate information from telemetry.csv
%==========================================================================
% rover frame is referenced from a Site frame:
%  +X: pointing north
%  +Y: east
%  +Z: nadir, pointing down 
% Note that elevation in the telemetry.csv is up positive

[rover_nav_coord] = get_rover_nav_coord(site_id,drive_id,pose_id);
rov_rot_mat = get_rot_mat(rover_nav_coord.ov_roll,rover_nav_coord.rov_pitch,rover_nav_coord.rov_yaw);

% % get telemetry data
% lbl_telemetry = pds3lblread(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS','telemetry.lbl'));
% telemetry = msl_telemetryCSVread(joinPath(mastcam_rootpath,'MSLPLC_1XXX/DATA/LOCALIZATIONS','telemetry.csv'),lbl_telemetry);
% site_telemetry = cat(1,telemetry.SITE);
% drive_telemetry = cat(1,telemetry.DRIVE);
% pose_telemetry = cat(1,telemetry.POSE);
% loc_telemetry = [site_telemetry drive_telemetry pose_telemetry];
% 
% % search matching (site,drive,pose) from the telemetry
% i_telemetry = find(all(loc_telemetry == [site_id, drive_id, pose_id],2));
% rover_nav_coord = telemetry(i_telemetry);
% 
% % get ROVER_NAV coordinate
% rov_plc_latitude = rover_nav_coord.PLANETOCENTRIC_LATITUDE;
% rov_longitude = rover_nav_coord.LONGITUDE;
% rov_northing = rover_nav_coord.NORTHING;
% rov_easting = rover_nav_coord.EASTING;
% rov_elevation = rover_nav_coord.ELEVATION;
% rov_roll = rover_nav_coord.ROLL;
% rov_pitch = rover_nav_coord.PITCH;
% rov_yaw = rover_nav_coord.YAW;
% 
% % compute rotation matrix
% rov_rot_mat = get_rot_mat(rov_roll,rov_pitch, rov_yaw);
% rov_rot_mat_inv = inv(rov_rot_mat);
% 
% % Get rover coordinate
% rover_nav = [];
% rover_nav.northing  = rov_northing; 
% rover_nav.easting   = rov_easting;
% rover_nav.elevation = rov_elevation;
% rover_nav.roll      = rov_roll;
% rover_nav.pitch     = rov_pitch;
% rover_nav.yaw       = rov_yaw;

%%
%==========================================================================
% convert to the camera xy to the coordinate of ROVER_NAV with no rotation
%==========================================================================
% rov0 indicates ROVER_NAV with no rotaion (roll, pitch ,yaw) = (0,0,0) to
% the reference coordinate Site.
cmmdl_A_rov0 = rov_rot_mat * cmmdl_A';
cmmdl_C_rov0 = rov_rot_mat * cmmdl_C';
imxy_rov0 = rov_rot_mat * imxy_rov;

%% 
%==========================================================================
% projection of ground reference coordinates to ROVER_NAV coordinate
%==========================================================================
geo_latitude_pc = DEdata.ddr.Latitude.img;
geo_longitude   = DEdata.ddr.Longitude.img;
geo_elevation   = DEdata.ddr.Elevation.img;
L_geo = TRRIFdata.hdr.lines; S_geo = TRRIFdata.hdr.samples;

% ----------------------------------------------------------
% planetocentric coordinate to northing easting coordinates
% ----------------------------------------------------------
Re = 3396190; % meters ellipsoid radius
geo_northing = Re .* (pi/180) .* geo_latitude_pc;
geo_easting  = Re .* (pi/180) .* geo_longitude;


% ------------------------------------------------------
% northing easting to rover coordinate with no rotation
% ------------------------------------------------------
% shift the coordinate
geo_rov0 = cat(3,geo_northing  - rov_northing,...
                 geo_easting   - rov_easting,...
                -geo_elevation - (-rov_elevation) );

geo_rov0 = permute(geo_rov0,[3,1,2]); % 3 x L_geo x S_geo

% ---------------------------------------
% rover coordinate (performing rotation)
% ---------------------------------------
geo_rov = mmx('mult',rov_rot_mat_inv,geo_rov0); % 3 x L_geo x S_geo
% geo_rov = reshape(geo_rov_2d,[3,L_geo,S_geo]);

%%
%==========================================================================
% EVALUATE camera field of view in the ground reference
%==========================================================================
right_dir = squeeze(mmx('mult',cmmdl_A,(geo_rov-cmmdl_C'))) > -1; % safeguard

geo_im = mmx('mult',[cmmdl_H;cmmdl_V],(geo_rov-cmmdl_C')) ./ mmx('mult',cmmdl_A,(geo_rov-cmmdl_C'));
geo_im = permute(geo_im,[2,3,1]);
geo_im_FOV = and(and(all(geo_im>-200,3),geo_im(:,:,1)<S_im+200),geo_im(:,:,2)<L_im+200);

geo_im_FOV_mask = and(right_dir,geo_im_FOV);
geo_im_FOV_mask_1nan = convertBoolTo1nan(geo_im_FOV_mask);

%%
% evaluation of DEM
basename_dem = 'MSL_Gale_DEM_Mosaic_1m_v3';
if ismac
    dpath_dem = '/Volumes/LaCie/data/';
elseif isunix
    dpath_dem = '/Volume2/yuki/mastcam/';
end
% hdr_dem = envihdrreadx(joinPath(dpath_dem,[basename_dem '.hdr']));
% geo_elevationl = lazyEnviReadl(joinPath(dpath_dem,[basename_dem '.img']),hdr_dem,1);

[geo_im_FOV_mask] = get_geo_im_FOV(basename_dem,dpath_dem,rover_nav_coord,cmmdl,im_size);


%%
%==========================================================================
% EVALUATE where each camera pointing are intersecting.
%==========================================================================
% effort on reducing computational burden
% only search within the image field of view. Safeguards are taken for the
% image field of view, so it will be sufficient to take the field of view
% defined in the last section.

for li_im = 1:L_im
    for si_im = 1:S_im
        % for each camera pointing vector, search the face that intersects
        % the geographical map.
        
        
    end
end



