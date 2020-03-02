%% get geographical information from image label.
global msl_env_vars
localrootDir = msl_env_vars.local_pds_msl_imaging_rootDir;
pds_msl_imaging_URL = msl_env_vars.pds_msl_imaging_URL;

mastcam_rootpath = joinPath(localrootDir,pds_msl_imaging_URL);
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

cmmdl = [];
cmmdl.C = cmmdl_C;
cmmdl.A = cmmdl_A;
cmmdl.H = cmmdl_H;
cmmdl.V = cmmdl_V;

S_im = lbl.OBJECT_IMAGE.LINE_SAMPLES; L_im = lbl.OBJECT_IMAGE.LINES;

%%
%==========================================================================
% Get camera pointing vector for each camera image pixel (x_im, y_im)
% using camera CAHV model
%==========================================================================
[imxy_direc_rov] = get_3d_pointing_from_CAHV([L_im,S_im],cmmdl); 
% 3 x L_im x S_im

    
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
rov_rot_mat = get_rot_mat(rover_nav_coord.ROLL,rover_nav_coord.PITCH,rover_nav_coord.YAW);
rov_rot_mat_inv = get_rot_mat_inv(rover_nav_coord.ROLL,rover_nav_coord.PITCH,rover_nav_coord.YAW);

%%
%==========================================================================
% convert to the camera xy to the coordinate of ROVER_NAV with no rotation
%==========================================================================
% rov0 indicates ROVER_NAV with no rotaion (roll, pitch ,yaw) = (0,0,0) to
% the reference coordinate Site.
cmmdl_A_rov0 = rov_rot_mat * cmmdl_A';
cmmdl_C_rov0 = rov_rot_mat * cmmdl_C';
imxy_direc_rov0 = mmx('mult', rov_rot_mat, imxy_direc_rov);
% rov_rot_mat = gpuArray(rov_rot_mat); imxy_direc_rov = gpuArray(imxy_direc_rov);
% mxy_direc_rov0 = pagefun(@mtimes, rov_rot_mat, imxy_direc_rov);
[rov_rot_mat,imxy_direc_rov,imxy_direc_rov0] = gather(rov_rot_mat,imxy_direc_rov,imxy_direc_rov0);

cmmdl_C_geo = cmmdl_C_rov0 + [rover_nav_coord.NORTHING; 
                              rover_nav_coord.EASTING;
                              -rover_nav_coord.ELEVATION];

%%
%==========================================================================
% evaluation of DEM
%==========================================================================
basename_dem = 'MSL_Gale_DEM_Mosaic_1m_v3';
if ismac
    % dpath_dem = '/Volumes/LaCie/data/';
    dpath_dem = '/Users/yukiitoh/src/matlab/mastcam';
elseif isunix
    dpath_dem = '/Volume2/yuki/mastcam/';
end
[geo_im_FOV_mask] = get_geo_im_FOV(basename_dem,dpath_dem,rover_nav_coord,cmmdl,[L_im,S_im]);

%%
%==========================================================================
% EVALUATE where each camera pointing are intersecting.
%==========================================================================
% effort on reducing computational burden
% only search within the image field of view. Safeguards are taken for the
% image field of view, so it will be sufficient to take the field of view
% defined in the last section.
imxy_direc_rov0_2d = reshape(imxy_direc_rov0,[3,L_im*S_im]);
tic; [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = get_intersect_geo(cmmdl_C_geo,...
    imxy_direc_rov0_2d,basename_dem,dpath_dem,geo_im_FOV_mask,'gpu',1); toc;




