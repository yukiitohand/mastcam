function [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = get_intersect_geo_forloop_v2(...
    cam_C_geo,im_size,basename_dem,dpath_dem,geo_im_FOV_mask,...
    basename_dem_imxy,dpath_dem_imxy,varargin)

is_gpu = false;
precision = 'double';
basename_dem = '';
dpath_dem    = '';
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

proc_page = batch_size>1;

L_im = im_size(1); S_im = im_size(2);

if is_gpu
    gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
end

if is_gpu
    cam_C_geo = gpuArray(cam_C_geo);
    geo_im_FOV_mask = gpuArray(geo_im_FOV_mask);
end

switch lower(precision)
    case 'double'
        cam_C_geo = double(cam_C_geo);
    case 'single'
        cam_C_geo = single(cam_C_geo);
end

% get image grid information
imx_im_1d = 1:S_im;
imy_im_1d = reshape(1:L_im,[],1);
[imx_im,imy_im] = meshgrid(imx_im_1d,imy_im_1d);
imxy_im_2d = permute(reshape(cat(2,imx_im,imy_im),[L_im*S_im,2]),[2,1,3]);

% camera center information
cam_C_geo = reshape(cam_C_geo,[],1);

% dem images
hdr_dem = envihdrreadx(joinPath(dpath_dem,[basename_dem '.hdr']));
demimg_path = joinPath(dpath_dem,[basename_dem '.img']);
L_geo = hdr_dem.lines; S_geo = hdr_dem.samples;
geo_northing   = reshape(hdr_dem.y,1,L_geo,1);
geol_easting    = reshape(hdr_dem.x,1,1,S_geo);

if is_gpu
    geo_northing = gpuArray(geo_northing);
    geol_easting = gpuArray(geol_easting);
end

switch lower(precision)
    case 'double'
        geo_northing = double(geo_northing);
        geol_easting = double(geol_easting);
    case 'single'
        geo_northing = single(geo_northing);
        geol_easting = single(geol_easting);
end


% dem imxy images
hdr_dem_imxy = envihdrreadx(joinPath(dpath_dem_imxy,[basename_dem_imxy '.hdr']));
dem_imxy_img_path = joinPath(dpath_dem_imxy,[basename_dem_imxy '.img']);


imxyz_geo = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_ref = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_range = inf(1,L_im*S_im,precision,gpu_varargin{:});

% if the whole line is outside of FOV, skip
% it is important to take the any in the dimension 1.
% last line is also removed because of the implementation we did.
valid_lines = any(geo_im_FOV_mask',1); valid_lines(L_geo) = false;
valid_lines = find(valid_lines);
if is_gpu
    valid_lines = gather(valid_lines);
end

geol = nan(3,2,S_geo,precision,gpu_varargin{:});
if is_gpu
    geol = gpuArray(geol);
end
geol(2,:,:) = repmat(geol_easting,[1,2,1]);

len_vl = length(valid_lines);
for li = 1:len_vl % l = 27270 
    tic;
    l = valid_lines(li);
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    % get geographical information from 
    % geol_northing   = repmat(geo_northing(:,[l l+1],:),[1,1,S_geo]);
    geol(1,:,:) = repmat(geo_northing(:,[l l+1],:),[1,1,S_geo]);
    geol_elevation  = lazyEnviReadl(demimg_path,hdr_dem,l);
    geolp1_elevation = lazyEnviReadl(demimg_path,hdr_dem,l+1);
    geol(3,:,:) = permute(cat(1, -geol_elevation, -geolp1_elevation),[3,1,2]);

    % geol = cat(1, geol_northing, repmat(geol_easting,[1,2,1]),...
    %     permute(cat(1, -geol_elevation, -geolp1_elevation),[3,1,2]) );
    % if is_gpu
    %     geol = gpuArray(geol);
    % end
    % switch lower(precision)
    %     case 'double'
    %         geol = double(geol);
    %     case 'single'
    %         geol = single(geol);
    % end
    
    % get supporting information from _imxy
    geol_imxy  = lazyEnviReadl(dem_imxy_img_path,hdr_dem_imxy,l-hdr_dem_imxy.line_offset);
    geolp1_imxy  = lazyEnviReadl(dem_imxy_img_path,hdr_dem_imxy,l-hdr_dem_imxy.line_offset+1);
    
    geol_imxy = permute( cat(3,geol_imxy,geolp1_imxy),[1,3,2]);

    valid_samples = find(geo_im_FOV_mask(l,:));
    if valid_samples(end) == S_geo
        valid_samples = valid_samples(1:end-1);
    end

    for s = valid_samples
        for j=1:2
            if j==1
                ppv1 = geol_imxy(:,1,s); % plane position vector in image space
                ppv2 = geol_imxy(:,1,s+1);
                ppv3 = geol_imxy(:,2,s);
                ppv1_geo = geol(:,1,s); % plane position vector
                ppv2_geo = geol(:,1,s+1);
                ppv3_geo = geol(:,2,s);
                % ppv_xyList = [geol_imxy(:,1,s) geol_imxy(:,1,s+1) geol_imxy(:,2,s)];
            elseif j==2
                ppv1 = geol_imxy(:,1,s+1);
                ppv2 = geol_imxy(:,2,s+1);
                ppv3 = geol_imxy(:,2,s);
                ppv1_geo = geol(:,1,s+1);
                ppv2_geo = geol(:,2,s+1);
                ppv3_geo = geol(:,2,s);
                % ppv_xyList = [geol_imxy(:,1,s+1) geol_imxy(:,2,s+1) geol_imxy(:,2,s)];
            end
            ppv_xyList = [ppv1 ppv2 ppv3];
            if any(isnan(ppv_xyList(:)))
                continue;
            end
            xy_min = ceil(min(ppv_xyList,[],2));
            xy_max = floor(max(ppv_xyList,[],2));

            xy_min = min(max(xy_min,[1;1]),[S_im;L_im]); xy_max = max(min(xy_max,[S_im;L_im]),[1;1]);
            x_list = xy_min(1):xy_max(1); y_list = xy_min(2):xy_max(2);
            idx_xy2d_list = L_im*(x_list-1) + y_list';
            idx_xy2d_list = idx_xy2d_list(:);

            % imxy_direc_geo_2d_focus = imxy_direc_geo_2d(:,idx_xy2d_list);

            % test line segment intersect with the plane determined by the
            % three points
            % [line_param,is_intersect] = line_plane_intersect_ldv(...
            %        cam_C_geo,imxy_direc_geo_2d_focus,ppv1,ppv2,ppv3,...
            %        is_gpu,proc_page);
            % if the line segment intersect with the plane, then test if the
            % intersection is within the triangle or not.
            % if is_intersect
            % pipv = cam_C_geo + imxy_direc_geo_2d_focus.*line_param; % plane intersection position vector
            % pipv = imxy_im_2d;
            [plane_param,is_in_face] = get_plane_param_coefficient(...
                ppv1,ppv2,ppv3,imxy_im_2d(:,idx_xy2d_list),precision,is_gpu,proc_page);
            
            imxyz_geo_s = ppv1_geo + ([ppv2_geo ppv3_geo] - ppv1_geo) *plane_param;
            imxyz_geo_range_s = sqrt(sum((imxyz_geo_s - cam_C_geo).^2,1));

            imls_update = find(and( is_in_face, imxyz_geo_range_s < imxyz_geo_range(:,idx_xy2d_list) ));
            if ~isempty(imls_update)
                idxes_update = idx_xy2d_list(imls_update);
                imxyz_geo(:,idxes_update) = imxyz_geo_s(:,imls_update);
                imxyz_geo_ref(:,idxes_update) = repmat([s;l;j],[1,length(imls_update)]);
                imxyz_geo_range(idxes_update) = imxyz_geo_range_s(1,imls_update);
            end
        end
    end
    toc;
end

if is_gpu
    [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = gather(imxyz_geo,imxyz_geo_ref,imxyz_geo_range);

end

imxyz_geo = permute(reshape(imxyz_geo,[3,L_im,S_im]),[2,3,1]);
imxyz_geo_ref = permute(reshape(imxyz_geo_ref,[3,L_im,S_im]),[2,3,1]);
imxyz_geo_range = reshape(imxyz_geo_range,[L_im,S_im]);

end
