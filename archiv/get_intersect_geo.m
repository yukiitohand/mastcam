function [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = get_intersect_geo(...
    cam_C_geo,imxy_direc_geo,basename_dem,dpath_dem,geo_im_FOV_mask,...
    basename_dem_imxy,dpath_dem_imxy,varargin)

is_gpu = false;
precision = 'double';
% proc_page = false;
batch_size = 10;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'GPU'
                is_gpu = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
            case 'BATCH_SIZE'
                batch_size = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

proc_page = batch_size>1;

if is_gpu
    gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
end

if is_gpu
    cam_C_geo = gpuArray(cam_C_geo);
    imxy_direc_geo = gpuArray(imxy_direc_geo);
    geo_im_FOV_mask = gpuArray(geo_im_FOV_mask);
end

switch lower(precision)
    case 'double'
        cam_C_geo = double(cam_C_geo);
        imxy_direc_geo = double(imxy_direc_geo);
        % geo_im_FOV_mask = double(geo_im_FOV_mask);
    case 'single'
        cam_C_geo = single(cam_C_geo);
        imxy_direc_geo = single(imxy_direc_geo);
        % geo_im_FOV_mask = single(geo_im_FOV_mask);
end

cam_C_geo = reshape(cam_C_geo,[],1);
[~,L_im,S_im] = size(imxy_direc_geo);
imxy_direc_geo_2d = reshape(imxy_direc_geo,[],L_im*S_im);

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

len_vl = length(valid_lines);
for li = 1:len_vl
    tic;
    l = valid_lines(li);
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    % get geographical information from 
    geol_northing   = repmat(geo_northing(:,[l l+1],:),[1,1,S_geo]);
    geol_elevation  = lazyEnviReadl(demimg_path,hdr_dem,l);
    geolp1_elevation = lazyEnviReadl(demimg_path,hdr_dem,l+1);

    geol = cat(1, geol_northing, repmat(geol_easting,[1,2,1]),...
        permute(cat(1, -geol_elevation, -geolp1_elevation),[3,1,2]) );
    if is_gpu
        geol = gpuArray(geol);
    end
    switch lower(precision)
        case 'double'
            geol = double(geol);
        case 'single'
            geol = single(geol);
    end
    
    % get supporting information from _imxy
    geol_imxy  = lazyEnviReadl(dem_imxy_img_path,hdr_dem_imxy,l-hdr_dem_imxy.line_offset);
    geolp1_imxy  = lazyEnviReadl(dem_imxy_img_path,hdr_dem_imxy,l-hdr_dem_imxy.line_offset+1);
    
    geol_imxy = reshape(permute(cat(3,geol_imxy,geolp1_imxy),[1,3,2]),[2,2*S_geo]);
    
    geol_2d = reshape(geol,[3,1,2*S_geo]);
    
    ppv1s = geol_2d(:,:,1:end-2);
    ppv2s = geol_2d(:,:,2:end-1);
    ppv3s = geol_2d(:,:,3:end);
% 
    valid_samples = find(geo_im_FOV_mask(l,:));
    % valid_samples = setdiff(valid_samples,[S_geo]); % remove the last index
    min_vs = max(valid_samples(1)-1,1); max_vs = valid_samples(end);
    if max_vs == valid_samples(end)==S_geo, max_vs = max_vs-1; end
    sidxes2 = (2*(min_vs)-1):(2*(max_vs));
    ppv1s = ppv1s(:,:,sidxes2);
    ppv2s = ppv2s(:,:,sidxes2);
    ppv3s = ppv3s(:,:,sidxes2);
    n_batch = ceil((length(sidxes2))/batch_size);

    for ni = 1:n_batch
        if ni~=n_batch
            sidx2 = (1+batch_size*(ni-1)):(batch_size*ni);
        else
            sidx2 = (1+batch_size*(ni-1)):length(sidxes2);
        end
        ppv1 = ppv1s(:,:,sidx2); ppv2 = ppv2s(:,:,sidx2); ppv3 = ppv3s(:,:,sidx2);
        ppv_imxy = geol_imxy(:,[sidxes2(sidx2) sidxes2(sidx2(end))+[1 2]]);
        % if sidx2(end)<(length(sidxes2)-1)
        %     ppv_imxy = geol_imxy(:,sidxes2([sidx2 sidx2(end)+[1 2]]));
        % elseif sidx2(end)==(length(sidxes2)-1)
        %     ppv_imxy = geol_imxy(:,[sidxes2([sidx2 sidx2(end)+1]) sidxes2(end)+1]);
        % elseif sidx2(end)==length(sidxes2)
        %     ppv_imxy = geol_imxy(:,[sidxes2(sidx2) sidxes2(end)+[1 2]]);
        % end
        % ppv_imxy = geol_imxy(:,sidxes2([sidx2 sidx2(end)+1 sidx2(end)+2]));
        
        xy_min = ceil(min(ppv_imxy,[],2));
        xy_max = floor(max(ppv_imxy,[],2));
        xy_min = min(max(xy_min,[1;1]),[S_im;L_im]); xy_max = max(min(xy_max,[S_im;L_im]),[1;1]);
        x_list = xy_min(1):xy_max(1); y_list = xy_min(2):xy_max(2);
        idx_xy2d_list = L_im*(x_list-1) + y_list';
        idx_xy2d_list = idx_xy2d_list(:);
        
        imxy_direc_geo_2d_focus = imxy_direc_geo_2d(:,idx_xy2d_list);
        
        % test line segment intersect with the plane determined by the
        % three points
        [line_param,is_intersect] = line_plane_intersect_ldv(...
                cam_C_geo,imxy_direc_geo_2d_focus,ppv1,ppv2,ppv3,...
                is_gpu,proc_page);
        % if the line segment intersect with the plane, then test if the
        % intersection is within the triangle or not.
        % if is_intersect
        pipv = cam_C_geo + imxy_direc_geo_2d_focus.*line_param; % plane intersection position vector
        [plane_param,is_in_face] = get_plane_param_coefficient(...
            ppv1,ppv2,ppv3,pipv,precision,is_gpu,proc_page);

        for j=1:length(sidx2)
            imls_update = find(and( is_in_face(:,:,j), line_param(:,:,j) < imxyz_geo_range(:,idx_xy2d_list) ));
            if ~isempty(imls_update)
                imxyz_geo(:,imls_update) = pipv(:,imls_update,j);
                s_ref = ceil(sidxes2(sidx2(j))/2); j_ref = sidxes2(sidx2(j)) - 2*(s_ref-1);
                imxyz_geo_ref(:,imls_update) = repmat([s_ref;l;j_ref],[1,length(imls_update)]);
                imxyz_geo_range(imls_update) = line_param(1,imls_update,j);
            end
        end
    end
    toc;
end

if is_gpu
    [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = gather(imxyz_geo,imxyz_geo_ref,imxyz_geo_range);

end


