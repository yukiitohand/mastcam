function [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = get_intersect_geo(...
    cam_C_geo,imxy_direc_geo_2d,basename_dem,dpath_dem,geo_im_FOV_mask,varargin)

is_gpu = false;
precision = 'double';
proc_page = false;
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
            case 'PAGE'
                proc_page = varargin{i+1};
            case 'BATCH_SIZE'
                batch_size = varargin{i+1};
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

if is_gpu
    cam_C_geo = gpuArray(cam_C_geo);
    imxy_direc_geo_2d = gpuArray(imxy_direc_geo_2d);
    geo_im_FOV_mask = gpuArray(geo_im_FOV_mask);
end

switch lower(precision)
    case 'double'
        cam_C_geo = double(cam_C_geo);
        imxy_direc_geo_2d = double(imxy_direc_geo_2d);
        % geo_im_FOV_mask = double(geo_im_FOV_mask);
    case 'single'
        cam_C_geo = single(cam_C_geo);
        imxy_direc_geo_2d = single(imxy_direc_geo_2d);
        % geo_im_FOV_mask = single(geo_im_FOV_mask);
end

cam_C_geo = reshape(cam_C_geo,[],1);
[~,LS_im] = size(imxy_direc_geo_2d);

hdr_dem = envihdrreadx(joinPath(dpath_dem,[basename_dem '.hdr']));
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

imxyz_geo = nan(3,LS_im,precision,gpu_varargin{:});
imxyz_geo_ref = nan(3,LS_im,precision,gpu_varargin{:});
imxyz_geo_range = inf(1,LS_im,precision,gpu_varargin{:});

% if the whole line is outside of FOV, skip
% it is important to take the any in the dimension 1.
% last line is also removed because of the implementation we did.
valid_lines = any(geo_im_FOV_mask',1); valid_lines(L_geo) = false;
valid_lines = find(valid_lines);
if is_gpu
    valid_lines = gather(valid_lines);
end

len_vl = length(valid_lines);
% counter = 0;
for li = 1:len_vl
    tic;
%     if floor(li*100/len_vl/5) > counter
%         fprintf('*');
%         counter = counter + 1;
%     end
    l = valid_lines(li);
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    % get geographical information from 
    geol_northing   = repmat(geo_northing(:,[l l+1],:),[1,1,S_geo]);
    geol_elevation  = lazyEnviReadl(joinPath(dpath_dem,[basename_dem '.img']),hdr_dem,l);
    geolp1_elevation = lazyEnviReadl(joinPath(dpath_dem,[basename_dem '.img']),hdr_dem,l+1);

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
    
    geol_2d = reshape(geol,[3,1,2*S_geo]);
    
    ppv1s = geol_2d(:,:,1:end-2);
    ppv2s = geol_2d(:,:,2:end-1);
    ppv3s = geol_2d(:,:,3:end);
% 
    valid_samples = find(geo_im_FOV_mask(l,:));
    % valid_samples = setdiff(valid_samples,[S_geo]); % remove the last index
    min_vs = valid_samples(1); max_vs = valid_samples(end);
    if max_vs == valid_samples(end)==S_geo, max_vs = max_vs-1; end
    sidxes2 = (2*min_vs-1):(2*max_vs);
    ppv1s = ppv1s(:,:,sidxes2);
    ppv2s = ppv2s(:,:,sidxes2);
    ppv3s = ppv3s(:,:,sidxes2);
    n_batch = ceil(length(sidxes2)/batch_size);
%     
%     [line_param,is_intersect] = line_plane_intersect_ldv(...
%                 cam_C_geo,imxy_direc_geo_2d,ppv1,ppv2,ppv3);
%     % if the line segment intersect with the plane, then test if the
%     % intersection is within the triangle or not.
%     % if is_intersect
%     pipv = cam_C_geo + imxy_direc_geo_2d.*line_param; % plane intersection position vector
%     [plane_param,is_in_face] = get_plane_param_coefficient(ppv1,ppv2,ppv3,pipv);
%     % else
%     %     is_in_face = false;
%     % end
%     line_param = line_param .* convertBoolTo1nan(is_in_face) .* convertBoolTo1nan(is_intersect);
%     [line_param_1d,intersect_s] = min(line_param,[],3);
%     intersect_s(isnan(line_param_1d)) = nan;
%     
%     imls_update = find(and(~isnan(intersect_s),line_param_1d<imxyz_geo_range));
%     
%     for imi_ls = imls_update
%         imxyz_geo(:,imi_ls) = pipv(:,imi_ls,intersect_s(imi_ls));
%         s_ref = ceil(intersect_s(imi_ls)/2); j_ref = intersect_s(imi_ls) - 2*(s_ref-1);
%         imxyz_geo_ref(:,imi_ls) = [s_ref,l,j_ref];
%         imxyz_geo_range(imi_ls) = line_param_1d(imi_ls);
%     end

    % valid_samples = find(geo_im_FOV_mask(l,:));
    % valid_samples = setdiff(valid_samples,[S_geo]); % remove the last index
    % min_vs = valid_samples(1); max_vs = valid_samples(end);
    % if max_vs == valid_samples(end)==S_geo, max_vs = max_vs-1; end
    % sidxes = min_vs:max_vs;
    for ni = 1:n_batch
        if ni~=n_batch
            sidx2 = (1+batch_size*(ni-1)):(batch_size*ni);
        else
            sidx2 = (1+batch_size*(ni-1)):length(sidxes2);
        end
        ppv1 = ppv1s(:,:,sidx2); ppv2 = ppv2s(:,:,sidx2); ppv3 = ppv3s(:,:,sidx2);
        % select three points
        % ppv11 = geol(:,1,s); % plane position vector
        % ppv12 = geol(:,1,s+1);
        % ppv13 = geol(:,2,s);
        % idx_vert = [l,s;l,s+1;l+1,s];
        % ppv21 = geol(:,1,s+1);
        % ppv22 = geol(:,2,s+1);
        % ppv23 = geol(:,2,s);
        % ppv1 = cat(3,ppv11,ppv21); ppv2 = cat(3,ppv12,ppv22); ppv3 = cat(3,ppv13,ppv23);
        % idx_vert = [l,s;l+1,s;l+1,s-1];
        % test line segment intersect with the plane determined by the
        % three points
        [line_param,is_intersect] = line_plane_intersect_ldv(...
                cam_C_geo,imxy_direc_geo_2d,ppv1,ppv2,ppv3,...
                is_gpu,proc_page);
        % if the line segment intersect with the plane, then test if the
        % intersection is within the triangle or not.
        % if is_intersect
        pipv = cam_C_geo + imxy_direc_geo_2d.*line_param; % plane intersection position vector
        [plane_param,is_in_face] = get_plane_param_coefficient(...
            ppv1,ppv2,ppv3,pipv,precision,is_gpu,proc_page);

        for j=1:length(sidx2)
            imls_update = find(and(is_intersect(:,:,j), ...
                and( is_in_face(:,:,j), line_param(:,:,j) < imxyz_geo_range )));
            imxyz_geo(:,imls_update) = pipv(:,imls_update,j);
            s_ref = ceil(sidxes2(sidx2(j))/2); j_ref = sidxes2(sidx2(j)) - 2*(s_ref-1);
            imxyz_geo_ref(:,imls_update) = repmat([s_ref;l;j_ref],[1,length(imls_update)]);
            imxyz_geo_range(imls_update) = line_param(imls_update,j);
        end
    end
    toc;
end

if is_gpu
    [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = gather(imxyz_geo,imxyz_geo_ref,imxyz_geo_range);

end


