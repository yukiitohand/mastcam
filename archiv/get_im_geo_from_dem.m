function [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = get_im_geo_from_dem(...
    cam_C_geo,im_size,basename_dem,dpath_dem,dem_imFOV_mask,...
    dem_imxy,hdr_dem_imxy,varargin)

is_gpu = false;
precision = 'double';
proc_page = 0;
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

% proc_page = batch_size>1;

L_im = im_size(1); S_im = im_size(2);

if is_gpu
    gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
end

if is_gpu
    cam_C_geo = gpuArray(cam_C_geo);
    dem_imFOV_mask = gpuArray(dem_imFOV_mask);
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
dem_northing   = reshape(hdr_dem.y,1,L_geo,1);
deml_easting    = reshape(hdr_dem.x,1,1,S_geo);

if is_gpu
    dem_northing = gpuArray(dem_northing);
    deml_easting = gpuArray(deml_easting);
end

switch lower(precision)
    case 'double'
        dem_northing = double(dem_northing);
        deml_easting = double(deml_easting);
    case 'single'
        dem_northing = single(dem_northing);
        deml_easting = single(deml_easting);
end

imxyz_geo = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_ref = nan(3,L_im*S_im,precision,gpu_varargin{:});
imxyz_geo_range = inf(1,L_im*S_im,precision,gpu_varargin{:});

% if the whole line is outside of FOV, skip
% it is important to take the any in the dimension 1.
% last line is also removed because of the implementation we did.
valid_lines = any(dem_imFOV_mask',1);
valid_lines = find(valid_lines);
valid_lines = valid_lines(1:end-1);
if is_gpu
    valid_lines = gather(valid_lines);
end

deml_geo = nan(3,2,S_geo,precision,gpu_varargin{:});
if is_gpu
    deml_geo = gpuArray(deml_geo);
end
deml_geo(2,:,:) = repmat(deml_easting,[1,2,1]);

len_vl = length(valid_lines);
for li = 1:len_vl % l = 27270 
    %tic;
    l = valid_lines(li);
    %==========================================================================
    % projection of ground reference coordinates to ROVER_NAV coordinate
    %==========================================================================
    % get geographical information from 
    % geol_northing   = repmat(geo_northing(:,[l l+1],:),[1,1,S_geo]);
    deml_geo(1,:,:) = repmat(dem_northing(:,[l l+1],:),[1,1,S_geo]);
    deml_elevation  = lazyEnviReadl(demimg_path,hdr_dem,l);
    demlp1_elevation = lazyEnviReadl(demimg_path,hdr_dem,l+1);
    deml_geo(3,:,:) = permute(cat(1, -deml_elevation, -demlp1_elevation),[3,1,2]);
    
    deml_imxy = permute( dem_imxy([l l+1]-hdr_dem_imxy.line_offset,:,:),[3,1,2]);

    valid_samples = find(dem_imFOV_mask(l,:));
    if valid_samples(end) == (hdr_dem_imxy.samples+hdr_dem_imxy.sample_offset)
        valid_samples = valid_samples(1:end-1);
    end

    for s = valid_samples
        s_dem_imxy = s - hdr_dem_imxy.sample_offset;
        for j=1:2
            if j==1
                ppv1 = deml_imxy(:,1,s_dem_imxy); % plane position vector in image space
                ppv2 = deml_imxy(:,1,s_dem_imxy+1);
                ppv3 = deml_imxy(:,2,s_dem_imxy);
                ppv1_geo = deml_geo(:,1,s); % plane position vector
                ppv2_geo = deml_geo(:,1,s+1);
                ppv3_geo = deml_geo(:,2,s);
            elseif j==2
                ppv1 = deml_imxy(:,1,s_dem_imxy+1);
                ppv2 = deml_imxy(:,2,s_dem_imxy+1);
                ppv3 = deml_imxy(:,2,s_dem_imxy);
                ppv1_geo = deml_geo(:,1,s+1);
                ppv2_geo = deml_geo(:,2,s+1);
                ppv3_geo = deml_geo(:,2,s);
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

            % test if the intersection is within the triangle or not.
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
    % toc;
end

if is_gpu
    [imxyz_geo,imxyz_geo_ref,imxyz_geo_range] = gather(imxyz_geo,imxyz_geo_ref,imxyz_geo_range);

end

imxyz_geo = permute(reshape(imxyz_geo,[3,L_im,S_im]),[2,3,1]);
imxyz_geo_ref = permute(reshape(imxyz_geo_ref,[3,L_im,S_im]),[2,3,1]);
imxyz_geo_range = reshape(imxyz_geo_range,[L_im,S_im]);

end
