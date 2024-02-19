% Code assumes follwoing folder structure:
%          home_dir
%             - life_life[patient num]...
%                - [scan num].[name]_[mode]_BH_[scan type]Ao
%                   - .tgz file
%                - [scan num].[name]_[mode]_BH_[scan type]Ao
%                   - .tgz file
%             ...
%             - results
%                - bssfp_masks
%                   - [patient num]-[scan num]_bssfp_[anatomy].mat
%                - pwv_masks
%                   - [patient num]-[pwv scan num]_bssfp_[anatomy].mat
%             ...
% Change home_dir to parent folder containing all patient folders
% Change folder to specfic patient folder
% Change result_dir to home_dir plus results (use double slashes cause
% sprintf is weird)
% Change num to last 2 digits of scan starting number and pwvnum to
% corresponding pwv scan, i.e., '00008.FIESTA...' -> num = 08 and 
% '00007.PWV...' -> pwvnum = 07
% Change anat to 'A' for ascending aorta or 'Abd' for abdominal

%% Loading bssfp images
home_dir = 'C:\Users\naren\Desktop\Aorta_Segmentation';
folder = 'life_life00016_11110_2020-03-09-h08';
result_dir = 'C:\\Users\\naren\\Desktop\\Aorta_Segmentation\\results\\';
num = '08';
pwvnum = '07';
anat = 'A';

scan = sprintf('000%s.FIESTA_Cine_BH_%sAo', num, anat);
patient = folder(13:14);
id = [patient '-' num];
if strcmp(anat, 'A')
    z = 2;
    asc_name = [result_dir 'bssfp_masks\\%s_bssfp_AscAo.mat'];
    AscAof = sprintf(asc_name, id);
    desc_name  = [result_dir 'bssfp_masks\\%s_bssfp_DescAo.mat'];
    DescAof = sprintf(desc_name, id);   
    if exist(AscAof, 'file') == 0
        AscAo = zeros(256,256,length(dirFiles));
    else
        load(AscAof);
    end
    if exist(DescAof, 'file') == 0
        DescAo = zeros(256,256,length(dirFiles));
    else
        load(DescAof);
    end
    % adjustable parameters
    crop = [155 410 101 356];
    radrange = [12 30];
    sens = 0.90;
else
    z = 1;
    abd_name  = [result_dir 'bssfp_masks\\%s_bssfp_AbdAo.mat'];
    AbdAof = sprintf(abd_name, id);
    if exist(AbdAof, 'file') == 0
        AscAo = zeros(256,256,length(dirFiles));
    else
        load(AbdAof);
    end
    % adjustable parameters
    crop = [105 360 121 376];
    radrange = [12 23];
    sens = 0.9;
end

img_name  = [result_dir 'bssfp_masks\\%s_bssfp_images.mat'];
bssfpf = sprintf(img_name, id);
bssfp = zeros(256,256,length(dirFiles));
lifeDir = [pwd '\' folder '\' scan];
cd(lifeDir);
untar([scan '.tgz'])
dirFiles = dir('*.dcm');
for i=1:length(dirFiles)
    image = dicomread(dirFiles(i).name);
    image = image(crop(1):crop(2),crop(3):crop(4));
    
    %% Find circles (Hough transform)
    image = imadjust(image);
    [centers,radii,~] = imfindcircles(image,radrange,'ObjectPolarity','bright','Sensitivity',sens);
    BW = false(size(image,1),size(image,2));
    [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
    for n = 1:z
        try
            BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
        catch error
            if strcmp(anat, 'A')
                save(AscAof, 'AscAo')    
                save(DescAof, 'DescAo')
            else    
                save(AbdAof, 'AscAo')
            end
            close all
            cd(home_dir)
            disp(i)
            rethrow(error)
        end
    end
    
    %% Morphological Dilation 
    BW = imdilate(BW,strel('disk',3)); %dilate so active contours pulls in segmentation
    
    %% Active Contours + Dilation
    BW = activecontour(image,BW,300,'Chan-Vese','SmoothFactor',5,'ContractionBias',0.2);
    BW = imdilate(BW,strel('disk',1));
    
    %% Freehand ROI conversion
    blocations = bwboundaries(BW,'noholes');
    f = figure;
    f.WindowState = 'maximized';
    imshow(image, []);
    for ind = 1:numel(blocations)
        % Convert to x,y order.
        pos = blocations{ind};
        pos = fliplr(pos);
        % Create a freehand ROI.
        h = drawfreehand('Position', pos);
        customWait(h);
    end
    
    hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
    editedMask = false(size(image));
    for ind = 1:numel(hfhs)
        % Accumulate the mask from each ROI
        editedMask = editedMask | hfhs(ind).createMask();
        % Include the boundary of the ROI in the final mask.
        boundaryLocation = round(hfhs(ind).Position);
        bInds = sub2ind(size(image), boundaryLocation(:,2), boundaryLocation(:,1));
        editedMask(bInds) = true;
    end
    
    
    %% Separate Ascending and Descending Aorta
    [r,c] = find(editedMask);
    [~,idx_min] = min(r);
    Asc = bwselect(editedMask,c(idx_min),r(idx_min),8);
    
    if strcmp(anat, 'A')
        [~,idx_max] = max(r);
        Desc = bwselect(editedMask,c(idx_max),r(idx_max),8);
    end
    
    %% Show and Save Images
    label = labeloverlay(image,editedMask);
    figure; imshow(label,[]);
    
    bssfp(:,:,i) = image;
    if strcmp(anat, 'A')
        AscAo(:,:,i) = Asc;
        DescAo(:,:,i) = Desc;
    else
        AscAo(:,:,i) = Asc; 
    end
end

save(bssfpf, 'bssfp')
if strcmp(anat, 'A')
    save(AscAof, 'AscAo')    
    save(DescAof, 'DescAo')
else    
    save(AbdAof, 'AscAo')
end
close all
cd(home_dir)

disp("bssfp masks saved.");

%% Opening PWV files
desiredSize = 512;
scan2 = sprintf('000%s.PWV_CartBH_%sAo', pwvnum, anat);
id2 = [patient '-' num2];
pwvDir = [pwd '\' folder '\' scan2];
cd(pwvDir);
untar([scan2 '.tgz'])
dirFiles2 = dir('*.dcm');
if strcmp(anat, 'A')
    masks = AscAo + DescAo;
    ascpwv_name  = [result_dir 'pwv_masks\\%s_pwv_AscAo.mat'];
    AscAopwvf = sprintf(ascpwv_name, id2);
    AscAopwv = zeros(desiredSize,desiredSize,40);
    descpwv_name  = [result_dir 'pwv_masks\\%s_pwv_DescAo.mat'];
    DescAopwvf = sprintf(descpwv_name, id2);
    DescAopwv = zeros(desiredSize,desiredSize,40);
else
    masks = AscAo;
    abdpwv_name  = [result_dir 'pwv_masks\\%s_pwv_AbdAo.mat'];
    AbdAof = sprintf(abdpwv_name, id2);
    AbdAopwv = zeros(desiredSize,desiredSize,40);
end

all = zeros(512,512,length(dirFiles2));
for i=1:length(dirFiles2)
    img = dicomread(dirFiles2(i).name);
    img = imresize(img,[desiredSize desiredSize]);
    all(:, :, i) = img;
end
pc = all(:,:,1:40);
mag = all(:,:,41:end);

%% interpolating bssfp masks to pwv
desiredFrames = size(mag,3);
interpMasks = zeros(size(masks,1),size(masks,2),desiredFrames);
for i=1:size(masks,1)
    for j=1:size(masks,2)
        pixelLine = squeeze(masks(i,j,:));
        x = 1:length(pixelLine);
        xi = linspace(1,length(pixelLine),desiredFrames);
        yi = interp1(x,pixelLine,xi,'spline');
        interpMasks(i,j,:) = yi;
    end 
end 
uncropped = padarray(interpMasks, [crop(1) crop(3)], 'pre');
uncropped = padarray(uncropped, [desiredSize-crop(2)-1 desiredSize-crop(4)-1], 'post');
final = round(uncropped);
areas = zeros(1,40);
flows = zeros(1,40);

%% separate ascending and descending aorta and calculate area and flow
for i=1:size(final,3)
    [r2,c2] = find(final(:,:,i));
    [~,idx_min2] = min(r2);
    Ascf = bwselect(final(:,:,i),c2(idx_min2),r2(idx_min2),8);
    if strcmp(anat, 'A')
        [~,idx_max2] = max(r2);
        Descf = bwselect(final(:,:,i),c2(idx_max2),r2(idx_max2),8);
    end
    if strcmp(anat, 'A')
        AscAopwv(:,:,i) = Ascf; % contains final ascending aorta masks
        DescAopwv(:,:,i) = Descf; % contains final descending aorta masks
        
        % TODO: calculate area and flow for each mask
    else
        AscAopwv(:,:,i) = Ascf; % contains final abdominal aorta masks
        
        % TODO: calculate area and flow for each mask
    end
end
 
if strcmp(anat, 'A')
    save(AscAopwvf, 'AscAopwv')    
    save(DescAopwvf, 'DescAopwv')
else    
    save(AbdAopwvf, 'AscAopwv')
end

disp("pwv masks saved.")


%% helper functions
function customWait(hROI)
    % Listen for mouse clicks on the ROI
    l = addlistener(hROI,'ROIClicked',@clickCallback);
    % Block program execution
    uiwait;
    % Remove listener
    delete(l);
    % Return the current position
    pos = hROI.Position;
end

function clickCallback(~,evt)
    if strcmp(evt.SelectionType,'double')
        uiresume;
    end
end