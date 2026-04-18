clear all; close all;

disp("Choose where to save analysis data")
save_loc = uigetdir(pwd, "Choose where to save analysis data");

plane = inputdlg("Type in the analysis region name");
plane = plane{1};

num_roi = inputdlg("Type in the number of ROIs in plane");
num_roi = str2num(num_roi{1});

list_options = {'kmeans','otsu','circle+contours','hough','hough+contours','edge+hough+contours'};
[choice, tf] = listdlg('ListString', list_options, 'InitialValue', 5, 'PromptString', "Choose a segmentation method", 'SelectionMode', 'single');
segType = list_options{choice};

% cardiac_answer = questdlg("Do cardiac output analysis?");
% switch cardiac_answer
%     case "Yes"
%         sys_flag = 1;
%     case "No"
%         sys_flag = 0;
%     case "Cancel"
%         sys_flag = 0;
% end

% plane = 'asc';
% segType = 'hough+contours';
% num_roi = 2;

date = string(datetime('now', 'Format', 'yyyy-MM-dd-HHmm'));

[~, save_name] = fileparts(save_loc);
if strcmp(save_name, '2DSegmentationAnalysis')
    result_dir = save_loc;
else
    result_dir = fullfile(save_loc, '2DSegmentationAnalysis');
    if ~exist(result_dir,'dir')
        mkdir(result_dir);
    end
end
% save_loc = 'D:\PWV\vol_Tarun';
analysis_dir = fullfile(result_dir, plane, strcat(segType, '_', date));
if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir)
end
disp("Created analysis folder")

%% Load 2DPC

disp("Loading 2D data...")
disp("Select 2D scan folder")
pc_folder = uigetdir(pwd,'Select 2D scan folder');
dirFiles = dir(fullfile(pc_folder));
dirFiles = dirFiles(~ismember({dirFiles.name},{'.','..'}));
pc_files = {dirFiles.name};

if any(strcmp(pc_files, 'pcvipr_header.txt'))
    data_type = 'radial';

elseif any(endsWith(pc_files, '.dcm'))
    data_type = 'dicom';
end

if strcmp(data_type, 'dicom')
    pcInfo = dicominfo(fullfile(dirFiles(1).folder, dirFiles(1).name)); %grab dicom info from 1st frame
    pcHeight = double(pcInfo.Height);
    pcWidth = double(pcInfo.Width);
    pcFrames = pcInfo.CardiacNumberOfImages;
    pcHR = pcInfo.HeartRate; % bpm
    pixelArea = pcInfo.PixelSpacing(1).*pcInfo.PixelSpacing(2);
    pcmr = zeros(pcHeight,pcWidth,pcFrames);
    for f=1:length(dirFiles)
        pcmr(:,:,f) = dicomread(fullfile(dirFiles(f).folder, dirFiles(f).name));
        temp = dicominfo(fullfile(dirFiles(f).folder, dirFiles(f).name));
        pcmrTimes(f) = temp.TriggerTime;
    end
    if pcInfo.Private_0019_10cc > 0 % check if venc set, if so its prob phase contrast image
        mag = pcmr(:,:,(pcFrames+1):end); 
    else
        mag = pcmr;
    end
    pcmrTimes = pcmrTimes(1:pcFrames);
    MAG = mean(mag, 3);
elseif strcmp(data_type, 'radial')
    fid = fopen(fullfile(dirFiles(1).folder,'pcvipr_header.txt'), 'r'); %open header
    dataArray = textscan(fid,'%s%s%[^\n\r]','Delimiter',' ', ...
        'MultipleDelimsAsOne',true,'ReturnOnError',false); %parse header info
    fclose(fid);
    dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:),'UniformOutput',false);
    pcviprHeader = cell2struct(dataArray{1,2}(:),dataArray{1,1}(:),1); %turn to structure
    fovx = pcviprHeader.fovx;
    fovy = pcviprHeader.fovy;
    resx = pcviprHeader.matrixx; %resolution in x
    resy = pcviprHeader.matrixy; %resolution in y
    nframes = pcviprHeader.frames; %number of cardiac frames
    timeres = pcviprHeader.timeres;
    MAG = load_dat(fullfile(dirFiles(1).folder,'MAG.dat'),[resx resy]); %Average magnitude
    CD = load_dat(fullfile(dirFiles(1).folder,'CD.dat'),[resx resy]); %Average complex difference

    % Initialize data time-resolved data arrays
    mag = zeros(resx,resy,nframes); %Time-resolved magnitude
    cd = zeros(resx,resy,nframes); %Time-resolved complex difference
    for j = 1:nframes  %velocity is placed in v3 for 2D (through-plane)
        mag(:,:,j) = load_dat(fullfile(dirFiles(1).folder,[filesep 'ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
        cd(:,:,j) = load_dat(fullfile(dirFiles(1).folder,[filesep 'ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
    end
  
    MAG = flipud(MAG);
    CD = flipud(CD);
    mag = flipud(mag);
    cd = flipud(cd); 
    
    pcHeight = resy;
    pcWidth = resx;
    pcFrames = nframes;
    pixelArea = (fovx/resx) * (fovy/resy);
    pcmrTimes = double(timeres.*(0:(nframes-1)));
    pcHR = 60/(timeres*nframes)*1000; %bpm
end



%% Initialize Masks
answer = questdlg("Load masks?");
switch answer
    case 'Yes'
        % select mask.mat file
        [mask_file, mask_path] = uigetfile({'*.mat'}, "Select mask file", fullfile(save_loc, 'masks.mat'));
        load(fullfile(mask_path, mask_file));
        load(fullfile(mask_path, 'crop.mat'))
        load(fullfile(mask_path, 'areas.mat'))
        if size(masks_saved,2) < pcFrames
            needsMask = 1;
            currFrame = size(masks_saved,2);
        else
            needsMask = 0;
            currFrame = pcFrames;
        end
        mask = zeros(num_roi, pcHeight, pcWidth, pcFrames);
        area = zeros(num_roi, pcFrames);
        for i=1:currFrame
            for j=1:num_roi
                mask(j,:,:,i) = masks_saved{j, i};
                area(j, i) = areas_saved{j, i};
            end
        end

    case 'No'
        % Zoom in Manually
        figure; imshow(MAG,[]); title('Select region to zoom in on');
        rect = drawrectangle;
        crop = round([rect.Position(2), rect.Position(2)+rect.Position(4), rect.Position(1), rect.Position(1)+rect.Position(3)]);
        save(fullfile(analysis_dir, 'crop.mat'), 'crop')
        close all;
        needsMask = 1;
        currFrame = 1;
        masks_saved = cell(num_roi, 1);
        areas_saved = cell(num_roi, 1);
        mask = zeros(num_roi, pcHeight, pcWidth, pcFrames);
        area = zeros(num_roi, pcFrames);

    case 'Cancel'
        return
end 

disp("Segmenting PC mag images...")
radrange = [10 30]; % radius range for Hough Transform
sens = 0.8; % sensitivity of circle find
%% Segment Each Frame Separately
if needsMask % if we haven't loaded in the mask...
    image_mag = MAG(crop(1):crop(2),crop(3):crop(4));
    for i=currFrame:pcFrames
        image = mag(crop(1):crop(2),crop(3):crop(4),i);
        % image = rescale(image);
        if i==currFrame
            switch segType
                case 'kmeans' % Kmeans Clustering
                    disp("ERROR! Not working at the moment. Please select a different method.")
                    return
                    % s = rng; rng('default');
                    % L = imsegkmeans(single(image),2,'NumAttempts',5);
                    % rng(s);
                    % labelOfInterest = L(round(size(image,1)/2),round(size(image,2)/2));
                    % BW = L == labelOfInterest;
                    % if sum(BW(:))==0
                    %     level = graythresh(image);
                    %     BW = imbinarize(image,level);
                    % end 
                case 'otsu' % Otsu's Automatic Threshold
                    disp("ERROR! Not working at the moment. Please select a different method.")
                    return
                    % level = graythresh(image);
                    % BW = imbinarize(image,level);
                case 'hough' % Hough Transform
                    [BW, centers, radii, radrange, sens] = circleFinder(image_mag, num_roi, radrange, sens);
                case 'hough+contours' % Hough Transform + Active Countours
                    [BW, centers, radii, radrange, sens] = circleFinder(image_mag, num_roi, radrange, sens);
                    try
                        BW2 = imdilate(BW,strel('disk',2)); %dilate so active contours pulls in segmentation
                        iters = 300;
                        smF = 3;
                        contrF = 0.1; %bias towards shrinking
                        method = 'Chan-Vese'; %or Chan-Vese (default)
                        BW2 = activecontour(image_mag,BW2,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                        BW2 = imdilate(BW2,strel('disk',1));
                    catch
                        disp('Active contours failed, continuing with circle mask');
                    end
                case 'edge+hough+contours' % Hough Transform + Active Countours
                    grad = rescale(imgradient(image_mag));
                    [BW, centers, radii, radrange, sens] = circleFinder(image_mag, num_roi, radrange, sens);
                    try
                        BW2 = imerode(BW,strel('disk',2)); %erode so active contours pushes out
                        iters = 300;
                        smF = 3;
                        contrF = -0.1; %bias towards growing
                        method = 'Chan-Vese'; %or Chan-Vese (default)
                        BW2 = activecontour(grad,BW2,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                        BW2 = imdilate(BW2,strel('disk',1));
                    catch
                        disp('Active contours failed, continuing with circle mask');
                    end
                case 'circle+contours'
                    f=figure; imshow(image_mag,[]); 
                    circle = drawcircle();
                    radius = circle.Radius; %get radius of circle
                    center = round(circle.Center); %get center coordinates
                    close(f); delete(f);
                    [X,Y] = ndgrid(1:size(image_mag,1),1:size(image_mag,2));
                    X = X-center(2); %shift coordinate grid
                    Y = Y-center(1);
                    BW = sqrt(X.^2 + Y.^2)<=radius; %anything outside radius is ignored
                    try
                        BW2 = imdilate(BW,strel('disk',2)); %dilate so active contours pulls in segmentation
                        iters = 100;
                        smF = 5;
                        contrF = 0.1; %bias towards shrinking
                        method = 'Chan-Vese'; %or Chan-Vese (default)
                        BW2 = activecontour(image_mag,BW2,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                        BW2 = imdilate(BW2,strel('disk',1));
                    catch
                        disp('Active contours failed, continuing with circle mask');
                    end
                otherwise
                    disp("Invalid segmentation choice!")
                    return
            end
        end

        % Freehand ROI conversion
        blocations = bwboundaries(BW2,'noholes');
        numBlobs = numel(blocations);
        if numBlobs ~= num_roi
            bwboundaries(BW,'noholes');
            numBlobs = numel(blocations);
        end
        f = figure;
        imshow(image, []);
        title(sprintf("Frame %d", i))
        f.WindowState = 'maximized';

        % Manual adjustment
        for ind = 1:numBlobs
            pos = blocations{ind}; %convert to x,y order.
            subx = int16(linspace(1,length(pos(:,1)),7)); %subsample positions (x)
            suby = int16(linspace(1,length(pos(:,2)),7)); %subsample positions (y)
            posSpline = interparc(90,pos(subx,1),pos(suby,2),'csape'); %closed spline
            posSpline = fliplr(posSpline);
            % Can 'ESC' to delete ROI, Right-click to Add Waypoint
            waypoints = zeros(90,1); %initialize index locations of waypoints
            waypoints(1:15:end) = 1; %set 8 evenly-spaced waypoints
            h = drawfreehand('Position', posSpline,'FaceAlpha',0.15,'LineWidth',1, ...
                'Multiclick',true,'Waypoints',logical(waypoints)); %create freehand ROI
            customWait(h); %see custom function below in 'helper functions'
        end

        hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
        editedMask = false(size(image));
        for ind = 1:numel(hfhs)
           editedMask = editedMask | hfhs(ind).createMask(); %accumulate mask from each ROI
           boundaryLocation = round(hfhs(ind).Position); %include ROI boundary
           bInds = sub2ind(size(image), boundaryLocation(:,2), boundaryLocation(:,1));
           editedMask(bInds) = true;
        end
        BW = editedMask;
        
        % sort circles by radii
        sortedCircles = sortrows([centers, radii], 3, 'descend');

        % Separate Vessel
        for j=1:num_roi
            vessel = bwselect(BW, sortedCircles(j, 1), sortedCircles(j, 2), 8);
            vessel = padarray(vessel, [crop(1) crop(3)], 'pre');
            vessel = padarray(vessel, [pcHeight-crop(2)-1 pcWidth-crop(4)-1], 'post');
            mask(j,:,:,i) = vessel;
            area(j, i) = sum(vessel(:))*pixelArea; % mm^2
            masks_saved{j, i} = mask(j,:,:,i);
            areas_saved{j, i} = area(j, i);
        end
        close all
        save(fullfile(analysis_dir, 'masks.mat'), 'masks_saved');
        save(fullfile(analysis_dir, 'areas.mat'), 'areas_saved');
    end
end


%% Save ROI images
g = figure; imshow(rescale(mag(:,:,1)), []);
frame_pc = getframe(g);
imwrite(frame2im(frame_pc), fullfile(analysis_dir, 'mag.png'));
close(g);
combined_mask = squeeze(sum(mask, 1));
h = figure; imshow(imoverlay(rescale(mag(:,:,1)), combined_mask(:,:,1)), []);
frame_pc_mask = getframe(h);
imwrite(frame2im(frame_pc_mask), fullfile(analysis_dir, 'mag_mask.png'));
close(h);


% Initialize array for Excel file
results = cell(num_roi+2, 6);
results{1, 1} = 'ROI';
results{2, 1} = num_roi;
results{1, 2} = 'Mean HR (bpm)';
results{2, 2} = pcHR;

results2 = cell(pcFrames+1,1+2*num_roi);
results2{1, 1} = 'Time (ms)';
for i=1:pcFrames
    results2{i+1, 1} = pcmrTimes(i);
end
for j = 1:num_roi
    results2{1, 2*j} = sprintf('ROI%d Area (mm^2)', j);
    results2{1, 2*j+1} = sprintf('ROI%d Diameter (mm)', j);
    for i = 1:pcFrames
        results2{i+1, 2*j} = area(j, i);
        results2{i+1, 2*j+1} = 2*sqrt(area(j, i)/pi);
    end
end

writecell(results, fullfile(analysis_dir, 'results.xlsx'))
writecell(results2, fullfile(analysis_dir, 'results.xlsx'), 'WriteMode', 'append')
fprintf("( つ ◕_◕ )つ Segmentation Data Saved to %s!\n", fullfile('2DSegmentationAnalysis', plane, strcat(segType, '_', date)))


%% Helper functions
function [BW, centers, radii, radrange, sens] = circleFinder(image, num_roi, radrange, sens)
    BW = false(size(image,1),size(image,2));
    [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
    flag = 0;
    mflag = 0;
    while ~flag
        g = figure;
        imshow(image, []);
        q = gcf;
        movegui(g, 'west');
        q.Position(3:4) = 5.*q.Position(3:4);
        g.WindowState = 'maximized';
        [centers,radii,~] = imfindcircles(image,radrange,'ObjectPolarity','bright','Sensitivity',sens);
        viscircles(centers, radii);
        options.WindowStyle = 'normal';
        answer = inputdlg({"Enter lower radius range:", "Enter upper radius range:", "Enter sensitivity:", "Draw manually? (enter 1)", "Accept? (enter 1)"}, ...
            "Circle Estimation Parameters", [1 30; 1 30; 1 30; 1 30; 1 30;], {num2str(radrange(1)), num2str(radrange(2)), ...
            num2str(sens), num2str(mflag), num2str(flag)}, options);
        close(g)
        radrange = [str2num(answer{1}) str2num(answer{2})];
        sens = str2double(answer{3});
        mflag = str2num(answer{4});
        flag = str2num(answer{5});
        if mflag
            centers = zeros(num_roi, 2);
            radii = zeros(num_roi, 1);
            k = figure; 
            imshow(image, []);
            k.WindowState = 'maximized';
            for z=1:num_roi
                % movegui(k, 'west');
                % q = gcf;
                % q.Position(3:4) = 5.*q.Position(3:4);
                circle = drawcircle();
                wait(circle)
                radii(z) = circle.Radius; %get radius of circle
                centers(z,:) = round(circle.Center); %get center coordinates
                % close(k);
            end
            close(k)
            flag = 1;
        end
        if flag
            if size(centers, 1) ~= num_roi
                fprintf("Error: Number of ROIs is greater than or less than %d.\n" + ...
                    "Please continue adjusting parameters\n", num_roi);
                flag = 0;
            end
        end
    end
    for n = 1:num_roi
        try
            BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
        catch error
            disp(error)
        end
    end
end

function circ = circleROI(x,y,r,imgSize)
    [xx,yy] = ndgrid((1:imgSize)-y,(1:imgSize)-x);
    circ = (xx.^2 + yy.^2)<r^2;
end

function customWait(hROI)
    % Listen for mouse clicks on the ROI
    l = addlistener(hROI,'ROIClicked',@clickCallback);
    uiwait; % Block program execution
    delete(l); % Remove listener
    pos = hROI.Position; % Return the current position
end

function clickCallback(~,evt)
    if strcmp(evt.SelectionType,'double')
        uiresume;
    end
end