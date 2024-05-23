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
if strcmp(save_name, '2DPC_QA_Analysis')
    result_dir = save_loc;
else
    result_dir = fullfile(save_loc, '2DPC_QA_Analysis');
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


%% Load BSSFP
bssfp_data_folder = fullfile(result_dir, plane, 'bssfp');
if ~exist(bssfp_data_folder,'dir')
    mkdir(bssfp_data_folder);
end
disp("Loading bSSFP data...")
if ~exist(fullfile(bssfp_data_folder, 'bssfp.mat'),'file')
    disp("Select bSSFP scan folder")
    bssfp_folder = uigetdir(pwd, 'Select bSSFP scan folder');
    % bssfp_folder = 'D:\PWV\vol_Tarun\standard\015_FIESTA_Cine_BH_AAo_NOPROSP_14VPS_30FRAMES';  
    dirFiles = dir(fullfile(bssfp_folder, '*.dcm'));
    if isempty(dirFiles) %check if dicoms are zipped
        zipFiles = dir(fullfile(bssfp_folder, '*.tgz'));
        if isempty(zipFiles)
            disp('NO READABLE DICOMS OR TGZ FILE FOUND, TRY ANOTHER FOLDER');
        else
            gunzip(fullfile(bssfp_folder, '*.tgz'), 'Dicoms'); %unzip first
            dd = dir(fullfile(bssfp_folder, 'Dicoms')); %get the name of the only file in the new dir
            untar(dd(3).name,'Dicoms'); %untar that file
            movefile(fullfile(bssfp_folder, 'Dicoms/*'), bssfp_folder); %move unzipped files back up
            rmdir(fullfile(bssfp_folder, 'Dicoms'),'s') %get rid of created dummy unzipping folder
        end 
    end     
    bssfpInfo = dicominfo(fullfile(dirFiles(1).folder, dirFiles(1).name)); %grab dicom info from 1st frame
    bssfpHR = bssfpInfo.HeartRate; % bpm
    bssfpHeight = double(bssfpInfo.Height);
    bssfpWidth = double(bssfpInfo.Width);
    bssfpFrames = bssfpInfo.CardiacNumberOfImages;
    pixelArea = bssfpInfo.PixelSpacing(1).*bssfpInfo.PixelSpacing(2); % mm^2
    bssfp = zeros(bssfpHeight,bssfpWidth,bssfpFrames);
    for f=1:length(dirFiles)
        bssfp(:,:,f) = dicomread(fullfile(dirFiles(f).folder, dirFiles(f).name));
        temp = dicominfo(fullfile(dirFiles(f).folder, dirFiles(f).name));
        bssfpTimes(f) = temp.TriggerTime;
    end 
    save(fullfile(bssfp_data_folder, 'bssfp.mat'), 'bssfp');
    save(fullfile(bssfp_data_folder, 'bssfpTimes.mat'), 'bssfpTimes');
    save(fullfile(bssfp_data_folder, 'bssfpInfo.mat'), 'bssfpInfo');
else
    load(fullfile(bssfp_data_folder, 'bssfp.mat'));
    load(fullfile(bssfp_data_folder, 'bssfpTimes.mat'));
    load(fullfile(bssfp_data_folder, 'bssfpInfo.mat'));
    bssfpHeight = double(bssfpInfo.Height);
    bssfpWidth = double(bssfpInfo.Width);
    bssfpFrames = bssfpInfo.CardiacNumberOfImages;
    bssfpHR = bssfpInfo.HeartRate; % bpm
    pixelArea = bssfpInfo.PixelSpacing(1).*bssfpInfo.PixelSpacing(2); % mm^2
end 
tavg_bssfp = mean(bssfp, 3);

%% Initialize Masks
answer = questdlg("Load masks?");
switch answer
    case 'Yes'
        % select mask.mat file
        [mask_file, mask_path] = uigetfile({'*.mat'}, "Select mask file", fullfile(bssfp_data_folder, 'masks.mat'));
        load(fullfile(mask_path, mask_file));
        load(fullfile(mask_path, 'crop.mat'))
        load(fullfile(mask_path, 'areas.mat'))
        if size(masks_saved,2) < bssfpFrames
            needsMask = 1;
            currFrame = size(masks_saved,2);
        else
            needsMask = 0;
            currFrame = bssfpFrames;
        end
        mask = zeros(num_roi, bssfpHeight, bssfpWidth, bssfpFrames);
        area = zeros(num_roi,bssfpFrames);
        for i=1:currFrame
            for j=1:num_roi
                mask(j,:,:,i) = masks_saved{j, i};
                area(j, i) = areas_saved{j, i};
            end
        end

    case 'No'
        % Zoom in Manually
        figure; imshow(tavg_bssfp,[]); title('Select region to zoom in on');
        rect = drawrectangle;
        crop = round([rect.Position(2), rect.Position(2)+rect.Position(4), rect.Position(1), rect.Position(1)+rect.Position(3)]);
        save(fullfile(bssfp_data_folder, 'crop.mat'), 'crop')
        close all;
        needsMask = 1;
        currFrame = 1;
        masks_saved = cell(num_roi, 1);
        areas_saved = cell(num_roi, 1);
        mask = zeros(num_roi, bssfpHeight, bssfpWidth, bssfpFrames);
        area = zeros(num_roi, bssfpFrames);

    case 'Cancel'
        return
end 

disp("Segmenting bSSFP images...")
radrange = [10 30]; % radius range for Hough Transform
sens = 0.8; % sensitivity of circle find
%% Segment Each Frame Separately
if needsMask % if we haven't loaded in the mask...
    for i=currFrame:bssfpFrames
        fprintf("Frame %d\n", i)
        image = bssfp(crop(1):crop(2),crop(3):crop(4),i);
        % image = rescale(image);
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
                [BW, centers, radii, radrange, sens] = circleFinder(image, num_roi, radrange, sens);
            case 'hough+contours' % Hough Transform + Active Countours
                [BW, centers, radii, radrange, sens] = circleFinder(image, num_roi, radrange, sens);
                try
                    BW = imdilate(BW,strel('disk',2)); %dilate so active contours pulls in segmentation
                    iters = 300;
                    smF = 5;
                    contrF = 0.2; %bias towards shrinking
                    method = 'Chan-Vese'; %or Chan-Vese (default)
                    BW = activecontour(image,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                    BW = imdilate(BW,strel('disk',1));
                catch
                    disp('Active contours failed, continuing with circle mask');
                end
            case 'edge+hough+contours' % Hough Transform + Active Countours
                grad = rescale(imgradient(image));
                [BW, centers, radii, radrange, sens] = circleFinder(image, num_roi, radrange, sens);
                try
                    BW = imerode(BW,strel('disk',2)); %erode so active contours pushes out
                    iters = 300;
                    smF = 5;
                    contrF = -0.1; %bias towards growing
                    method = 'Chan-Vese'; %or Chan-Vese (default)
                    BW = activecontour(grad,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                    BW = imdilate(BW,strel('disk',1));
                catch
                    disp('Active contours failed, continuing with circle mask');
                end
            case 'circle+contours'
                f=figure; imshow(image,[]); 
                circle = drawcircle();
                radius = circle.Radius; %get radius of circle
                center = round(circle.Center); %get center coordinates
                close(f); delete(f);
                [X,Y] = ndgrid(1:size(image,1),1:size(image,2));
                X = X-center(2); %shift coordinate grid
                Y = Y-center(1);
                BW = sqrt(X.^2 + Y.^2)<=radius; %anything outside radius is ignored
                try
                    BW = imdilate(BW,strel('disk',2)); %dilate so active contours pulls in segmentation
                    iters = 100;
                    smF = 5;
                    contrF = 0.1; %bias towards shrinking
                    method = 'Chan-Vese'; %or Chan-Vese (default)
                    BW = activecontour(image,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                    BW = imdilate(BW,strel('disk',1));
                catch
                    disp('Active contours failed, continuing with circle mask');
                end
            otherwise
                disp("Invalid segmentation choice!")
                return
        end

        % Freehand ROI conversion
        blocations = bwboundaries(BW,'noholes');
        numBlobs = numel(blocations);
        f = figure;
        imshow(image, []);
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
            vessel = padarray(vessel, [bssfpHeight-crop(2)-1 bssfpWidth-crop(4)-1], 'post');
            mask(j,:,:,i) = vessel;
            area(j, i) = sum(vessel(:))*pixelArea; % mm^2
            masks_saved{j, i} = mask(j,:,:,i);
            areas_saved{j, i} = area(j, i);
        end
        close all
        save(fullfile(bssfp_data_folder, 'masks.mat'), 'masks_saved');
        save(fullfile(bssfp_data_folder, 'areas.mat'), 'areas_saved');
    end
    % save(fullfile(bssfp_data_folder, 'masks.mat'), 'mask');
    % save(fullfile(bssfp_data_folder, 'areas.mat'), 'area');
    disp('Mask and Area Data Saved');
end

%% Load 2DPC

pc_data_folder = fullfile(result_dir, plane, 'pc');
if ~exist(pc_data_folder,'dir')
    mkdir(pc_data_folder);
end
disp("Loading 2DPC data...")
if ~exist(fullfile(pc_data_folder, 'pc.mat'),'file')
    disp("Select 2DPC scan folder")
    pc_folder = uigetdir(pwd,'Select 2DPC scan folder');
    % pc_folder = 'D:\PWV\vol_Tarun\standard\008_PWV_CartBH_AAo_NOPROSP_4VPS_40FRAMES';
    dirFiles = dir(fullfile(pc_folder, '*.dcm'));
    if isempty(dirFiles) %check if dicoms are zipped
        zipFiles = dir(fullfile(pc_folder, '*.tgz'));
        if isempty(zipFiles)
            disp('NO READABLE DICOMS OR TGZ FILE FOUND, TRY ANOTHER FOLDER');
        else
            gunzip(fullfile(pc_folder, '*.tgz'), 'Dicoms'); %unzip first
            dd = dir(fullfile(pc_folder, 'Dicoms')); %get the name of the only file in the new dir
            untar(dd(3).name,'Dicoms'); %untar that file
            movefile(fullfile(bssfp_folder, 'Dicoms/*'), bssfp_folder); %move unzipped files back up
            rmdir(fullfile(bssfp_folder, 'Dicoms'),'s') %get rid of created dummy unzipping folder
        end 
    end 
    
    pcInfo = dicominfo(fullfile(dirFiles(1).folder, dirFiles(1).name)); %grab dicom info from 1st frame
    pcHeight = double(pcInfo.Height);
    pcWidth = double(pcInfo.Width);
    pcFrames = pcInfo.CardiacNumberOfImages;
    pcHR = pcInfo.HeartRate; % bpm
    pcmr = zeros(pcHeight,pcWidth,pcFrames);
    for f=1:length(dirFiles)
        pcmr(:,:,f) = dicomread(fullfile(dirFiles(f).folder, dirFiles(f).name));
        temp = dicominfo(fullfile(dirFiles(f).folder, dirFiles(f).name));
        pcmrTimes(f) = temp.TriggerTime;
    end
    pc = pcmr(:,:,1:pcFrames);
    mag = pcmr(:,:,(pcFrames+1):end);
    pcmrTimes = pcmrTimes(1:pcFrames);
    save(fullfile(pc_data_folder, 'pc.mat'), 'pc');
    save(fullfile(pc_data_folder, 'mag.mat'), 'mag');
    save(fullfile(pc_data_folder, 'pcmrTimes.mat'), 'pcmrTimes');
    save(fullfile(pc_data_folder, 'pcInfo.mat'), 'pcInfo');
else
    load(fullfile(pc_data_folder, 'pc.mat'));
    load(fullfile(pc_data_folder, 'mag.mat'));
    load(fullfile(pc_data_folder, 'pcmrTimes.mat'))
    load(fullfile(pc_data_folder, 'pcInfo.mat'));
    pcHeight = double(pcInfo.Height);
    pcWidth = double(pcInfo.Width);
    pcFrames = pcInfo.CardiacNumberOfImages;
    pcHR = pcInfo.HeartRate; % bpm
end


%% Make temporal resolutions equivalent
% desiredFrames = max(bssfpFrames, pcFrames);
% pcmrTimes_int = interp1((1:pcFrames),pcmrTimes,linspace(1,pcFrames,desiredFrames),'linear');
% bssfpTimes_int = interp1((1:bssfpFrames),bssfpTimes,linspace(1,bssfpFrames,desiredFrames),'linear');

% bssfpTimes_int = interp1(bssfpTimes, 1:1:maxTime, 'linear', 'extrap');
% pcmrTimes_int = interp1(1:pcmrFrames, pcmrTimes, timesInterp, 'linear', 'extrap');

% mag_int = imresize3(mag,[desiredHeight desiredWidth desiredFrames]);

% interp_quest = questdlg("Load interpolated data?");
% switch interp_quest
%     case "Yes"
        % minTime = min(min(pcmrTimes),min(bssfpTimes));
        % if minTime < 1
        %     minTime = 1;
        % end
        % maxTime = max(max(pcmrTimes),max(bssfpTimes));
        % timesInterp = 1:1:maxTime;
        % desiredFrames = length(timesInterp);
        % desiredHeight = max(bssfpHeight, pcHeight);
        % desiredWidth = max(bssfpWidth, pcWidth);
        % bssfpTimes_int = timesInterp;
        % pcmrTimes_int = timesInterp;
        % 
        % load(fullfile(pc_data_folder, 'pc_interp.mat'), 'pc_int');
        % load(fullfile(bssfp_data_folder, 'bssfp_interp.mat'), 'bssfp_int');
        % load(fullfile(bssfp_data_folder, 'mask_interp.mat'), 'mask_int');
        % load(fullfile(bssfp_data_folder, 'area_interp.mat'), 'area_int');

    % case "No"

% add duplicated beginning and end points to both ends to make interpolation cylical
disp("Interpolating data...")
desiredHeight = max(bssfpHeight, pcHeight);
desiredWidth = max(bssfpWidth, pcWidth);

minTime = min(min(pcmrTimes),min(bssfpTimes));
if minTime < 1
    minTime = 1;
end
maxTime = max(max(pcmrTimes),max(bssfpTimes));
timesInterp = 1:1:maxTime;
desiredFrames = length(timesInterp);

bssfp_adj = cat(3, bssfp(:,:,end), bssfp, bssfp(:,:,1));
bssfp_adj = imresize3(bssfp_adj,[desiredHeight desiredWidth bssfpFrames+2]);
bssfpTimes_adj = cat(2, minTime-1, bssfpTimes, maxTime+1);

pc_adj = cat(3, pc(:,:,end), pc, pc(:,:,1));
pc_adj = imresize3(pc_adj,[desiredHeight desiredWidth pcFrames+2]);
pcmrTimes_adj = cat(2, minTime-1, pcmrTimes, maxTime+1);

mask_adj = cat(4, mask(:,:,:,end), mask, mask(:,:,:,1));
area_adj = cat(2, area(:,end), area, area(:,1));

bssfpTimes_int = timesInterp;
[borig_x, borig_y, borig_t] = meshgrid(1:desiredHeight, 1:desiredWidth, bssfpTimes_adj);
[bnew_x, bnew_y, bnew_t] = meshgrid(1:desiredHeight, 1:desiredWidth, bssfpTimes_int);
bssfp_int = interp3(borig_x, borig_y, borig_t, bssfp_adj, bnew_x, bnew_y, bnew_t, 'linear');

pcmrTimes_int = timesInterp;
[porig_x, porig_y, porig_t] = meshgrid(1:desiredHeight, 1:desiredWidth, pcmrTimes_adj);
[pnew_x, pnew_y, pnew_t] = meshgrid(1:desiredHeight, 1:desiredWidth, pcmrTimes_int);
pc_int = interp3(porig_x, porig_y, porig_t, pc_adj, pnew_x, pnew_y, pnew_t, 'linear');

mask_int = zeros(num_roi, desiredHeight, desiredWidth, desiredFrames);
area_int = zeros(num_roi, desiredFrames);
for j=1:num_roi
    % mask_int(j, :,:,:) = imresize3(squeeze(mask(j, :,:,:)), [desiredHeight desiredWidth desiredFrames]);
    mask_int(j, :,:,:) = interp3(borig_x, borig_y, borig_t, squeeze(mask_adj(j, :,:,:)), bnew_x, bnew_y, bnew_t, 'linear');
    % area_int(j, :) = interp1((1:length(area),squeeze(area(j,:)),linspace(1,length(area),desiredFrames), 'linear');
    area_int(j, :) = interp1(bssfpTimes_adj, squeeze(area_adj(j, :)), bssfpTimes_int, 'linear');
end
        
%         save(fullfile(pc_data_folder, 'pc_interp.mat'), 'pc_int');
%         save(fullfile(bssfp_data_folder, 'bssfp_interp.mat'), 'bssfp_int');
%         save(fullfile(bssfp_data_folder, 'mask_interp.mat'), 'mask_int');
%         save(fullfile(bssfp_data_folder, 'area_interp.mat'), 'area_int');
%         disp("Interpolated data saved")
% 
%     case "Cancel"
%         return
% end


%% Save ROI images
e = figure; imshow(bssfp(:,:,1), []);
frame_bssfp = getframe(e);
imwrite(frame2im(frame_bssfp), fullfile(analysis_dir, 'bssfp.png'));
close(e);
f = figure; imshow(imoverlay(rescale(bssfp(:,:,1)), squeeze(sum(sum(mask, 1), 4))));
frame_bssfp_mask = getframe(f);
imwrite(frame2im(frame_bssfp_mask), fullfile(analysis_dir, 'bssfp_mask.png'));
close(f);
g = figure; imshow(rescale(pc_int(:,:,1)), []);
frame_pc = getframe(g);
imwrite(frame2im(frame_pc), fullfile(analysis_dir, 'pc.png'));
close(g);
combined_mask = squeeze(sum(mask_int, 1));
h = figure; imshow(imoverlay(rescale(pc_int(:,:,1)), combined_mask(:,:,1)), []);
frame_pc_mask = getframe(h);
imwrite(frame2im(frame_pc_mask), fullfile(analysis_dir, 'pc_mask.png'));
close(h);


%% Calculating flow
% flow_quest = questdlg("Load flow data?");
% switch flow_quest
%     case "Yes"
%         load(fullfile(pc_data_dir, 'flows.mat'), 'flow');
%     case "No"
        
disp("Calculating flow...")
flow = zeros(num_roi, desiredFrames);
for j=1:num_roi
    for i=1:desiredFrames
        ROI = imbinarize(mask_int(j,:,:,i));
        vTemp = pc_int(:,:,i); %single frame velocities (mm/s)
        ROIindex = double(vTemp(ROI)); %indexed velocities within mask
        vMean = mean(ROIindex); %mean velocity in frame i (mm/s)
        flow(j, i) = area_int(j, i).*vMean.*0.001; %flow in frame i (mL/s)
    end
end
save(fullfile(pc_data_folder, 'flows.mat'), 'flow');
disp('Flow data saved')
%     case "Cancel"
%         return
% end


%% Shift the data due to periph. gating lag + prosp. gating lags
disp("Shifting velocity data...")
flag = 0;
shift = 0;
options.WindowStyle = 'normal';
while ~flag
    f = figure;
    if mean(flow(1,:))<0
        plot(circshift(-1*flow(1, :), shift)); 
    else
        plot(circshift(flow(1, :), shift));
    end
    title("Shift flow curve")
    xlabel("Time (ms)")
    ylabel("Flow (mL/s)")
    % movegui(f, 'west');
    answer = inputdlg({"Enter shift amount (ms):", "Accept? (enter 1)"}, "Velocity curve shift", ...
        [1 30; 1 30], {num2str(shift), num2str(flag)}, options);
    shift = str2double(answer{1});
    flag = str2double(answer{2});
    close(f)
end
bssfp_int = circshift(bssfp_int,shift,3);
pc_int = circshift(pc_int,shift,3);
% mag_int = circshift(mag_int,shift,3);
for i=1:num_roi
    mask_int(i, :,:,:) = circshift(mask_int(i, :,:,:), shift, 4);
    area_int(i, :) = circshift(area_int(i, :), shift, 2);
    flow(i, :) = circshift(flow(i, :), shift, 2);
end


%% Compute PWV
disp("Calculating PWV...")
PWV = zeros(num_roi, 1);
SV = zeros(num_roi, 1);

flow_calc = flow;
for j=1:num_roi
    if mean(flow_calc(j,:))<0
        flow_calc(j, :) = -1*flow_calc(j, :);
    end

    %Define Early Systolic Region (flow vs. time)
    figure; plot(flow_calc(j, :)); 
    title(sprintf('%s: ROI %d - Early Systole', plane, j)); xlabel('Time (ms)'); ylabel('Flow (mL/s)');
    disp("Draw around the early systole (upslope) portion of the flow curve")
    flow_plot = gcf;
    free = drawfreehand;
    earlySystolePts_int = find(inpolygon(linspace(1,length(flow_calc(j, :)),length(flow_calc(j, :))),flow_calc(j, :),free.Position(:,1),free.Position(:,2)));
    earlySystoleTimes_int = pcmrTimes_int(earlySystolePts_int);
    earlySystolePts = ismember(earlySystoleTimes_int, pcmrTimes);
    earlySystoleTimes = earlySystoleTimes_int(earlySystolePts);
    delete(free);
    saveas(flow_plot, fullfile(analysis_dir, sprintf('flowplot_ROI_%d', j)));
    close(flow_plot); clear flow_plot;
    SV(j) = trapz(flow_calc(j, systolePts));

    % if sys_flag
    %     %Define Systolic Region
    %     figure; plot(flow_calc(j, :)); 
    %     title(sprintf('%s: ROI %d - Systole', plane, j)); xlabel('Time (ms)'); ylabel('Flow (mL/s)');
    %     disp("Draw around the systole portion of the flow curve")
    %     flow_plot = gcf;
    %     free = drawfreehand;
    %     systolePts = find(inpolygon(linspace(1,length(flow_calc(j, :)),length(flow_calc(j, :))),flow_calc(j, :),free.Position(:,1),free.Position(:,2)));
    %     systoleTimes = pcmrTimes_int(systolePts);
    %     delete(free); close(flow_plot); clear flow_plot;
    %     ESV(j) = trapz(systoleTimes*0.001, flow_calc(j, systolePts));
    % 
    %     %Define Diastolic Region
    %     figure; plot(flow_calc(j, :)); 
    %     title(sprintf('%s: ROI %d - Diastole', plane, j)); xlabel('Time (ms)'); ylabel('Flow (mL/s)');
    %     disp("Draw around the diastole portion of the flow curve")
    %     flow_plot = gcf;
    %     free = drawfreehand;
    %     diastolePts = find(inpolygon(linspace(1,length(flow_calc(j, :)),length(flow_calc(j, :))),flow_calc(j, :),free.Position(:,1),free.Position(:,2)));
    %     diastoleTimes = pcmrTimes_int(diastolePts);
    %     delete(free); close(flow_plot); clear flow_plot;
    %     EDV(j) = trapz(diastoleTimes*0.001, flow_calc(j, diastolePts));
    % end
    
    %Define Linear QA Region (flow vs. area)
    % x = area_int(j, earlySystolePts_int);
    % y = flow_calc(j, earlySystolePts_int);
    x = area_int(j, earlySystolePts);
    y = flow_calc(j, earlySystolePts);
    figure; scatter(x,y); 
    title(sprintf('%s: ROI %d - QA Plot', plane, j)); xlabel('Area (mm^2)'); ylabel('Flow (mL/s)');
    earlySystolePts2 = [x; y]';
    [coef,stats] = polyfit(earlySystolePts2(:,1),earlySystolePts2(:,2),1);
    
    %Linear regression
    minArea = min(earlySystolePts2(:,1));
    maxArea = max(earlySystolePts2(:,1));
    xq = linspace(minArea-0.1,maxArea+0.1,100);
    yq = coef(1)*xq + coef(2); %y = mx+b
    hold on; 
    scatter(earlySystolePts2(:,1),earlySystolePts2(:,2),72,'k','x');
    plot(xq,yq); 
    str = [];
    text(min(xq),max(yq)-0.1*max(yq),sprintf('Slope = %.4f \n Intcpt = %.4f', coef(1), coef(2)));
    hold off;
    QAplot = gcf; 
    PWV(j) = coef(1); % m/s
    fprintf('    ROI %d PWV_QA = %.4f m/s\n', j, PWV(j));
    fprintf('    ROI %d Stroke Volume = %.4f mL\n', j, SV(j));
    fprintf('    ROI %d Cardiac Output = %.4f (mL/min)\n', j, SV(j)*pcHR);
    saveas(QAplot,fullfile(analysis_dir, sprintf('QAplot_ROI_%d', j)));
    save(fullfile(analysis_dir, sprintf('coef_ROI_%d.mat', j)), 'coef');
    save(fullfile(analysis_dir, sprintf('stats_ROI_%d.mat', j)), 'stats');
end

% Initialize array for Excel file
results = cell(num_roi+2, 6);
results{1, 1} = 'ROI';
results{1, 2} = 'PWV (m/s)';
if sys_flag
    results{1, 3} = 'Stroke Volume (mL)';
    results{1, 4} = 'Cardiac Output (mL/min)';
    results{1, 5} = 'Mean HR (bpm)';
    results{2, 5} = pcHR;
end

for j = 1:num_roi
    results{j+1, 1} = j;
    results{j+1, 2} = PWV(j);
    if sys_flag
        results{j+1, 3} = SV(j);
        results{j+1, 4} = SV(j) * pcHR;
    end
end

results2 = cell(desiredFrames+1,1+2*num_roi);
results2{1, 1} = 'Time (ms)';
for i=1:desiredFrames
    results2{i+1, 1} = timesInterp(i);
end
for j = 1:num_roi
    results2{1, 2*j} = sprintf('ROI%d Area (mm^2)', j);
    results2{1, 2*j+1} = sprintf('ROI%d Flow (mL/s)', j);
    for i = 1:desiredFrames
        results2{i+1, 2*j} = area_int(j, i);
        results2{i+1, 2*j+1} = flow(j, i);
    end
end

writecell(results, fullfile(analysis_dir, 'results.xlsx'))
writecell(results2, fullfile(analysis_dir, 'results.xlsx'), 'WriteMode', 'append')
disp("( つ ◕_◕ )つ PWV-QA Data Saved! ( つ ◕_◕ )つ")


%% Helper functions
function [BW, centers, radii, radrange, sens] = circleFinder(image, num_roi, radrange, sens)
    BW = false(size(image,1),size(image,2));
    [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
    flag = 0;
    mflag = 0;
    while ~flag
        g = figure;
        imshow(image, []);
        movegui(g, 'onscreen');
        [centers,radii,~] = imfindcircles(image,radrange,'ObjectPolarity','bright','Sensitivity',sens);
        for n = 1:num_roi
            try
                BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
            catch error
                
            end
        end
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
            for z=1:num_roi
                k = figure; 
                imshow(image, []);
                circle = drawcircle();
                radii(z) = circle.Radius; %get radius of circle
                centers(z,:) = round(circle.Center); %get center coordinates
                close(k);
                BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
            end
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