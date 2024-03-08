clear all; close all; clc;

options = {'kmeans','otsu','circle+contours','hough','hough+contours','edge+hough+contours'};
[choice, tf] = listdlg('ListString', options, 'InitialValue', 5, 'PromptString', "Choose a segmentation method", 'SelectionMode', 'single');
segType = options{choice};

plane = inputdlg("Type in image plane name");
plane = plane{1};

num_roi = inputdlg("Type in number of ROIs in plane");
num_roi = str2num(num_roi{1});

% plane = 'asc';
% segType = 'hough+contours';
% num_roi = 2;

date = string(datetime('now', 'Format', 'yyyy-MM-dd-HHmm'));
save_loc = uigetdir(pwd, "Choose where to save analysis data");
% save_loc = 'D:\PWV\vol_Tarun';
disp("Created analysis folder")
cd(save_loc)
result_dir = fullfile(save_loc, '2DPC_QA_Analysis');
if ~exist(result_dir,'dir')
    mkdir(result_dir);
end
analysis_dir = fullfile(result_dir, plane, strcat(segType, '_', date));

%% Load BSSFP
bssfp_folder = uigetdir(pwd, 'Select bSSFP Scan');
% bssfp_folder = 'D:\PWV\vol_Tarun\standard\015_FIESTA_Cine_BH_AAo_NOPROSP_14VPS_30FRAMES';
disp("Loading bSSFP data...")
bssfp_data_folder = fullfile(result_dir, plane, 'bssfp');
if ~exist(bssfp_data_folder,'dir')
    mkdir(bssfp_data_folder);
end

if ~exist(fullfile(bssfp_data_folder, 'bssfp.mat'),'file')
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
    bssfpSize = double(bssfpInfo.Height); %assumes matrix size is equal in x and y
    bssfpFrames = bssfpInfo.CardiacNumberOfImages;
    bssfp = zeros(bssfpSize,bssfpSize,bssfpFrames);
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
end 
bssfpSize = double(bssfpInfo.Height); %assumes matrix size is equal in x and y
bssfpFrames = bssfpInfo.CardiacNumberOfImages;
pixelArea = bssfpInfo.PixelSpacing(1).*bssfpInfo.PixelSpacing(2)*0.01; %cm^2

tavg_bssfp = mean(bssfp, 3);

%% Initialize Masks
if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir)
end
answer = questdlg("Load masks?");
switch answer
    case 'Yes'
        % select mask.mat file
        [mask_file, mask_path] = uigetfile({'*.mat'}, "Select mask file");
        load(fullfile(mask_path, mask_file));
        load(fullfile(mask_path, 'crop.mat'))
        load(fullfile(mask_path, 'areas.mat'))
        needsMask = 0;

    case 'No'
        % Zoom in Manually
        figure; imshow(tavg_bssfp,[]); title('ZOOM');
        rect = drawrectangle;
        crop = round([rect.Position(2), rect.Position(2)+rect.Position(4), rect.Position(1), rect.Position(1)+rect.Position(3)]);
        save(fullfile(bssfp_data_folder, 'crop.mat'), 'crop')
        close all;

        mask = zeros(num_roi, bssfpSize, bssfpSize, bssfpFrames);
        area = zeros(num_roi, bssfpFrames);
        needsMask = 1;

    case 'Cancel'
        return
end 

radrange = [10 30]; % radius range for Hough Transform
sens = 0.8; % sensitivity of circle find
%% Segment Each Frame Separately
if needsMask % if we haven't loaded in the mask...
    for i=1:bssfpFrames
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
                BW = imdilate(BW,strel('disk',2)); %dilate so active contours pulls in segmentation
                iters = 300;
                smF = 5;
                contrF = 0.2; %bias towards shrinking
                method = 'Chan-Vese'; %or Chan-Vese (default)
                BW = activecontour(image,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                BW = imdilate(BW,strel('disk',1));
            case 'edge+hough+contours' % Hough Transform + Active Countours
                grad = rescale(imgradient(image));
                [BW, centers, radii, radrange, sens] = circleFinder(image, num_roi, radrange, sens);
                BW = imerode(BW,strel('disk',2)); %erode so active contours pushes out
                iters = 300;
                smF = 5;
                contrF = -0.1; %bias towards growing
                method = 'Chan-Vese'; %or Chan-Vese (default)
                BW = activecontour(grad,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                BW = imdilate(BW,strel('disk',1));
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
                BW = imdilate(BW,strel('disk',2)); %dilate so active contours pulls in segmentation
                iters = 100;
                smF = 5;
                contrF = 0.1; %bias towards shrinking
                method = 'Chan-Vese'; %or Chan-Vese (default)
                BW = activecontour(image,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
                BW = imdilate(BW,strel('disk',1));
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
            vessel = padarray(vessel, [bssfpSize-crop(2)-1 bssfpSize-crop(4)-1], 'post');
            mask(j,:,:,i) = vessel;
            area(j, i) = sum(vessel(:))*pixelArea; %cm^2
        end
        close all
    end
    save(fullfile(bssfp_data_folder, 'masks.mat'), 'mask');
    save(fullfile(bssfp_data_folder, 'areas.mat'), 'area');
    disp('Mask and Area Data Saved');
end

%% Load 2DPC
pc_folder = uigetdir(pwd,'Select folder 2DPC scan');
% pc_folder = 'D:\PWV\vol_Tarun\standard\008_PWV_CartBH_AAo_NOPROSP_4VPS_40FRAMES';
disp("Loading 2DPC data...")
pc_data_folder = fullfile(result_dir, plane, 'pc');
if ~exist(pc_data_folder,'dir')
    mkdir(pc_data_folder);
end
if ~exist(fullfile(pc_data_folder, 'pc.mat'),'file')
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
    pcSize = double(pcInfo.Height); %assumes matrix size is equal in x and y
    pcFrames = pcInfo.CardiacNumberOfImages;
    pcmr = zeros(pcSize,pcSize,pcFrames);
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
end 
pcSize = double(pcInfo.Height); %assumes matrix size is equal in x and y
pcFrames = pcInfo.CardiacNumberOfImages;


%% Make temporal resolutions equivalent
desiredSize = max(bssfpSize, pcSize);
desiredFrames = max(bssfpFrames, pcFrames);

bssfp = imresize3(bssfp,[desiredSize desiredSize desiredFrames]);
pc = imresize3(pc,[desiredSize desiredSize desiredFrames]);
mag = imresize3(mag,[desiredSize desiredSize desiredFrames]);
temp_mask = zeros(num_roi, desiredSize, desiredSize, desiredFrames);
temp_area = zeros(num_roi, desiredFrames);
for i=1:num_roi
    temp_mask(i, :,:,:) = imresize3(squeeze(mask(i, :,:,:)), [desiredSize desiredSize desiredFrames]);
    temp_area(i, :) = interp1((1:bssfpFrames),squeeze(area(i,:)),linspace(1,bssfpFrames,desiredFrames), 'linear');
end
mask = temp_mask;
area = temp_area;

pcmrTimes = interp1((1:pcFrames),pcmrTimes,linspace(1,pcFrames,desiredFrames),'linear');
bssfpTimes = interp1((1:bssfpFrames),bssfpTimes,linspace(1,bssfpFrames,desiredFrames),'linear');


%% Shift the data due to periph. gating lag + prosp. gating lags
%bssfp = circshift(bssfp,15,3);
%pc = circshift(pc,15,3);
%mag = circshift(mag,15,3);


%% Calculating Velocity and Flow
for j=1:num_roi
    for i=1:desiredFrames
        ROI = imbinarize(mask(j,:,:,i));
        vTemp = pc(:,:,i)*0.1; %single frame velocities (cm/s)
        ROIindex = double(vTemp(ROI)); %indexed velocities within mask
        meanV = mean(ROIindex); %mean velocity in frame i (cm/s)
        flow(j, i) = area(j, i).*meanV; %flow in frame i (cm^3/s)
    end
end
stroke_volume = sum(flow, 2);
save(fullfile(analysis_dir, 'flows.mat'), 'flow');
disp('Flow Data Saved')

for j=1:num_roi
    %% Compute PWV
    if mean(flow(j,:))<0
        flow(j, :) = -1*flow(j, :);
    end      
    maxTime = max(max(pcmrTimes),max(bssfpTimes));
    area_int = interp1(bssfpTimes,area(j, :),linspace(1,maxTime,maxTime));
    
    %Define Systolic Region (flow vs. time)
    figure; plot(flow(j, :)); 
    title(sprintf('%s: ROI %d - Flow', plane, j)); xlabel('Time Frame'); ylabel('Flow (cm^3/s)');
    disp("Draw around the upslope portion of the curve")
    flow_plot = gcf;
    free = drawfreehand;
    systolePts = find(inpolygon(linspace(1,length(flow(j, :)),length(flow(j, :))),flow(j, :),free.Position(:,1),free.Position(:,2)));
    systoleTimes = pcmrTimes(systolePts);
    saveas(flow_plot, fullfile(analysis_dir, sprintf('flowplot_ROI_%d', j)));
    delete(free); close(flow_plot); clear flow_plot;
    
    %Define Linear QA Region (flow vs. area)
    x = area_int(systoleTimes); %time adjustment (fixing delays)
    y = flow(j, systolePts);
    figure; scatter(x,y); 
    title(sprintf('%s: ROI %d - QA Plot', plane, j)); xlabel('Area (cm^2)'); ylabel('Flow (cm^3/s)');
    free2 = drawfreehand;
    in = inpolygon(x,y,free2.Position(:,1),free2.Position(:,2));
    systolePts2 = [x(in); y(in)]';
    [coef,stats] = polyfit(systolePts2(:,1),systolePts2(:,2),1);
    
    %Linear regression
    minArea = min(systolePts2(:,1));
    maxArea = max(systolePts2(:,1));
    xq = linspace(minArea-0.1,maxArea+0.1,100);
    yq = coef(1)*xq + coef(2); %y = mx+b
    delete(free2)
    hold on; 
    scatter(systolePts2(:,1),systolePts2(:,2),72,'k','x');
    plot(xq,yq); 
    str = ['Slope = ' num2str(coef(1)) '\newlineIntcpt = ' num2str(coef(2))];
    text(min(xq),max(yq)-0.1*max(yq),str);
    hold off;
    QAplot = gcf; 
    DescAo_PWV = coef(1);
    fprintf('ROI %d PWV_QA = %.4f m/s\n', j, coef(1)*0.01);
    fprintf('ROI %d Stroke Volume = %.4f mL\n', j, stroke_volume(j)*0.1)
    saveas(QAplot,fullfile(analysis_dir, sprintf('QAplot_ROI_%d', j)));
    save(fullfile(analysis_dir, sprintf('coef_ROI_%d.mat', j)), 'coef');
    save(fullfile(analysis_dir, sprintf('stats_ROI_%d.mat', j)), 'stats');
end
disp('PWV-QA Data Saved!')
disp("(•_•)    ( •_•)>⌐■-■     (⌐■_■)")


%% Helper functions
function [BW, centers, radii, radrange, sens] = circleFinder(image, num_roi, radrange, sens)
    BW = false(size(image,1),size(image,2));
    [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
    flag = 0;
    while ~flag
        g = figure;
        imshow(image, []);
        [centers,radii,~] = imfindcircles(image,radrange,'ObjectPolarity','bright','Sensitivity',sens);
        for n = 1:num_roi
            try
                BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
            catch error
                
            end
        end
        viscircles(centers, radii);
        options.WindowStyle = 'normal';
        answer = inputdlg({"Enter lower radius range:", "Enter upper radius range:", "Enter sensitivity:", "Accept? (enter 1)"}, ...
            "Circle Estimation Parameters", [1 30; 1 30; 1 30; 1 30;], {num2str(radrange(1)), num2str(radrange(2)), ...
            num2str(sens), num2str(flag)}, options);
        close(g)
        radrange = [str2num(answer{1}) str2num(answer{2})];
        sens = str2double(answer{3});
        flag = str2num(answer{4});
        if flag
            if size(centers, 1) ~= num_roi
                fprintf("Error: Number of ROIs is greater than or less than %d.\n" + ...
                    "Please continue adjusting parameters", num_roi);
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