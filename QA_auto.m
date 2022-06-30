clear all; close all; clc;

%home_dir = 'C:\Users\naren\Desktop\Aorta_Segmentation\';
roi = 'AbdAo';
%options: 'kmeans','otsu','hough','hough+contours','edge+hough+contours,'circle+contours'
segType = 'hough+contours'; 
%options: 0=automatic; 1=manual fine-tune
manualFlag = 1; 

home_dir = 'D:\PWV\PWV_data\Volunteers\Grant_dQdA_Testing_110521';
cd(home_dir)
analysis_dir = [home_dir '\2DPC_QA_Analysis\' roi '_' segType];
if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir);
end 

%% Load BSSFP
bssfp_folder = uigetdir(home_dir,'Select AAo or AbdAo FIESTA Scan');
cd(bssfp_folder)

if ~exist('bssfp.mat','file')
    dirFiles = dir('*.dcm');
    if isempty(dirFiles) %check if dicoms are zipped
        zipFiles = dir('*.tgz');
        if isempty(zipFiles)
            disp('NO READABLE DICOMS OR TGZ FILE FOUND, TRY ANOTHER FOLDER');
        else
            gunzip('*.tgz', 'Dicoms'); %unzip first
            cd('Dicoms') %move to unzipped folder
            dd = dir(); %get the name of the only file in the new dir
            untar(dd(3).name,'Dicoms'); %untar that file
            movefile('Dicoms/*','..'); %move unzipped files back up
            cd('..') %move up a directory
            rmdir('Dicoms','s') %get rid of created dummy unzipping folder
        end 
    end     
    bssfpInfo = dicominfo(dirFiles(1).name); %grab dicom info from 1st frame
    bssfpSize = double(bssfpInfo.Height); %assumes matrix size is equal in x and y
    bssfpFrames = bssfpInfo.CardiacNumberOfImages;
    bssfp = zeros(bssfpSize,bssfpSize,bssfpFrames);
    for f=1:length(dirFiles)
        bssfp(:,:,f) = dicomread(dirFiles(f).name);
        temp = dicominfo(dirFiles(f).name);
        bssfpTimes(f) = temp.TriggerTime;
    end 
    save('bssfp.mat','bssfp');
    save('bssfpTimes.mat','bssfpTimes');
    save('bssfpInfo.mat','bssfpInfo');
else
    load('bssfp.mat');
    load('bssfpTimes.mat');
    load('bssfpInfo.mat');
end 
bssfpSize = double(bssfpInfo.Height); %assumes matrix size is equal in x and y
bssfpFrames = bssfpInfo.CardiacNumberOfImages;
pixelArea = bssfpInfo.PixelSpacing(1).*bssfpInfo.PixelSpacing(2)*0.01; %cm^2


%% Load 2DPC
pc_folder = uigetdir(home_dir,'Select folder 2DPC scan');
cd(pc_folder)

if ~exist('pc.mat','file')
    dirFiles = dir('*.dcm');
    if isempty(dirFiles) %check if dicoms are zipped
        zipFiles = dir('*.tgz');
        if isempty(zipFiles)
            disp('NO READABLE DICOMS OR TGZ FILE FOUND, TRY ANOTHER FOLDER');
        else
            gunzip('*.tgz', 'Dicoms'); %unzip first
            cd('Dicoms') %move to unzipped folder
            dd = dir(); %get the name of the only file in the new dir
            untar(dd(3).name,'Dicoms'); %untar that file
            movefile('Dicoms/*','..'); %move unzipped files back up
            cd('..') %move up a directory
            rmdir('Dicoms','s') %get rid of created dummy unzipping folder
        end 
    end 
    
    pcInfo = dicominfo(dirFiles(1).name); %grab dicom info from 1st frame
    pcSize = double(pcInfo.Height); %assumes matrix size is equal in x and y
    pcFrames = pcInfo.CardiacNumberOfImages;
    pcmr = zeros(pcSize,pcSize,pcFrames);
    for f=1:length(dirFiles)
        pcmr(:,:,f) = dicomread(dirFiles(f).name);
        temp = dicominfo(dirFiles(f).name);
        pcmrTimes(f) = temp.TriggerTime;
    end
    pc = pcmr(:,:,1:pcFrames);
    mag = pcmr(:,:,(pcFrames+1):end);
    pcmrTimes = pcmrTimes(1:pcFrames);
    save('pc.mat','pc');
    save('mag.mat','mag');
    save('pcmrTimes.mat','pcmrTimes');
    save('pcInfo.mat','pcInfo');
else
    load('pc.mat');
    load('mag.mat','mag');
    load('pcmrTimes.mat')
    load('pcInfo.mat');
end 
pcSize = double(pcInfo.Height); %assumes matrix size is equal in x and y
pcFrames = pcInfo.CardiacNumberOfImages;


%% Make temporal resolutions equivalent
desiredSize = max(bssfpSize,pcSize);
desiredFrames = max(bssfpFrames,pcFrames);

bssfp = imresize3(bssfp,[desiredSize desiredSize desiredFrames]);
pc  = imresize3(pc,[desiredSize desiredSize desiredFrames]);
mag = imresize3(mag,[desiredSize desiredSize desiredFrames]);
pcmrTimes = interp1((1:pcFrames),pcmrTimes,linspace(1,pcFrames,desiredFrames),'linear');
bssfpTimes = interp1((1:bssfpFrames),bssfpTimes,linspace(1,bssfpFrames,desiredFrames),'linear');

%% Shift the data due to periph. gating lag + prosp. gating lags
%bssfp = circshift(bssfp,15,3);
%pc = circshift(pc,15,3);
%mag = circshift(mag,15,3);
MAG = mean(mag,3);

%% Initialize Masks
cd(analysis_dir)
if ~exist('mask','file')
    mask = zeros(desiredSize,desiredSize,desiredFrames);
    needsMask = 1;
else
    load('mask.mat');
    needsMask = 0;
end 

%% Zoom in Manually
figure; imshow(MAG,[]); title('ZOOM');
rect = drawrectangle;
crop = round([rect.Position(2), rect.Position(2)+rect.Position(4), rect.Position(1), rect.Position(1)+rect.Position(3)]);
close all;

%% Segment Each Frame Separately
if needsMask %if we haven't loaded in the mask...
    for i=1:desiredFrames
        image = bssfp(crop(1):crop(2),crop(3):crop(4),i);
        image = rescale(image);
        switch segType
            case 'kmeans' % Kmeans Clustering
                s = rng; rng('default');
                L = imsegkmeans(single(image),2,'NumAttempts',5);
                rng(s);
                labelOfInterest = L(round(size(image,1)/2),round(size(image,2)/2));
                BW = L == labelOfInterest;
                if sum(BW(:))==0
                    level = graythresh(image);
                    BW = imbinarize(image,level);
                end 
            case 'otsu' % Otsu's Automatic Threshold
                level = graythresh(image);
                BW = imbinarize(image,level);
            case 'hough' % Hough Transform + Active Countours
                radrange = [10 100];
                sens = 0.9;
                [centers,radii,~] = imfindcircles(image,[10 100],'ObjectPolarity','bright','Sensitivity',0.9);
                BW = false(size(image,1),size(image,2));
                [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
                BW = BW | (hypot(Xgrid-centers(1,1),Ygrid-centers(1,2)) <= radii(1));   
            case 'hough+contours' % Hough Transform + Active Countours
                radrange = [10 100];
                sens = 0.9;
                [centers,radii,~] = imfindcircles(image,[10 100],'ObjectPolarity','bright','Sensitivity',0.9);
                BW = false(size(image,1),size(image,2));
                [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
                BW = BW | (hypot(Xgrid-centers(1,1),Ygrid-centers(1,2)) <= radii(1));   
                BW = imdilate(BW,strel('disk',2)); %dilate so active contours pulls in segmentation
                iters = 100;
                smF = 5;
                contrF = 0.1; %bias towards shrinking
                method = 'Chan-Vese'; %or Chan-Vese (default)
                BW = activecontour(image,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
            case 'edge+hough+contours' % Hough Transform + Active Countours
                grad = rescale(imgradient(image));
                radrange = [10 100];
                sens = 0.9;
                [centers,radii,~] = imfindcircles(grad,[10 100],'ObjectPolarity','bright','Sensitivity',0.9);
                BW = false(size(image,1),size(image,2));
                [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
                BW = BW | (hypot(Xgrid-centers(1,1),Ygrid-centers(1,2)) <= radii(1));   
                BW = imerode(BW,strel('disk',2)); %erode so active contours pushes out
                iters = 300;
                smF = 5;
                contrF = -0.1; %bias towards growing
                method = 'Chan-Vese'; %or Chan-Vese (default)
                BW = activecontour(grad,BW,iters,method,'SmoothFactor',smF,'ContractionBias',contrF);
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
        end 
        BW = bwselect(BW,round(size(image,1)/2),round(size(image,2)/2),8);
        BW = imdilate(BW,strel('disk',1));
        
        if manualFlag
            % Freehand ROI conversion
            blocations = bwboundaries(BW,'noholes');
            pos = blocations{1}; %convert to x,y order.
            subx = int16(linspace(1,length(pos(:,1)),15)); %subsample positions (x)
            suby = int16(linspace(1,length(pos(:,2)),15)); %subsample positions (y)
            upscaleSize = 45;
            posSpline = interparc(upscaleSize,pos(subx,1),pos(suby,2),'csape'); %closed spline

            posSpline = fliplr(posSpline);
            % Can 'ESC' to delete ROI, Right-click to Add Waypoint
            waypoints = zeros(upscaleSize,1); %initialize index locations of waypoints
            waypoints(1:5:end) = 1; %set 8 evenly-spaced waypoints
            f = figure; f.WindowState = 'maximized';
            imshow(image,[]);
            h = drawfreehand('Position', posSpline,'FaceAlpha',0.15,'LineWidth',1, ...
               'Multiclick',true,'Waypoints',logical(waypoints)); %create freehand ROI
            customWait(h); %see custom function below in 'helper functions'

            hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
            editedMask = false(size(image));
            for ind = 1:numel(hfhs)
               editedMask = editedMask | hfhs(ind).createMask(); %accumulate mask from each ROI
               boundaryLocation = round(hfhs(ind).Position); %include ROI boundary
               bInds = sub2ind(size(image), boundaryLocation(:,2), boundaryLocation(:,1));
               editedMask(bInds) = true;
            end
        else
            f = figure; f.WindowState = 'maximized';
            imshow(image,[]);
            hold on; visboundaries(BW);
        end 
                
        % Separate Vessel
        BW = padarray(BW, [crop(1) crop(3)], 'pre');
        BW = padarray(BW, [desiredSize-crop(2)-1 desiredSize-crop(4)-1], 'post');
        mask(:,:,i) = BW;
        area(i) = sum(BW(:))*pixelArea; %cm^2
    end
    save('mask.mat','mask');
    save('area.mat','area');
    disp('Mask Data Saved');
    close all
end 

%% Calculating Velocity and Flow
for i=1:desiredFrames
    ROI = imbinarize(mask(:,:,i));
    vTemp = pc(:,:,i)*0.1; %single frame velocities (cm/s)
    ROIindex = double(vTemp(ROI)); %indexed velocities within mask
    meanV = mean(ROIindex); %mean velocity in frame i (cm/s)
    flow(i) = area(i).*meanV; %flow in frame i (cm^3/s)
end
save('flow.mat','flow');
disp('Flow and Area Data Saved')


%% Compute PWV
if mean(flow)<0
    flow = -1*flow;
end      
maxTime = max(max(pcmrTimes),max(bssfpTimes));
area_int = interp1(bssfpTimes,area,linspace(1,maxTime,maxTime));

%Define Systolic Region (flow vs. time)
figure; plot(flow); 
title([roi ' - Flow']); xlabel('Frame'); ylabel('Flow (cm^3/s)');
flow_plot = gcf;
free = drawfreehand;
systolePts = find(inpolygon(linspace(1,length(flow),length(flow)),flow,free.Position(:,1),free.Position(:,2)));
systoleTimes = pcmrTimes(systolePts);
saveas(flow_plot,'flowplot');
delete(free); close(flow_plot); clear flow_plot;

%Define Linear QA Region (flow vs. area)
x = area_int(systoleTimes); %time adjustment (fixing delays)
y = flow(systolePts);
figure; scatter(x,y); 
title([roi ' - QA Plot']); xlabel('Area (cm^2)'); ylabel('Flow (cm^3/s)');
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
disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
saveas(QAplot,'QAplot');
save('coef.mat','coef');
save('stats.mat','stats');
disp('QA-PWV Data Saved')


%% Helper functions
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