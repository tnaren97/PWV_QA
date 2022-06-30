clear all; close all; clc;

roi = 'AscAo';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Tarun_dQdA_Testing_110521\011_PWV_CartBH_AAo_PROSP_4VPS_47FRAMES_3p0ACC';
home_dir = 'D:\PWV\PWV_data\Volunteers\Tarun_dQdA_Testing_110521\012_PWV_CartBH_AAo_PROSP_2VPS_81FRAMES_3p0ACC';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Tarun_dQdA_Testing_110521\013_PWV_CartBH_AAo_NOPROSP_1VPS_70FRAMES_3p5ACC';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Grant_dQdA_Testing_110521\008_PWV_CartBH_AAo_PROSP_4VPS_51FRAMES';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Grant_dQdA_Testing_110521\009_PWV_CartBH_AAo_NOPROSP_1VPS_80FRAMES';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Grant_v2_dQdA_Testing_110821\005_PWV_CartBH_AAo_4VPS_33FRAMES_ECG';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Grant_v2_dQdA_Testing_110821\006_PWV_CartBH_AAo_2VPS_66FRAMES_ECG';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Grant_v2_dQdA_Testing_110821\007_PWV_CartBH_AAo_1VPS_88FRAMES_ECG';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Oliver_dQdA_Testing_110821\007_PWV_CartBH_AAo_4VPS_57FRAMES_ECG';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Oliver_dQdA_Testing_110821\009_PWV_CartBH_AAo_2VPS_87FRAMES_ECG';
%home_dir = 'D:\PWV\PWV_data\Volunteers\Oliver_dQdA_Testing_110821\008_PWV_CartBH_AAo_1VPS_109FRAMES_ECG';

cd(home_dir)
analysis_dir = [home_dir '\2DPC_QA_Analysis\' roi];
if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir);
end 

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
pixelArea = pcInfo.PixelSpacing(1).*pcInfo.PixelSpacing(2)*0.01; %cm^2
VPS = pcInfo.Private_0021_1035;
TR = pcInfo.RepetitionTime;
acquiredTR = VPS*TR*2; %2-point encoding
disp(['Acquired TR = ' num2str(acquiredTR)]);
MAG = mean(mag,3);

%% Initialize Masks
cd(analysis_dir)
if ~exist('masked.mat','file')
    masked = zeros(pcSize,pcSize,pcFrames);
    needsMask = 1;
else
    load('masked.mat');
    needsMask = 0;
end 

%% Zoom in Manually
figure; imshow(MAG,[]); title('ZOOM');
rect = drawrectangle;
crop = round([rect.Position(2), rect.Position(2)+rect.Position(4), rect.Position(1), rect.Position(1)+rect.Position(3)]);
close all;

%% Draw Rough ROI to find systole
figure; imshow(MAG,[0 0.6*max(MAG(:))]); title('DRAW ROUGH ROI');
circle = drawcircle();
radius = circle.Radius; %get radius of circle
center = round(circle.Center); %get center coordinates

[X,Y] = ndgrid(1:size(pc,1),1:size(pc,2));
X = X-center(2); %shift coordinate grid
Y = Y-center(1);
roiMask = sqrt(X.^2 + Y.^2)<=radius; %anything outside radius is ignored
AREA = sum(roiMask(:))*pixelArea; %ROI area (cm^2)
for i=1:size(pc,3)
    vTemp = pc(:,:,i)*0.1; %through-plane velocity in frame i
    roiDataRaw(:,i) = double(vTemp(roiMask)); %indexed velocities within mask
    meanROI(i) = mean(roiDataRaw(:,i)); %mean velocity in frame i (mm/s)
    flowROI(i) = AREA.*meanROI(i); %flow in frame i (cm^3/s = mL/s)
end 
close;

if mean(flowROI)<0
    flowROI = -1*flowROI;
end 
figure; plot(flowROI);
pts = ginput(1);
lateSystole = round(pts(1,1));
close;

shift = 2;
allPoints = 1:pcFrames;
shiftPoints = circshift(allPoints,shift);
numPoints = lateSystole + shift*2;

%% Segment Each Frame Separately
area = zeros(1,pcFrames);
if needsMask %if we haven't loaded in the mask...
    for i=1:40
        image = mag(crop(1):crop(2),crop(3):crop(4),i);
        image = rescale(image);
        f = figure; f.WindowState = 'maximized';
        imshow(image, []);
        
        poly = drawpolygon;
        pos = poly.Position;
        upscaleSize = 30;
        posSpline = interparc(upscaleSize,[pos(:,1); pos(1,1)],[pos(:,2); pos(1,2)],'csape'); %closed spline
        %posSpline = fliplr(posSpline);
        % Can 'ESC' to delete ROI, Right-click to Add Waypoint
        waypoints = zeros(upscaleSize,1); %initialize index locations of waypoints
        waypoints(1:5:end) = 1; %set 8 evenly-spaced waypoints
        delete(poly); clear poly;
        h = drawfreehand('Position', posSpline,'FaceAlpha',0.15,'LineWidth',1, ...
            'Multiclick',true,'Waypoints',logical(waypoints)); %create freehand ROI
        customWait(h); %see custom function below in 'helper functions'
        blocations = posSpline;

        hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
        editedMask = false(size(image));
        for ind = 1:numel(hfhs)
           editedMask = editedMask | hfhs(ind).createMask(); %accumulate mask from each ROI
           boundaryLocation = round(hfhs(ind).Position); %include ROI boundary
           bInds = sub2ind(size(image), boundaryLocation(:,2), boundaryLocation(:,1));
           editedMask(bInds) = true;
        end
        area1 = polyarea(posSpline(:,1),posSpline(:,2));
        area(:,i) = area1.*pixelArea;
        BW = padarray(editedMask, [crop(1) crop(3)], 'pre');
        BW = padarray(BW, [pcSize-crop(2)-1 pcSize-crop(4)-1], 'post');
        masked(:,:,i) = BW;
        
    end
    save('masked.mat','masked');
    save('area.mat','area');
    disp('Mask Data Saved');
    close all
end 

%% Calculating Velocity and Flow
for i=1:pcFrames
    ROI = imbinarize(masked(:,:,i));
    vTemp = pc(:,:,i)*0.1; %single frame velocities (cm/s)
    ROIindex = double(vTemp(ROI)); %indexed velocities within mask
    meanV = mean(ROIindex); %mean velocity in frame i (cm/s)
    flow(i) = area(i).*meanV; %flow in frame i (cm^3/s)
end
save('flow.mat','flow');
disp('Flow and Area Data Saved')


%% Compute PWV
if nanmean(flow)<0
    flow = -1*flow;
end      

%Define Systolic Region (flow vs. time)
% figure; plot(flow); 
% title([roi ' - Flow']); xlabel('Frame'); ylabel('Flow (cm^3/s)');
% flow_plot = gcf;
% free = drawfreehand;
% systolePts = find(inpolygon(linspace(1,length(flow),length(flow)),flow,free.Position(:,1),free.Position(:,2)));
% saveas(flow_plot,'flowplot');
% delete(free); close(flow_plot); clear flow_plot;

%Define Linear QA Region (flow vs. area)
pointsPadded = shiftPoints(1:numPoints);
%x = circshift(area(pointsPadded),5);
x = area(pointsPadded);
y = flow(pointsPadded);

figure; scatter(x,y,[],1:length(x)); colorbar; 
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