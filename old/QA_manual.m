%clear all; close all; clc;

%home_dir = 'C:\Users\naren\Desktop\Aorta_Segmentation\';
home_dir = 'D:\PWV\PWV_data\life_life00016_11110_2020-03-09-h08\dicoms';
cd(home_dir)

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
bssfp = circshift(bssfp,15,3);
pc = circshift(pc,15,3);
mag = circshift(mag,15,3);

%% Initialize Masks
if contains(pc_folder,'AAo','IgnoreCase',true)
    numROIs = 2; 
    if ~exist('AscAo_mask.mat','file')
        needsMask = 1;
        AscAo_mask = zeros(desiredSize,desiredSize,desiredFrames);
    else
        needsMask = 0;    
        load('AscAo_mask.mat');
    end
    if ~exist('DescAo_mask.mat','file')
        DescAo_mask = zeros(desiredSize,desiredSize,desiredFrames);
    else
        load('DescAo_mask.mat');
    end
else
    numROIs = 1;
    if ~exist('AbdAo_mask.mat','file')
        needsMask = 1;
        AbdAo_mask = zeros(desiredSize,desiredSize,desiredFrames);
    else
        needsMask = 0;
        load('AbdAo_mask.mat');
    end
end

%% Draw Rough ROI to find systole
MAG = mean(mag,3);
% figure; imshow(MAG,[0 0.6*max(MAG(:))]); title('DRAW ROUGH ROI');
% circle = drawcircle();
% radius = circle.Radius; %get radius of circle
% center = round(circle.Center); %get center coordinates
% 
% [X,Y] = ndgrid(1:size(pc,1),1:size(pc,2));
% X = X-center(2); %shift coordinate grid
% Y = Y-center(1);
% roiMask = sqrt(X.^2 + Y.^2)<=radius; %anything outside radius is ignored
% area = sum(roiMask(:))*pixelArea; %ROI area (cm^2)
% for i=1:size(pc,3)
%     vTemp = pc(:,:,i)*0.1; %through-plane velocity in frame i
%     roiDataRaw(:,i) = double(vTemp(roiMask)); %indexed velocities within mask
%     meanROI(i) = mean(roiDataRaw(:,i)); %mean velocity in frame i (mm/s)
%     flowROI(i) = area.*meanROI(i); %flow in frame i (cm^3/s = mL/s)
% end 
% close;
% 
% if mean(flowROI)<0
%     flowROI = -1*flowROI;
% end 
% figure; plot(flowROI);
% pts = ginput(2);
% pts = round(pts);
% earlySystole = pts(1,1);
% lateSystole = pts(2,1);
% close;

%% Zoom in Manually
f=figure; imshow(MAG,[]); title('ZOOM');
rect = drawrectangle;
crop = round([rect.Position(2), rect.Position(2)+rect.Position(4), rect.Position(1), rect.Position(1)+rect.Position(3)]);
close all;

%% Segment Each Frame Separately
if needsMask %if we haven't loaded in the mask...
    for i=1:desiredFrames
        image = bssfp(crop(1):crop(2),crop(3):crop(4),i);
        f = figure; f.WindowState = 'maximized';
        imshow(image, []);

        for ind = 1:numROIs
            poly = drawpolygon;
            pos = poly.Position;
            upscaleSize = 35;
            posSpline = interparc(upscaleSize,[pos(:,1); pos(1,1)],[pos(:,2); pos(1,2)],'csape'); %closed spline
            %posSpline = fliplr(posSpline);
            % Can 'ESC' to delete ROI, Right-click to Add Waypoint
            waypoints = zeros(upscaleSize,1); %initialize index locations of waypoints
            waypoints(1:7:end) = 1; %set 8 evenly-spaced waypoints
            delete(poly); clear poly;
            h = drawfreehand('Position', posSpline,'FaceAlpha',0.15,'LineWidth',1, ...
                'Multiclick',true,'Waypoints',logical(waypoints)); %create freehand ROI
            customWait(h); %see custom function below in 'helper functions'
            blocations(:,:,ind) = posSpline;
        end

        hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
        editedMask = false(size(image));
        for ind = 1:numel(hfhs)
            editedMask = editedMask | hfhs(ind).createMask(); %accumulate mask from each ROI
            boundaryLocation = round(hfhs(ind).Position); %include ROI boundary
            bInds = sub2ind(size(image), boundaryLocation(:,2), boundaryLocation(:,1));
            editedMask(bInds) = true;
        end

        %% Separate Ascending and Descending Aorta
        [r,c] = find(editedMask);
        [~,idx_min] = min(r);
        aorta1 = bwselect(editedMask,c(idx_min),r(idx_min),8);
        [~,idx_roi] = min(min(blocations(:,2,:))); %roi with min rows
        area1 = polyarea(blocations(:,1,idx_roi),blocations(:,2,idx_roi));
        aorta1 = padarray(aorta1, [crop(1) crop(3)], 'pre');
        aorta1 = padarray(aorta1, [desiredSize-crop(2)-1 desiredSize-crop(4)-1], 'post');
        if numROIs==2
            [~,idx_max] = max(r);
            aorta2 = bwselect(editedMask,c(idx_max),r(idx_max),8);
            [~,idx_roi] = max(max(blocations(:,2,:))); %roi with min rows
            area2 = polyarea(blocations(:,1,idx_roi),blocations(:,2,idx_roi));
            aorta2 = padarray(aorta2, [crop(1) crop(3)], 'pre');
            aorta2 = padarray(aorta2, [desiredSize-crop(2)-1 desiredSize-crop(4)-1], 'post');
        end

        if numROIs==2
            AscAo_mask(:,:,i) = aorta1;
            AscAo_area(i) = area1*pixelArea; %cm^2
            DescAo_mask(:,:,i) = aorta2;
            DescAo_area(i) = area2*pixelArea; %cm^2
        else
            AbdAo_mask(:,:,i) = aorta1; 
            AbdAo_area(i) = area1*pixelArea; %cm^2
        end
    end
end 

if numROIs==2
    save('AscAo_mask.mat','AscAo_mask')    
    save('DescAo_mask.mat','DescAo_mask')
else    
    save('AbdAo_mask.mat','AbdAo_mask')
end
close all

disp("Masks saved.");

%% Calculating Velocity and Flow

% Separate ascending and descending aorta and calculate area and flow
for i=1:desiredFrames
    if numROIs==2     
        % Calculate Ascending Aorta flow
        Asc = imbinarize(AscAo_mask(:,:,i));
        vTemp = pc(:,:,i)*0.1; %single frame velocities (cm/s)
        AscAo_ROI = double(vTemp(Asc)); %indexed velocities within mask
        AscAo_meanV = mean(AscAo_ROI); %mean velocity in ROI in frame i (cm/s)
        AscAo_flow(i) = AscAo_area(i).*AscAo_meanV; %flow in frame i (cm^3/s)
        
        % Calculate Descending Aorta flow
        Desc = imbinarize(DescAo_mask(:,:,i));
        Desc_ROI = double(vTemp(Desc)); %indexed velocities within mask
        Desc_meanV = mean(Desc_ROI); %mean velocity in frame i (cm/s)
        DescAo_flow(i) = DescAo_area(i).*Desc_meanV; %flow in frame i (cm^3/s)
    else
        % Calculate Abdominal Aorta flow
        Abd = imbinarize(AbdAo_mask(:,:,i));
        vTemp = pc(:,:,i)*0.1; %single frame velocities (cm/s)
        AbdAo_ROI = double(vTemp(Abd)); %indexed velocities within mask
        AbdAo_meanV = mean(AbdAo_ROI); %mean velocity in frame i (cm/s)
        AbdAo_flow(i) = AbdAo_area(i).*AbdAo_meanV; %flow in frame i (cm^3/s)
    end
end

if ~exist('2DPC_QA_Analysis','dir')
    mkdir('2DPC_QA_Analysis');
end 

if numROIs==2
    movefile('AscAo_mask.mat','2DPC_QA_Analysis');
    movefile('DescAo_mask.mat','2DPC_QA_Analysis');
    cd('2DPC_QA_Analysis');
    save('AscAo_flow.mat','AscAo_flow');
    save('AscAo_area.mat','AscAo_area');
    save('DescAo_flow.mat','DescAo_flow');
    save('DescAo_area.mat','DescAo_area');
else   
    movefile('AbdAo_mask.mat','2DPC_QA_Analysis');
    cd('2DPC_QA_Analysis');
    save('AbdAo_flow.mat','AbdAo_flow');
    save('AbdAo_area.mat','AbdAo_area');
end

disp("Flow and Area data saved.")


%% Compute PWV


if numROIs==2
    %Ascending Aorta
    if mean(AscAo_flow)<0
        AscAo_flow = -1*AscAo_flow;
    end 
    AscAo_area = interp1(bssfpTimes,AscAo_area,(1:maxTime),'linear');
    AscAo_flow = interp1(pcmrTimes,AscAo_flow,(1:maxTime),'linear');
    
    figure; 
    scatter(AscAo_area,AscAo_flow,[],(1:length(AscAo_area)));
    c = colorbar;
    c.Label.String = 'Cardiac Frame #';
    title('Ascending Aorta - QA'); 
    xlabel('Area (cm^2)'); 
    ylabel('Flow (cm^3/s)');
    free = drawfreehand;
    in = inpolygon(AscAo_area,AscAo_flow,free.Position(:,1),free.Position(:,2));
    systolePts = [AscAo_area(in); AscAo_flow(in)]';
    [coef,stats] = polyfit(systolePts(:,1),systolePts(:,2),1);
    minArea = min(systolePts(:,1));
    maxArea = max(systolePts(:,1));
    xq = linspace(minArea-0.1,maxArea+0.1,100);
    yq = coef(1)*xq + coef(2); %y = mx+b
    delete(free)
    hold on; 
    scatter(systolePts(:,1),systolePts(:,2),72,'k','x');
    plot(xq,yq); 
    str = ['Slope = ' num2str(coef(1)) '\newlineIntcpt = ' num2str(coef(2))];
    text(max(AscAo_area)+0.1,max(AscAo_flow)+0.1,str);
    hold off;
    AscAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    
    mkdir('AscAo');
    cd('AscAo');
    saveas(gcf,'AscAo_QAplot');
    save('AscAo_PWV.mat','AscAo_PWV');
    save('coef.mat','coef');
    save('stats.mat','stats');
    cd ..
    
    %Descending Aorta
    if mean(DescAo_flow)<0
        DescAo_flow = -1*DescAo_flow;
    end     
    DescAo_area = interp1(bssfpTimes,DescAo_area,(1:maxTime),'linear');
    DescAo_flow = interp1(pcmrTimes,DescAo_flow,(1:maxTime),'linear');
    
    figure; 
    scatter(DescAo_area,DescAo_flow,[],(1:length(DescAo_area)));
    c = colorbar;
    c.Label.String = 'Cardiac Frame #';
    title('Descending Aorta - QA'); 
    xlabel('Area (cm^2)'); 
    ylabel('Flow (cm^3/s)');
    free = drawfreehand;
    in = inpolygon(DescAo_area,DescAo_flow,free.Position(:,1),free.Position(:,2));
    systolePts = [DescAo_area(in); DescAo_flow(in)]';
    [coef,stats] = polyfit(systolePts(:,1),systolePts(:,2),1);
    minArea = min(systolePts(:,1));
    maxArea = max(systolePts(:,1));
    xq = linspace(minArea-0.1,maxArea+0.1,100);
    yq = coef(1)*xq + coef(2); %y = mx+b
    delete(free)
    hold on; 
    scatter(systolePts(:,1),systolePts(:,2),72,'k','x');
    plot(xq,yq); 
    str = ['Slope = ' num2str(coef(1)) '\newlineIntcpt = ' num2str(coef(2))];
    text(max(DescAo_area)+0.1,max(DescAo_flow)+0.1,str);
    hold off;
    DescAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    
    mkdir('DescAo');
    cd('DescAo');
    saveas(gcf,'DescAo_QAplot');
    save('DescAo_PWV.mat','DescAo_PWV');
    save('coef.mat','coef');
    save('stats.mat','stats');
    cd ..
else
    %Abdominal Aorta
    if mean(AbdAo_flow)<0
        AbdAo_flow = -1*AbdAo_flow;
    end      
    AbdAo_area_int = interp1(bssfpTimes,AbdAo_area,linspace(1,maxTime,maxTime));
    
    figure; 
    plot(AbdAo_flow); 
    title('Abdominal Aorta - Flow'); 
    xlabel('Frame'); 
    ylabel('Flow (cm^3/s)');
    flow_plot = gcf;
    free = drawfreehand;
    systolePts = find(inpolygon(linspace(1,length(AbdAo_flow),length(AbdAo_flow)),AbdAo_flow,free.Position(:,1),free.Position(:,2)));
    systoleTimes = pcmrTimes(systolePts);
    figure; 
    scatter(AbdAo_area_int(systoleTimes),AbdAo_flow(systolePts),[],1:length(systolePts)); 
    colormap jet; colorbar
    title('Abdominal Aorta - QA'); 
    xlabel('Delta Area (cm^2)'); 
    ylabel('Delta Flow (cm^3/s)');
    x = AbdAo_area_int(systoleTimes);
    y = AbdAo_flow(systolePts);
    [coef,stats] = polyfit(x,y,1);
    xq = linspace(min(x)-0.1,max(x)+0.1,100);
    yq = coef(1)*xq + coef(2); %y = mx+b
    delete(free)
    hold on; 
    plot(xq,yq); 
    str = ['Slope = ' num2str(coef(1)) '\newlineIntcpt = ' num2str(coef(2))];
    text(min(xq),max(yq)-0.1*max(yq),str);
    hold off;
    QAplot = gcf; 
    
    AbdAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    
    mkdir('AbdAo');
    cd('AbdAo');
    saveas(flow_plot,'AbdAo_flowplot');
    saveas(QAplot,'AbdAo_QAplot');
    save('AbdAo_PWV.mat','AbdAo_PWV');
    save('coef.mat','coef');
    save('stats.mat','stats');
end


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