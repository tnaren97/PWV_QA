% Code assumes follwoing folder structure:
%          home_dir
%             - life_life[patient num]...
%                - [scan num].[name]_[mode]_BH_[scan type]AAo
%                   - .tgz file or DICOMS
%                - [scan num].[name]_[mode]_BH_[scan type]AbdAo
%                   - .tgz file or DICOMS
% Select patient folder directory
clear all
home_dir = 'C:\Users\naren\Desktop\Aorta_Segmentation\';
% home_dir = 'D:\PWV\PWV_data\life_life00030_11178_2020-09-14-h08\dicoms';
cd(home_dir)

% result_dir = [home_dir 'results\'];
bssfp_folder = uigetdir(home_dir,'Select AAo or AbdAo FIESTA Scan');
cd(bssfp_folder)

if ~(exist('bssfp.mat','file')&&exist('bssfp_times.mat','file') )
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
        dirFiles = dir('*.dcm');
    end 
    
    bssfp_info = dicominfo(dirFiles(1).name); %grab dicom info from 1st frame
    bssfpSize = double(bssfp_info.Height); %assumes matrix size is equal in x and y
    bssfpFrames = bssfp_info.CardiacNumberOfImages;
    bssfp = zeros(bssfpSize,bssfpSize,bssfpFrames);
    for f=1:length(dirFiles)
        bssfp(:,:,f) = dicomread(dirFiles(f).name);
        temp = dicominfo(dirFiles(f).name);
        bssfp_times(f) = temp.TriggerTime;
    end
else
    load('bssfp.mat');
    load('bssfp_info.mat');
    load('bssfp_times.mat');
end
%     disp(bssfp_times)
bssfpFrames = bssfp_info.CardiacNumberOfImages;
bssfpSize = double(bssfp_info.Height); %assumes matrix size is equal in x and y
cropSize = 256; %BASED ON CROP SIZES BELOW (currently 256x256 to better visualize aorta)
pixelArea = bssfp_info.PixelSpacing(1).*bssfp_info.PixelSpacing(2)*0.01; %cm^2

if contains(bssfp_folder,'AAo','IgnoreCase',true)
    numROIs = 2; 
    if ~exist('AscAo_mask.mat','file')
        needsMask = 1;
        AscAo_mask = zeros(cropSize,cropSize,bssfpFrames);
    else
        needsMask = 0;    
        load('AscAo_mask.mat');
    end
    if ~exist('DescAo_mask.mat','file')
        DescAo_mask = zeros(cropSize,cropSize,bssfpFrames);
    else
        load('DescAo_mask.mat');
    end
    % adjustable parameters
%     crop = [145 145+cropSize-1 121 121+cropSize-1];
    crop = [105 105+cropSize-1 131 131+cropSize-1];%crop dims (1:2 = row width, 3:4 = col width)
%     crop = [180 180+cropSize-1 190 190+cropSize-1];
    radrange = [6 25]; %radius range for Hough Transform
%     sens = 0.90; %sensitivity of circle find
    sens = 0.95;
else
    numROIs = 1;
    if ~exist('AbdAo_mask.mat','file')
        needsMask = 1;
        AbdAo_mask = zeros(cropSize,cropSize,bssfpFrames);
    else
        needsMask = 0;
        load('AbdAo_mask.mat');
    end
    % adjustable parameters
    crop = [105 360 121 376]; %crop dims (1:2 = row width, 3:4 = col width)
    radrange = [12 23]; %radius range for Hough Transform
    sens = 0.9; %sensitivity of circle find
end

%% Segment Each Frame Separately
if needsMask %if we haven't loaded in the mask...
    for i=1:bssfpFrames
        image = bssfp(crop(1):crop(2),crop(3):crop(4),i);

        %% Find circles (Hough transform)
        image = rescale(image); %normalizes image
%         [centers,radii,~] = imfindcircles(image,radrange,'ObjectPolarity','bright','Sensitivity',sens);
        BW = false(size(image,1),size(image,2));
%         [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
%         for n = 1:numROIs
%             try
%                 BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
%             catch error
%                 if numROIs==2
%                     save('AscAo_mask.mat', 'AscAo_mask')    
%                     save('DescAo_mask.mat', 'DescAo_mask')
%                 else    
%                     save('AbdAo_mask.mat', 'AbdAo_mask')
%                 end
%                 close all
%                 cd(home_dir)
%                 disp(['Frame: ' num2str(i)])
%                 rethrow(error)
%             end
%         end
% 
%         %% Morphological Dilation 
%         BW = imdilate(BW,strel('disk',3)); %dilate so active contours pulls in segmentation
% 
%         %% Active Contours + Dilation
%         BW = activecontour(image,BW,300,'Chan-Vese','SmoothFactor',5,'ContractionBias',0.2);
%         BW = imdilate(BW,strel('disk',2));

%         x1=55;
%         y1=40;
%         radius1 = 22;
%         [xx1,yy1] = ndgrid((1:cropSize)-y1,(1:cropSize)-x1);
%         circ1 = (xx1.^2 + yy1.^2)<radius1^2;
%         x2=95;
%         y2=100;
%         radius2 = 18;
%         [xx2,yy2] = ndgrid((1:cropSize)-y2,(1:cropSize)-x2);
%         circ2 = (xx2.^2 + yy2.^2)<radius2^2;
%         BW = BW + circ1 + circ2;

        x1=90;
        y1=100;
        radius1 = 40;
        [xx1,yy1] = ndgrid((1:cropSize)-y1,(1:cropSize)-x1);
        circ1 = (xx1.^2 + yy1.^2)<radius1^2;
        x2=140;
        y2=180;
        radius2 = 30;
        [xx2,yy2] = ndgrid((1:cropSize)-y2,(1:cropSize)-x2);
        circ2 = (xx2.^2 + yy2.^2)<radius2^2;
        BW = BW + circ1 + circ2;
        
%         %% Freehand ROI conversion
%         blocations = bwboundaries(BW,'noholes');
%         numBlobs = numel(blocations);
%         f = figure;
%         imshow(image, []);
%         f.WindowState = 'maximized';
% 
%         for ind = 1:numBlobs
%             pos = blocations{ind}; %convert to x,y order.
%             subx = int16(linspace(1,length(pos(:,1)),7)); %subsample positions (x)
%             suby = int16(linspace(1,length(pos(:,2)),7)); %subsample positions (y)
%             posSpline = interparc(90,pos(subx,1),pos(suby,2),'csape'); %closed spline
%             posSpline = fliplr(posSpline);
%             % Can 'ESC' to delete ROI, Right-click to Add Waypoint
%             waypoints = zeros(90,1); %initialize index locations of waypoints
%             waypoints(1:15:end) = 1; %set 8 evenly-spaced waypoints
%             h = drawfreehand('Position', posSpline,'FaceAlpha',0.15,'LineWidth',1, ...
%                 'Multiclick',true,'Waypoints',logical(waypoints)); %create freehand ROI
%             customWait(h); %see custom function below in 'helper functions'
%         end
% 
%         hfhs = findobj(gca, 'Type', 'images.roi.Freehand');
%         editedMask = false(size(image));
%         for ind = 1:numel(hfhs)
%             editedMask = editedMask | hfhs(ind).createMask(); %accumulate mask from each ROI
%             boundaryLocation = round(hfhs(ind).Position); %include ROI boundary
%             bInds = sub2ind(size(image), boundaryLocation(:,2), boundaryLocation(:,1));
%             editedMask(bInds) = true;
%         end
        editedMask = BW;
        
        %% Separate Ascending and Descending Aorta
        [r,c] = find(editedMask);
        [~,idx_min] = min(r);
        aorta1 = bwselect(editedMask,c(idx_min),r(idx_min),8);

        if numROIs==2
            [~,idx_max] = max(r);
            aorta2 = bwselect(editedMask,c(idx_max),r(idx_max),8);
        end

        %% Show and Save Images
%         label = labeloverlay(image,editedMask);
%         figure; imshow(label,[]);

        if numROIs==2
            AscAo_mask(:,:,i) = aorta1;
            DescAo_mask(:,:,i) = aorta2;
        else
            AbdAo_mask(:,:,i) = aorta1; 
        end
    end
end 

save('bssfp.mat','bssfp');
save('bssfp_info.mat','bssfp_info');
save('bssfp_times.mat','bssfp_times');
if numROIs==2
    save('AscAo_mask_manual.mat','AscAo_mask')    
    save('DescAo_mask_manual.mat','DescAo_mask')
else    
    save('AbdAo_mask.mat','AbdAo_mask')
end

close all
cd(home_dir)

disp("bssfp masks saved.");

%% Calculating Velocity and Flow
pc_folder = uigetdir(home_dir,'Select AAo or AbdAo 2DPC Scan');
cd(pc_folder)

dirFiles2 = dir('*.dcm');
if isempty(dirFiles2) %check if dicoms are zipped
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
    dirFiles2 = dir('*.dcm');
end 
pc_info = dicominfo(dirFiles2(1).name);
%pixelArea = pc_info.PixelSpacing(1).*pc_info.PixelSpacing(2);

% Uncrop (zoom out to original FOV) and, if 2 ROIs, combine masks
if numROIs==2
    uncropped = AscAo_mask + DescAo_mask;
else
    uncropped = AbdAo_mask;
end
uncropped = padarray(uncropped, [crop(1)-1 crop(3)-1], 'pre');
uncropped = padarray(uncropped, [bssfpSize-crop(2) bssfpSize-crop(4)], 'post');
% masks = imbinarize(uncropped);
masks = uncropped;

pcFrames = round(length(dirFiles2)/2);
pcmr_info = dicominfo(dirFiles2(1).name);
pcSize = double(pcmr_info.Height); %assumes matrix size is equal in x and y
rrInterval = pcmr_info.NominalInterval; %avg RR int. (ms) 
tempRes = rrInterval/pcFrames;

pcmr = zeros(pcSize,pcSize,length(dirFiles2));
for f=1:length(dirFiles2)
    pcmr(:,:,f) = dicomread(dirFiles2(f).name);
    temp = dicominfo(dirFiles2(f).name);
    pcmr_times(f) = temp.TriggerTime;
end
pc = pcmr(:,:,1:pcFrames);
mag = pcmr(:,:,(pcFrames+1):end);
pcmr_times = pcmr_times(1:pcFrames);
% disp(pcmr_times)
% disp(bssfp_times)
save('pc.mat','pc');
save('mag.mat','mag');
save('pcmr_times.mat','pcmr_times');

% Make spatial resolutions equivalent
desiredSize = max(bssfpSize,pcSize);
desiredFrames = max(bssfpFrames,pcFrames);
bssfp = imresize3(bssfp,[desiredSize desiredSize desiredFrames]);
pc  = imresize3(pc,[desiredSize desiredSize desiredFrames]);
mag = imresize3(mag,[desiredSize desiredSize desiredFrames]);
masks = imresize3(masks,[desiredSize desiredSize desiredFrames]);
pcmr_times = round(interp1((1:pcFrames),pcmr_times,linspace(1,pcFrames,desiredFrames),'linear'));
bssfp_times = interp1((1:bssfpFrames),bssfp_times,linspace(1,bssfpFrames,desiredFrames),'linear');

save('interp_masks_manual.mat','masks');
disp("Calculating flow data.")

% Separate ascending and descending aorta and calculate area and flow
for i=1:desiredFrames
    if numROIs==2
        [r2,c2] = find(masks(:,:,i));
        [~,idx_min2] = min(r2);
        [~,idx_max2] = max(r2);
        
        % Calculate Ascending Aorta flow
        Asc = bwselect(masks(:,:,i),c2(idx_min2),r2(idx_min2),8);
        area = sum(Asc(:)).*pixelArea; %ROI area (cm^2)
        vTemp = pc(:,:,i); %single frame velocities
        AscAo_ROI = double(vTemp(Asc)); %indexed velocities within mask
        AscAo_meanV = mean(AscAo_ROI)*0.1; %mean velocity in ROI in frame i (cm/s)
        AscAo_flow(i) = -area.*AscAo_meanV; %flow in frame i (cm^3/s = mL/s)
        AscAo_area(i) = area; %cm^2
        
        % Calculate Descending Aorta flow
        Desc = bwselect(masks(:,:,i),c2(idx_max2),r2(idx_max2),8);
        area = sum(Desc(:)).*pixelArea; %ROI area (cm^2)
        Desc_ROI = double(vTemp(Desc)); %indexed velocities within mask
        Desc_meanV = mean(Desc_ROI)*0.1; %mean velocity in frame i (cm/s)
        DescAo_flow(i) = area.*Desc_meanV; %flow in frame i (cm^3/s = mL/s)
        DescAo_area(i) = area; %cm^2
    else
        % Calculate Abdominal Aorta flow
        [r2,c2] = find(masks(:,:,i));
        [~,idx_min2] = min(r2);
        [~,idx_max2] = max(r2);
        
        Abd = bwselect(masks(:,:,i),c2(idx_min2),r2(idx_min2),8);
        area = sum(Abd(:)).*pixelArea; %ROI area (cm^2)
        vTemp = pc(:,:,i); %single frame velocities
        AbdAo_ROI = double(vTemp(Abd)); %indexed velocities within mask
        AbdAo_meanV = mean(AbdAo_ROI)*0.1; %mean velocity in frame i (cm/s)
        AbdAo_flow(i) = area.*AbdAo_meanV; %flow in frame i (cm^3/s = mL/s)
        AbdAo_area(i) = area; %cm^2
    end
end


% save_folder = uigetdir(home_dir,'Select folder to save PWV data');
save_folder = '..';
cd(save_folder)
if ~exist('PWV_QA_Analysis','dir')
    mkdir('PWV_QA_Analysis');
end 
cd('PWV_QA_Analysis');
AscAo_flow = circshift(AscAo_flow,5);
AscAo_area = circshift(AscAo_area,5);
if numROIs==2
    save('AscAo_flow.mat','AscAo_flow');
    save('DescAo_flow.mat','DescAo_flow');
    save('AscAo_area.mat','AscAo_area');
    save('DescAo_area.mat','DescAo_area');
else    
    save('AbdAo_flow.mat','AbdAo_flow');
    save('AbdAo_area.mat','AbdAo_area');
end

disp("Flow and area data saved.")

%% Compute PWV

%z = 5/2;
maxTime = max(max(pcmr_times),max(bssfp_times));
if numROIs==2
    % Ascending Aorta
    %interpolate areas to 1 ms resolution 
    temp = interp1(bssfp_times,AscAo_area,linspace(1,maxTime,maxTime),'linear','extrap'); 
    AscAo_area_int = temp(pcmr_times); %get areas only at flow trigger times
    figure;
%     plot(circshift(AscAo_flow,round(length(AscAo_flow)/z)))
    plot(AscAo_flow)
    title('Ascending Aorta - Flow'); 
    xlabel('Frame'); 
    ylabel('Flow (cm^3/s)');
    flow_plot = gcf;
    free1 = drawfreehand;
    asystole = inpolygon(linspace(1,length(AscAo_flow),length(AscAo_flow)), AscAo_flow,free1.Position(:,1),free1.Position(:,2));
%     asystole = inpolygon(circshift(linspace(1,length(AscAo_flow),length(AscAo_flow)),round(length(AscAo_flow)/z)), circshift(AscAo_flow,round(length(AscAo_flow)/z)),free1.Position(:,1),free1.Position(:,2));
%     free2 = drawfreehand;
%     early_sys = inpolygon(AscAo_area(systole),AscAo_flow(systole),free2.Position(:,1),free2.Position(:,2));
    figure; 
%     scatter([AscAo_area_int(3*end/4+1:end) AscAo_area_int(1:end/4)], [AscAo_flow(3*end/4+1:end) AscAo_flow(1:end/4)],[], linspace(1,length(AscAo_area)/2,length(AscAo_area)/2))
    scatter(AscAo_area_int(1:end/2), AscAo_flow(1:end/2),[], linspace(1,length(AscAo_area)/2,length(AscAo_area)/2))
    colormap jet
    % colorbar
    title('Ascending Aorta - QA'); 
    xlabel('Area (cm^2)'); 
    ylabel('Flow (cm^3/s)');
    systolePts = [AscAo_area_int(asystole); AscAo_flow(asystole)]';
    [coef,stats] = polyfit(systolePts(:,1),systolePts(:,2),1);
    minArea = min(systolePts(:,1));
    maxArea = max(systolePts(:,1));
    xq = linspace(minArea-0.1,maxArea+0.1,100);
    yq = coef(1)*xq + coef(2); %y = mx+b
    delete(free1)
    hold on; 
    scatter(systolePts(:,1),systolePts(:,2),72,'k','x');
    plot(xq,yq); 
    str = ['Slope = ' num2str(coef(1)) '\newlineIntcpt = ' num2str(coef(2))];
    text(min(xq),max(yq)-0.1*max(yq),str);
    hold off;
    AscAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    mkdir('AscAo');
    cd('AscAo');
    saveas(flow_plot, 'AscAo_flowplot');
    saveas(gcf,'AscAo_QAplot');
    save('systole_mask.mat', 'asystole');
    save('AscAo_PWV.mat','AscAo_PWV');
    save('coef.mat','coef');
    save('stats.mat','stats');
    cd ..
    movefile('AscAo_area.mat','AscAo')
    movefile('AscAo_flow.mat','AscAo')
    
    %Descending Aorta
    temp = interp1(bssfp_times,DescAo_area,linspace(1,maxTime,maxTime),'linear','extrap'); 
    DescAo_area_int = temp(pcmr_times); %get areas only at flow trigger times
    figure;
    plot(DescAo_flow)
%     plot(circshift(DescAo_flow,round(length(DescAo_flow)/z)))
    title('Descending Aorta - Flow'); 
    xlabel('Frame'); 
    ylabel('Flow (cm^3/s)');
    flow_plot = gcf;
    free2 = drawfreehand;
%     dsystole = inpolygon(circshift(linspace(1,length(DescAo_flow),length(DescAo_flow)),round(length(DescAo_flow)/z)),circshift(DescAo_flow,round(length(DescAo_flow)/z)),free2.Position(:,1),free2.Position(:,2));
    dsystole = inpolygon(linspace(1,length(DescAo_flow),length(DescAo_flow)),DescAo_flow,free2.Position(:,1),free2.Position(:,2));
%     free2 = drawfreehand;
%     early_sys = inpolygon(AscAo_area(systole),AscAo_flow(systole),free2.Position(:,1),free2.Position(:,2));
    figure; 
%     scatter([DescAo_area_int(3*end/4+1:end) DescAo_area_int(1:end/4)], [DescAo_flow(3*end/4+1:end) DescAo_flow(1:end/4)],[], linspace(1,length(DescAo_area)/2,length(DescAo_area)/2))
    scatter(DescAo_area_int(1:end/2), DescAo_flow(1:end/2),[], linspace(1,length(DescAo_area)/2,length(DescAo_area)/2))
    colormap jet
    % colorbar
    title('Descending Aorta - QA'); 
    xlabel('Delta Area (cm^2)'); 
    ylabel('Delta Flow (cm^3/s)');
%     free2 = drawfreehand;
%     in2 = inpolygon(DescAo_area(in),DescAo_flow(in),free2.Position(:,1),free2.Position(:,2));
    systolePts = [DescAo_area_int(dsystole); DescAo_flow(dsystole)]';
    [coef,stats] = polyfit(systolePts(:,1),systolePts(:,2),1);
    minArea = min(systolePts(:,1));
    maxArea = max(systolePts(:,1));
    xq = linspace(minArea-0.1,maxArea+0.1,100);
    yq = coef(1)*xq + coef(2); %y = mx+b
    delete(free2)
    hold on; 
    scatter(systolePts(:,1),systolePts(:,2),72,'k','x');
    plot(xq,yq); 
    str = ['Slope = ' num2str(coef(1)) '\newlineIntcpt = ' num2str(coef(2))];
    text(max(DescAo_area_int)+0.1,max(DescAo_flow)+0.1,str);
    hold off;
    DescAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    mkdir('DescAo');
    cd('DescAo');
    saveas(flow_plot, 'DescAo_flowplot');
    saveas(gcf,'DescAo_QAplot');
    save('systole_mask.mat', 'dsystole');
    save('DescAo_PWV.mat','DescAo_PWV');
    save('coef.mat','coef');
    save('stats.mat','stats');
    cd ..
    movefile('DescAo_area.mat','DescAo')
    movefile('DescAo_flow.mat','DescAo')
else
    %Abdominal Aorta
    temp = interp1(bssfp_times,AbdAo_area,linspace(1,maxTime,maxTime),'linear','extrap'); 
    AbdAo_area_int = temp(pcmr_times); %get areas only at flow trigger times
    figure;
    plot(AbdAo_flow)
%     scatter(linspace(1,length(AbdAo_flow),length(AbdAo_flow)), circshift(AbdAo_flow,round(length(AbdAo_flow)/z)))
    title('Abdominal Aorta - Flow'); 
    xlabel('Frame'); 
    ylabel('Flow (cm^3/s)');
    flow_plot = gcf;
    free = drawfreehand;
    systole = inpolygon(linspace(1,length(AbdAo_flow),length(AbdAo_flow)),AbdAo_flow,free.Position(:,1),free.Position(:,2));
%     systole = inpolygon(circshift(linspace(1,length(AbdAo_flow),length(AbdAo_flow)),round(length(AbdAo_flow)/z)),circshift(AbdAo_flow,round(length(AbdAo_flow)/z)),free.Position(:,1),free.Position(:,2));
%     free2 = drawfreehand;
%     early_sys = inpolygon(AscAo_area(systole),AscAo_flow(systole),free2.Position(:,1),free2.Position(:,2));
    figure; 
    scatter([AbdAo_area_int(3*end/4+1:end) AbdAo_area_int(1:end/4)], [AbdAo_flow(3*end/4+1:end) AbdAo_flow(1:end/4)],[], linspace(1,length(AbdAo_area)/2,length(AbdAo_area)/2))
    colormap jet
    % colorbar
    title('Abdominal Aorta - QA'); 
    xlabel('Delta Area (cm^2)'); 
    ylabel('Delta Flow (cm^3/s)');
%     free2 = drawfreehand;
%     in2 = inpolygon(AbdAo_area(in),AbdAo_flow(in),free2.Position(:,1),free2.Position(:,2));
    systolePts = [AbdAo_area_int(systole); AbdAo_flow(systole)]';
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
    text(max(AbdAo_area_int)+0.1,max(AbdAo_flow)+0.1,str);
    hold off;
    AbdAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    mkdir('AbdAo');
    cd('AbdAo');
    saveas(flow_plot, 'AbdAo_flowplot');
    saveas(gcf,'AbdAo_QAplot');
    save('systole_mask.mat', 'systole');
    save('AbdAo_PWV.mat','AbdAo_PWV');
    save('coef.mat','coef');
    save('stats.mat','stats');
    cd ..
    movefile('AbdAo_area.mat','AbdAo')
    movefile('AbdAo_flow.mat','AbdAo')
end 
close all
cd(home_dir)

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