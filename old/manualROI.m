% Code assumes follwoing folder structure:
%          home_dir
%             - life_life[patient num]...
%                - [scan num].[name]_[mode]_BH_[scan type]AAo
%                   - .tgz file or DICOMS
%                - [scan num].[name]_[mode]_BH_[scan type]AbdAo
%                   - .tgz file or DICOMS
% Select patient folder directory

home_dir = 'D:\PWV\PWV_data\';
cd(home_dir)

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
    
    bssfp_info = dicominfo(dirFiles(1).name); %grab dicom info from 1st frame
    bssfpSize = double(bssfp_info.Height); %assumes matrix size is equal in x and y
    bssfpFrames = bssfp_info.CardiacNumberOfImages;
    bssfp = zeros(bssfpSize,bssfpSize,bssfpFrames);
    for f=1:length(dirFiles)
        bssfp(:,:,f) = dicomread(dirFiles(f).name);
    end 
else
    load('bssfp.mat');
    load('bssfp_info.mat');
end 

bssfpFrames = bssfp_info.CardiacNumberOfImages;
bssfpSize = double(bssfp_info.Height); %assumes matrix size is equal in x and y
cropSize = 256; %BASED ON CROP SIZES BELOW (currently 256x256 to better visualize aorta)
pixelArea = bssfp_info.PixelSpacing(1).*bssfp_info.PixelSpacing(2); %mm^2

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
    crop = [155 410 101 356]; %crop dims (1:2 = row width, 3:4 = col width)
    radrange = [13 40]; %radius range for Hough Transform
    sens = 0.90; %sensitivity of circle find
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
        f = figure;
        imshow(image, []);
        f.WindowState = 'maximized';

        for ind = 1:numROIs
            poly = drawpolygon;
            pos = poly.Position;
            posSpline = interparc(90,[pos(:,1); pos(1,1)],[pos(:,2); pos(1,2)],'csape'); %closed spline
            %posSpline = fliplr(posSpline);
            % Can 'ESC' to delete ROI, Right-click to Add Waypoint
            waypoints = zeros(90,1); %initialize index locations of waypoints
            waypoints(1:15:end) = 1; %set 8 evenly-spaced waypoints
            delete('poly');
            h = drawfreehand('Position', posSpline,'FaceAlpha',0.15,'LineWidth',1, ...
                'Multiclick',true,'Waypoints',logical(waypoints)); %create freehand ROI
            customWait(h); %see custom function below in 'helper functions'
            blocations{ind} = pos;
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

        if numROIs==2
            [~,idx_max] = max(r);
            aorta2 = bwselect(editedMask,c(idx_max),r(idx_max),8);
        end

        %% Show and Save Images
        label = labeloverlay(image,editedMask);
        figure; imshow(label,[]);

        if numROIs==2
            AscAo_mask(:,:,i) = aorta1;
            DescAo_mask(:,:,i) = aorta2;
        else
            AbdAo_mask(:,:,i) = aorta1; 
        end
    end
end 

if numROIs==2
    save('bssfp.mat','bssfp');
    save('bssfp_info.mat','bssfp_info');
    save('AscAo_mask.mat','AscAo_mask')    
    save('DescAo_mask.mat','DescAo_mask')
else    
    save('bssfp.mat','bssfp');
    save('bssfp_info.mat','bssfp_info');
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
end 

% Uncrop (zoom out to original FOV) and, if 2 ROIs, combine masks
if numROIs==2
    uncropped = AscAo_mask + DescAo_mask;
else
    uncropped = AbdAo_mask;
end
uncropped = padarray(uncropped, [crop(1) crop(3)], 'pre');
uncropped = padarray(uncropped, [bssfpSize-crop(2)-1 bssfpSize-crop(4)-1], 'post');
masks = imbinarize(uncropped);

pcFrames = round(length(dirFiles2)/2);
pcmr_info = dicominfo(dirFiles2(1).name);
pcSize = double(pcmr_info.Height); %assumes matrix size is equal in x and y
rrInterval = pcmr_info.NominalInterval; %avg RR int. (ms) 
tempRes = rrInterval/pcFrames;

pcmr = zeros(pcSize,pcSize,length(dirFiles2));
for i=1:length(dirFiles2)
    pcmr(:,:,i) = dicomread(dirFiles2(i).name);
end
pc = pcmr(:,:,1:pcFrames);
mag = pcmr(:,:,(pcFrames+1):end);
save('pc.mat','pc');
save('mag.mat','mag');

% Make spatial resolutions equivalent
desiredSize = max(bssfpSize,pcSize);
if pcSize<bssfpSize
    pc =  imresize(pc,[desiredSize desiredSize]);
    mag = imresize(pc,[desiredSize desiredSize]);
elseif bssfpSize<pcSize
    masks = imresize(masks,[desiredSize desiredSize]);
    bssfp = imresize(bssfp,[desiredSize desiredSize]);
end 

% Make temporal resolutions equivalent
desiredFrames = max(bssfpFrames,pcFrames);
if bssfpFrames<pcFrames
    masks = imresize3(masks,[desiredSize desiredSize desiredFrames]);
    bssfp = imresize3(bssfp,[desiredSize desiredSize desiredFrames]);
elseif pcFrames<bssfpFrames
    pc  = imresize3(pc,[desiredSize desiredSize desiredFrames]);
    mag = imresize3(mag,[desiredSize desiredSize desiredFrames]);
end 

% Separate ascending and descending aorta and calculate area and flow
for i=1:desiredFrames
    if numROIs==2
        [r2,c2] = find(masks(:,:,i));
        [~,idx_min2] = min(r2);
        [~,idx_max2] = max(r2);
        
        % Calculate Ascending Aorta flow
        Asc = bwselect(masks(:,:,i),c2(idx_min2),r2(idx_min2),8);
        area = sum(Asc(:)).*pixelArea; %ROI area (mm^2)
        vTemp = pc(:,:,i); %single frame velocities
        AscAo_ROI = double(vTemp(Asc)); %indexed velocities within mask
        AbdAo_meanV = mean(AscAo_ROI); %mean velocity in ROI in frame i (mm/s)
        AscAo_flow(i) = area.*AbdAo_meanV.*0.001; %flow in frame i (mm^3/s*0.001 = mL/s)
        AscAo_area(i) = area*0.01; %cm^2
        
        % Calculate Descending Aorta flow
        Desc = bwselect(masks(:,:,i),c2(idx_max2),r2(idx_max2),8);
        area = sum(Desc(:)).*pixelArea; %ROI area (mm^2)
        Desc_ROI = double(vTemp(Desc)); %indexed velocities within mask
        Desc_meanV = mean(Desc_ROI); %mean velocity in frame i (mm/s)
        DescAo_flow(i) = area.*Desc_meanV.*0.001; %flow in frame i (mm^3/s*0.001 = mL/s)
        DescAo_area(i) = area*0.01; %cm^2
    else
        % Calculate Abdominal Aorta flow
        [r2,c2] = find(masks(:,:,i));
        [~,idx_min2] = min(r2);
        [~,idx_max2] = max(r2);
        
        Abd = bwselect(masks(:,:,i),c2(idx_min2),r2(idx_min2),8);
        area = sum(Abd(:)).*pixelArea;%ROI area (mm^2)
        vTemp = pc(:,:,i); %single frame velocities
        AbdAo_ROI = double(vTemp(Abd)); %indexed velocities within mask
        AbdAo_meanV = mean(AbdAo_ROI); %mean velocity in frame i (mm/s)
        AbdAo_flow(i) = area.*AbdAo_meanV.*0.001; %flow in frame i (mm^3/s*0.001 = mL/s)
        AbdAo_area(i) = area*0.01; %cm^2
    end
end

save_folder = uigetdir(home_dir,'Select folder to save PWV data');
cd(save_folder)
if ~exist('PWV_QA_Analysis','dir')
    mkdir('PWV_QA_Analysis');
end 
cd('PWV_QA_Analysis');

% Save image overlays with ROIs
if numROIs==2
    figure; imshowpair(bssfp(:,:,1),pc(:,:,1)); 
    title('BSSFP-PC Overlay Frame 1');
    hold on; visboundaries(masks(:,:,1),'LineWidth',1); hold off;
    saveas(gcf,'AAo_overlay_frame1');
    figure; imshowpair(bssfp(:,:,desiredFrames/2),pc(:,:,desiredFrames/2)); 
    title(['BSSFP-PC Overlay Frame ' num2str(desiredFrames/2)]);
    hold on; visboundaries(masks(:,:,desiredFrames/2),'LineWidth',1); hold off;
    saveas(gcf,['AAo_overlay_frame'  num2str(desiredFrames/2)]);
    figure; imshowpair(bssfp(:,:,desiredFrames),pc(:,:,desiredFrames)); 
    title(['BSSFP-PC Overlay Frame ' num2str(desiredFrames)]);
    hold on; visboundaries(masks(:,:,desiredFrames),'LineWidth',1); hold off;
    saveas(gcf,['AAo_overlay_frame'  num2str(desiredFrames)]);
else
    figure; imshowpair(bssfp(:,:,1),pc(:,:,1)); 
    title('BSSFP-PC Overlay Frame 1');
    hold on; visboundaries(masks(:,:,1),'LineWidth',1); hold off;
    saveas(gcf,'AbdAo_overlay_frame1');
    figure; imshowpair(bssfp(:,:,desiredFrames/2),pc(:,:,desiredFrames/2)); 
    title(['BSSFP-PC Overlay Frame ' num2str(desiredFrames/2)]);
    hold on; visboundaries(masks(:,:,desiredFrames/2),'LineWidth',1); hold off;
    saveas(gcf,['AbdAo_overlay_frame'  num2str(desiredFrames/2)]);
    figure; imshowpair(bssfp(:,:,desiredFrames),pc(:,:,desiredFrames)); 
    title(['BSSFP-PC Overlay Frame ' num2str(desiredFrames)]);
    hold on; visboundaries(masks(:,:,desiredFrames),'LineWidth',1); hold off;
    saveas(gcf,['AbdAo_overlay_frame'  num2str(desiredFrames)]);
end 

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
if numROIs==2
    %Ascending Aorta
    if mean(AscAo_flow)<0
        AscAo_flow = -1*AscAo_flow;
    end 
    figure; 
    scatter(AscAo_area,AscAo_flow,[],(1:40));
    c = colorbar;
    c.Label.String = 'Cardiac Frame #';
    title('Ascending Aorta - QA'); 
    xlabel('\Delta Area (cm^2)'); 
    ylabel('\Delta Flow (cm^3/s)');
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
    movefile('AscAo_area.mat','AscAo')
    movefile('AscAo_flow.mat','AscAo')
    
    %Descending Aorta
    if mean(DescAo_flow)<0
        DescAo_flow = -1*DescAo_flow;
    end     
    figure; 
    scatter(DescAo_area,DescAo_flow,[],(1:40));
    c = colorbar;
    c.Label.String = 'Cardiac Frame #';
    title('Descending Aorta - QA'); 
    xlabel('\Delta Area (cm^2)'); 
    ylabel('\Delta Flow (cm^3/s)');
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
    movefile('DescAo_area.mat','DescAo')
    movefile('DescAo_flow.mat','DescAo')
else
    %Abdominal Aorta
    if mean(AbdAo_flow)<0
        AbdAo_flow = -1*AbdAo_flow;
    end      
    figure; 
    scatter(AbdAo_area,AbdAo_flow,[],(1:40));
    c = colorbar;
    c.Label.String = 'Cardiac Frame #';
    title('Abdominal Aorta - QA'); 
    xlabel('\Delta Area (cm^2)'); 
    ylabel('\Delta Flow (cm^3/s)');
    free = drawfreehand;
    in = inpolygon(AbdAo_area,AbdAo_flow,free.Position(:,1),free.Position(:,2));
    systolePts = [AbdAo_area(in); AbdAo_flow(in)]';
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
    text(max(AbdAo_area)+0.1,max(AbdAo_flow)+0.1,str);
    hold off;
    AbdAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    mkdir('AbdAo');
    cd('AbdAo');
    saveas(gcf,'AbdAo_QAplot');
    save('AbdAo_PWV.mat','AbdAo_PWV');
    save('coef.mat','coef');
    save('stats.mat','stats');
    cd ..
    movefile('AbdAo_area.mat','AbdAo')
    movefile('AbdAo_flow.mat','AbdAo')
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