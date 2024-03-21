pc_folder = 'C:\Users\naren\Desktop\Aorta_Segmentation\vol_Grant\high_res\009_PWV_CartBH_AAo_NOPROSP_1VPS_80FRAMES';
dicoms = read_dicoms(pc_folder);
save_folder = 'C:\Users\naren\Desktop\Aorta_Segmentation\vol_Grant\high_res';
% load('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00016_11110_2020-03-09-h08\00007.PWV_CartBH_AAo\AscAo_mask2.mat')
% load('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00016_11110_2020-03-09-h08\00007.PWV_CartBH_AAo\DescAo_mask2.mat')
% for i=41:size(dicoms, 3)
% 
%         %% Find circles (Hough transform)
%         image = dicoms(:,:,i);
%         image = rescale(image); %normalizes image
% %         [centers,radii,~] = imfindcircles(image,radrange,'ObjectPolarity','bright','Sensitivity',sens);
%         BW = zeros(256,256);
% %         [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
% %         for n = 1:numROIs
% %             try
% %                 BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
% %             catch error
% %                 if numROIs==2
% %                     save('AscAo_mask.mat', 'AscAo_mask')    
% %                     save('DescAo_mask.mat', 'DescAo_mask')
% %                 else    
% %                     save('AbdAo_mask.mat', 'AbdAo_mask')
% %                 end
% %                 close all
% %                 cd(home_dir)
% %                 disp(['Frame: ' num2str(i)])
% %                 rethrow(error)
% %             end
% %         end
% % 
% %         %% Morphological Dilation 
% %         BW = imdilate(BW,strel('disk',3)); %dilate so active contours pulls in segmentation
% % 
% %         %% Active Contours + Dilation
% %         BW = activecontour(image,BW,300,'Chan-Vese','SmoothFactor',5,'ContractionBias',0.2);
% %         BW = imdilate(BW,strel('disk',3));
%         
%         x1=107;
%         y1=122;
%         radius1 = 12;
%         [xx1,yy1] = ndgrid((1:256)-y1,(1:256)-x1);
%         circ1 = (xx1.^2 + yy1.^2)<radius1^2;
%         x2=135;
%         y2=162;
%         radius2 = 9;
%         [xx2,yy2] = ndgrid((1:256)-y2,(1:256)-x2);
%         circ2 = (xx2.^2 + yy2.^2)<radius2^2;
%         BW = BW + circ1 + circ2;
% 
%         
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
% 
%         %% Separate Ascending and Descending Aorta
%         [r,c] = find(editedMask);
%         [~,idx_min] = min(r);
%         aorta1 = bwselect(editedMask,c(idx_min),r(idx_min),8);
% 
%         if numROIs==2
%             [~,idx_max] = max(r);
%             aorta2 = bwselect(editedMask,c(idx_max),r(idx_max),8);
%         end
% 
%         %% Show and Save Images
%         label = labeloverlay(image,editedMask);
%         figure; imshow(label,[]);
% 
%         if numROIs==2
%             AscAo_mask(:,:,i) = aorta1;
%             DescAo_mask(:,:,i) = aorta2;
%         else
%             AbdAo_mask(:,:,i) = aorta1; 
%         end
% end
%
% save('AscAo_mask2.mat','AscAo_mask')    
% save('DescAo_mask2.mat','DescAo_mask')
% numROIs=2;
% load('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00013_11099_2020-03-05-h09\00007.PWV_CartBH_AAo\manual_Asc.mat')
% load('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00013_11099_2020-03-05-h09\00007.PWV_CartBH_AAo\manual_Desc.mat')


%% Calculating Velocity and Flow

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

% Uncrop (zoom out to original FOV) and, if 2 ROIs, combine masks
if numROIs==2
    uncropped = AscAo_mask(:,:,41:80) + DescAo_mask(:,:,41:80);
else
    uncropped = AbdAo_mask;
end
% uncropped = padarray(uncropped, [crop(1)-1 crop(3)-1], 'pre');
% uncropped = padarray(uncropped, [bssfpSize-crop(2) bssfpSize-crop(4)], 'post');
% masks = imbinarize(uncropped);
pcmasks = uncropped;

pcFrames = round(length(dirFiles2)/2);
pcmr_info = dicominfo(dirFiles2(1).name);
pcSize = double(pcmr_info.Height); %assumes matrix size is equal in x and y
rrInterval = pcmr_info.NominalInterval; %avg RR int. (ms) 
tempRes = rrInterval/pcFrames;
pixelArea2 = pcmr_info.PixelSpacing(1).*pcmr_info.PixelSpacing(2);

pcmr = zeros(pcSize,pcSize,length(dirFiles2));
for i=1:length(dirFiles2)
    pcmr(:,:,i) = dicomread(dirFiles2(i).name);
end
pc = pcmr(:,:,1:pcFrames);
mag = pcmr(:,:,(pcFrames+1):end);
save('pc2.mat','pc');
save('mag2.mat','mag');

% % Make spatial resolutions equivalent
% desiredSize = max(bssfpSize,pcSize);
% if pcSize<bssfpSize
%     pc =  imresize(pc,[desiredSize desiredSize]);
%     mag = imresize(pc,[desiredSize desiredSize]);
% elseif bssfpSize<pcSize
%     masks = imresize(masks,[desiredSize desiredSize]);
% end
% 
% % Make temporal resolutions equivalent
% desiredFrames = max(bssfpFrames,pcFrames);
% if bssfpFrames<pcFrames
%     masks = imresize3(masks,[desiredSize desiredSize desiredFrames]);
% elseif pcFrames<bssfpFrames
%     pc  = imresize3(pc,[desiredSize desiredSize desiredFrames]);
%     mag = imresize3(mag,[desiredSize desiredSize desiredFrames]);
% end 
% save('interp_masks2.mat','masks');
% 
% disp("Calculating flow data.")

% Separate ascending and descending aorta and calculate area and flow
for i=1:40
    if numROIs==2
        [r2,c2] = find(pcmasks(:,:,i));
        [~,idx_min2] = min(r2);
        [~,idx_max2] = max(r2);
        
        % Calculate Ascending Aorta flow
        Asc = bwselect(pcmasks(:,:,i),c2(idx_min2),r2(idx_min2),8);
        area = sum(Asc(:)).*pixelArea2; %ROI area (mm^2)
        vTemp = pc(:,:,i); %single frame velocities
        AscAo_ROI = double(vTemp(Asc)); %indexed velocities within mask
        AscAo_meanV = mean(AscAo_ROI); %mean velocity in ROI in frame i (mm/s)
        AscAo_flow(i) = area.*AscAo_meanV.*0.001; %flow in frame i (mm^3/s*0.001 = mL/s)
        AscAo_area(i) = area*0.01; %cm^2
        
        % Calculate Descending Aorta flow
        Desc = bwselect(pcmasks(:,:,i),c2(idx_max2),r2(idx_max2),8);
        area = sum(Desc(:)).*pixelArea2; %ROI area (mm^2)
        Desc_ROI = double(vTemp(Desc)); %indexed velocities within mask
        Desc_meanV = mean(Desc_ROI); %mean velocity in frame i (mm/s)
        DescAo_flow(i) = area.*Desc_meanV.*0.001; %flow in frame i (mm^3/s*0.001 = mL/s)
        DescAo_area(i) = area*0.01; %cm^2
    else
        % Calculate Abdominal Aorta flow
        [r2,c2] = find(pcmasks(:,:,i));
        [~,idx_min2] = min(r2);
        [~,idx_max2] = max(r2);
        
        Abd = bwselect(pcmasks(:,:,i),c2(idx_min2),r2(idx_min2),8);
        area = sum(Abd(:)).*pixelArea2;%ROI area (mm^2)
        vTemp = pc(:,:,i); %single frame velocities
        AbdAo_ROI = double(vTemp(Abd)); %indexed velocities within mask
        AbdAo_meanV = mean(AbdAo_ROI); %mean velocity in frame i (mm/s)
        AbdAo_flow(i) = area.*AbdAo_meanV.*0.001; %flow in frame i (mm^3/s*0.001 = mL/s)
        AbdAo_area(i) = area*0.01; %cm^2
    end
end



cd(save_folder)
if ~exist('PWV_QA_Analysis','dir')
    mkdir('PWV_QA_Analysis');
end 
cd('PWV_QA_Analysis');

if numROIs==2
    save('AscAo_flow2.mat','AscAo_flow');
    save('DescAo_flow2.mat','DescAo_flow');
    save('AscAo_area2.mat','AscAo_area');
    save('DescAo_area2.mat','DescAo_area');
else    
    save('AbdAo_flow.mat','AbdAo_flow');
    save('AbdAo_area.mat','AbdAo_area');
end

disp("Flow and area data saved.")

AscAo_area = [AscAo_area(1:15) AscAo_area(36:end)];
AscAo_flow = [-AscAo_flow(1:15) -AscAo_flow(36:end)];
DescAo_area = [DescAo_area(1:15) DescAo_area(36:end)];
DescAo_flow = [DescAo_flow(1:15) DescAo_flow(36:end)];
%% Compute PWV
if numROIs==2
    %Ascending Aorta
%     AscAo_area = interp1((1:length(AscAo_area)),(AscAo_area),(1:0.5:length(AscAo_area)),'pchip');
%     AscAo_flow = interp1((1:length(AscAo_flow)),(AscAo_flow),(1:0.5:length(AscAo_flow)),'pchip');
    figure;
%     plot(circshift(AscAo_flow,round(length(AscAo_flow)/z)))
    plot(AscAo_flow)
    title('Ascending Aorta - Flow2'); 
    xlabel('Frame'); 
    ylabel('Flow (cm^3/s)');
    flow_plot = gcf;
    free1 = drawfreehand;
    asystole = inpolygon(linspace(1,length(AscAo_flow),length(AscAo_flow)), AscAo_flow,free1.Position(:,1),free1.Position(:,2));
%     asystole = inpolygon(circshift(linspace(1,length(AscAo_flow),length(AscAo_flow)),round(length(AscAo_flow)/z)), circshift(AscAo_flow,round(length(AscAo_flow)/z)),free1.Position(:,1),free1.Position(:,2));
%     free2 = drawfreehand;
%     early_sys = inpolygon(AscAo_area(systole),AscAo_flow(systole),free2.Position(:,1),free2.Position(:,2));
    figure; 
    scatter(AscAo_area, AscAo_flow,[], linspace(1,length(AscAo_area),length(AscAo_area)))
    colormap jet
    % colorbar
    title('Ascending Aorta - QA2'); 
    xlabel('Delta Area (cm^2)'); 
    ylabel('Delta Flow (cm^3/s)');
    systolePts = [AscAo_area(asystole); AscAo_flow(asystole)]';
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
    text(max(AscAo_area)+0.1,max(AscAo_flow)+0.1,str);
    hold off;
    AscAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    mkdir('AscAo');
    cd('AscAo');
    saveas(gcf,'AscAo_QAplot2');
    saveas(flow_plot, 'AscAo_flowplot2');
    save('systole_mask2.mat', 'asystole');
    save('AscAo_PWV2.mat','AscAo_PWV');
    save('coef2.mat','coef');
    save('stats2.mat','stats');
    cd ..
    movefile('AscAo_area2.mat','AscAo')
    movefile('AscAo_flow2.mat','AscAo')
    
    %Descending Aorta
%     DescAo_area = interp1((1:length(DescAo_area)),(DescAo_area),(1:0.5:length(DescAo_area)),'pchip');
%     DescAo_flow = interp1((1:length(DescAo_flow)),(DescAo_flow),(1:0.5:length(DescAo_flow)),'pchip');
    figure;
    plot(DescAo_flow)
%     plot(circshift(DescAo_flow,round(length(DescAo_flow)/z)))
    title('Descending Aorta - Flow2'); 
    xlabel('Frame'); 
    ylabel('Flow (cm^3/s)');
    flow_plot = gcf;
    free2 = drawfreehand;
%     dsystole = inpolygon(circshift(linspace(1,length(DescAo_flow),length(DescAo_flow)),round(length(DescAo_flow)/z)),circshift(DescAo_flow,round(length(DescAo_flow)/z)),free2.Position(:,1),free2.Position(:,2));
    dsystole = inpolygon(linspace(1,length(DescAo_flow),length(DescAo_flow)),DescAo_flow,free2.Position(:,1),free2.Position(:,2));
%     free2 = drawfreehand;
%     early_sys = inpolygon(AscAo_area(systole),AscAo_flow(systole),free2.Position(:,1),free2.Position(:,2));
    figure; 
    scatter(DescAo_area, DescAo_flow, [], linspace(1,length(DescAo_area),length(DescAo_area)))
    colormap jet
    % colorbar
    title('Descending Aorta - QA2'); 
    xlabel('Delta Area (cm^2)'); 
    ylabel('Delta Flow (cm^3/s)');
%     free2 = drawfreehand;
%     in2 = inpolygon(DescAo_area(in),DescAo_flow(in),free2.Position(:,1),free2.Position(:,2));
    systolePts = [DescAo_area(dsystole); DescAo_flow(dsystole)]';
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
    text(max(DescAo_area)+0.1,max(DescAo_flow)+0.1,str);
    hold off;
    DescAo_PWV = coef(1);
    disp(['PWV_QA = ' num2str(coef(1)*0.01) ' m/s']);
    mkdir('DescAo');
    cd('DescAo');
    saveas(gcf,'DescAo_QAplot2');
    saveas(flow_plot, 'DescAo_flowplot2');
    save('systole_mask2.mat', 'dsystole');
    save('DescAo_PWV2.mat','DescAo_PWV');
    save('coef2.mat','coef');
    save('stats2.mat','stats');
    cd ..
    movefile('DescAo_area2.mat','DescAo')
    movefile('DescAo_flow2.mat','DescAo')
else
    %Abdominal Aorta
    figure;
    plot(AbdAo_flow)
%     scatter(linspace(1,length(AbdAo_flow),length(AbdAo_flow)), circshift(AbdAo_flow,round(length(AbdAo_flow)/z)))
    title('Abdominal Aorta - Flow2'); 
    xlabel('Frame'); 
    ylabel('Flow (cm^3/s)');
    flow_plot = gcf;
    free = drawfreehand;
    systole = inpolygon(linspace(1,length(AbdAo_flow),length(AbdAo_flow)),AbdAo_flow,free.Position(:,1),free.Position(:,2));
%     systole = inpolygon(circshift(linspace(1,length(AbdAo_flow),length(AbdAo_flow)),round(length(AbdAo_flow)/z)),circshift(AbdAo_flow,round(length(AbdAo_flow)/z)),free.Position(:,1),free.Position(:,2));
%     free2 = drawfreehand;
%     early_sys = inpolygon(AscAo_area(systole),AscAo_flow(systole),free2.Position(:,1),free2.Position(:,2));
    figure; 
    scatter(AbdAo_area, AbdAo_flow,[], linspace(1,length(AbdAo_area),length(AbdAo_area)))
    colormap jet
    % colorbar
    title('Abdominal Aorta - QA2'); 
    xlabel('Delta Area (cm^2)'); 
    ylabel('Delta Flow (cm^3/s)');
%     free2 = drawfreehand;
%     in2 = inpolygon(AbdAo_area(in),AbdAo_flow(in),free2.Position(:,1),free2.Position(:,2));
    systolePts = [AbdAo_area(systole); AbdAo_flow(systole)]';
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
    saveas(gcf,'AbdAo_QAplot2');
    saveas(flow_plot, 'AbdAo_flowplot2');
    save('AbdAo_PWV2.mat','AbdAo_PWV');
    save('systole_mask2.mat', 'systole');
    save('coef2.mat','coef');
    save('stats2.mat','stats');
    cd ..
    movefile('AbdAo_area2.mat','AbdAo')
    movefile('AbdAo_flow2.mat','AbdAo')
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