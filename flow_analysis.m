clear all; close all;

disp("Choose where to save analysis data")
save_loc = uigetdir(pwd, "Choose where to save analysis data");

plane = inputdlg("Type in the analysis region name");
plane = plane{1};

num_roi = inputdlg("Type in the number of ROIs in plane");
num_roi = str2num(num_roi{1});

list_options = {'kmeans','otsu','hough','hough+contours','edge+hough+contours', 'circle+contours', 'manual'};
[choice, tf] = listdlg('ListString', list_options, 'InitialValue', 4, 'PromptString', "Choose a segmentation method", 'SelectionMode', 'single');
segType = list_options{choice};

data_options = {'dicom', 'dat', 'hdf5'};
[data_choice, data_tf] = listdlg('ListString', data_options, 'InitialValue', 1, 'PromptString', "Choose a input data type", 'SelectionMode', 'single');
dataType = data_options{data_choice};

% plane = 'asc';
% segType = 'hough+contours';
% num_roi = 2;

date = string(datetime('now', 'Format', 'yyyy-MM-dd-HHmm'));

[~, save_name] = fileparts(save_loc);
if strcmp(save_name, '2D_Flow_Analysis')
    result_dir = save_loc;
else
    result_dir = fullfile(save_loc, '2D_Flow_Analysis');
    if ~exist(result_dir,'dir')
        mkdir(result_dir);
    end
end

% save_loc = 'D:\PWV\vol_Tarun';
plane_dir = fullfile(result_dir, plane);
modal_dir = fullfile(plane_dir, dataType);
analysis_dir = fullfile(modal_dir, strcat(segType, '_', date));
if ~exist(analysis_dir,'dir')
    mkdir(analysis_dir)
end
disp("Created analysis folder")
disp("Loading data...")


data_dir = fullfile(modal_dir, "data");
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

switch dataType
    case 'dicom'
        disp("Select DICOM scan folder")
        dicom_folder = uigetdir(pwd, 'Select 2DPC DICOM folder');
        dirFiles = dir(fullfile(dicom_folder, '*.dcm'));
        dicomInfo = dicominfo(fullfile(dirFiles(1).folder, dirFiles(1).name)); %grab dicom info from 1st frame
        height = double(dicomInfo.Height);
        width = double(dicomInfo.Width);
        frames = dicomInfo.CardiacNumberOfImages;
        pixelArea = dicomInfo.PixelSpacing(1).*dicomInfo.PixelSpacing(2); % mm^2
        images = zeros(height,width,frames*2);
        for f=1:length(dirFiles)
            images(:,:,f) = dicomread(fullfile(dirFiles(f).folder, dirFiles(f).name));
            temp = dicominfo(fullfile(dirFiles(f).folder, dirFiles(f).name));
            pcmrTimes(f) = temp.TriggerTime;
        end
        pcmrTimes = pcmrTimes(1:frames);
        timeRes = frames/max(pcmrTimes);
        pc = images(:,:,1:frames);
        mag = images(:,:,frames+1:end);

    case 'dat'
        disp("Select folder with dat files")
        binary_data = loadPCVIPR();
        header = binary_data.header;
        width = header.matrixx;
        height = header.matrixy;
        frames = header.frames;
        timeRes = header.timeres;
        pcmrTimes = double(timeRes.*(0:(frames-1)));
        venc = header.VENC;
        pixelArea = header.matrixx/header.fovx;
        mag = flip(squeeze(binary_data.mag));
        pc = flip(squeeze(binary_data.vel3));
        
    case 'hdf5'
        %% Load H5
        disp("Select HDF5 file")
        [h5_file, h5_dir] = uigetfile('*.h5', pwd, 'Select .h5 file');
        mag = flip(squeeze(h5read(fullfile(h5_dir, h5_file), "/MAG")));
        pc = flip(squeeze(h5read(fullfile(h5_dir, h5_file), "/VZ")));
        try
            width = double(h5readatt(fullfile(h5_dir, h5_file), "/HEADER", "matrixx"));
            height = double(h5readatt(fullfile(h5_dir, h5_file), "/HEADER", "matrixy"));
            frames = h5readatt(fullfile(h5_dir, h5_file), "/HEADER", "frames");
            venc = h5readatt(fullfile(h5_dir, h5_file), "/HEADER", "venc");
            fovx = h5readatt(fullfile(h5_dir, h5_file), "/HEADER", "fovx");
            timeRes = h5readatt(fullfile(h5_dir, h5_file), "/HEADER", "time_res");
        catch 
            disp("ERROR: reading h5 header failed, attributes manually set. Please adjust values in code.")
            width = size(mag, 1);
            height = size(mag, 2);
            frames = size(mag, 3);
            venc = 80;
            fovx = 320;
            timeRes = 30;
        end
        pixelArea = fovx/width;
        pcmrTimes = double(timeRes.*(0:(frames-1)));
end
tavg_mag = mean(mag, 3);

%% Initialize Masks
answer = questdlg("Load masks?");
switch answer
    case 'Yes'
        % select mask.mat file
        [mask_file, mask_path] = uigetfile({'*.mat'}, "Select mask file", fullfile(data_dir, 'masks.mat'));
        load(fullfile(mask_path, mask_file));
        load(fullfile(mask_path, 'crop.mat'))
        load(fullfile(mask_path, 'areas.mat'))
        if size(masks_saved,2) < frames
            needsMask = 1;
            currFrame = size(masks_saved,2);
        else
            needsMask = 0;
            currFrame = frames;
        end
        mask = zeros(num_roi, height, width, frames);
        area = zeros(num_roi,frames);
        for i=1:currFrame
            for j=1:num_roi
                mask(j,:,:,i) = masks_saved{j, i};
                area(j, i) = areas_saved{j, i};
            end
        end

    case 'No'
        % Zoom in Manually
        figure; imshow(tavg_mag,[]); title('Select region to zoom in on');
        rect = drawrectangle;
        crop = round([rect.Position(2), rect.Position(2)+rect.Position(4), rect.Position(1), rect.Position(1)+rect.Position(3)]);
        save(fullfile(data_dir, 'crop.mat'), 'crop')
        close all;
        needsMask = 1;
        currFrame = 1;
        masks_saved = cell(num_roi, 1);
        areas_saved = cell(num_roi, 1);
        mask = zeros(num_roi, height, width, frames);
        area = zeros(num_roi, frames);

    case 'Cancel'
        return
end 

disp("Segmenting magnitude images...")
radrange = [10 30]; % radius range for Hough Transform
sens = 0.8; % sensitivity of circle find
%% Segment Each Frame Separately
if needsMask % if we haven't loaded in the mask...
    for i=currFrame:frames
        fprintf("Frame %d\n", i)
        image = mag(crop(1):crop(2),crop(3):crop(4),i);
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
                    if sum(BW) == 0
                        error('');
                    end
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
                    if sum(BW) == 0
                        error('');
                    end
                catch
                    disp('Active contours failed, continuing with circle mask');
                end
            case 'circle+contours'
                f=figure; imshow(image,[]); 
                circle = drawcircle();
                radius = circle.Radius; %get radius of circle
                center = round(circle.Center); %get center coordinates
                centers = center;
                radii = radius;
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
                    if sum(BW) == 0
                        error('Active contours failed, continuing with circle mask');
                    end
                catch
                    disp('Active contours failed, continuing with circle mask');
                end
            case 'manual'
                f=figure; imshow(image,[]); 
                circle = drawcircle();
                radius = circle.Radius; %get radius of circle
                center = round(circle.Center); %get center coordinates
                centers = center;
                radii = radius;
                close(f); delete(f);
                [X,Y] = ndgrid(1:size(image,1),1:size(image,2));
                X = X-center(2); %shift coordinate grid
                Y = Y-center(1);
                BW = sqrt(X.^2 + Y.^2)<=radius; %anything outside radius is ignored
            otherwise
                disp("Invalid segmentation choice!")
                return
        end

        % Freehand ROI conversion
        blocations = bwboundaries(BW,'noholes');
        numBlobs = numel(blocations);
        f = figure;
        imshow(image, []);
        set(f,'WindowStyle','normal')
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
            vessel = padarray(vessel, [height-crop(2)-1 width-crop(4)-1], 'post');
            mask(j,:,:,i) = vessel;
            area(j, i) = sum(vessel(:))*pixelArea; % mm^2
            masks_saved{j, i} = mask(j,:,:,i);
            areas_saved{j, i} = area(j, i);
        end
        close all
        save(fullfile(data_dir, 'masks.mat'), 'masks_saved');
        save(fullfile(data_dir, 'areas.mat'), 'areas_saved');
    end
    % save(fullfile(data_dir, 'masks.mat'), 'mask');
    % save(fullfile(data_dir, 'areas.mat'), 'area');
    disp('Mask and Area Data Saved');
end

%% Save ROI images
e = figure; imshow(mag(:,:,1), []);
frame_mag = getframe(e);
imwrite(frame2im(frame_mag), fullfile(analysis_dir, 'mag.png'));
close(e);
f = figure; imshow(imoverlay(rescale(mag(:,:,1)), squeeze(sum(sum(mask, 1), 4))));
frame_mag_mask = getframe(f);
imwrite(frame2im(frame_mag_mask), fullfile(analysis_dir, 'mag_mask.png'));
close(f);
g = figure; imshow(rescale(pc(:,:,1)), []);
frame_pc = getframe(g);
imwrite(frame2im(frame_pc), fullfile(analysis_dir, 'pc.png'));
close(g);
combined_mask = squeeze(sum(mask, 1));
h = figure; imshow(imoverlay(rescale(pc(:,:,1)), combined_mask(:,:,1)), []);
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
flow = zeros(num_roi, frames);
vel = zeros(num_roi, frames, 3); %mean, min, max
for j=1:num_roi
    for i=1:frames
        ROI = imbinarize(mask(j,:,:,i));
        vTemp = pc(:,:,i); %single frame velocities (mm/s)
        ROIindex = double(vTemp(ROI)); %indexed velocities within mask
        vMean = mean(ROIindex); %mean velocity
        vMin = min(ROIindex); %min velocity
        vMax = max(ROIindex); %max velcoity
        vel(j, i, :) = [vMean vMin vMax];
        flow(j, i) = area(j, i).*vMean.*0.001; %flow in frame i (mL/s)
    end
end
save(fullfile(data_dir, 'flows.mat'), 'flow');
save(fullfile(data_dir, 'vel_stats.mat'), 'vel');
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
mag = circshift(mag,shift,3);
pc = circshift(pc,shift,3);
area_calc = area;
flow_calc = flow;
vel_calc = vel;
for i=1:num_roi
    mask(i, :,:,:) = circshift(mask(i, :,:,:), shift, 4);
    area_calc(i, :) = circshift(area_calc(i, :), shift, 2);
    flow_calc(i, :) = circshift(flow_calc(i, :), shift, 2);
    vel_calc(i, :, :) = circshift(vel_calc(i, :, :), shift, 2);
end


%% Compute PWV
disp("Calculating PWV...")
PWV = zeros(num_roi, 1);

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
    earlySystolePts = find(inpolygon(linspace(1,length(flow_calc(j, :)),length(flow_calc(j, :))),flow_calc(j, :),free.Position(:,1),free.Position(:,2)));
    earlySystoleTimes = pcmrTimes(earlySystolePts);
    earlySystolePts = ismember(earlySystoleTimes, pcmrTimes);
    earlySystoleTimes = earlySystoleTimes(earlySystolePts);
    delete(free);
    saveas(flow_plot, fullfile(analysis_dir, sprintf('flowplot_ROI_%d', j)));
    close(flow_plot); clear flow_plot;

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
    x = area(j, earlySystolePts);
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
    saveas(QAplot,fullfile(analysis_dir, sprintf('QAplot_ROI_%d', j)));
    save(fullfile(analysis_dir, sprintf('coef_ROI_%d.mat', j)), 'coef');
    save(fullfile(analysis_dir, sprintf('stats_ROI_%d.mat', j)), 'stats');
end

% Initialize array for Excel file
results = cell(num_roi+2, 6);
results{1, 1} = 'ROI';
results{1, 2} = 'PWV (m/s)';
results{1, 3} = 'Mean Flow (mL/s)';
results{1, 4} = 'Min Flow (mL/s)';
results{1, 5} = 'Max Flow (mL/s)';
results{1, 6} = 'Peak Velocity (mm/s)';

for j = 1:num_roi
    results{j+1, 1} = j;
    results{j+1, 2} = PWV(j);
    results{j+1, 3} = mean(flow_calc(j,:)); 
    results{j+1, 4} = min(flow_calc(j,:));
    results{j+1, 5} = max(flow_calc(j,:));
    results{j+1, 6} = max(vel_calc(j,:,3));
end

results2 = cell(frames+1,1+2*num_roi);
results2{1, 1} = 'Time (ms)';
for i=0:frames - 1
    results2{i+2, 1} = i*timeRes;
end
for j = 1:num_roi
    results2{1, 2*j} = sprintf('ROI%d Area (mm^2)', j);
    results2{1, 2*j+1} = sprintf('ROI%d Flow (mL/s)', j);
    for i = 1:frames
        results2{i+1, 2*j} = area_calc(j, i);
        results2{i+1, 2*j+1} = flow_calc(j, i);
    end
end

writecell(results, fullfile(analysis_dir, 'results.xlsx'))
writecell(results2, fullfile(analysis_dir, 'results.xlsx'), 'WriteMode', 'append')
disp("( つ ◕_◕ )つ 2D Flow Data Saved! ( つ ◕_◕ )つ")


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