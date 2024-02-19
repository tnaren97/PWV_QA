folder = 'life_life00016_11110_2020-03-09-h08';
num = '07';
num2 = '08';
mode = 'PWV';
type = 'CartBH';
anat = 'A';

desiredSize = 512;

scan = sprintf('000%s.%s_%s_%sAo', num, mode, type, anat);
patient = folder(13:14);
id = [patient '-' num2];

if strcmp(anat, 'A')
    AscAof = sprintf('C:\\Users\\naren\\Desktop\\Aorta_Segmentation\\results\\%s_AscAo.mat', id);
    DescAof = sprintf('C:\\Users\\naren\\Desktop\\Aorta_Segmentation\\results\\%s_DescAo.mat', id);
    load(AscAof);
    load(DescAof);
    crop = [155 410 101 356];
    masks = AscAo + DescAo;
else
    AbdAof = sprintf('C:\\Users\\naren\\Desktop\\Aorta_Segmentation\\results\\%s_AbdAo.mat', id);
    load(AbdAof);
    crop = [105 360 121 376];
    masks = AscAo;
end

lifeDir = [pwd '\' folder '\' scan];
cd(lifeDir);
% untar([scan '.tgz'])
dirFiles = dir('*.dcm');
for i=1:length(dirFiles)
    img = dicomread(dirFiles(i).name);
    img = imresize(img,[desiredSize desiredSize]);
    all(:, :, i) = img;
end
pc = all(:,:,1:40);
mag = all(:,:,41:end);


desiredFrames = size(mag,3);
interpMasks = zeros(size(masks,1),size(masks,2),desiredFrames);
for i=1:size(masks,1)
    for j=1:size(masks,2)
        pixelLine = squeeze(masks(i,j,:));
        x = 1:length(pixelLine);
        xi = linspace(1,length(pixelLine),desiredFrames);
        yi = interp1(x,pixelLine,xi,'spline');
        interpMasks(i,j,:) = yi;
    end 
end 

uncropped = padarray(interpMasks, [crop(1) crop(3)], 'pre');
uncropped = padarray(uncropped, [desiredSize-crop(2)-1 desiredSize-crop(4)-1], 'post');
final = round(uncropped);

for i=1:size(mag, 3)
    %% Show and Save Images
%     label = labeloverlay(mag(:,:,i), final(:,:,i));
%     figure; imshow(label,[]);
    figure; imshowpair(pc(:,:,i),final(:,:,i));
    
    %% To do: flow calculations using masks
end
close all






