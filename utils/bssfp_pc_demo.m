PCdir = 'C:\Users\naren\Desktop\Aorta_Segmentation\life_life00013_11099_2020-03-05-h09\00007.PWV_CartBH_AAo';
cd(PCdir)
dirFiles = dir('*.dcm');
for i=1:length(dirFiles)
    image = dicomread(dirFiles(i).name);
    temp(:,:,i) = image;
end
pc = temp(:,:,1:40);
mag = temp(:,:,41:end);
desiredSize = 512; %upsample to match resolution of bssfp
pc = imresize(pc,[desiredSize desiredSize]);
mag = imresize(mag,[desiredSize desiredSize]);
clear temp PCdir desiredSize

ANATdir = 'C:\Users\naren\Desktop\Aorta_Segmentation\life_life00013_11099_2020-03-05-h09\00008.FIESTA_Cine_BH_AAo';
cd(ANATdir);

dirFiles = dir('*.dcm');
for i=1:length(dirFiles)
    image = dicomread(dirFiles(i).name);
    temp(:,:,i) = image;
end
temp = double(temp);

desiredFrames = size(mag,3);
bssfp = zeros(size(temp,1),size(temp,2),desiredFrames);
for i=1:size(temp,1)
    for j=1:size(temp,2)
        pixelLine = squeeze(temp(i,j,:));
        x = 1:length(pixelLine);
        xi = linspace(1,length(pixelLine),desiredFrames);
        yi = interp1(x,pixelLine,xi,'spline');
        bssfp(i,j,:) = yi;
    end 
end 
bssfp = int16(bssfp);
clear ANATdir ans desiredFrames pixelLine temp x xi yi i j

for f=1:size(bssfp,3)
    figure; imshowpair(mag(:,:,f),bssfp(:,:,f));
end 

% close all
