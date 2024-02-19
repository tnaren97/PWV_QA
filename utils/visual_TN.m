% load('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00003_10781_2019-11-08-h07\00009.FIESTA_Cine_BH_AAo\AscAo_mask.mat');
% load('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00003_10781_2019-11-08-h07\00009.FIESTA_Cine_BH_AAo\DescAo_mask.mat');
% bssfp = read_dicoms('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00003_10781_2019-11-08-h07\00009.FIESTA_Cine_BH_AAo');
% pwv = read_dicoms('C:\Users\naren\Desktop\Aorta_Segmentation\life_life00003_10781_2019-11-08-h07\00007.PWV_CartBH_AAo');
crop = [155 410 101 356];
% masks = AscAo_mask + DescAo_mask;
% load('C:\Users\naren\Desktop\Aorta_Segmentation\results\16-11_AbdAo.mat');
% lifeDir = 'C:\Users\naren\Desktop\Aorta_Segmentation\life_life00016_11110_2020-03-09-h08\00011.FIESTA_Cine_BH_AbdAo';
% crop = [105 360 121 376];
% masks = AscAo;
% bssfp = imresize3(bssfp, [256 256 40]);
% cd(lifeDir);
load('vol_Dan\high_res\048_PWV_CartBH_AAo_PROSP_1VPS_98FRAMES\interp_masks.mat')
% dirFiles = dir('*.dcm');
% bssfp = read_dicoms('C:\Users\naren\Desktop\Aorta_Segmentation\vol_Dan\high_res\038_FIESTA_Cine_BH_AAo_PROSP_2VPS_283FRAMES');
load('vol_Dan\high_res\048_PWV_CartBH_AAo_PROSP_1VPS_98FRAMES\mag.mat')
mag = imresize3(mag,[desiredSize desiredSize desiredFrames]);
for i=1:49
%      image = dicomread(dirFiles(i+49).name);
% %     x = dicominfo(dirFiles(i).name);
% %     image = image(crop(1):crop(2),crop(3):crop(4));
% %     image = imresize(image, [512,512]);
%      image = imadjust(image);
     im1 = masks(:,:,i);
     im2 = mag(:,:,i);
     figure; imshowpair(imtranslate(imresize(im1,0.99), [29, -4]), im2);
%      imtranslate(im1, [25, -11])
%      figure; imshow(imasjuim1);
    
%     disp(x.PixelSpacing(1).*x.PixelSpacing(2));
    %% Show and Save Images
%     label = labeloverlay(image, masks(:,:,i));
%     figure; imshow(label,[]);
%      figure; imshowpair(image, masks(:,:,i));
%      figure; imshowpair(image, pcmasks(:,:,i));
%      if i == 10
%         break
%      end
%      f.WindowState = 'maximized';
end
close all