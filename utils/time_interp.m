% change pcDir
pcDir = 'C:\Users\naren\Desktop\Aorta_Segmentation\life_life00047_11454_2021-02-25-h08\00009.PWV_CartBH_AAo';
cd(pcDir);
d = dir('*.dcm');
numPCframes = 40;
PCtimes = zeros(1,numPCframes);
for i=1:length(PCtimes)
temp = dicominfo(d(i).name);
PCtimes(i) = temp.TriggerTime;
end



% change bssfpDir
bssfpDir = 'C:\Users\naren\Desktop\Aorta_Segmentation\life_life00047_11454_2021-02-25-h08\00010.FIESTA_Cine_BH_AAo';
cd(bssfpDir);
d = dir('*.dcm');
numBSSFPframes = 30;
BSSFPtimes = zeros(1,numBSSFPframes);
for i=1:length(BSSFPtimes)
temp = dicominfo(d(i).name);
BSSFPtimes(i) = temp.TriggerTime;
end



% interpolate the times like we would for the images (30-->40)
BSSFP_interp = interp1((1:numBSSFPframes),BSSFPtimes,linspace(1,numBSSFPframes,numPCframes),'linear');
figure; plot(PCtimes);
hold on; plot(BSSFP_interp);
xlabel('Frames'); ylabel('Time (ms)');
legend('PC', 'BSSFP')
