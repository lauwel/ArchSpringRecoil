% Assess the stiffness and arch height from the pictures in the PRFC study


% points required for measurement are posterior heel, toe, met head, arch
% height at half of foot length, and floor

base_path = 'C:\Users\Lauren\Documents\School\PhD\Research\AllStudies\3_PRFC\Data\';
list_subj_files = dir([base_path '*PRFC*']);

%  dataStruct = struct('subject','','AHI',[],'RAD',[]);

for i = 1:length(list_subj_files)
sit_string = '*it*';
stand_string = '*tand*';

assess_path = fullfile(base_path,list_subj_files(i).name,'Initial Assessment',filesep);

path_sit = fullfile(base_path,list_subj_files(i).name ,sit_string);
path_stand = fullfile(base_path,list_subj_files(i).name ,stand_string);

list_sit = dir([assess_path,sit_string]);
list_stand = dir([assess_path,stand_string]);

subj_str = list_subj_files(i).name;
% subj_ind = str2double(list_subj_files(i).name(end-2:end));
imgSit = imread(fullfile(assess_path,list_sit(1).name));
imgStand = imread(fullfile(assess_path,list_stand(1).name));

[footLengthSit, truncfootLengthSit, arch_heightSit] = calcArchParams(imgSit);
% [footLengthStand, truncfootLengthStand, arch_heightStand] = calcArchParams(imgStand);

AHI = arch_heightSit/truncfootLengthSit;
% RAD = (arch_heightSit-arch_heightStand)/arch_heightSit;

dataStruct(i).subject = subj_str;
dataStruct(i).AHI = AHI;
% dataStruct(i).RAD = RAD;
% dataStruct(i).footLength = [footLengthSit footLengthStand];
% dataStruct(i).truncFootLength = [truncfootLengthSit truncfootLengthStand];
% dataStruct(i).archHeight = [arch_heightSit arch_heightStand];

end

function [footLength, truncfootLength, arch_height] = calcArchParams(img)
[ht, w,~] = size(img);
rat = ht/w;
hf = figure;
set(hf,'units','normalized','color',[ 1 1 1],'Position',[1 0.4 rat 1])
ha_text = axes(gcf,'Position',[0.1 0.95 0.8 0.05]);
ha_img = axes(gcf);%,'Position',[0.1 0.1 0.5 0.5*ht/w]);
imagesc(ha_img,img);
axis image

hold on;
messages = {'Select heel.', 'Select toe.', 'Select met head.',...
    'Select arch height on line.', 'Select floor on line.',...
    'Select measurement point 1.', 'Select measurement point 2.' ,...
    'Input mm distance of selected points.'};

% select the four corners of the plot
for i = 1:7
      
    if ismember(i,[4,5])
        x_half = sum(x_lims(1:2))/2;
%         y_half = sum(y_lims(1:2))/2;
           plot(ha_img,[x_half x_half],[1 ht-1],'color','k','linestyle','--')
    end
    % display the messages
    h = text(ha_text,0,0.5,messages{i},'fontsize',20);
    set(ha_text,'visible','off') % turn off the axes
    
    % get the selected points
    [x_lims(i),y_lims(i)] = ginput(1);
    
    % get rid of the text
    delete(h)
    
    % draw the selected points
    scatter(ha_img,x_lims(i),y_lims(i),'wx')
    drawnow
  
end

mmLengthCell = inputdlg(messages{8},'Scale size');
mmLength = str2double(mmLengthCell{1});
scale_PX = norm(diff([x_lims(6:7); y_lims(6:7)],1,2));
scale_PXtoMM = mmLength/scale_PX;

x_mm = x_lims * scale_PXtoMM;
y_mm = y_lims * scale_PXtoMM;

footLength = norm(diff([x_mm(1:2); y_mm(1:2)],1,2));
truncfootLength = norm(diff([x_mm([1,3]); y_mm([1,3])],1,2));
arch_height = norm(diff([x_mm([4,5]); y_mm([4,5])],1,2));
end




