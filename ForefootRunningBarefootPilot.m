% analyse Chris' mocap from forefoot running in HMRL Sept 2018
close all
clear
clc
cd('C:\Users\Lauren\Documents\School\PhD\Research\AllStudies\3B_BarefootRunning\BarefootRunning\Data\')

% before this, ran MocapDataTransform
bw = 702;
for i = 3%,3]
    switch i
        case 1
            load('BFWp001_norm_bf_ffs_fore0001_processedmotion.mat')
            f1(i) = 263;
            f2(i) = 292;
        case 2
            load('BFWp001_norm_bf_ffs_full0003_processedmotion.mat')

            f1 = 540;
            f2 = 620;
        case 3
%             load('BFWp001_norm_bf_ffs_toes0002_processedmotion.mat')

            load('BFWp001_norm_bf_ffs_fore0003_processedmotion.mat')
            f1(i) = 387-50;
            f2(i) = 417+50;
    end
    
    nfrM = marker_data.nFrames;
    
    % zero the force plates and resample over mocap
%     force_temp = resample(force_back_zero',125,round(3/4*125))';
    force_back_nonzero = resample(force_data(3).Force',125,625)';
    force_front_nonzero = resample(force_data(2).Force',125,625)';
%     force_back_nonzero =  force_data(4).Force;
%     force_front_nonzero =force_data(3).Force;
    
    nfrF = length(force_back_nonzero); % force nframes
    
    force_back_zero = force_back_nonzero - repmat(mean(force_back_nonzero(1:3,1:5),2),1,nfrF);
    force_front_zero = force_front_nonzero - repmat(mean(force_front_nonzero(1:3,1:5),2),1,nfrF);
    
    force_back = force_back_zero;%normalise(force_back_zero',nfrM-1)';
    force_front = force_front_zero;%normalise(force_front_zero',nfrM-1)';
    
    % create the foot model
    
    data_struct(i) = createFootModel(marker_data);
    force_struct(i).force_back = force_back/bw; 
    
    force_struct(i).force_front = force_front/bw;
    if i == 2
    force_struct(i).force_front = [];
    end

% figure
subplot(4,1,1);hold on;
    plot_data_1 =  data_struct(i).marker_data.MH1(3,f1(i):f2(i)) ;%data_struct(i).F2Ps(f1(i):f2(i));
%     nans = isnan(plot_data_1);
%     plot_data_1(nans) = 0;
plot(normaliseNaN(plot_data_1,2))
title('Toe angle')
xlim([0,100])

subplot(4,1,2);hold on;
    plot_data_2 = data_struct(i).marker_data.CA_(3,f1(i):f2(i)); %data_struct(i).sagittal_arch.flex_ankle(f1(i):f2(i));
%     nans = isnan(plot_data_2);
%     plot_data_2(nans) = 0;
plot(normaliseNaN(plot_data_2,2))

title('Ankle Angle')
xlim([0,100])
subplot(4,1,3);hold on;
plot(normalise(force_struct(i).force_back(3,f1(i):f2(i))' )')
xlim([0,100])
if i ~= 2
subplot(4,1,4);hold on;
plot(normalise(force_struct(i).force_front(3,f1(i):f2(i))' )')
title('Force (Front)')
xlim([0,100])
end
end




% 
% figure;
% 
% for fr = f1:5:f2
% hold on
%     plotPointsAndCoordSys1([],data_struct(i).pose.foot(:,:,fr),100,'k')
%     hold on
% plotPointsAndCoordSys1([],data_struct(i).pose.shank(:,:,fr),100,'k')
% end
%%
% % close all
% figure;
% 
% for fr = f1:5:f2
% hold on
%     plotPointsAndCoordSys1(data_struct(i).marker_data.ME_(:,fr),data_struct(i).pose.foot(:,:,fr),100,'k');
%     hold on
% plotPointsAndCoordSys1(data_struct(i).marker_data.LE_(:,fr),data_struct(i).pose.shank(:,:,fr),100,'k');
% hold on
% plot3quick([data_struct(i).marker_data.MH1(:,fr),data_struct(i).marker_data.MH5(:,fr),data_struct(i).marker_data.MB1(:,fr),data_struct(i).marker_data.MB5(:,fr),data_struct(i).marker_data.CA_(:,fr),data_struct(i).marker_data.MM_(:,fr),data_struct(i).marker_data.LM_(:,fr)],'k','o');
% 
% end
% axis equal
%% Using Bruening's Data


load('/home/lauren/Documents/VirtualDesktopTemp/PhDResearch/BarefootRunning/ffs_dors.mat')
load('/home/lauren/Documents/VirtualDesktopTemp/PhDResearch/BarefootRunning/ffs_toe.mat')
load('/home/lauren/Documents/VirtualDesktopTemp/PhDResearch/BarefootRunning/ffs_force.mat')
all3 = [dors_angle;force_val;toe_angle];
plot3quick(all3,'k','x');
xlabel('Ankle angle (dors +)')
ylabel('Force val (BW)')
zlabel('toe angle (dors +)')
figure;

toe_sort = sortrows(all3',3);
%  plot(toe_sort(:,3),toe_sort(:,2),'x')
 ankle_sort = sortrows(all3',1);
 
 cmap = colormap('jet');
 c =  interp1(linspace(1,100,64),cmap,0:100);
 mdl = fitlm(linspace(-30,30,100),1:100);
 
xval = @(x) round(mdl.Coefficients{1,1} + x * mdl.Coefficients{2,1});


hold on

 for i = 1:101
     c_ind = xval(toe_sort(i,1));
 h = scatter(toe_sort(i,3),toe_sort(i,2),'Marker','square','SizeData',400,'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceColor',c(c_ind,:));
 end
 colorbar
 caxis([min(toe_sort(:,1)) max(toe_sort(:,1))])
 
xlabel('Toe Angle (Dorsiflexion +)')
ylabel('Force to the forefoot (BW)' )
title('Ankle angle (Dors +) ')


figure;
% subplot(3,1,1)
 plot(0:100,toe_angle)
xlabel('% Stance')
ylabel('Toe Angle (Dorsiflexion +)')
% subplot(3,1,2)

figure;
plot(0:100,force_val);
xlabel('% Stance')
ylabel('Force to the forefoot (BW)' )
% subplot(3,1,3)

figure;

 hold on
 for i = 1:101
 c_ind = xval(dors_angle(i));
scatter(i,dors_angle(i),'Marker','square','SizeData',100,'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceColor',c(c_ind,:));
 end
xlabel('% Stance')
ylabel('Ankle angle (Dors +) ')