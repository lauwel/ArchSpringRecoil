close all
clear
load('SOL001_qtm_norm_jog0001_processedmotion.mat')
m{1}(1,:) = marker_data.r_lowback_S(2,1:end-20);
m{1}(2,:) = marker_data.l_lowback_S(2,1:end-20);
load('SOL001_qtm_norm_jog0002_processedmotion.mat')
m{2}(1,:) = marker_data.r_lowback_S(2,1:end-20);
m{2}(2,:) = marker_data.l_lowback_S(2,1:end-20);
load('SOL001_qtm_norm_jog0004_processedmotion.mat')
m{3}(1,:) = marker_data.r_lowback_S(2,1:end-20);
m{3}(2,:) = marker_data.l_lowback_S(2,1:end-20);
figure;
hold on;
for i = 1:3
    
vel{i} = diff(m{i}')*250;
mean_vel(i,:)= mean(vel{i})/1000;
std_vel(i,:) = std(vel{i})/1000;

plot(vel{i})
end

vel = [];
for i = 1:3
    
vel = [vel;diff(m{i}')*250];

end
mean_vel= mean(vel(:))/1000;
std_vel = std(vel(:))/1000;



