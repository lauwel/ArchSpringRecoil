close all
clear
clc


load('P:\Personnel\LaurenWelte\Research\Lauren_FootModel\nwalk0001.mat')
load('P:\Personnel\LaurenWelte\Research\Lauren_FootModel\hmrl061_qtm_norm_swalk_0004.mat')

force_half1 = nwalk0001.Force(2).Force;
force_half2 = nwalk0001.Force(1).Force;

force_full = hmrl061_qtm_norm_swalk_0004.Force(1).Force;

zero_1 = mean(force_half1(1:3,1:10),2);
zero_2 = mean(force_half2(1:3,1:10),2);
zero_full = mean(force_full(1:3,1:10),2);

force_half1 = force_half1 -zero_1;
force_half2 = force_half2 -zero_2;
force_full = force_full -zero_full;

ind_start1 = find(force_half1(3,:) > 10);
ind_start2 = find(force_half2(3,:) > 10);

force_half1 = force_half1(1:3,ind_start1(1):ind_start2(end));
force_half2 = force_half2(1:3,ind_start1(1):ind_start2(end));

ind_start = find(force_full(3,:) > 10);
force_full = force_full(1:3,ind_start(1):ind_start(end));
force_full(1:2,:) = force_full(1:2,:);

force_fullhalf = force_half1 + force_half2;
force_full_norm = normaliseNaN(force_full',1,100);
force_split_norm = normaliseNaN(force_fullhalf',1,100);
for i = 1:100
        force_full_mag(i) =  norm(force_full_norm(i,1:3));
        force_split_mag(i) =  norm(force_split_norm(i,1:3));
end
figure(1)
plot(force_half1')

figure(2)
plot(force_half2')

figure(3)
plot(force_full')

figure(4)
plot(force_full_norm)
hold on;
plot(force_split_norm)
legend('x','y','z','xsplit','ysplit','zsplit')


figure(5)
plot(force_full_mag)
hold on
plot(force_split_mag)
legend('full','split')






%% 
close all 
clear 
clc
data1 = importdata('P:\Personnel\LaurenWelte\ForcePlateValidation\GRF_plotData.txt');


% first column is subject, second is fx, fy, fz, third column is reg, target, split
tag_matrix = nan(350,3);

for i = 2:350
    if regexp(data1.textdata{i,1},'lw_')
        tag_matrix(i,1) = 3;
    elseif regexp(data1.textdata{i,1},'lh_')
        tag_matrix(i,1) = 4;
    elseif regexp(data1.textdata{i,1},'ao_')
        tag_matrix(i,1) = 1;
    elseif regexp(data1.textdata{i,1},'ad_')
        tag_matrix(i,1) = 2;
    else
        tag_matrix(i,1) = NaN;
    end

    if regexp(data1.textdata{i,4},'FX')
        tag_matrix(i,2) = 1;
    elseif regexp(data1.textdata{i,4},'FY')
        tag_matrix(i,2) = 2;
    elseif regexp(data1.textdata{i,4},'FZ')
        tag_matrix(i,2) = 3;
    else
        tag_matrix(i,2) = NaN;
    end

    if regexp(data1.textdata{i,4},'REG')
        tag_matrix(i,3) = 1;
    elseif regexp(data1.textdata{i,4},'target')
        tag_matrix(i,3) = 2;
    elseif regexp(data1.textdata{i,4},'split')
        tag_matrix(i,3) = 3;
    else
        tag_matrix(i,3) = NaN;
    end

end


% [B,ind] = sortrows(tag_matrix,2)
figure(1);
hold on;
cols = colormap('jet');
c = cols(1:4:end,:);
figure(2);
hold on
figure(3)
hold on
dat_org = cell(12,3);
dat_avg = cell(12,3);
dat_std = cell(12,3);

for i = 2:350
    ind = (tag_matrix(i,1)-1)*3 + tag_matrix(i,3);
    switch tag_matrix(i,2)
        case 1
            figure(1)
            h = plot(data1.data(i,:));
            h_save(ind,1) = h;
            
            dat_org{ind,1}(end+1,:) = data1.data(i,:);
       
        case 2
            figure(2)
            h = plot(data1.data(i,:));
            h_save(ind,2) = h;
            
            dat_org{ind,2}(end+1,:) = data1.data(i,:);
 
        case 3
            figure(3)
            h = plot(data1.data(i,:));
            h_save(ind,3) = h;
            
            dat_org{ind,3}(end+1,:) = data1.data(i,:);
     
    end
    
    set(h,'color',c(ind,:))
    
end

for i = 1:3
    figure(i)
    legend(h_save(:,i),{'s1 - normal','s1- target','s1 - split','s2 - normal','s2- target','s2 - split','s3 - normal','s3- target','s3 - split','s4 - normal','s4- target','s4 - split'})
end

%% get the mean difference and std dev
for i = 1:12
    for j = 1:3
    dat_avg{i,j} = mean(dat_org{i,j},1);
    dat_std{i,j} = std(dat_org{i,j},[],1);
    end
end

j = 2;

figure(j+3);
hold on
for i = 1:12
    PrettyStdDevGraphs(0:100,dat_avg{i,j},dat_std{i,j},c(i,:),1)
    hold on
end

% figure(j+6)
norm = [1,4,7,10];
target = [2,5,8,11];
split = [3,6,9,12];
yval = 120;
for i = 1:4
    figure(1)
    
    subplot(4,1,i)
    hold on
    plot(0:100,dat_avg{target(i),j}-dat_avg{norm(i),j},'color',c(target(i),:))
    PrettyStdDevGraphs(0:100,zeros(1,101),dat_std{norm(i),j},c(target(i),:),3)
    ylim([-yval  yval ])
    xlabel('% stance')
    ylabel('Force [N]')
    title('Fy mean differences, targeted')
    figure(2)
    
    subplot(4,1,i)
    hold on
    plot(0:100,dat_avg{split(i),j}-dat_avg{norm(i),j},'color',c(split(i),:))
    PrettyStdDevGraphs(0:100,zeros(1,101),dat_std{norm(i),j},c(split(i),:),3)
    ylim([-yval  yval ])
    xlabel('% stance')
    ylabel('Force [N]')
    title('Fy mean differences, split force plates')
end





