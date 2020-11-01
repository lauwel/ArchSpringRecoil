% write XMA bead positions to transforms

% bead positions are contained in a file exported from the
% sphereFitToBeads.m

close all
clear
clc
trial_name = 'T0034_SOL001_srun_pref_barefoot';
subj_dir = 'E:\SOL001_VISIT2\';

filter_flag = 0;
autoscoper_animation_flag = 1; % make animations of the autoscoped data

fc = [20 60]; % cut-off frequency if filtering
wc = [20 60];
xmadir = [subj_dir trial_name '\XMA_csv\'];
xmafile = dir([xmadir '*3D*nterp.csv']);
xmafile2D = dir([xmadir '*MH1*2D*.csv']);
subj_name = 'SOL001B';

% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\'; % if saving animation to
% server
ivdir = [subj_dir 'Models\IV\'];
trialdir = [subj_dir trial_name '\'];

%% Interpolate the 2D xma points if necessary
interp_flag = 0;
if interp_flag == 1
    interpolateXMA2DPoints([xmadir xmafile2D(1).name])
    
end




%% get the bead locations and save into a structure
ivDir = 'E:\SOL001_VISIT2\Models\IV\Beads\';

datBeadLocs = csvread('E:\SOL001_VISIT2\Models\bead_positions.txt',1,1);

dat = readtable('E:\SOL001_VISIT2\Models\bead_positions.txt');

beadStructTemp = table2struct(dat);
nbeads = length(beadStructTemp);
for b = 1:nbeads
    
    fieldname = beadStructTemp(b).bead;
    beadnum = str2double(fieldname(4));
    beadCT.(fieldname(1:3))(beadnum,:) = datBeadLocs(b,:);
    
end

 

%% load the XMA file containing the bead co-ordinates
i = 1;

[beadPos,names] = xlsread([xmadir xmafile(i).name]); % load the file
beadPos = beadPos(:,2:end)'; % get rid of the first column and transpose
names = names(1,2:end)'; % get rid of first column and transpose
nfr = size(beadPos,2);
% save all the beads into the appropriate structure location
for b = 1:nbeads*3
    beadname = names{b}(2:4);
    beadnum = str2double(names{b}(5));
    if filter_flag == 1
        if ismember(b,87:89)
            pflag = 1;
        else 
            pflag = 0;
        end
        if sum(~isnan(beadPos(b,:))) < 5
            
            tempBead = beadPos(b,:);
        else
            tempBead = adaptiveLowPassButterworth(beadPos(b,:),fc,250,pflag);
        end
    else
        tempBead = beadPos(b,:);
    end
    switch names{b}(7)
        case 'X'
            beadXMA.(beadname)(beadnum,1,:) = tempBead;
        case 'Y'
            beadXMA.(beadname)(beadnum,2,:) = tempBead;
        case 'Z'
            beadXMA.(beadname)(beadnum,3,:) = tempBead;
            
        otherwise
            error('Mislabelled bead.')
    end
end

%% match the bead positions to the CT to get the transform
T_save = [];
bonesCell = fields(beadXMA);

nbones = length(bonesCell);
% tracked_frames = nan(nbones,nfr);
for bn = 1:nbones % for every bone, compute the transform between the beads and xma points
    bonebeadsCT = beadCT.(bonesCell{bn});
    ind = 1;
    for fr = 1:nfr
        % get the transform from CT to x-ray space
        T.(bonesCell{bn})(:,:,fr) = eye(4,4);
        Tanim.(bonesCell{bn})(:,:,fr) = eye(4,4);
        bonebeadsXMA = beadXMA.(bonesCell{bn})(:,:,fr);
        ind_nonan = ~isnan(bonebeadsXMA(:,1));
        
        if sum(ind_nonan) >= 3 % if there are at least 3 beads with data
            
            [R,d,rms] = soder(bonebeadsCT(ind_nonan,:), bonebeadsXMA(ind_nonan,:));
            
            T.(bonesCell{bn})(1:3,1:3,fr) = R;
            T.(bonesCell{bn})(1:3,4,fr) = d;
            
            Tanim.(bonesCell{bn})(1:3,1:3,fr) = R;
            Tanim.(bonesCell{bn})(1:3,4,fr) = d;
            
            rms_save(bn,fr) = rms;
            
            tracked_frames(bn,ind) = fr;
            ind = ind+1;
        else
            T.(bonesCell{bn})(:,:,fr) = nan(4,4);
            Tanim.(bonesCell{bn})(:,:,fr) = ones(4,4);
        end
        
    end
    
    % write the file with the bone transforms
%     q = convertRotation(T.(bonesCell{bn}),'4x4xn','quaternion');
%     q_filt = adaptiveLowPassButterworth(q',wc,250)';
%     Tstacked = convertRotation(q_filt,'quaternion','autoscoper');
    Tstacked.(bonesCell{bn}) = convertRotation(T.(bonesCell{bn}),'4x4xn','autoscoper');
    
    TstackedR.(bonesCell{bn}) = convertRotation(T.(bonesCell{bn}),'4x4xn','autoscoperRows');
%     
    
    autoscoper_dir = [trialdir 'Autoscoper\'];
    if exist(autoscoper_dir,'dir') == 0
        mkdir(autoscoper_dir);
    end
    if filter_flag == 1
        auto_file = [autoscoper_dir trial_name '_' (bonesCell{bn}) '_filt.tra'];
    else
        auto_file = [autoscoper_dir trial_name '_' (bonesCell{bn}) '_unfilt.tra'];
    end
    dlmwrite(auto_file,Tstacked.(bonesCell{bn}))
    
    T_save = [T_save TstackedR.(bonesCell{bn})];
end



% Save the ALL bones file - in rows
if filter_flag == 1
    auto_file = [autoscoper_dir 'allBones_' trial_name '_rows_filt.tra'];
else
    auto_file = [autoscoper_dir 'allBones_' trial_name '_rows_unfilt.tra'];
end
dlmwrite(auto_file,T_save)


%% Animate the trial

if filter_flag == 1
anim_folder = strcat(trialdir,'POS\Filtered\');
else
    
anim_folder = strcat(trialdir,'POS\Unfiltered\');
end
rigidiv_folder = fullfile(anim_folder,'rigidiv',filesep);

if exist(anim_folder,'dir')==0;     mkdir(anim_folder);     end
if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end


first_fr = nanmin(nanmin(tracked_frames(:,:)));
if first_fr == 0
    first_fr = 1;
end
end_fr   = nanmax(nanmax(tracked_frames(:,:)));

for bn = 1:nbones
    
    ivstring = createInventorHeader();
    % make the linked iv file
    ivstring = [ivstring createInventorLink([ivdir subj_name '_'  bonesCell{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
    
    fid = fopen(fullfile(rigidiv_folder,[bonesCell{bn} '.iv']),'w');
    fprintf(fid,ivstring);
    fclose(fid);
end


for bn = 1:nbones
    write_RTp(bonesCell{bn} , Tanim.(bonesCell{bn})(:,:,first_fr:end_fr) , anim_folder)
end

pos_text = write_pos(bonesCell,anim_folder,trial_name);

    filename = fullfile(anim_folder, [trial_name '.pos']);



fid = fopen(filename,'w'); % open the file to write
fprintf(fid,pos_text);
fclose(fid);

%% plots the RMS error for each frame for each bone
figure;
cmap = colormap('jet');

plot(rms_save(:,tracked_frames(3,:)),'.','Markersize',15)
hold on;
for i = 1:size(tracked_frames,2)
    
    for bn = 1:10
        if tracked_frames(bn,i) == 0
            continue
        end
        text(bn,rms_save(bn,tracked_frames(bn,i)),num2str(tracked_frames(bn,i)))
    end
end
title('RMS values at tracked frames')
ylabel('RMS Error')

hold on;
plot([0 10.5], [0.15 0.15],'k--')
xlim([0 10.5])
set(gca,'Xtick',1:10)
set(gca,'Xticklabel',[bonesCell])
% legend(num2str(tracked_frames(1,:)'))

