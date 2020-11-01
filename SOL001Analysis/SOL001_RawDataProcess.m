
% Process all the raw data and output the bone structure and data structure
close all
clear 
clc

subj_name = 'SOL001B';
subj_dir = 'E:\SOL001_VISIT2\';
% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);


% get the trials that we want
list_files = dir([subj_dir '*run*barefoot*']);
trialname_list = {list_files(:).name};
ntrials = length(trialname_list);


% total number of frames for force and xray
nFrX = 400;
nFrF = 1600;

% framerates for xray and force
frRateX = 250;
frRateF = 1000;
% filtering frame rate
% wc = [10,30];
fcF = [40 80];

% coregistration matrix
coreg = csvread(fullfile(subj_dir,'Calibration','Set 1','end_pylon','end_pylon_COREG.csv'));


%% Re-calculate all the raw data
% load all the motion capture data

load([subj_dir  '\Calibration\Sync\syncFrames.mat'])
struct_data = struct('filename',[],'marker_data',[],'force_data',[],'analog_data',[]);
% load the force data
for i = 1:ntrials
    mocap_find_file = ls(fullfile(subj_dir,trialname_list{i},'Mocap','*_tracked*'));
    if isempty(mocap_find_file)
        fprintf('No Mocap file for %s.\n',trialname_list{i})
        struct_data(i).filename = trialname_list{i};
        continue
    else
    mocap_file = fullfile(subj_dir,trialname_list{i},'Mocap',mocap_find_file);
    struct_data(i) = MocapDataTransform(mocap_file,'filter','adaptive','forceCutOff',fcF,'saveProc','off','lab','SOL','resample','force');
    
        struct_data(i).filename = trialname_list{i};
    end
%     struct_data(i) = MocapDataTransform(mocap_file,'filter','off','saveProc','off','lab','SOL','resample','force');
end




% load all the x-ray data
first_flag = 0;
for i = 1:ntrials
    bone_transform_file = ls(fullfile(subj_dir,trialname_list{i},'BoneTransforms','*transforms2DFILT*'));
    if isempty(bone_transform_file)
        fprintf('No Bone transform file for %s.\n',trialname_list{i})
        continue
    end
    temp_bonesT = load(fullfile(subj_dir,trialname_list{i},'BoneTransforms',bone_transform_file));
    
    struct_data(i).T = temp_bonesT.T; % add the transforms to the data structure
    if first_flag == 0
        bone_list = fields(struct_data(i).T);
        nBones = length(bone_list);
        first_flag = 1;
    end
end


% load all the bones, build the anatomical co-ordinate systems for 
%  cal, mt1
struct_bone = struct('bone_name',bone_list,'pts',[],'cns',[]);
for bn = 1:nBones
    bone_file = fullfile(subj_dir,'Models','IV', ls(fullfile(iv_dir,['*' bone_list{bn} '*aligned*'])) );
    [pts,cns] = read_vrml_fast( bone_file);
    cns = cns+1;
    struct_bone(bn).pts = pts;
    struct_bone(bn).cns = cns;
    
    [cent,~,~,CoM_ev123,CoM_eigenvectors,I1,I2,I_CoM,I_origin,patches] = mass_properties(bone_file);
    struct_bone(bn).cent = cent';
    % build the anatomical pose
    T_init = eye(4);
    T_init(1:3,1:3) = CoM_eigenvectors;
    T_init(1:3,4) = cent';
    
    switch bn
        case 1 % cal
            
            T_new = rotateCoordSys(T_init,90,1);
            T_new = rotateCoordSys(T_new,90,3);
        case 5
            T_new = rotateCoordSys(T_init,-90,1);
            T_new = rotateCoordSys(T_new,90,3);
        case 8
            T_new = rotateCoordSys(T_init,90,3);
            T_new = rotateCoordSys(T_new,180,2);
            
        case 9
            T_new = rotateCoordSys(T_init,-90,3);
            
    end
    struct_bone(bn).T_ACS = T_new;
end    






% load the bone transforms + filter the data
for i = 1:ntrials
    
    if isempty(struct_data(i).T)
       continue
    end
    for bn = 1:nBones
        Tfilt = struct_data(i).T.(bone_list{bn});
        if all(isnan(Tfilt),'all')
            continue
        end
        
        for fr = 1:nFrX
            struct_data(i).Tm.(bone_list{bn})(:,:,fr) = coreg * Tfilt(:,:,fr);
            struct_data(i).Tm.(bone_list{bn})(1:3,4,fr) = struct_data(i).Tm.(bone_list{bn})(1:3,4,fr)/1000;
            % determine the linear position of the centre of mass of the bone
            
            p_com(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr), struct_bone(bn).cent/1000);
                %         p_com = adaptiveLowPassButterworth( p_com,wc,frRateX);
            struct_data(i).pcom.(bone_list{bn})  = p_com;% save it
           
  
        end
        
        % calculate velocity
        struct_data(i).vcom.(bone_list{bn})  = calculateVelocity(p_com,frRateX);
        
        
        % determine the angular velocity
        struct_data(i).w.(bone_list{bn}) = calculateRotMatAngularVelocity(struct_data(i).Tm.(bone_list{bn})(1:3,1:3,:),frRateX,'rad');
        %         struct_data(i).w.(bone_list{bn})= adaptiveLowPassButterworth(struct_data(i).w.(bone_list{bn}),wc,250);
        
        
        
        
    end
    
    struct_data(i).cropFrs = syncFrs.(trialname_list{i}); % add the offset between force and x-ray data
    
    force1 = struct_data(i).force_data(1).Force;
    force2 = struct_data(i).force_data(2).Force;
    normForce = norm3d(force1+force2);
    %
    %     ind = find((force1(3,1:1000) > 100) | (force2(3,1:1000) > 100));
    %     s1 = ind(1) -20; e1 = ind(end) + 20;
    %     ind = find((force1(3,s1:e1) > 20) | (force2(3,s1:e1) > 20));
    %     s2 = s1 + ind(1)-2; e2 = s1 + ind(end)+1;
    
    ind = find(normForce(1:1000) > 200) ;%| (force2(3,1:1000) > 100));
    s1 = ind(1) -20; e1 = ind(end) + 20;
    ind = find(normForce(s1:e1) > 20);% | (force2(3,s1:e1) > 20));
    s2 = s1 + ind(1)-2; e2 = s1 + ind(end)+1;
    struct_data(i).cropFrsForce = [s2 e2];
    
    struct_data(i).cropFrsXray = [s2 e2] - syncFrs.(trialname_list{i}) -250; % add the offset between force and x-ray data
    
    figure
%     plot(1:500,force2(3,1:500)')
    hold on;
%     plot(1:500,force1(3,1:500)')
    plot(1:500,normForce(1:500)')
    plot([s2 e2],normForce([s2 e2]),'o')
    plot([s2 e2],normForce([s2 e2]),'o')
    
    fprintf('Trial %s added to the data structure. \n',trialname_list{i})
    
end

save([subj_dir 'struct_bone.mat'],'struct_bone');
save([subj_dir 'struct_data.mat'],'struct_data');
