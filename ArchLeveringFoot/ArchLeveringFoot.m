% Arch Levering Foot
close all
clear all
clc
dir_analy = 'E:\ArchSpring_Gearing\';%'C:\Users\Lauren\OneDrive - Queen''s University\Research\Projects\2019 Beaded Mechanisms\Arch Levering\';
bone_list = {'cal','mt1','tal','tib'};
nBones = length(bone_list);

dir_ref = 'E:\ReferenceBones\';% fullfile(dir_analy,'References',filesep);


data_filename = '_filt10_30';
Fc = [10 30];
%% Load/Create data structures as required - load trialStruct and create subjStruct


demo.weight = [60.0,61.7,76.0,77.5,76.3,63.3,72.3,83];
demo.height = [1.73, 1.66, 1.68, 1.735, 1.815, 1.6, 1.8, 1.75];
demo.age = [22, 21, 29, 23, 25, 23, 20, 51];
demo.gender = {'F','F','F','M','M','F','M','M'};


load([dir_analy 'trialStruct' data_filename '.mat'])



if exist([dir_analy 'subjStruct' data_filename '.mat'],'file')
    load([dir_analy 'subjStruct' data_filename '.mat']);
    nsubj = length(subjStruct);
     ntrials = length(trialStruct);
else
    


    subj_list = unique({trialStruct(:).subject})';
    nsubj = length(subj_list);
    ntrials = length(trialStruct);
    
    
    
    % this gives you the corresponding subject information so if there are
    % multiple files in the trialStruct, it doesn't require repeated bone info
    for s = 1:nsubj
        subjStruct(s).subject = subj_list{s};
        subjStruct(s).trialNums = findInStruct(trialStruct,'subject',subj_list{s}); % indices of trials in TrialStruct
    end
    
    % Load the bone data structure, create distance fields if needed
    
    for s = 1:nsubj
        
        % take the first trial in the to get the necessary info in
        % trialStruct
        t1 = subjStruct(s).trialNums(1);
        subj_dir = trialStruct(t1).metadata.subj_dir;
        
        load(fullfile(subj_dir,'Models','IV','3Aligned Reduced','bonestruct.mat'))
        
        subjStruct(s).bones = bonestruct;
        
        % create or load distance fields
        for b = 3:4 % tal and tib in bone_list
            
            if ~isfield(subjStruct(s).bones.(bone_list{b}),'dfield') % if the dfield is not in the structure
                dfield_fileInfo = dir(fullfile(dir_analy,'dfields',[subjStruct(s).subject '*' (bone_list{b}) '*']));
                
                
                if ~isempty(dfield_fileInfo) % check if a dfield has already been made
                    dfield_filepath = fullfile(dfield_fileInfo.folder, dfield_fileInfo.name);
                    load(dfield_filepath);
                    subjStruct(s).bones.(bone_list{b}).dfield = Dfield;
                else % otherwise make a new distance field
                    fprintf('Creating dfield based on : %s \n',subjStruct(s).bones.(bone_list{b}).metadata.orig_file)
                    dfield = createDfield(subjStruct(s).bones.(bone_list{b}).metadata.orig_file, 0.25, fullfile(dir_analy,'dfields'));
                    subjStruct(s).bones.(bone_list{b}).dfield = dfield;
                    
                end
            end
        end
    end
    
    save([dir_analy 'subjStruct' data_filename '.mat'],'subjStruct') % save the subjStruct
    
end



% for each trial, set all the values that are the same to be NaN's
for t = 1:ntrials
    bones_tr = fields(trialStruct(t).Tx_mm);
    for b = 1:length(bones_tr)
        %             (diff(squeeze(trialStruct(t).Tm_mm.cal(1,1,[1:end,end])))~=0) & (~isnan(squeeze(trialStruct(t).Tm_mm.cal(1,1,:))));
        
        ind_starts = find(diff(squeeze(trialStruct(t).Tx_mm.(bones_tr{b})(1,4,:))) ~= 0);
        Is = [ind_starts(1), ind_starts(end)+1];
%         figure;
%         plot(squeeze(trialStruct(t).Tx_mm.(bones_tr{b})(1,4,:)))
%         hold on
%         plot(Is,squeeze(trialStruct(t).Tx_mm.(bones_tr{b})(1,4,Is)),'o')
        trialStruct(t).Tx_mm.(bones_tr{b})(:,:,[1:Is(1)-1,Is(2)+1:end]) = NaN;
        trialStruct(t).Tx_m.(bones_tr{b})(:,:,[1:Is(1)-1,Is(2)+1:end]) = NaN;
        
        trialStruct(t).Tm_mm.(bones_tr{b})(:,:,[1:Is(1)-1,Is(2)+1:end]) = NaN;
        trialStruct(t).Tm_m.(bones_tr{b})(:,:,[1:Is(1)-1,Is(2)+1:end]) = NaN;
%         plot(squeeze(trialStruct(t).Tx_mm.(bones_tr{b})(1,4,:)))
    end
end

%% associate a static mocap file with each trial
close all
pelvis_marker_names = {'r_ASIS','r_PSIS','r_pelvis','l_ASIS','l_PSIS','l_pelvis'};

for t = 1:ntrials
    subj_dir = trialStruct(t).metadata.subj_dir;
    
    static_info = dir([subj_dir '\**\Mocap\*static*mat']);
    static_folder = [static_info.folder '\'];
    static_file = static_info.name;
    
    frs = 10:20; % average 10 frames for averaging of segment
    
    struct_data = MocapDataTransform(fullfile(static_folder,static_file),'filter','adaptive','forceCutOff',Fc,'mocapCutOff',Fc,'saveProc','off','lab','SOL','resample','force');
    md_stat = struct_data.marker_data;
    
    
    marker_name = fields(md_stat);
    marker_name(1:2) = []; % get rid of the framerate and nframes variables
    % Average all the markers over the selected frames
    for m = 1:length(marker_name)
        
        md_stat.(marker_name{m}) = mean(md_stat.(marker_name{m})(:,frs),2);
        
    end
    
    [T_pelvis,T_pelvis_R,T_pelvis_L,hipC_R,hipC_L] = makePelvisCoordSys(md_stat);
    struct_data(t).full_body.static.markers = md_stat;
    struct_data(t).full_body.static.T_pelvis = T_pelvis;
    struct_data(t).full_body.static.leg_length = hipC_R;
    
    nM = length(pelvis_marker_names);
    frsM = trialStruct(t).keyFrames.cropFrsForce;
    T_pelvis_mov = nan(4,4,frsM(2));rms_pelvis = nan(1,frsM(2));com_vel = nan(3,1,frsM(2));
    for fr = frsM(1):frsM(2)
        for m = 1:nM
            markers_anat(m,:) = md_stat.(pelvis_marker_names{m})';
            if ~isfield(trialStruct(t).marker_data,pelvis_marker_names{m})
                 markers_mov(m,:) = nan(3,1);
            else
                markers_mov(m,:) = trialStruct(t).marker_data.(pelvis_marker_names{m})(:,fr)';
            end
        end
        
        Inan = ~isnan(markers_mov(:,1));
        
        if sum(Inan) >= 3 % there are at least three markers
            
            [R,d,rms] = soder(markers_anat(Inan,:),markers_mov(Inan,:));
            
            T(:,:) = eye(4,4);
            T(1:3,1:3) = R;
            T(1:3,4) = d;
            T_pelvis_mov(:,:,fr) = T*T_pelvis;
            rms_pelvis(fr) = rms;
        end
    end
    struct_data(t).full_body.moving.T_pelvis = T_pelvis_mov;
    com = adaptiveLowPassButterworth(squeeze(T_pelvis_mov(1:3,4,frsM(1):frsM(2))),Fc,trialStruct(t).marker_data.FrameRate) ;
    figure; subplot(3,1,1);plot(squeeze(T_pelvis_mov(1:3,4,frsM(1):frsM(2)))')
    
    com_vel(:,frsM(1):frsM(2)) = calculateVelocity(com,trialStruct(t).marker_data.FrameRate);
    
   subplot(3,1,2);plot(com_vel(:,frsM(1):frsM(2))')
    legend('medial lateral','ant-post','vertical')
    
   subplot(3,1,3);plot(rms_pelvis(:,frsM(1):frsM(2))')
   
   vel_ant(t) = nanmean(com_vel(2,:));
   Fr(t) = vel_ant(t)^2/(9.81*hipC_R(3));
   leg_l(t) = hipC_R(3);
   
   md=trialStruct(t).marker_data;
   figure
%    hold on;
   plot(md.l_heel(3,:)); hold on
   plot(md.r_heel(3,:));
   plot([frsM(1) frsM(1)],[0 .5],'k')
   plot([frsM(2) frsM(2)],[0 .5],'k')
end
cmap = parula(5);
figure;
makeStatsDotPlot(repmat(1:2,1,nsubj)',vel_ant',{'walk','run'},cmap([1,3],:),'o')
title('anterior velocity')

makeNicePlotsFunction

figure;
[mdata,sdata] = makeStatsDotPlot(repmat(1:2,1,nsubj)',Fr',{'walk','run'},cmap([1,3],:),'o')
title('Froude number')

makeNicePlotsFunction
figure; 
subplot(2,1,1)
plot(leg_l(1:2:end),vel_ant(1:2:end),'o');
xlabel('leg length')
ylabel('anterior velocity')
title('walk')
subplot(2,1,2); plot(leg_l(2:2:end),vel_ant(2:2:end),'o');
xlabel('leg length')
ylabel('anterior velocity')
title('run')
makeNicePlotsFunction

figure; 
subplot(2,1,1)
plot(leg_l(1:2:end),Fr(1:2:end),'o');
xlabel('leg length')
ylabel('Froude #')
title('walk')
subplot(2,1,2); plot(leg_l(2:2:end),Fr(2:2:end),'o');
xlabel('leg length')
ylabel('Froude #')
title('run')

makeNicePlotsFunction
%% set up a colororder based on the subjects
subj_list_unique = subjStruct(:).subject;
figure;axes();hold on;
col_subj = jet(nsubj);%get(gca,'ColorOrder');
col_order = [6 3 7 1 5 2 4];% HA to low arch
col_subj(col_order,:) = col_subj;

marker_list.srun = 'v';
marker_list.nrun = 'x';
marker_list.frun = '^';
marker_list.walk = 'o';
marker_list.unknown = 'x';

marker_list.rfs = 'filled';
marker_list.ffs = 'none';

line_list.minimal = '-';
line_list.barefoot = ':';


for t = 1:ntrials
    
    trialStruct(t).plot.col = col_subj(findInStruct(subjStruct,'subject',trialStruct(t).subject),:); % get the subject's specific colour
    
    % initialise the conditions:
    trialStruct(t).condition.task = ' '; 
    trialStruct(t).condition.speed = ' ';
    trialStruct(t).condition.footwear = ' ';
    trialStruct(t).condition.strike = ' ';
    
    % determine if it's a walk or run
    
    if contains(trialStruct(t).trial,'walk')
        trialStruct(t).condition.task = 'walk';
        trialStruct(t).plot.marker_type = marker_list.walk;
        trialStruct(t).plot.marker_fill = 'none';
    elseif contains(trialStruct(t).trial,'run')
        
        trialStruct(t).condition.task = 'run';
        
        iR = strfind(trialStruct(t).trial,'run'); % index of the r
        
        lett = trialStruct(t).trial(iR-1); % the type of run, either nrun,frun, srun
        if ~contains('snf',lett)
            trialStruct(t).plot.marker_type = marker_list.unknown;
        else
            trialStruct(t).plot.marker_type = marker_list.([lett 'run']);
            
            trialStruct(t).condition.speed = [lett 'run'];
        end
        
        if contains(trialStruct(t).trial,'rfs') || contains(trialStruct(t).trial,'pref') 
            trialStruct(t).plot.marker_fill = trialStruct(t).plot.col;
            
            trialStruct(t).condition.strike = 'rfs';
        elseif contains(trialStruct(t).trial,'ffs')
            trialStruct(t).plot.marker_fill = marker_list.ffs;
            
        trialStruct(t).condition.strike = 'ffs';
        else
            
            trialStruct(t).plot.marker_fill = marker_list.ffs;
        end
        
    end
    
    % now check for minimal or barefoot
    
    if contains(trialStruct(t).trial,'barefoot')
        trialStruct(t).condition.footwear = 'barefoot';
        trialStruct(t).plot.line = line_list.barefoot;
        
    elseif contains(trialStruct(t).trial,'minimal')
        trialStruct(t).condition.footwear = 'minimal';
        trialStruct(t).plot.line = line_list.minimal;
    end
    
    c = trialStruct(t).condition;
    p = trialStruct(t).plot;
    h = plot([1 2],[t t],'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    text(2.5,t,[trialStruct(t).subject ' - ' c.task ' - ' c.speed  ' - ' c.footwear  ' - ' c.strike])
        
end

ylim([0 t+1]); xlim([0 5])



%% ANGLES : Find the frame to reference to for each trial: when the heel comes off the back plate, or the frame of max dorsiflexion of the tibia relative to calc

% also find the frame where the tibia is most upright

close all;
t_start = 1;
for t = t_start:1:ntrials
    nfrs_tot = size(trialStruct(t).Tm_mm.cal,3);
    fr_track = (diff(squeeze(trialStruct(t).Tm_mm.cal(1,1,[1:end,end])))~=0) & (~isnan(squeeze(trialStruct(t).Tm_mm.cal(1,1,:))));
    
    frs = findGroupsOfLogicals(fr_track')+1;
    nfrs_tr = diff(frs)+1;
    mid_pt = floor(nfrs_tr/2)+frs(1);
    
    Fs = trialStruct(t).marker_data.FrameRate;
      % Forces -----------------------
     
    frsX = trialStruct(t).keyFrames.cropFrsXray;
    frsF = trialStruct(t).keyFrames.cropFrsForce;
    nfrsX = abs(diff(frsX))+1;
    if contains(trialStruct(t).condition.task,'walk')
        frs0 = find(trialStruct(t).force_data(2).Force(3,frsF(1)+15:frsF(2)) < 50);
        figure(302); plot(trialStruct(t).force_data(2).Force(3,frsF(1)+15:frsF(2)))
        trialStruct(t).keyFrames.heelupForce = frsF(1)+15+frs0(1)-1;
        
    elseif contains(trialStruct(t).condition.task,'run')
        force_z{t} =  trialStruct(t).force_data(1).Force(3,frsF(1):frsF(2)) + trialStruct(t).force_data(2).Force(3,frsF(1):frsF(2));
        [~,Ifzmax] = max(force_z{t});
        trialStruct(t).keyFrames.heelupForce = frsF(1)+Ifzmax-1;
%         
%         figure(99)
%         hold on; plot(force_z{t})
%         plot(Ifzmax,force_z{t}(Ifzmax),'o')
    end
    
    
    
    Imax_hu = trialStruct(t).keyFrames.heelupForce - frsF(1)+1;
    
    trialStruct(t).keyFrames.heelupX = frsX(1)+Imax_hu-1;%Imax_hu -1;
     
     
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    ang_DP_tibtal = nan(nfrs_tot,1); align_z = nan(nfrs_tot,1);align_z_rot = nan(nfrs_tot,1);ang_DP_mtp = nan(nfrs_tot,1); ang_DP_calmet = nan(nfrs_tot,1);
    ang_DP_cal = nan(nfrs_tot,1);ang_DP_ph1 = nan(nfrs_tot,1); ang_DP_mt1 = nan(nfrs_tot,1); ang_DP_tibcal = nan(nfrs_tot,1);
    for fr = frsX(1):frsX(2)
        [ang_DP_tibtal(fr),ang_SP_tibtal(fr),ang_AD_tibtal(fr)] = eulerYZX(   trialStruct(t).Tm_mm.tib(:,:,fr),trialStruct(t).Tm_mm.tal(:,:,fr) ,...
            subjStruct(s_ind).bones.tib.T_ACS.T_TC,  subjStruct(s_ind).bones.tal.T_ACS.T_TC);
        
        
        [ang_DP_tibcal(fr),ang_SP_tibcal(fr),ang_AD_tibcal(fr)]  = eulerYZX(  trialStruct(t).Tm_mm.tib(:,:,fr),trialStruct(t).Tm_mm.cal(:,:,fr) ,...
            subjStruct(s_ind).bones.tib.T_ACS.T_TC,  subjStruct(s_ind).bones.cal.T_Aligned);
        
        [ang_DP_calmet(fr),ang_SP_calmet(fr),ang_AD_calmet(fr)]  = eulerYZX(  trialStruct(t).Tm_mm.cal(:,:,fr),trialStruct(t).Tm_mm.mt1(:,:,fr) ,...
            subjStruct(s_ind).bones.cal.T_Aligned,  subjStruct(s_ind).bones.mt1.T_Aligned);
        
        
        [ang_DP_mtp(fr),ang_SP_mtp(fr),ang_AD_mtp(fr)]  = eulerYZX(  trialStruct(t).Tm_mm.mt1(:,:,fr),trialStruct(t).Tm_mm.ph1(:,:,fr) ,...
            subjStruct(s_ind).bones.mt1.T_Aligned,  subjStruct(s_ind).bones.ph1.T_Aligned);
        
        
        [ang_DP_cal(fr),~,~]  = eulerYZX(  eye(4,4),trialStruct(t).Tm_mm.cal(:,:,fr) ,...
            eye(4,4),  eye(4,4));
        [ang_DP_ph1(fr),~,~]  = eulerYZX(  eye(4,4),trialStruct(t).Tm_mm.ph1(:,:,fr) ,...
            eye(4,4),  eye(4,4));
        [ang_DP_mt1(fr),~,~]  = eulerYZX(  eye(4,4),trialStruct(t).Tm_mm.mt1(:,:,fr) ,...
            eye(4,4),  eye(4,4));
        % find when the tibia z axis is most aligned with the global z
        T_acs_tib_x(:,:,fr) =  trialStruct(t).Tm_mm.tib(:,:,fr)* subjStruct(s_ind).bones.tib.T_ACS.T_TC;
        
        align_z(fr) = dot( T_acs_tib_x(1:3,3,fr),[0 0 1]');
        [align_z_rot(fr),~,~] = eulerYZX(eye(4,4),T_acs_tib_x(:,:,fr),eye(4,4),eye(4,4)); % the global axis rotation
        % get the angular velocity of the bones but put it in their
        % anatomical cosys
        Tw.tibACS(:,:,fr) = invTranspose(subjStruct(s_ind).bones.tib.T_ACS.T_TC)*trialStruct(t).Tm_mm.tib(:,:,fr);
        Tw.calACS(:,:,fr) = invTranspose(subjStruct(s_ind).bones.cal.T_Aligned)*trialStruct(t).Tm_mm.cal(:,:,fr);
        Tw.talACS(:,:,fr) = invTranspose(subjStruct(s_ind).bones.tal.T_Aligned)*trialStruct(t).Tm_mm.tal(:,:,fr);
        Tw.mt1ACS(:,:,fr) = invTranspose(subjStruct(s_ind).bones.mt1.T_Aligned)*trialStruct(t).Tm_mm.mt1(:,:,fr);
        Tw.ph1ACS(:,:,fr) = invTranspose(subjStruct(s_ind).bones.ph1.T_Aligned)*trialStruct(t).Tm_mm.ph1(:,:,fr);
        
        Tw.mt1cal(:,:,fr) = invTranspose(trialStruct(t).Tm_mm.cal(:,:,fr)) *trialStruct(t).Tm_mm.mt1(:,:,fr) ;
        Tw.taltib(:,:,fr) = invTranspose(trialStruct(t).Tm_mm.tib(:,:,fr)) *trialStruct(t).Tm_mm.tal(:,:,fr) ;
    end
    
    w(t).tib(:,frsX(1):frsX(2)) = calculateRotMatAngularVelocity( Tw.tibACS(1:3,1:3,frsX(1):frsX(2)),Fs,'deg');
    w(t).cal(:,frsX(1):frsX(2)) = calculateRotMatAngularVelocity( Tw.calACS(1:3,1:3,frsX(1):frsX(2)),Fs,'deg');
    w(t).tal(:,frsX(1):frsX(2)) = calculateRotMatAngularVelocity( Tw.talACS(1:3,1:3,frsX(1):frsX(2)),Fs,'deg');
    w(t).mt1(:,frsX(1):frsX(2)) = calculateRotMatAngularVelocity( Tw.mt1ACS(1:3,1:3,frsX(1):frsX(2)),Fs,'deg');
    w(t).ph1(:,frsX(1):frsX(2)) = calculateRotMatAngularVelocity( Tw.ph1ACS(1:3,1:3,frsX(1):frsX(2)),Fs,'deg');
     
    
    w(t).mt1cal(:,frsX(1):frsX(2)) = calculateRotMatAngularVelocity( Tw.mt1cal(1:3,1:3,frsX(1):frsX(2)),Fs,'deg');
    w(t).taltib(:,frsX(1):frsX(2)) = calculateRotMatAngularVelocity( Tw.taltib(1:3,1:3,frsX(1):frsX(2)),Fs,'deg');
    
      trialStruct(t).kinemat.ang_DP_calmet = ang_DP_calmet;
      trialStruct(t).kinemat.ang_DP_tibtal = ang_DP_tibtal;
    
    % Event detection
    [v_maxmtp(t),fr_mtp] = nanmax(ang_DP_mtp(mid_pt:end));
    trialStruct(t).keyFrames.max_mtp = fr_mtp + mid_pt-1;
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    
    [~,trialStruct(t).keyFrames.max_upright] = max(align_z);
    
    if rem(t,2) == 0 % run
        [v_maxad(t),trialStruct(t).keyFrames.max_dors] = max(ang_DP_tibtal(mid_pt:end));
        
        trialStruct(t).keyFrames.max_dors = trialStruct(t).keyFrames.max_dors + mid_pt;
    else % walk
        
        [v_maxad(t),trialStruct(t).keyFrames.max_dors] = max(ang_DP_tibtal(mid_pt+10:end));
        
        trialStruct(t).keyFrames.max_dors = trialStruct(t).keyFrames.max_dors + mid_pt + 10;
      
    end
    
    
    [v_minmtp(t),trialStruct(t).keyFrames.min_mtp] = nanmin(ang_DP_mtp(mid_pt:end));
    trialStruct(t).keyFrames.min_mtp = trialStruct(t).keyFrames.min_mtp + mid_pt;
    [v_maxarch(t),trialStruct(t).keyFrames.max_arch] = nanmax(ang_DP_calmet(1:end));
    trialStruct(t).keyFrames.max_arch = trialStruct(t).keyFrames.max_arch;
    
    F_AP = trialStruct(t).force_data(1).globForce(2,frsF(1):frsF(2))+trialStruct(t).force_data(2).globForce(2,frsF(1):frsF(2)); % sum of both FP AP forces
    F_SI = trialStruct(t).force_data(1).globForce(3,frsF(1):frsF(2))+trialStruct(t).force_data(2).globForce(3,frsF(1):frsF(2)); % sum of both FP SI forces
    
    Iprop = find(F_AP(10:end) > 0,1,'first') + 9;
    [maxSIForce(t),ImaxGRF] = max(F_SI(25:end));

    ImaxGRF = ImaxGRF + 24;
   
    Imax_mtp = trialStruct(t).keyFrames.max_mtp - frsX(1) + 1 ;
    Imin_mtp = trialStruct(t).keyFrames.min_mtp - frsX(1) + 1 ;
    Imax_adors = trialStruct(t).keyFrames.max_dors - frsX(1) + 1;
    Imax_arch = trialStruct(t).keyFrames.max_arch - frsX(1) + 1;
    Imax_start = trialStruct(t).keyFrames.heelupX-frsX(1) + 1;%trialStruct(t).keyFrames.max_dors - frsX(1);

     force_at_recoil(t) = F_SI(Imax_arch);
     trialStruct(t).keyFrames.Fprop = Iprop + frsF(1)-1;
     trialStruct(t).keyFrames.prop = Iprop + frsX(1)-1;
    
 
     % CALCULATE POWER
      
   
    
    bone = 'tal';
    T_UD = trialStruct(t).Tm_m.(bone)(:,:,frsX(1):frsX(2));
    cm_UD = subjStruct(s_ind).bones.(bone).centroid'/1000;
    F_GRF = trialStruct(t).force_data(1).globForce(:,frsF(1):frsF(2))/demo.weight(s_ind);
    M_free = trialStruct(t).force_data(1).globFreeMoment(:,frsF(1):frsF(2))/demo.weight(s_ind);
    cop = trialStruct(t).force_data(1).globCOP(:,frsF(1):frsF(2));
   
    
    
     [P_UD,w_UD,v_cm_UD] = UDDeformablePower(T_UD,cm_UD,F_GRF,M_free,cop,Fs );
    [Pmax,I_pow] = max(P_UD);
    E_UD = nancumtrapz(0:1/Fs:(nfrsX-1)/Fs,P_UD);
    data_mat_P(t,:) = [Pmax E_UD(I_pow)];
    
    p = trialStruct(t).plot;
    x_vals = linspace(0,100,nfrsX);
    figure(5); hold on; plot(x_vals,P_UD,'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill); ylabel('Power [W/kg]')

%     calculate impulse
impulse = [];
    for d = 1:3
    impulse(d,:) = nancumtrapz(0:1/Fs:(nfrsX-1)/Fs,F_GRF(d,:));
    end
    impulse_norm = norm3d(impulse);
  
    %     figure(6); hold on; plot(x_vals,E_UD/demo.weight(s_ind),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill); ylabel('Work [J/kg]')
%     figure(7); hold on; plot(E_UD(I_pow),Pmax,'Color',p.col,'marker',p.marker_type)
%     figure(8); hold on; plot(ROM_ankle(t),E_UD(I_pow),'Color',p.col,'marker',p.marker_type)
%     figure(9); hold on; plot(ROM_ankle(t),Pmax,'Color',p.col,'marker',p.marker_type)
%     figure(10); hold on; plot(v_maxarch(t),E_UD(I_pow),'Color',p.col,'marker',p.marker_type)
%     figure(11); hold on; plot(v_maxarch(t),Pmax,'Color',p.col,'marker',p.marker_type)
%      
%      
     
     
     % Calculate the timing - -------------------------------------
     
     cont_time(t) = diff(frsX)/Fs;%(fr_mtp_s'-fr_start_s');
     prop_time_Force(t) =  (Imax_mtp - Iprop + 1)/Fs; 
     arch_rec_time(t) = (Imax_mtp-Imax_arch+1)/Fs;
%      fr_diff_beg = (fr_tib_glob_align'-Imax_arch);
%      early_t = cont_time-fr_diff_beg;
%      perc = early_t./cont_time(t)*100;
%      
     
     
     p = trialStruct(t).plot;
    
    stancePts = linspace(0,100,frsX(2)-frsX(1)+1);

    I_labels = {'Maximum MTP dorsiflexion','Maximum ankle dorsiflexion','Maximum arch dorsiflexion','Peak vertical force','Minimum MTP dorsiflexion','Propulsion','max Talus power'};
    I_frs_save(t,:) = stancePts([Imax_mtp,Imax_adors,Imax_arch,ImaxGRF,Imin_mtp, Iprop,I_pow]);
    
    
    
    % measure tibial lean -------------------------------------
    align_z_s(t) = align_z(trialStruct(t).keyFrames.max_mtp);
    align_z_R_s(1,t) = align_z_rot(trialStruct(t).keyFrames.prop); % alignment at propulsion
    align_z_R_s(2,t) = align_z_rot(trialStruct(t).keyFrames.max_arch); % alignment at max archdors
    align_z_R_s(3,t) = align_z_rot(trialStruct(t).keyFrames.max_dors); % alignment at max dors
    align_z_R_s(4,t) = align_z_rot(trialStruct(t).keyFrames.max_mtp); % alignment at Max mtp
    
    arch_ang_s(1,t) = ang_DP_calmet(trialStruct(t).keyFrames.prop); % alignment at propulsion
    arch_ang_s(2,t) = ang_DP_calmet(trialStruct(t).keyFrames.max_arch); % alignment at max archdors
    arch_ang_s(3,t) = ang_DP_calmet(trialStruct(t).keyFrames.max_dors); % alignment at max dors
    arch_ang_s(4,t) = ang_DP_calmet(trialStruct(t).keyFrames.max_mtp); % alignment at Max mtp
    
    ankle_ang_s(1,t) = ang_DP_tibtal(trialStruct(t).keyFrames.prop); % alignment at propulsion
    ankle_ang_s(2,t) = ang_DP_tibtal(trialStruct(t).keyFrames.max_arch); % alignment at max archdors
    ankle_ang_s(3,t) = ang_DP_tibtal(trialStruct(t).keyFrames.max_dors); % alignment at max dors
    ankle_ang_s(4,t) = ang_DP_tibtal(trialStruct(t).keyFrames.max_mtp); % alignment at Max mtp

    figure(145); hold on;
%     subplot(2,1,1); hold on; grid on;
%     plot(stancePts,align_z(frsX(1):frsX(2))','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     ylabel('dot product of long tib axis with glob')
%     
%     subplot(2,1,2)
    hold on; grid on;
    plot(stancePts,-align_z_rot(frsX(1):frsX(2))','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    ylabel('rotational tibia alignment')
    makeNicePlotsFunction
    trialStruct(t).kinemat.tib_lean.dot_mov = align_z;
    trialStruct(t).kinemat.tib_lean.rot_mov = align_z_rot;
    
    figure(146); hold on

    plot([1 2 3 4], align_z_R_s(:,t)','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    xlim([0 4])
    set(gca,'XTick',[1 2 3],'XTickLabel',{'propulsion','max ankle dors','max MTP'})
    makeNicePlotsFunction
    
    
    figure(196+t_start)
    subplot(3,1,1); hold on
    plot(stancePts,norm3d(w(t).mt1cal(:,frsX(1):frsX(2))'),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    ylabel('mt1 dors ang vel [deg]')
    makeNicePlotsFunction
    
    subplot(3,1,2); hold on
    plot(stancePts,norm3d(w(t).taltib(:,frsX(1):frsX(2))'),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    ylabel('tib dors ang vel [deg]')
    makeNicePlotsFunction    
    
    subplot(3,1,3); hold on
    plot(stancePts,w(t).tal(1,frsX(1):frsX(2))','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    ylabel('tal dors ang vel [deg]')
    makeNicePlotsFunction
    
%     subplot(4,1,4); hold on
%     plot(stancePts([Imax_mtp,Imax_adors,Imax_arch,Imin_mtp,ImaxGRF,Iprop]),1:6,'Color',p.col,'marker',p.marker_type)  
%     set(gca,'YTick',1:5,'YTickLabel',{'Max MTP','Max ank dors','Max arch dors','Min MTP','Force SI','Force AP'})
%     xlim([0 100])
%     ylim([0 7])
%     
    
    
    figure(198+t_start)
    subplot(3,1,3); hold on
    plot(stancePts,ang_DP_calmet(frsX(1):frsX(2))','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    
    plot(stancePts(Imax_mtp),ang_DP_calmet(trialStruct(t).keyFrames.max_mtp),'Color','k','marker',p.marker_type,'MarkerFaceColor','k');
    ylabel('arch dorsiflexion [deg]')
    makeNicePlotsFunction

    subplot(3,1,2); hold on
    plot(stancePts,ang_DP_tibtal(frsX(1):frsX(2))','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    plot(stancePts(Imax_mtp),ang_DP_tibtal(trialStruct(t).keyFrames.max_mtp),'Color','k','marker',p.marker_type,'MarkerFaceColor','k');
    ylabel('tib-tal dors [deg]')
    makeNicePlotsFunction    
    
    subplot(3,1,1); hold on
    plot(stancePts,ang_DP_mtp(frsX(1):frsX(2))','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
    plot(stancePts(Imax_mtp),ang_DP_mtp(trialStruct(t).keyFrames.max_mtp),'Color','k','marker',p.marker_type,'MarkerFaceColor','k');
    ylabel('MTP dors [deg]')
    makeNicePlotsFunction
    
    figure(201);hold on;
     plot(ang_DP_calmet( (frsX(1) + Imax_arch - 1):fr_mtp)',ang_DP_tibtal((frsX(1)+Imax_arch-1):fr_mtp)','Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
   xlabel('arch angle')
   ylabel('tib-tal angle')
   makeNicePlotsFunction
   
   
    trialStruct(t).kinemat.mtp = ang_DP_mtp;
%     trialStruct(t).kinemat.mtp_ROM = diff(ang_DP_mtp([trialStruct(t).keyFrames.max_dors,trialStruct(t).keyFrames.max_mtp]));
    
    ROM_vals = ang_DP_calmet([trialStruct(t).keyFrames.Fprop-frsF(1)+1:trialStruct(t).keyFrames.max_mtp]);
    ROM_arch(t) = diff([max(ROM_vals),min(ROM_vals)]);
    ROM_vals = ang_DP_mtp([trialStruct(t).keyFrames.Fprop-frsF(1)+1:trialStruct(t).keyFrames.max_mtp]);
    ROM_mtp(t) = diff([max(ROM_vals),min(ROM_vals)]);
    ROM_vals = ang_DP_tibtal([trialStruct(t).keyFrames.Fprop-frsF(1)+1:trialStruct(t).keyFrames.max_mtp]);
    ROM_ankle(t) = diff([max(ROM_vals),min(ROM_vals)]);
    
    trialStruct(t).kinemat.ROM_mtp = ROM_mtp(t);
    trialStruct(t).kinemat.ROM_arch = ROM_arch(t);
    trialStruct(t).kinemat.ROM_ankle = ROM_ankle(t);
    
    if isnan(trialStruct(t).Tm_mm.tib(:,:,trialStruct(t).keyFrames.max_mtp))
        warning(sprintf('Not enough tibia frames tracked %s %s',trialStruct(t).subject,trialStruct(t).trial));
    end
    
    
% impulse vs ROM    
       figure(6); hold on; plot(ROM_arch(t),impulse_norm(end),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
       ylabel('Impusle [Ns/kg]')
    xlabel('Arch ROM')
%   ROM PLOTS AGAINST EACH OTHER------------------------------------------     
%     figure(208)
%     hold on;
%     plot(ROM_mtp(t),ROM_arch(t),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     xlabel('ROM MTP')
%     ylabel('ROM arch')
%     %     bar(s,ROM_arch(s),'Facecolor', plot_cols{s,1})
%     title('ROM of arch to MTP ROM between heel up and max mtp')
%     
%     figure(209)
%     hold on;
%     plot(ROM_ankle(t),ROM_arch(t),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     xlabel('ROM ankle')
%     ylabel('ROM arch')
%     title('ROM of arch to ROM of ankle between heel up and max mtp')

    
    fr_tracked = find(diff(squeeze(trialStruct(t).Tm_mm.tal(1,4,:)))~=0);
    last_fr_tracked = fr_tracked(end);
    trialStruct(t).keyFrames.last_fr_tracked = last_fr_tracked;
    
  
end


makeGroupedStatsDotPlot(align_z_R_s',repmat(1:2,1,ntrials/2)',{'walk','run'},{'propulsion','max arch dors','max ankle dors','max MTP'},parula(5),'d')
ylabel('tibia alignment with global axes [^o]')
makeNicePlotsFunction

makeGroupedStatsDotPlot(arch_ang_s'-arch_ang_s(2,:)',repmat(1:2,1,ntrials/2)',{'walk','run'},{'propulsion','max arch dors','max ankle dors','max MTP'},parula(5),'d')
ylabel('arch angle [^o]')
makeNicePlotsFunction

makeGroupedStatsDotPlot(ankle_ang_s'-ankle_ang_s(3,:)',repmat(1:2,1,ntrials/2)',{'walk','run'},{'propulsion','max arch dors','max ankle dors','max MTP'},parula(5),'d')
ylabel('ankle angle [^o]')
makeNicePlotsFunction



% Timing --------------------------------



res.timing.cont_time = cont_time;
res.timing.prop_time_force = prop_time_Force;
res.timing.arch_rec_time = arch_rec_time;


res.tib_lean.rot_mov = align_z_R_s;
res.tib_lean.dot_mov = align_z_s;


cmap = turbo(14);
% showColorMap(cmap)

     
     figure(120); hold on;
     h(1) = plot(ones(nsubj,1)-0.1,cont_time(1:2:end),'o','color',cmap(14,:));
     plot(ones(nsubj,1)*2-0.1,cont_time(2:2:end),'x','color',cmap(14,:));
    
     h(2) =plot(ones(nsubj,1),prop_time_Force(1:2:end),'o','color',cmap(5,:));
     plot(ones(nsubj,1)*2,prop_time_Force(2:2:end),'x','color',cmap(5,:));
     
     h(3) =plot(ones(nsubj,1)+0.1,arch_rec_time(1:2:end),'o','color',cmap(2,:));
     plot(ones(nsubj,1)*2+0.1,arch_rec_time(2:2:end),'x','color',cmap(2,:));
     xlim([0 3])
     ylabel('Time [s]')
     set(gca,'XTick',[1 2],'XTickLabel',{'walk','run'})
     legend(h,{'contact time','propulsion time (force)','arch dors peak to mtp'})


% ALL PEAK COMPARISONS OF ANGLES
% 
% 
% figure; 
% hold on;
% for t = 1:2:ntrials
%     
%      p = trialStruct(t).plot;
%     
%  plot([ROM_arch(t);ROM_arch(t+1)],'o','color',p.col,'Linestyle','-')
% 
% end
% xlim([0 3])
% title('arch ROM')
% set(gca,'XTick',1:2,'XTickLabel',{'walk','run'})
% makeNicePlotsFunction
% 
% 
% figure; 
% hold on;
% for t = 1:2:ntrials
%     
%      p = trialStruct(t).plot;
%     
%  plot([ROM_ankle(t);ROM_ankle(t+1)],'o','color',p.col,'Linestyle','-')
% 
% end
% xlim([0 3])
% title('ankle ROM')
% set(gca,'XTick',1:2,'XTickLabel',{'walk','run'})
% makeNicePlotsFunction
% 
% 
% figure; 
% hold on;
% for t = 1:2:ntrials
%     
%      p = trialStruct(t).plot;
%     
%  plot([v_maxarch(t);v_maxarch(t+1)],'o','color',p.col,'Linestyle','-')
% 
% end
% xlim([0 3])
% title('arch peak')
% set(gca,'XTick',1:2,'XTickLabel',{'walk','run'})
% makeNicePlotsFunction
% 
% figure; 
% hold on;
% for t = 1:2:ntrials
%     
%      p = trialStruct(t).plot;
%     
%  plot([v_maxad(t);v_maxad(t+1)],'o','color',p.col,'Linestyle','-')
% %  plot([v_maxarch(1:2:ntrials);v_maxarch(2:2:ntrials)],'k')
% % [mar,sar] =makeStatsDotPlot(1:2,abs([v_maxarch(1:2:ntrials);v_maxarch(2:2:ntrials)]),{'walk','run'},[0 0.25 0.5;0 0.25 0.5],'o')
% 
% end
% xlim([0 3])
% title('ankle peak')
% set(gca,'XTick',1:2,'XTickLabel',{'walk','run'})
% makeNicePlotsFunction

%% Abstract plot showing timing
figure;
group_names = I_labels;
trial_names = {'walk','run'};
cmap = jet(20);
% showColorMap(cmap)
colW = cmap(15,:);
colR = cmap(1,:);
for group_num = 1:7 % each index that's saved
    
    meanW = mean(I_frs_save(1:2:ntrials,group_num));
    stdW = std(I_frs_save(1:2:ntrials,group_num));
    h(1) = scatter(I_frs_save(1:2:ntrials,group_num),repmat(group_num-0.1,1,nsubj),'markerfacecolor',colW,'markeredgecolor',colW,'marker','o');hold on;
    errorbar(meanW,group_num-0.1,stdW,'horizontal','marker','square','capsize',20,'markersize',10,'color',colW)
    
    meanR = mean(I_frs_save(2:2:ntrials,group_num));
    stdR = std(I_frs_save(2:2:ntrials,group_num));
    h(2) = scatter(I_frs_save(2:2:ntrials,group_num),repmat(group_num+0.1,1,nsubj),'markerfacecolor',colR,'markeredgecolor',colR,'marker','o');
    errorbar(meanR,group_num+0.1,stdR,'horizontal','marker','square','capsize',20,'markersize',10,'color',colR)
    
    fprintf('%s W: %0.0f %s %0.0f, R: %0.0f %s %0.0f \n',I_labels{group_num},meanW,0177,stdW/sqrt(nsubj),meanR,0177,stdR/sqrt(nsubj))
end

set(gca,'YTick',1:length(I_labels),'YTickLabel',I_labels,'XGrid','on','YGrid','on')
% grid on
xlim([30 100])
ylim([0 length(I_labels)+1])
legend(h(1:2),{'walk','run'})
xlabel('% stance')
makeNicePlotsFunction

%% Load all the reference meshes for CPD


addpath(genpath('C:\Users\welte\Documents\Code\coherent_point_drift_cuda-master\'))


dir_ref_bn{1} = fullfile(dir_ref,'Calcaneus',filesep);
dir_ref_bn{3} = fullfile(dir_ref,'Talus',filesep);
dir_ref_bn{4} = fullfile(dir_ref,'Tibia',filesep);

for bn = [1,3,4]%nBones
    bone_file = ls(fullfile(dir_ref_bn{bn},['*' bone_list{bn} '*.iv']));
    
    load(fullfile(dir_ref_bn{bn},'bonestruct.mat'))
    T_ACS = bonestruct.(bone_list{bn}).T_Aligned;
    
    [ref_pts, cns_temp] = read_vrml_fast(fullfile(dir_ref_bn{bn},bone_file));
    % load the specifically created ref files
    
    cpd_ref.(bone_list{bn}).cns = cns_temp(:,1:3)+1;
    cpd_ref.(bone_list{bn}).pts = transformPoints(T_ACS,ref_pts,-1);
    cpd_ref.(bone_list{bn}).T_ACS = T_ACS;
end
%% Orient all the bones in the anatomical frame and run CPD

col_map = colormap('summer');
ind = round(linspace(1,64,nsubj));
col_map = col_map(ind,:);

bone_cpd_list = fields(cpd_ref);
for bn = [1,3:4] % 1:nBones
    ivstring = createInventorHeader();
    
    for s  = 1:nsubj%6%[1:ref_ind-1,ref_ind+1:nsubj]
        T_ACS = subjStruct(s).bones.(bone_list{bn}).T_Aligned; % gives the aligned inertial axes
        bone_pts = subjStruct(s).bones.(bone_list{bn}).pts;
        pts_anat_pre = transformPoints(T_ACS,bone_pts,-1);
       
        
        
        maxIter = 100;
        
% MATLAB implementation to do the rigid alignment first
        opt.method = 'rigid'; opt.corresp = 0; opt.max_it = maxIter; opt.viz = 0;
        
        Tra = cpd_register(cpd_ref.(bone_list{bn}).pts,pts_anat_pre,opt);
        
         pts_anat_rot_sc = Tra.Y;
          Trot = eye(4,4);
            Trot(1:3,:) = [Tra.R Tra.t]; % this is the rotation that transforms pts_anat_pre to pts_anat_rot
            
            
              omega = 0.1; beta = 2; lambda = 3; maxIter = 3000; tol = 1e-5;
                 %         Parameter lambda represents the trade off between data fitting and smoothness regularization.
        %       Parameter beta reflects the strength of interaction between points. Small values of produce locally smooth transformation, while large
        % values of beta correspond to nearly pure translation transformation
            [pts_anat_def_sc,C] = cpd_cuda(pts_anat_rot_sc,cpd_ref.(bone_list{bn}).pts,omega, beta, lambda, maxIter, tol);
            
            pts_anat_def = 1/Tra.s * pts_anat_def_sc; % anatomically located deformed points
            pts_anat_corr = 1/Tra.s * pts_anat_rot_sc(C,:); % anatomically located, corresponding points (from original mesh)
            
            pts_CT_corr = bone_pts(C,:);
                pts_temp = transformPoints(Trot,pts_anat_def,-1); % add the inital rotation back
            pts_CT_def = transformPoints(T_ACS,pts_temp); % now return it to CT space
            
        save(fullfile(dir_analy,'CPD',[ subjStruct(s).subject '_' bone_list{bn} '_CPD.mat']),'pts_anat_corr','pts_anat_def','pts_CT_def','pts_CT_corr','C');
    
%         save(fullfile(dir_analy,'CPD',[ subjStruct(s).subject '_' bone_list{bn} '_CPD.mat']),'pts_corr','pts_anat','pts_CT','C');

        
        % %
        patch2iv( pts_CT_corr, cpd_ref.(bone_list{bn}).cns ,fullfile(dir_analy,'CPD','Viz',[ subjStruct(s).subject '_' bone_list{bn} '_CPD.iv']),[0.3 0 0.8])
        patch2iv( pts_anat_pre, subjStruct(s).bones.(bone_list{bn}).cnt(:,1:3) ,fullfile(dir_analy,'CPD','Viz',[ subjStruct(s).subject '_' bone_list{bn} '_Anat.iv']))
        patch2iv( pts_anat_def, cpd_ref.(bone_list{bn}).cns ,fullfile(dir_analy,'CPD','Viz',[ subjStruct(s).subject '_' bone_list{bn} '_CPDAnatDef.iv']))
        patch2iv( pts_anat_corr, cpd_ref.(bone_list{bn}).cns ,fullfile(dir_analy,'CPD','Viz',[ subjStruct(s).subject '_' bone_list{bn} '_CPDAnatCorr.iv']))
        
        ivstring = [ivstring createInventorLink(fullfile(dir_analy,'CPD','Viz',[ subjStruct(s).subject '_' bone_list{bn} '_Anat.iv']),eye(3,3),[0 0 0],col_map(s,:),0.5)];
        
        ivstring = [ivstring createInventorLink(fullfile(dir_analy,'CPD','Viz',[ subjStruct(s).subject '_' bone_list{bn} '_Anat.iv']),eye(3,3),[s*50,0,100],col_map(s,:),0)];
        ivstring = [ivstring createInventorLink(fullfile(dir_analy,'CPD','Viz',[ subjStruct(s).subject '_' bone_list{bn} '_CPDAnatDef.iv']),eye(3,3),[s*50,0,200],col_map(s,:),0)];
    end
    fid = fopen(fullfile(dir_analy,'CPD','Viz',['ALL_ANAT_' bone_list{bn} '.iv']),'w');
    fprintf(fid,ivstring);
    fclose(fid);
end

%% segment the talar surface and tibia surface


[tal_dome.pts,tal_dome.cns] = read_vrml_fast(fullfile(dir_ref,'Talus','ta_dome_large.iv'));
tal_dome.cns = tal_dome.cns(:,1:3) + 1;
tal_dome.pts = transformPoints(cpd_ref.tal.T_ACS,tal_dome.pts,-1);
[tib_dome.pts,tib_dome.cns] = read_vrml_fast(fullfile(dir_ref,'Tibia','talocrural_surf.iv'));
tib_dome.cns = tib_dome.cns(:,1:3) + 1;

tib_dome.pts = transformPoints(cpd_ref.tib.T_ACS,tib_dome.pts,-1);
% calculateSegIVindices(ref_iv,crop_iv)

[~,~,tal_dome.iRef] = intersect(round(tal_dome.pts,2),round(cpd_ref.tal.pts,2),'rows');
[~,~,tib_dome.iRef] = intersect(round(tib_dome.pts,2),round(cpd_ref.tib.pts,2),'rows');

% [tal_dome.iRef2,Itd] = sort(tal_dome.iRef);
% [tib_dome.iRef2,Itbd] = sort(tib_dome.iRef);
% [new_tal_pts,tal_dome.cns] = rewriteCns(cpd_ref.tal.pts,cpd_ref.tal.cns,tal_dome.iRef2);
% [~,tib_dome.cns] = rewriteCns(cpd_ref.tib.pts,cpd_ref.tib.cns,tib_dome.iRef2);

[~,tal_dome.cns] = rewriteCns(cpd_ref.tal.pts,cpd_ref.tal.cns,tal_dome.iRef);
[~,tib_dome.cns] = rewriteCns(cpd_ref.tib.pts,cpd_ref.tib.cns,tib_dome.iRef);


bn = 3;% tal
for s = 1:nsubj
    
    % load the raw points
    pts_raw = subjStruct(s).bones.(bone_list{bn}).pts;
    % load the CPD
    load(fullfile(dir_analy,'CPD',[ subjStruct(s).subject '_' bone_list{bn} '_CPD.mat']))
    
    % get the CPD points
    
    [pts_crop,~] = rewriteCns(pts_CT_def,cpd_ref.tal.cns,tal_dome.iRef);
   
    % save the points of the dome
    subjStruct(s).tal_dome.pts = pts_crop;
    subjStruct(s).tal_dome.cns = tal_dome.cns;
    
    if s ==1
        figure;
        patch('faces',tal_dome.cns,'vertices',subjStruct(s).tal_dome.pts,'facealpha',0.5); hold on
        patch('faces',subjStruct(s).bones.(bone_list{bn}).cnt(:,1:3),'vertices',pts_raw,'facealpha',0.1,'facecolor','m');
    end
end

bn = 4; % tib
for s = 1:nsubj
    
    % load the raw points
    pts_raw = subjStruct(s).bones.(bone_list{bn}).pts;
    % load the CPD
    load(fullfile(dir_analy,'CPD',[ subjStruct(s).subject '_' bone_list{bn} '_CPD.mat']))
   
    % get the CPD points
    
    [pts_crop,~] = rewriteCns(pts_CT_def,cpd_ref.tib.cns,tib_dome.iRef);
   
    % save the points
    subjStruct(s).tib_dome.pts = pts_crop;
    subjStruct(s).tib_dome.cns = tib_dome.cns;
    
    
    if s ==1
        figure;
        patch('faces',tib_dome.cns,'vertices',subjStruct(s).tib_dome.pts,'facealpha',0.5); hold on
    %     patch('faces',trialStruct(s).bones.(bone_list{bn}).cnt(:,1:3),'vertices',pts_raw,'facealpha',0.1,'facecolor','m');
    %     figure;
    %     pcshow(pointCloud(trialStruct(s).tib_dome.pts))
    %     pcshow(pointCloud(trialStruct(s).tal_dome.pts))
    end
end
% see how far the tibia is away in CT space
figure
for s = 1:nsubj
    subplot(3,3,s)
    D = lookUpDfieldPts(subjStruct(s).bones.tal.dfield,subjStruct(s).tib_dome.pts,eye(3,3),[0 0 0]);
    
    hold on;
    hist(D)
    title(subjStruct(s).subject)
    surface_dist(s) = prctile(D,95);
    plot(repmat(surface_dist(s),2,1),[0 40])
end


%% measure the fixed arch transforms
close all

for t = 1:ntrials-2

    Tm = trialStruct(t).Tm_mm;
    nfr_tot = size(Tm.cal,3);
    
    frsX = trialStruct(t).keyFrames.cropFrsXray;
    frsF = trialStruct(t).keyFrames.cropFrsForce;
    nfrs_stance = abs(diff(frsX))+1;
    % important gait events
    fr_start = trialStruct(t).keyFrames.max_arch;%
%     fr_start_s(t) - fr_star;
    fr_mtp = trialStruct(t).keyFrames.max_mtp;%fr_mtp_s(t);
    % I frames depict that they are from the first frame of cropped data
    Ifr_hu = trialStruct(t).keyFrames.heelupForce - frsF(1)+1;
    Ifr_min_mtp = trialStruct(t).keyFrames.min_mtp - frsX(1) + 1 ;
    Ifr_mxad = trialStruct(t).keyFrames.max_dors - frsX(1) + 1;
    Ifr_mxarch = trialStruct(t).keyFrames.max_arch - frsX(1) + 1;


    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    
     % RIGID------------------------------------------------------------------------------
    
    T_mt1_fix = Tm.mt1(:,:,fr_start);
    T_mt1 = Tm.mt1;
    clearvars('T_fix')
    bone_listM = fields(Tm);
    nBones = length(bone_listM);
    for bn = 1:nBones
        T_bone_fix = Tm.(bone_listM{bn})(:,:,fr_start);
        if strcmp(bone_listM{bn},'tib') % give it tibio-talar kinematics
            for fr = frsX(1):frsX(2)
                % give all the fixed transforms with the mt1 relative to the selected max dors position
                T_fix.(bone_listM{bn})(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * Tm.tal(:,:,fr_start) *  invTranspose(Tm.tal(:,:,fr)) * Tm.tib(:,:,fr);
                % lock it in the plantarflexed position
                T_bone_fix = Tm.(bone_listM{bn})(:,:,fr_mtp);
                T_fix.tib_plant(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * Tm.tal(:,:,fr_start) *  invTranspose(Tm.tal(:,:,fr_mtp)) * Tm.tib(:,:,fr_mtp);%T_mt1(:,:,fr) * invTranspose(Tx.mt1(:,:,fr_mtp)) * T_bone_fix;
            end
        elseif  strcmp(bone_listM{bn},'ph1')% give it MTPJ kinematics
            
            for fr = frsX(1):frsX(2)
                % give all the fixed transforms with the mt1 relative to the selected max dors position
                T_fix.(bone_listM{bn})(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * Tm.mt1(:,:,fr_start) *  invTranspose(Tm.mt1(:,:,fr)) * Tm.ph1(:,:,fr);
                
                % THIS IS REDUNDANT BECAUSE IT COULD JUST BE Tx.ph1 -
                % because we're fixing it relative to the bone we're moving
                % everything with
            end
        else
            for fr = frsX(1):frsX(2)
                % give all the fixed transforms with the mt1 relative to the selected max dors position
                T_fix.(bone_listM{bn})(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * T_bone_fix;
            end
        end
    end
    
    trialStruct(t).T_fix    = T_fix;
    
    
    
    % -------------------- push off frame if tibia is in max plantarflexion----
    T_ACS_tib = subjStruct(s_ind).bones.tib.T_ACS.T_TC;
    tib_pts = subjStruct(s_ind).bones.tib.pts;
    % find the tibia's position
    align_z = [];
    T_act_plant = Tm.tib(:,:,fr_mtp) * T_ACS_tib;
    for fr = fr_start:fr_mtp
        
        T_fix_plant = T_fix.tib_plant(:,:,fr) * T_ACS_tib;
        
        
        %            figure(101)
        %             pcshow(pointCloud(transformPoints(T_fix.tib_plant(:,:,fr),tib_pts)));hold on;
        %             pcshow(pointCloud(transformPoints(Tx.tib(:,:,fr_mtp),tib_pts)));
        %             drawnow
        %             figure(100)
        %             hold on;
        %             plotvector3([0 0 0]',T_fix_plant(1:3,3),'r')
        %             plotvector3([0 0 0]',T_act_plant(1:3,3),'b')
        %             drawnow
        
        
        align_z(fr) = dot(T_fix_plant(1:3,3),T_act_plant(1:3,3));
    end
    
    [~,fr_tib_glob_align(t)] = max(align_z);
    
    trialStruct(t).keyFrames.fr_tg = fr_tib_glob_align(t); %
    
    
    
    % finally, take the tib @ max mtp in the fixed position, and rotate it
    % to align with the global pos of the tib @ max mtp, see where it is on
    % the talus
    T_ACS_tal = subjStruct(s_ind).bones.tal.T_ACS.T_TC;
    
    Tm_fr_mtp = Tm.tib(:,:,fr_mtp);
    T_fix_mtp =  T_fix.tib(:,:,fr_mtp);
     T_helical = invTranspose(T_fix_mtp) * Tm_fr_mtp ;
    hel_params = convertRotation(T_helical,'4x4xn','helical');
%     s_xr = transformPoints(T_fix_mtp,hel_params(6:8));
%     n_xr = transformVectors(T_fix_mtp,hel_params(2:4));
%     T_helical_tib = invTranspose(Tm.tib(:,:,fr_mtp))
%     [phi,n,L,s] = helicalInstantaneous(Tm.tal(:,:,[fr_start,fr_mtp]),Tm.tib(:,:,[fr_start,fr_mtp]));
    
    % use the axis from the full range of push-off
    trans = 0;
    q = hel_params(6:8);%[0 0 0];
%     n = T_ACS_tal(1:3,1)';%hel_params(2:4);
%     phi = hel_params(1);
%     
%     [R,T] = Helical_To_RT(phi, n, trans, q);
%     Tt = eye(4);
%     Tt(1:3,1:3) = R;
%     Tt(1:3,4) = T;
%     
%     
%     T_tib_fix_act_mtp =  T_fix.tib(:,:,fr_mtp)* T_ACS_tal *  Tt * invTranspose(T_ACS_tal);
%     
%     
%     
%     figure; visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(:,1:3),Tm_fr_mtp)
%     hold on; 
%     visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(:,1:3),T_fix_mtp,[0.5 0.5 0.5])
%     visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(:,1:3), T_tib_fix_act_mtp,[0.5 0 0])
%     
%      view([82 13])
%      title(trialStruct(t).trial)
%      axis off;
%     
    n = hel_params(2:4);
    phi = hel_params(1);
    
    [R,T] = Helical_To_RT(phi, n, trans, q);
    Tt = eye(4);
    Tt(1:3,1:3) = R;
    Tt(1:3,4) = T;
    
    
    T_tib_fix_act_mtp =  T_fix.tib(:,:,fr_mtp)*  Tt;
    figure; hold on
    visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(:,1:3), T_tib_fix_act_mtp,[0 0.5 0])
    
%     legend('real tib @ mtp','fixed tib','moved with TC axis','moved with helical axis')

%     visualizeBone(subjStruct(s_ind).bones.tal.pts,subjStruct(s_ind).bones.tal.cnt(:,1:3),T_fix.tal(:,:,fr_mtp),[0.5 0.5 0.5])
%------------for tal points in tib dfield---------------
%     tal_pts = subjStruct(s_ind).tal_dome.pts;
%     tal_cns = subjStruct(s_ind).tal_dome.cns(:,1:3);
%     tal_pts_x = transformPoints(T_fix.tal(:,:,fr_mtp),tal_pts); 
%     [centrs,norms] = calculateMeshCentroidsAndNormals(tal_pts_x,tal_cns);
%     
%      dfield = subjStruct(s_ind).bones.tib.dfield;
%     % this will provide the distance of all the triangle centroids from the
%     % tibia
%     dd = lookUpDfieldPts(dfield,centrs,T_tib_fix_act_mtp(1:3,1:3),T_tib_fix_act_mtp(1:3,4)');
%     
%     t_new = optimizeCongruence(centrs,dfield,T_tib_fix_act_mtp);
%     
%     ---------------for tib points in tal dfield
    
    tib_pts = subjStruct(s_ind).tib_dome.pts;
    tib_cns = subjStruct(s_ind).tib_dome.cns(:,1:3);
%     tib_pts_x = transformPoints( T_tib_fix_act_mtp,tib_pts); 
%     [centrs,norms] = calculateMeshCentroidsAndNormals(tib_pts_x,tib_cns);
    [centrs,norms] = calculateMeshCentroidsAndNormals(tib_pts,tib_cns);
    
     dfield = subjStruct(s_ind).bones.tal.dfield;
     
     
    visualizeBone(tib_pts,tib_cns, T_tib_fix_act_mtp,[0.5 0 0])
    hold on; visualizeBone(subjStruct(s_ind).bones.tal.pts,subjStruct(s_ind).bones.tal.cnt(:,1:3),T_fix.tal(:,:,fr_mtp))
    % this will provide the distance of all the triangle centroids from the
    % tibia
    centrs_x = transformPoints(T_tib_fix_act_mtp,centrs);
    dd = lookUpDfieldPts(dfield,centrs_x,T_fix.tal(1:3,1:3,fr_mtp),T_fix.tal(1:3,4,fr_mtp)');
    
    t_new = optimizeCongruence(centrs,dfield,T_tib_fix_act_mtp,T_fix.tal(:,:,fr_mtp));
    
     T_tib_fix_act_mtp(1:3,4) = t_new';
    trialStruct(t).kinemat.T_tib_fix_act_mtp = T_tib_fix_act_mtp;
%     
    visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(:,1:3),T_tib_fix_act_mtp)
    hold on; visualizeBone(subjStruct(s_ind).bones.tal.pts,subjStruct(s_ind).bones.tal.cnt(:,1:3),T_fix.tal(:,:,fr_mtp))
     view([76 0])
     title(trialStruct(t).trial)
     axis off;
     drawnow
%     figure; visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(1:10:end,1:3),invTranspose(Tm.tal(:,:,fr_mtp)) * Tm.tib(:,:,fr_mtp))
%     hold on; visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(1:10:end,1:3),invTranspose(Tm.tal(:,:,fr_start)) * Tm.tib(:,:,fr_start))
%     visualizeBone(subjStruct(s_ind).bones.tal.pts,subjStruct(s_ind).bones.tal.cnt(1:10:end,1:3),invTranspose(Tm.tal(:,:,fr_mtp)) * Tm.tal(:,:,fr_mtp))
%     h = plotvector3(s,n*100,'r');
%     h.LineWidth = 5;
end
%% Calculate the plane of the achilles on the calcaneus


cmap = get(gca,'ColorOrder');
[ach_surf.pts,ach_surf.cns] = read_vrml_fast(fullfile(dir_ref,'Calcaneus','Achilles','achilles_surf.iv'));
ach_surf.pts = transformPoints(cpd_ref.cal.T_ACS,ach_surf.pts,-1); % put the achilles in the same cosys
[~,~,iRef] = intersect(ach_surf.pts,cpd_ref.cal.pts,'rows');
% figure; plot3quick_scatter(ach_surf.pts');
% hold on; plot3quick_scatter(cpd_ref.cal.pts')
bn = 1; % cal


for t = 1:ntrials
    
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    
    % load the raw points
    pts_raw = subjStruct(s_ind).bones.(bone_list{bn}).pts;
    % load the CPD
    load(fullfile(dir_analy,'CPD',[ subjStruct(s_ind).subject '_' bone_list{bn} '_CPD.mat']))
    % get the CPD points
    pts_cpd_CT = pts_CT_corr;
    
    % segment the achilles surface
    pts_ach = pts_cpd_CT(iRef,:);
    
    pts_tib_CT = subjStruct(s_ind).bones.tib.pts;
    
    
    % save the points of the achilles
    subjStruct(s_ind).achilles.pts = pts_ach;
    
    
    
    % move the tibia and calc/achilles insertion to be at the position at
    % max upright
    fr_mu = trialStruct(t).keyFrames.max_upright	;
    
    T_align_tib = trialStruct(t).Tm_mm.tib(:,:,fr_mu);
    T_align_calc = trialStruct(t).Tm_mm.cal(:,:,fr_mu);
    
    % move the achilles points to max dorsiflexion
    pts_ach_md = transformPoints(T_align_calc,pts_ach);
    pts_tib_x = transformPoints(T_align_tib,pts_tib_CT);
    
    % fit a plane to the points
    [plane_normal,~,plane_point] = affine_fit(pts_ach_md);
    
    mean_ach_pt = mean(pts_ach_md);
    
    d =  - plane_point * plane_normal;
    pModel = planeModel([plane_normal;d]);
    
    point_plane = []; dist_save = [];
    % find the closest point on the tibia to the plane
    for p= 1:length(pts_tib_CT)
        point_plane(p,:) = closestPointonPlanealongVector(pts_tib_x(p,:),plane_normal',plane_point,plane_normal') ;
        dist_saveAch(p) = distBtwPoints3xN(point_plane(p,:), mean_ach_pt);%  pts_tib_x(p,:));
        dist_saveTib(p) = distBtwPoints3xN(point_plane(p,:),pts_tib_x(p,:));
    end
    
    %     mean(point_plane)
    [~,piA] = min(dist_saveAch);
    
    [~,piT] = min(dist_saveTib);
    
    figure; % verify we've got the achilles surface
    pcshow(pointCloud(transformPoints(T_align_calc,pts_raw)));hold on;
    h = pcshow(pointCloud(transformPoints(T_align_calc,pts_ach)),'markersize',50);% the achilles points
    pcshow(pointCloud(transformPoints(T_align_tib,pts_tib_CT)))% the tibia @ max dors
    %     h = pModel.plot; h.FaceAlpha = 0.9;
    h = plot3quick(point_plane(piA,:)',[1 1 1 ],'o','none');
    h = plot3quick(point_plane(piT,:),cmap(2,:),'o','none');
    h = plot3quick(mean(point_plane),[1 1 1],'o');
    h = plot3quick(mean_ach_pt,[1 1 1],'o');
    view([-90 0])
    
    
    trialStruct(t).achilles.cal_pt_CT = transformPoints(T_align_calc, mean_ach_pt,-1);
    trialStruct(t).achilles.prox_pt_CT  =  transformPoints(T_align_tib, point_plane(piA,:),-1);
end





%% Identify the contact of the tibia and the talus

% out = JointContact_Dfield(Dfield_ref,pts,conns,threshold,ptsconns,visualize,RT_ref,RT_pts)
% % OUT = JOINTCONTACT_DFIELD(DFIELD_REF,PTS,CONNS,THRESHOLD,PTSCONNS,VISUALIZE,RT_REF,RT_PTS)
%
%   inputs:
%       Dfield_ref - reference distance field
%       pts - test bone points
%       conns - test bone connections
%       threshold - distance threshold for contact patch (e.g. [surf1 surf2])
%       visualize - 1 = display results, 0 = do not display results
%       ptsconn - 1 = output pts and connections for the contact patch, 0 =
%           do not output pts and connections for the contact patch
%       RT_ref (optional) - 4x4 transformation matrix for the bone
%       associated with the distance field
%       RT_pts (optional) - 4x4 transformation matrix for the test points

for t = 1:ntrials-2
    
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    fprintf('Working on subject %s...',trialStruct(t).subject,trialStruct(t).trial)
    
    fr_start = trialStruct(t).keyFrames.Fprop - trialStruct(t).keyFrames.frOffset ;%trialStruct(t).keyFrames.max_arch;% 
%     min([trialStruct(t).keyFrames.heelupForce-trialStruct(t).keyFrames.frOffset, ...
%         trialStruct(t).keyFrames.min_mtp]); % minimum ofmtp angle or heel up force value
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    
    fr_start_s(t) =trialStruct(t).keyFrames.max_arch;% fr_start;
    fr_mtp_s(t)=  fr_mtp;

    dfield = subjStruct(s_ind).bones.tal.dfield;
    pts = subjStruct(s_ind).tib_dome.pts;%trialStruct(s).bones.tib.pts;
    conns = subjStruct(s_ind).tib_dome.cns;% trialStruct(s).bones.tib.cnt(:,1:3);
    thresh = surface_dist(s_ind);
    viz = 0;
    
    
    for fr = fr_start:fr_mtp
        RT_ref = trialStruct(t).Tm_mm.tal(:,:,fr);
        RT_pts = trialStruct(t).Tm_mm.tib(:,:,fr);
        if any(isnan([RT_ref(:);RT_pts(:)]))
            continue
        end
        %         if rem(fr,5) == 0
        %             out = JointContact_Dfield(dfield,pts,conns,thresh,1,viz,RT_ref,RT_pts);
        %         else
        
        % visualise the movement 
%         figure(401);
%       hold off;
%       cla
%         visualizeBone(subjStruct(s_ind).bones.tal.pts,subjStruct(s_ind).bones.tal.cnt(:,1:3),RT_ref);  hold on;
%         visualizeBone(pts,conns,RT_pts)
%         view([-104,5]); drawnow;
%         pause(0.1)
%          drawnow
        out = JointContact_Dfield(dfield,pts,conns,thresh,1,viz,RT_ref,RT_pts);
       
        trialStruct(t).contact.tib(fr).ContactCenter = out.ContactCenter;
        trialStruct(t).contact.tib(fr).ContactCenterDistance = out.ContactCenterDistance;
        trialStruct(t).contact.tib(fr).ContactArea = out.ContactArea;
        trialStruct(t).contact.tib(fr).MinDistance = out.MinDistance;
        trialStruct(t).contact.tib(fr).MinDistanceCentroid = out.MinDistanceCentroid;
        trialStruct(t).contact.tib(fr).AverageDistance = out.AverageDistance;
        trialStruct(t).contact.tib(fr).Penetration = out.Penetration;
        trialStruct(t).contact.tib(fr).Patch = out.Patch;
        
        cc_CT = out.ContactCenter;
        cc_xr = transformPoints(trialStruct(t).Tm_mm.tib(:,:,fr),out.ContactCenter);
        cN_CT = nanmean(out.Patch.wtdist.* out.Patch.norm);
        cN_xr = unit(transformPoints(trialStruct(t).Tm_mm.tib(:,:,fr), cN_CT));
        
        trialStruct(t).contact.tib(fr).cc_CT = cc_CT;
        trialStruct(t).contact.tib(fr).cc_xr = cc_xr;
        trialStruct(t).contact.tib(fr).cN_CT = cN_CT;
        trialStruct(t).contact.tib(fr).cN_xr= cN_xr;
        if out.Penetration == 1
            fprintf(' in frame %i.\n',fr)
        end
    end
   


    % keep everything fixed relative to the talus
    
%      
%     dfield = subjStruct(s_ind).bones.tib.dfield;
%     pts = subjStruct(s_ind).bones.tal.pts;%subjStruct(s_ind).tal_dome.pts;
%     conns =  subjStruct(s_ind).bones.tal.cnt(:,1:3);%subjStruct(s_ind).tal_dome.cns;%
%     thresh = 5;
%      cmap_jet = flipud(jet(thresh*10));
%     viz = 0;
%     clearvars('contact_dist_save','contact_dist_save2')
%     for fr = fr_start:fr_mtp
%        disp(fr)
%       
%         RT_pts = eye(4,4);
%         RT_ref =  invTranspose(trialStruct(t).Tm_mm.tal(:,:,fr)) * trialStruct(t).Tm_mm.tib(:,:,fr);
%         if any(isnan([RT_ref(:);RT_pts(:)]))
%             continue
%         end
%         % visualise the movement 
% %         figure(401);
% %       hold off;
% %       cla
% % %         visualizeBone(subjStruct(s_ind).bones.tal.pts,subjStruct(s_ind).bones.tal.cnt(:,1:3),RT_ref);  hold on;
% %         visualizeBone(subjStruct(s_ind).bones.tib.pts,subjStruct(s_ind).bones.tib.cnt(:,1:3),RT_ref);  hold on;
% %         visualizeBone(pts,conns,RT_pts)
% %         view([-126,-10])
% %         
% %         drawnow
% %    
% %         
%         
% %         figure(402);
%         out = JointContact_Dfield(dfield,pts,conns,thresh,1,viz,RT_ref,RT_pts);
% %         drawnow
%         
% %         pause(0.5)
%     
%     
%         trialStruct(t).contact.tal_fixed_tib(fr).ContactCenter = out.ContactCenter;
%         trialStruct(t).contact.tal_fixed_tib(fr).ContactCenterDistance = out.ContactCenterDistance;
%         trialStruct(t).contact.tal_fixed_tib(fr).ContactArea = out.ContactArea;
%         trialStruct(t).contact.tal_fixed_tib(fr).MinDistance = out.MinDistance;
%         trialStruct(t).contact.tal_fixed_tib(fr).MinDistanceCentroid = out.MinDistanceCentroid;
%         trialStruct(t).contact.tal_fixed_tib(fr).AverageDistance = out.AverageDistance;
%         trialStruct(t).contact.tal_fixed_tib(fr).Penetration = out.Penetration;
%         trialStruct(t).contact.tal_fixed_tib(fr).Patch = out.Patch;
%         
%         cc_xr = transformPoints(RT_pts,out.ContactCenter);
%         cN_CT = nanmean(out.Patch.wtdist.* out.Patch.norm);
%         cN_xr = unit(transformPoints(RT_pts, cN_CT));
%         
%         trialStruct(t).contact.tal_fixed_tib(fr).cc_xr = cc_xr;
%         trialStruct(t).contact.tal_fixed_tib(fr).cN_CT = cN_CT;
%         trialStruct(t).contact.tal_fixed_tib(fr).cN_xr= cN_xr;
%         
% %         contact_dist_save(:,fr) = out.Patch.wtdist;
%         
%         pts_xr = transformPoints(RT_pts,pts);
%         D = lookUpDfieldPts(dfield,pts_xr,RT_ref(1:3,1:3),RT_ref(1:3,4));
%         D(D>thresh) = thresh;
%         D(D<0.1) = 0.1;
%         contact_dist_save2(:,fr) = D;
%         cmap_pts = cmap_jet(round(D*10),:);
%         figure(403);
%         
%           visualizeBone(subjStruct(s_ind).bones.tal.pts,subjStruct(s_ind).bones.tal.cnt(:,1:3),RT_pts); hold on;
%         scatter3(pts_xr(:,1),pts_xr(:,2),pts_xr(:,3),5,cmap_pts); hold off;
%           view([120,30])
%           drawnow
%           pause(0.2)
%         if out.Penetration == 1
%             fprintf(' in frame %i.\n',fr)
%         end
%     end
%     
%     figure(150+t);
%    meanD =  rms(contact_dist_save2(:,fr_start:fr_mtp)');
%     cmap_pts = cmap_jet(round(meanD*10),:);
% %         figure(403);
%         scatter3(pts_xr(:,1),pts_xr(:,2),pts_xr(:,3),5,cmap_pts)
%           view([89 3])
%           drawnow
%           title(trialStruct(t).trial)
%           axis equal
    
    
%     
    
    
    
    dfield = subjStruct(s_ind).bones.tib.dfield;
    pts = subjStruct(s_ind).tal_dome.pts;% subjStruct(s_ind).bones.tal.pts;%
    conns = subjStruct(s_ind).tal_dome.cns;% subjStruct(s_ind).bones.tal.cnt(:,1:3);%
    thresh = surface_dist(s_ind);
    viz = 0;
    
    
    for fr = fr_start:fr_mtp
        RT_ref = trialStruct(t).Tm_mm.tib(:,:,fr);
        RT_pts = trialStruct(t).Tm_mm.tal(:,:,fr);
         if any(isnan([RT_ref(:);RT_pts(:)]))
            continue
        end
        %         figure(1); patch('faces',conns,'vertices',transformPoints(RT_ )
        %         if rem(fr,5) == 0
        %             out = JointContact_Dfield(dfield,pts,conns,thresh,1,viz,RT_ref,RT_pts);
        %         else
        
        out = JointContact_Dfield(dfield,pts,conns,thresh,1,0,RT_ref,RT_pts);
        
        trialStruct(t).contact.tal(fr).ContactCenter = out.ContactCenter;
        trialStruct(t).contact.tal(fr).ContactCenterDistance = out.ContactCenterDistance;
        trialStruct(t).contact.tal(fr).ContactArea = out.ContactArea;
        trialStruct(t).contact.tal(fr).MinDistance = out.MinDistance;
        trialStruct(t).contact.tal(fr).MinDistanceCentroid = out.MinDistanceCentroid;
        trialStruct(t).contact.tal(fr).AverageDistance = out.AverageDistance;
        trialStruct(t).contact.tal(fr).Penetration = out.Penetration;
        trialStruct(t).contact.tal(fr).Patch = out.Patch;
        
        % get the normal at the contact centroid surface - weight all the
        % normals by their distance
        cc_CT = out.ContactCenter;
        cc_xr = transformPoints(trialStruct(t).Tm_mm.tal(:,:,fr),out.ContactCenter);
        cN_CT = nanmean(out.Patch.wtdist.* out.Patch.norm);
        cN_xr = unit(transformPoints(trialStruct(t).Tm_mm.tal(:,:,fr), cN_CT));
        
        trialStruct(t).contact.tal(fr).cc_CT = cc_CT;
        trialStruct(t).contact.tal(fr).cc_xr = cc_xr;
        trialStruct(t).contact.tal(fr).cN_CT= cN_CT;
        trialStruct(t).contact.tal(fr).cN_xr = cN_xr;
    end
    
    
    
    % Now see where the tibia is on the talus if we rotate it backwards
    % around the talocrural axis (i.e. SIMULATED TIBIA
    
    dfield = subjStruct(s_ind).bones.tib.dfield;
    pts = subjStruct(s_ind).tal_dome.pts;% subjStruct(s_ind).bones.tal.pts;%
    conns = subjStruct(s_ind).tal_dome.cns;% subjStruct(s_ind).bones.tal.cnt(:,1:3);%
    RT_ref= trialStruct(t).kinemat.T_tib_fix_act_mtp;
    RT_pts = trialStruct(t).T_fix.tal(:,:,fr_mtp);
    %         thresh = 8;
    
    if any(isnan([RT_ref(:);RT_pts(:)]))
        continue
    end
    %         figure(1); patch('faces',conns,'vertices',transformPoints(RT_ )
    %         if rem(fr,5) == 0
    %             out = JointContact_Dfield(dfield,pts,conns,thresh,1,viz,RT_ref,RT_pts);
    %         else
    dd = lookUpDfieldPts(dfield,transformPoints(RT_pts,pts),RT_ref(1:3,1:3),RT_ref(1:3,4));
    thresh = median(dd);
    out = JointContact_Dfield(dfield,pts,conns,thresh,1,0,RT_ref,RT_pts);
    
    trialStruct(t).contact.tib_sim(fr_mtp).ContactCenter = out.ContactCenter;
    trialStruct(t).contact.tib_sim(fr_mtp).ContactCenterDistance = out.ContactCenterDistance;
    trialStruct(t).contact.tib_sim(fr_mtp).ContactArea = out.ContactArea;
    trialStruct(t).contact.tib_sim(fr_mtp).MinDistance = out.MinDistance;
    trialStruct(t).contact.tib_sim(fr_mtp).MinDistanceCentroid = out.MinDistanceCentroid;
    trialStruct(t).contact.tib_sim(fr_mtp).AverageDistance = out.AverageDistance;
    trialStruct(t).contact.tib_sim(fr_mtp).Penetration = out.Penetration;
    trialStruct(t).contact.tib_sim(fr_mtp).Patch = out.Patch;
    
    % get the normal at the contact centroid surface - weight all the
    % normals by their distance
    cc_CT = out.ContactCenter;
    cc_xr = transformPoints(RT_pts,out.ContactCenter);
    cN_CT = nanmean(out.Patch.wtdist.* out.Patch.norm);
    cN_xr = unit(transformPoints(RT_pts, cN_CT));
    
    trialStruct(t).contact.tib_sim(fr_mtp).cc_CT = cc_CT;
    trialStruct(t).contact.tib_sim(fr_mtp).cc_xr = cc_xr;
    trialStruct(t).contact.tib_sim(fr_mtp).cN_CT= cN_CT;
    trialStruct(t).contact.tib_sim(fr_mtp).cN_xr = cN_xr;
    
    
    
    fprintf('... Complete!\n')
end

%% Wrist Viz : Visualize the contact at any frame but show ALL subjs in one, in tal co sys - make wrist viz

ivstring = createInventorHeader();

fr_name = 'Max_MTP'; % name this frame
contact_type = 'tib_sim'; %'tal' % typically 'tal' but if you wanted to see the frame of tibia simulation - 'tib_sim'

offset_plot = 0;
for t = 1:ntrials-2%[[1:2:ntrials-2], [ 2:2:ntrials-2]]
    
    fr_look = trialStruct(t).keyFrames.max_mtp; %trialStruct(t).keyFrames.max_dors; % the frame to look at
    
    
    offset_plot = offset_plot + 1; % to plot them all in order of walking and then running
    if rem(t,2) == 0 % put the runs forward
        forward = -300;
        col = [0 0.2 0.2]; % green blue
        run_flag = 1;
    else
        forward = 0;
        col = [0.4 0 0]; % red
        run_flag = 0;
    end
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    Tm = trialStruct(t).Tm_mm;
    nfr = size(Tm.tib,3);
    bone_listM = fields(Tm);
    nBones = length(bone_listM);
    
    
    %     fr_tg = trialStruct(t).keyFrames.fr_tg ;
    
    anim_dir = fullfile(dir_analy, 'Animation',['ALL ' contact_type ' Contact_' fr_name],filesep);
    if ~exist(anim_dir,'dir')
        mkdir(anim_dir)
    end
    
    
    
    %------ make the actual tibia/rigid foot iv file -------------------
    pat = trialStruct(t).contact.(contact_type)(fr_look).Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,[trialStruct(t).trial '_surfaceSimNormal.iv']))
    
    if strcmp(contact_type,'tib_sim')
        T_tal = trialStruct(t).T_fix.tal(:,:,fr_look);
    else
        T_tal = Tm.tal(:,:,fr_look);
    end
    
    T_anat = subjStruct(s_ind).bones.tal.T_ACS.T_TC;
    T_anatInv = T_tal;%invTranspose(T_anat);
    cc_anat = transformPoints(T_anatInv,trialStruct(t).contact.(contact_type)(fr_look).cc_CT);
    %     ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),T_tal(1:3,1:3),T_tal(1:3,4)')];
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,[trialStruct(t).trial '_surfaceSimNormal.iv']),T_anatInv(1:3,1:3),T_anatInv(1:3,4)'+[50*offset_plot forward 0])];
    ivstring = [ivstring createInventorSphere(cc_anat,2,col,0.2)];
    ivstring = [ivstring createInventorSphere(cc_anat+[50*offset_plot forward 0],2,col,0.2)];
    %
    %     if run_flag == 1 % plot the walking centre on the running one
    %         fr_lookA = trialStruct(t-1).keyFrames.max_dors;
    %         cc_anat = transformPoints(T_anatInv,trialStruct(t-1).contact.(contact_type)(fr_lookA).cc_CT); % plots the same subject's walk
    %         ivstring = [ivstring createInventorSphere(cc_anat+[100*(offset_plot) forward 0],2,[0.8 0 0],0.2)];
    %     end
    
    for bn = 1:nBones
        if strcmp(contact_type,'tib_sim')
            Tf = trialStruct(t).T_fix.(bone_listM{bn})(:,:,fr_look);
        else
            Tf =   Tm.(bone_listM{bn})(:,:,fr_look); % T of this frame
        end
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            if strcmp(contact_type,'tib_sim')
                Tf = trialStruct(t).kinemat.T_tib_fix_act_mtp;
            end
            
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)'+[50*offset_plot forward 0],[0.7 0.7 0.7],0.8)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)'+[50*offset_plot forward 0],[0.7 0.7 0.7],0.2)];
        end
        
    end
    
    
    
    
    
end
fid = fopen(fullfile(anim_dir,['TalusContactMoving' fr_name '.iv']),'w');
fprintf(fid,ivstring);
fclose(fid);
%% Wrist Viz : Visualize the contact at any frame for each subject


fr_name = 'Max_MTP'; % name this frame
contact_type = 'tib_sim'; %'tal' % typically 'tal' but if you wanted to see the frame of tibia simulation - 'tib_sim'

offset_plot = 0;
for t = 1:ntrials-2%[[1:2:ntrials-2], [ 2:2:ntrials-2]]
    
ivstring = createInventorHeader();
    fr_look = trialStruct(t).keyFrames.max_mtp; %trialStruct(t).keyFrames.max_dors; % the frame to look at
    
    
    offset_plot = offset_plot + 1; % to plot them all in order of walking and then running
    if rem(t,2) == 0 % put the runs forward
        forward = -300;
        col = [0 0.2 0.2]; % green blue
        run_flag = 1;
    else
        forward = 0;
        col = [0.4 0 0]; % red
        run_flag = 0;
    end
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    Tm = trialStruct(t).Tm_mm;
    nfr = size(Tm.tib,3);
    bone_listM = fields(Tm);
    nBones = length(bone_listM);
    
    
    %     fr_tg = trialStruct(t).keyFrames.fr_tg ;
    
    anim_dir = fullfile(dir_analy, 'Animation',['INDIVIDUAL ' contact_type ' Contact_' fr_name],filesep);
    if ~exist(anim_dir,'dir')
        mkdir(anim_dir)
    end
    
    
    
    %------ make the actual tibia/rigid foot iv file -------------------
    pat = trialStruct(t).contact.(contact_type)(fr_look).Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,[trialStruct(t).trial '_surfaceSimNormal.iv']))
    
    if strcmp(contact_type,'tib_sim')
        T_tal = trialStruct(t).T_fix.tal(:,:,fr_look);
    else
        T_tal = Tm.tal(:,:,fr_look);
    end
    
    T_anat = subjStruct(s_ind).bones.tal.T_ACS.T_TC;
    T_anatInv = T_tal;%invTranspose(T_anat);
    cc_anat = transformPoints(T_anatInv,trialStruct(t).contact.(contact_type)(fr_look).cc_CT);
    %     ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),T_tal(1:3,1:3),T_tal(1:3,4)')];
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,[trialStruct(t).trial '_surfaceSimNormal.iv']),T_anatInv(1:3,1:3),T_anatInv(1:3,4)')];
    ivstring = [ivstring createInventorSphere(cc_anat,2,col,0.2)];
    ivstring = [ivstring createInventorSphere(cc_anat+[100*offset_plot forward 0],2,col,0.2)];
    %
    %     if run_flag == 1 % plot the walking centre on the running one
    %         fr_lookA = trialStruct(t-1).keyFrames.max_dors;
    %         cc_anat = transformPoints(T_anatInv,trialStruct(t-1).contact.(contact_type)(fr_lookA).cc_CT); % plots the same subject's walk
    %         ivstring = [ivstring createInventorSphere(cc_anat+[100*(offset_plot) forward 0],2,[0.8 0 0],0.2)];
    %     end
    
    for bn = 1:nBones
        if strcmp(contact_type,'tib_sim')
            Tf = trialStruct(t).T_fix.(bone_listM{bn})(:,:,fr_look);
        else
            Tf =   Tm.(bone_listM{bn})(:,:,fr_look); % T of this frame
        end
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            if strcmp(contact_type,'tib_sim')
                Tf = trialStruct(t).kinemat.T_tib_fix_act_mtp;
            end
            
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
    end
    
    
    if ~exist(fullfile(anim_dir,trialStruct(t).trial),'dir')
        mkdir(fullfile(anim_dir,trialStruct(t).trial));
    end
fid = fopen(fullfile(anim_dir,trialStruct(t).trial,['TalusContactMoving' fr_name '.iv']),'w');
fprintf(fid,ivstring);
fclose(fid);
    
    
end
%% Measure the contact area for each trial
close all
plotI = 0;
for t = 1:ntrials-2%[1:2:12,2:2:12]%
    plotI = plotI +1;
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    fr_start = trialStruct(t).keyFrames.Fprop - trialStruct(t).keyFrames.frOffset;% fr_start_s(t);%trialStruct(t).keyFrames.heelupForce-trialStruct(t).keyFrames.frOffset;
    cmapP = parula(fr_mtp-fr_start+1);
    
    %     T_anat = rotateCoordSys(subjStruct(s_ind).bones.tal.T_ACS.T_TC,-45,1); % talocrural axis
    %     T_anat = rotateCoordSys(T_anat,45,3); % talocrural axis
    pts_dome_anat = transformPoints(subjStruct(s_ind).bones.tal.T_ACS.T_TC,subjStruct(s_ind).tal_dome.pts,-1);
    coeff = pca(pts_dome_anat);
    
    T_pca = eye(4);
    T_pca(1:3,1:3) = coeff;
    pts_dome = transformPoints(T_pca, pts_dome_anat,-1);
    
    T_total =  invTranspose(T_pca) * invTranspose(subjStruct(s_ind).bones.tal.T_ACS.T_TC);
    
    pts_dome2 = transformPoints(T_total,subjStruct(s_ind).tal_dome.pts);
    
    
    figure;
    plot3quick_scatter(pts_dome'); hold on; % plot the talar dome in the rotated co-sys
    plot3quick_scatter(pts_dome2'); hold on; axis equal; % plot talar dome transformed with whole system
    
    [xv,yv,curve] = fitCurveToTalarDome(pts_dome,1000); % the fit on each subjects dome
    
    if curve(-100) > curve(100)% the curve is flipped, the PCA has put it backwards, switch it and recalculate it such that anterior is always +x
        T_pca = rotateCoordSys(T_pca,180,3); % 180 deg flip
        pts_dome = transformPoints(T_pca,pts_dome_anat,-1);
        T_total =  invTranspose(T_pca) * invTranspose(subjStruct(s_ind).bones.tal.T_ACS.T_TC);
        [xv,yv,curve] = fitCurveToTalarDome(pts_dome,1000); % recompute the fit on each subjects dome
        
        pts_dome2 = transformPoints(T_total,subjStruct(s_ind).tal_dome.pts);
    end
    
    cla
    
    plot3quick_scatter(pts_dome'); hold on; % plot the talar dome in the rotated co-sys
    plot3quick_scatter(pts_dome2'); hold on; axis equal; % plot talar dome transformed with whole system
    % find all points that are within a threshold value of the centre plane
    % (sagittal plane)
    pts_plane = [];dis = [];
    for p = 1:size(pts_dome,1)
        pts_plane(p,:) = closestPointonPlanealongVector(pts_dome(p,:),[0 1 0],[0 0 0], [0 1 0]);  % talar dome in the sagittal plane
        dis(p) = norm(pts_dome(p,:)-pts_plane(p,:)); % distance of point from the plane
    end
    Iclose = dis < 4; % where distance from the plane is within +- 4mm
    
    [xv,yv,curve] = fitCurveToTalarDome(pts_dome(Iclose,:),1000);
    
    totalROM(t) = lineLength([xv,yv]); % the entire length of the curve
    
    
    
    
    %     figure(999);
    %     subplot(4,3,t)
    %     plot3quick_scatter(pts_dome_anat'); hold on; % plot the talar dome in the rotated co-sys
    %     plot3quick_scatter(pts_dome(Iclose,:)'); hold on; axis equal;
    
    tt = 1;
    all_pts = [];
    
    for fr = fr_start:fr_mtp-1
        %         plot3quick(trialStruct(t).contact.tal(fr).Patch.pts,[0.5 0 0],'.','none'); hold on
        points_init = transformPoints(T_total,trialStruct(t).contact.tal(fr).Patch.pts); % patch pts are in CT
        points_used = intersect(round(points_init,2),round(pts_dome(Iclose,:),2),'rows');
        %       plot3quick(points_init,cmapP(tt,:),'.','none'); hold on;
        %         plot3quick(points_used,cmapP(tt,:),'.','none'); hold on;
        concat_pts = [all_pts;points_used];
        [~,Iinc]= unique(concat_pts(:,1));
        all_pts = concat_pts(Iinc,:);
        %         plot3quick(all_pts,cmapP(tt,:),'o','none'); hold on;
        %         drawnow
        %         view([0,0])
        tt = tt+1;
    end
    %    figure;
    %    plot3quick(all_pts,cmapP(1,:),'.','none'); hold on;
    %     T_Al = subjStruct(s_ind).bones.tal.T_Aligned;
    
    
    
    %    plot3quick(pts_dome,cmapP(1,:),'.','none'); hold on;
    
    pts_ROM = all_pts; % now project all the points into the plane and find the range of x-values
    pts_ROM_plane =[];
    dis = [];
    for p = 1:size(pts_ROM,1)
        pts_ROM_plane(p,:) = closestPointonPlanealongVector(pts_ROM(p,:),[0 1 0],[0 0 0], [0 1 0]);  % talar dome in the sagittal plane
        %         dis(p) = norm(pts_ROM(p,:)-pts_ROM_plane(p,:));
    end
    
    
    %     figure; plot( pts_ROM_plane(:,2), pts_ROM_plane(:,3),'.');hold on;
    
    x_range = linspace(min(pts_ROM_plane(:,1)),max(pts_ROM_plane(:,1)),200)';
    y_range = curve(x_range);
    rom_length(t) = lineLength([x_range,y_range]);
    %    plot(x_range,y_range)
    
    
    first_val = linspace(xv(1),x_range(1),200)';
    first_y = curve(first_val);
    first_seg = lineLength([first_val,first_y]);
    start_perc(t) = first_seg/totalROM(t);
    
    last_val = linspace(x_range(end),xv(end),200)';
    last_y = curve(last_val);
    last_seg = lineLength([last_val,last_y]);
    end_perc(t) = (totalROM(t)-last_seg)/totalROM(t);
    %
    
    
    %     ----------- see the farthest posterior point in the simulated
    %     positions-----------------------------------
    tt = 1;
    all_pts = [];
    
    for fr = fr_mtp
        
        points_init = transformPoints(T_total,trialStruct(t).contact.tib_sim(fr).Patch.pts); % patch pts are in CT
        points_used = intersect(round(points_init,2),round(pts_dome(Iclose,:),2),'rows');
        
        plot3quick(points_init,cmapP(tt,:),'.','none'); hold on;
        plot3quick(points_used,cmapP(tt+10,:),'.','none'); hold on;
        concat_pts = [all_pts;points_used];
        [~,Iinc]= unique(concat_pts(:,1));
        all_pts = concat_pts(Iinc,:);
        %         plot3quick(all_pts,cmapP(tt,:),'o','none'); hold on;
        %         drawnow
        %         view([0,0])
        tt = tt+1;
    end
    
    pts_ROM = all_pts; % now project all the points into the plane and find the x-value
    pts_ROM_plane =[];
    dis = [];
    for p = 1:size(pts_ROM,1)
        pts_ROM_plane(p,:) = closestPointonPlanealongVector(pts_ROM(p,:),[0 1 0],[0 0 0], [0 1 0]);  % talar dome in the sagittal plane
        dis(p) = norm(pts_ROM(p,:)-pts_ROM_plane(p,:));
    end
    x_lim = linspace(min(pts_ROM_plane(:,1)),x_range(1),200)';%min(pts_ROM_plane(:,1));
    y_lim = curve(x_lim);
    
    totalROM_lost(t) = lineLength([x_lim,y_lim]);
    
    x2lost = linspace(xv(1),x_lim(1),200)';
    y2lost = curve(x2lost);
    lost_seg = lineLength([x2lost,y2lost]);
    lost_seg_perc(t) = lost_seg/totalROM(t);
    
    first_val = linspace(xv(1),x_range(1),200)';
    first_y = curve(first_val);
    first_seg = lineLength([first_val,first_y]);
    start_perc(t) = first_seg/totalROM(t);
    
    y_p = @(c,x) 5*c.p1 * x.^4 + 4*c.p2*x.^3  + 3* c.p3*x.^2  + 2*c.p4*x+ c.p5;
    
    y_pp = @(c,x) 20*c.p1*x.^3 + 12*c.p2*x.^2 +6 * c.p3 * x + 2*c.p4;
    
    
    
    xx = -30:0.1:30;
    s = [xx;curve(xx)'];
    ds_dt = [repmat(0.1,1,length(xx));y_p(curve,xx)];
    T = ds_dt./sqrt(ds_dt(1,:).^2 + ds_dt(2,:).^2);
    
    %     figure; hold on; plot(s(1,:),s(2,:)); quiver(s(1,1:10:end),s(2,1:10:end),T(1,1:10:end),T(2,1:10:end)); axis equal
    %     dT_dt = ([diff(T(1,:));diff(T(2,:))]/0.1);
    %     kur = sqrt(dT_dt(1,:).^2 + dT_dt(2,:).^2)./sqrt(ds_dt(1,1:end-1).^2 + ds_dt(2,1:end-1).^2)
    %     curvatu = deriv(curve,xx);
    %     figure; yyaxis right; plot(xx,curvatu); ylim([-1 1]);hold on;plot(xx([1,end]),[0 0 ],'k:'); plot(xx(1:end-1),kur)
    % hold on;
    % yyaxis left; plot(xx,curve(xx));
    % make the curved dome plots with ROM
    p = trialStruct(t).plot;
    figure(1001);
    %     plot(first_val,first_y+plotI*6,'k.');
    
    %     plot(last_val,last_y+plotI*6,'k.'); hold on;
    
    pI = ceil(t/2);
    subplot(2,3,pI)
    plot(xv,yv+pI*0,'k.');hold on; % the black line
    
    if rem(t,2) == 1 % walk
        ex = 3;
        h(1) = plot(x_range,y_range+pI*0+ex,'-.','Color',p.col); hold on;
        
        plot(x_range([1,end]),y_range([1,end])+pI*0+ex,'Color',p.col,'marker',p.marker_type,'linestyle','none');% plot the ends of the curve
        plot(x_lim,y_lim+ex,':','Color',p.col)
    else % run
        ex = -3;
        
        h(2) = plot(x_range,y_range+pI*0+ex,'-','Color',p.col); hold on;
        plot(x_range([1,end]),y_range([1,end])+pI*0+ex,'Color',p.col,'marker',p.marker_type,'linestyle','none'); hold on;
        h(3) = plot(x_lim,y_lim+ex,':','Color',p.col);
    end
    
    % plot(x_lim,y_lim,'o')
    %     title('Visual representation of Joint surface use')
    axis equal; axis off
    makeNicePlotsFunction
    figure(1000)
    hold on;
    if rem(t,2) == 1
        plot([start_perc(t),end_perc(t)]*100,[plotI,plotI],'Color',p.col,'marker',p.marker_type,'linestyle','-.','MarkerFaceColor',p.marker_fill);
    else
        plot([start_perc(t),end_perc(t)]*100,[plotI,plotI],'Color',p.col,'marker',p.marker_type,'linestyle','-','MarkerFaceColor',p.marker_fill);
    end
    plot([lost_seg_perc(t) ,start_perc(t)]*100,[plotI,plotI],'Color',p.col,'marker',p.marker_type,'linestyle',':','MarkerFaceColor',p.marker_fill);
    xlabel('% of ROM')
    ylabel('Subject')
    makeNicePlotsFunction
    xlim([0 100])
end
legend(h(1:3),{ 'Walk','Run','Required ROM if arch is rigid'})

figure; plot(repmat(1:2,1,6)',totalROM_lost,'o')%,{'walk','run'})
xlim([0 3])
ylabel('Additional required ankle ROM with rigid arch')

%% Talar centroid analysis

for t = 1:ntrials-2
    T_fix = trialStruct(t).T_fix;
    Tm = trialStruct(t).Tm_mm;
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    fr_start = fr_start_s(t);
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    %------------- TALAR HEIGHT ANALYSIS---------------------------------
    trialStruct(t).kinemat.tal_cent_fix = transformPoints(T_fix.tal,subjStruct(s_ind).bones.tal.centroid);
    trialStruct(t).kinemat.tal_cent_short = transformPoints(Tm.tal,subjStruct(s_ind).bones.tal.centroid);
    
    T_ACS_ph1{t} = Tm.ph1(:,:,fr_start) * subjStruct(s_ind).bones.ph1.T_Aligned;
    
    tal_ht_init(t,:) = transformPoints(T_ACS_ph1{t} , trialStruct(t).kinemat.tal_cent_fix(fr_start,:),-1);
    tal_ht_act(t,:) = transformPoints(T_ACS_ph1{t} , trialStruct(t).kinemat.tal_cent_short(fr_mtp,:),-1);
    tal_ht_fix(t,:) = transformPoints(T_ACS_ph1{t} , trialStruct(t).kinemat.tal_cent_fix(fr_mtp,:),-1);
    
    
    tal_range(t,:) = tal_ht_act(t,:)- tal_ht_init(t,:);
    tal_diff(t,:) = (tal_ht_fix(t,:) - tal_ht_act(t,:));
    tal_perc(t,:) = tal_diff(t,:)./tal_range(t,:) * 100;
    
    
    
    mtp_vals =  trialStruct(t).kinemat.mtp([fr_start,fr_tib_glob_align(t),fr_mtp]);
    mtp_perc_tg(t) = (mtp_vals(2)-mtp_vals(1))/ (mtp_vals(3)-mtp_vals(1))*100;
    %       tal_ht_act_tg(s,:) = transformPoints(T_ACS_ph1{s} , trialStruct(s).kinemat.tal_cent_short(fr_tib_glob_align(s),:),-1);
    tal_ht_fix_tg(t,:) = transformPoints(T_ACS_ph1{t} , trialStruct(t).kinemat.tal_cent_fix(fr_tib_glob_align(t),:),-1); % look at take off of rigid foot
    
    tal_diff_tg(t,:) = (tal_ht_act(t,:) - tal_ht_fix_tg(t,:))./(tal_ht_act(t,:) -tal_ht_init(t,:))*100;
    
    
    
    % -----------determine the global tibia lean of the fixed arch-------
    
    align_z_rot_mov = trialStruct(t).kinemat.tib_lean.rot_mov;
   T_acs_tib_x=[];align_z_rot_fix=nan(nfr_tot,1);
    for fr = fr_start:fr_mtp
        T_acs_tib_x(:,:,fr) = T_fix.tib(:,:,fr)* subjStruct(s_ind).bones.tib.T_ACS.T_TC;
        [align_z_rot_fix(fr),~,~] = eulerYZX(eye(4,4),T_acs_tib_x(:,:,fr),eye(4,4),eye(4,4)); % the global axis rotation
        
    end
    [~,fr_align_s(t)] = min(abs(align_z_rot_mov(fr_mtp) - align_z_rot_fix));
    
    trialStruct(t).kinemat.tib_lean.rot_fix = align_z_rot_fix;
    
end
%% calculate all the axes, moment arms etc
% glb is global, identity


% clearvars -except trialStruct subjStruct dir_analy cmap_crop cmap nsubj ntrials dir_ref fr_start_s fr_mtp_s demo res
clear stancePts

load('colorblind_map.mat') % loads the cmap_cb variable

close all
bone_refs = {'tib','tal';'glb','mt1';'cal','mt1'};

for t = 1:ntrials-2
    
    Tm = trialStruct(t).Tm_mm;
    nfr_tot = size(Tm.cal,3);
    
    frsX = trialStruct(t).keyFrames.cropFrsXray;
    frsF = trialStruct(t).keyFrames.cropFrsForce;
    nfrs_stance = abs(diff(frsX))+1;
    % important gait events
    fr_start = fr_start_s(t);
    fr_mtp = fr_mtp_s(t);
    % I frames depict that they are from the first frame of cropped data
    Ifr_hu = trialStruct(t).keyFrames.heelupForce - frsF(1)+1;
    Ifr_min_mtp = trialStruct(t).keyFrames.min_mtp - frsX(1) + 1 ;
    Ifr_mxad = trialStruct(t).keyFrames.max_dors - frsX(1) + 1;
    Ifr_mxarch = trialStruct(t).keyFrames.max_arch - frsX(1) + 1;
%     Ifr_prop = trialStruct(t).keyFrames.Fprop - frsF(1) + 1;
    
    
    frs_interest = [fr_start-frsX(1)+1,Ifr_hu,Ifr_min_mtp,Ifr_mxad,Ifr_mxarch,fr_mtp-frsX(1)];
    frs_interest_s(t,:) = frs_interest;
    
    stancePts{t} = linspace(0,100,nfrs_stance);

    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    for br = 2;%1:size(bone_refs,1)
        % total number of frames
        ntot = size(Tm.(bone_refs{br,2}),3);
        
        n_ct = nan(ntot,3); n_xr = nan(ntot,3);
        s_ct = nan(ntot,3); s_xr = nan(ntot,3);
        phi = nan(1,ntot); L = nan(1,ntot);
        s_xrC = nan(ntot,3); cent_xr = nan(ntot,3);
        
        
        
        % take every  frame between max ankle dorsiflexion and max mtp dorsiflexion
        
         % !!! added 5 frames to stabilize the helical axis and then removed
            % them after
            fr_add = 0;
            T = Tm.(bone_refs{br,2})(:,:,fr_start-fr_add:fr_mtp);
            nfrs = fr_mtp-fr_start+fr_add+1;
            TR = repmat(nan(4,4),1,1,nfrs);

            if strcmp(bone_refs{br,1},'glb')
%                 for fr = fr_start-fr_add:fr_mtp
%                     TR(:,:,1:nfrs) = repmat(Tm.(bone_refs{br,2})(:,:,fr_mtp),1,1,nfrs);% * subjStruct(s_ind).bones.(bone_refs{br,2}).T_Aligned; % reference to global, i.e. identity matrix
%                 end
                TR = repmat(eye(4),1,1,nfrs);
            else
                TR = Tm.(bone_refs{br,1})(:,:,fr_start-fr_add:fr_mtp);
            end
%         if contains(trialStruct(t).trial,'run')

           T = filterTransforms(T,[14 20],trialStruct(t).marker_data.FrameRate,0);
%            TR = filterTransforms(TR,[10,30],trialStruct(t).marker_data.FrameRate,0);
%             
%             if length(fr_start-fr_add:fr_mtp) < 9
%                  [phi(fr_start-fr_add:fr_mtp),~,L(fr_start-fr_add:fr_mtp),~] = helicalInstantaneousGapFill(TR,T,3);
%             else
%             [phi(fr_start-fr_add:fr_mtp),~,L(fr_start-fr_add:fr_mtp),~] = stabilizeHelicalAxis(TR,T,5,7);
%             
%             end
%             
%             
%             mp1_2 = find(phi(fr_start-fr_add:fr_mtp) > 1,1,'first');% get the midpoint frames to average the HA between - sufficient rotation for stability
%             mp3_4 = nfrs;% round(3*(nfrs)/4);
%             
%             disp(mp1_2-mp3_4)
% %             [~,n_xr_mp,~,s_xr_mp] = helicalInstantaneous(TR(:,:,[mp1_2,mp3_4]) ,T(:,:,[mp1_2,mp3_4]));
%             [~,n_xr(:,fr_start+mp1_2-1:fr_mtp),~,s_xr(:,fr_start+mp1_2-1:fr_mtp)] = stabilizeHelicalAxis(TR(:,:,mp1_2:mp3_4),T(:,:,mp1_2:mp3_4),3,5);
% 
%         
            
            
%         n_xr_mp = n_xr_mp';
%         s_xr_mp = s_xr_mp';
% 
%         n_ct_mp = transformVectors(T(:,:,mp1_2),n_xr_mp,-1);
%         s_ct_mp = transformPoints(T(:,:,mp1_2),s_xr_mp,-1);
%         
%         n_xr(fr_start:fr_mtp,:) = transformVectors(T,n_ct_mp);
%         s_xr(fr_start:fr_mtp,:) = transformPoints(T,s_ct_mp);
        
        
           
% 
%         n_xr = n_xr';
%         s_xr = s_xr';
%         n_xr(fr_start:fr_start+mp1_2-1,:) = repmat(n_xr(fr_start+mp1_2,:),mp1_2,1);
%         s_xr(fr_start:fr_start+mp1_2-1,:) = repmat(s_xr(fr_start+mp1_2,:),mp1_2,1);
%     
%         n_ct(fr_start:fr_mtp,:)  = transformVectors(T(:,:,[repmat(mp1_2,1,mp1_2-1),mp1_2:mp3_4]),n_xr(fr_start:fr_mtp,:),-1);
%         s_ct(fr_start:fr_mtp,:)  = transformPoints(T(:,:,[repmat(mp1_2,1,mp1_2-1),mp1_2:mp3_4]),s_xr(fr_start:fr_mtp,:),-1);
%         
%         n_xr(fr_start:fr_mtp,:) = transformVectors(T,n_ct(fr_start:fr_mtp,:) );
%         s_xr(fr_start:fr_mtp,:) = transformPoints(T,s_ct(fr_start:fr_mtp,:) );
%         
        frs =(fr_start+1:fr_mtp);
         [phi_temp,n_temp,L_temp,s_temp] = helicalInstantaneous(TR ,T);
         phi(frs) = phi_temp;
         n_xr(frs,:) = n_temp';
         L(frs) = L_temp;
         s_xr(frs,:)=s_temp';
         
         
        n_ct(frs,:) = transformVectors(T,n_xr(frs,:),-1);
        s_ct(frs,:) = transformPoints(T,s_xr(frs,:),-1);
        
        cent = subjStruct(s_ind).bones.(bone_refs{br,2}).centroid;
        cent_xr(fr_start:fr_mtp,:) = transformPoints(T,cent);
      
%         cnt = 1;
        for fr = frs% fr_start:fr_mtp
            s_xrC(fr,:) =  closestPointonVector( cent_xr(fr,:),s_xr(fr,:),n_xr(fr,:),0);
%             figure(1000);hold on;
%             
%             visualizeBone(subjStruct(s_ind).bones.cal.pts,subjStruct(s_ind).bones.cal.cnt(:,1:3),Tm.cal(:,:,fr),[0 0 0.3])
%             visualizeBone(subjStruct(s_ind).bones.mt1.pts,subjStruct(s_ind).bones.mt1.cnt(:,1:3),T(:,:,cnt),[0 0 0.3])
%             plotvector3(s_ct(fr,:),n_ct(fr,:)*100)
%             cnt = cnt+1;
        end
        
        
        
        trialStruct(t).helical(br).ref_bone = bone_refs{br,1};
        trialStruct(t).helical(br).mv_bone = bone_refs{br,2};
        trialStruct(t).helical(br).phi = phi;
        trialStruct(t).helical(br).phi_sum = nancumtrapz(phi);
        trialStruct(t).helical(br).n_xr = n_xr;
        trialStruct(t).helical(br).s_xr = s_xr;
        trialStruct(t).helical(br).n_ct = n_ct;
        trialStruct(t).helical(br).s_ct = s_ct;
        trialStruct(t).helical(br).cent_xr = cent_xr;
        trialStruct(t).helical(br).s_xrC = s_xrC;
        trialStruct(t).helical(br).L = L;
        
        
    end
%     ind_nonan = find(~isnan(trialStruct(t).helical(3).phi_sum));
%     phi_ROM_arch(t) = trialStruct(t).helical(3).phi_sum(ind_nonan(end));
%     
    
    
    
    % calculate the moment arms for the moving and the rigid foot
    
    % MOVING-----------------------------------------------------------
    T_tib = Tm.tib;
    tib_ach_p = transformPoints(T_tib, trialStruct(t).achilles.prox_pt_CT);
    T_cal =Tm.cal;
    cal_ach_d = transformPoints(T_cal, trialStruct(t).achilles.cal_pt_CT	);
    ach_vec = tib_ach_p - cal_ach_d;
    
    MA_ank = nan(ntot,3);
    MA_ach =  nan(ntot,3);
    for fr = fr_start:fr_mtp
        
        if fr > length(trialStruct(t).contact.tal)
            continue
        end
         % use the talus/ or tibia
        TC_pt =  transformPoints(Tm.tal(:,:,fr),trialStruct(t).contact.tal(fr).cc_CT);
        TC_vec =   transformVectors(Tm.tal(:,:,fr),trialStruct(t).contact.tal(fr).cN_CT);
       
        % ankle moment arm
        [s1,s2] = closestPointsBtw2Lines(TC_pt  , trialStruct(t).helical(2).s_xrC(fr,:),...
            -TC_vec, trialStruct(t).helical(2).n_xr(fr,:));
        MA_ank(fr,:) = s2-s1;
        
        % achilles moment arm
        [s1,s2] = closestPointsBtw2Lines(cal_ach_d(fr,:), trialStruct(t).helical(2).s_xrC(fr,:),...
            ach_vec(fr,:), trialStruct(t).helical(2).n_xr(fr,:));
        MA_ach(fr,:) =s2-s1;
        
        
    end
    
    trialStruct(t).achilles.tib_ach_xr = tib_ach_p;
    trialStruct(t).achilles.cal_ach_xr = cal_ach_d;
    trialStruct(t).achilles.ach_vec = ach_vec;
    trialStruct(t).MA.ank = MA_ank;
    trialStruct(t).MA.ach = MA_ach;
    
   T_fix = trialStruct(t).T_fix;
    
    
    T_tibF = T_fix.tib;
    tib_ach_pF = transformPoints(T_tibF, trialStruct(t).achilles.prox_pt_CT);
    T_calF = T_fix.cal;
    cal_ach_dF = transformPoints(T_calF, trialStruct(t).achilles.cal_pt_CT	);
    ach_vecF = tib_ach_pF - cal_ach_dF;
    
    
    
    MA_ankF = nan(ntot,3);
    MA_achF =  nan(ntot,3);
    for fr = fr_start:fr_mtp
        
          if fr > length(trialStruct(t).contact.tal)
            continue
        end
        % use the talus/ or tibia
        TC_ptF =  transformPoints(T_fix.tal(:,:,fr),trialStruct(t).contact.tal(fr).cc_CT);
        TC_vecF =   transformVectors(T_fix.tal(:,:,fr),trialStruct(t).contact.tal(fr).cN_CT);
        %     TC_ptF =  transformPoints(T_fix.tib(:,:,fr),trialStruct(s).contact.tib(fr).ContactCenter);
        %    	TC_vecF =   transformVectors(T_fix.tib(:,:,fr),-trialStruct(s).contact.tib(fr).cN_CT);
        
        
        % ankle moment arm
        [s1,s2] = closestPointsBtw2Lines(TC_ptF, trialStruct(t).helical(2).s_xrC(fr,:),...
            -TC_vecF, trialStruct(t).helical(2).n_xr(fr,:));
        MA_ankF(fr,:) = s2-s1;
        
        % achilles moment arm
        [s1,s2] = closestPointsBtw2Lines(cal_ach_dF(fr,:), trialStruct(t).helical(2).s_xrC(fr,:),...
            ach_vecF(fr,:), trialStruct(t).helical(2).n_xr(fr,:));
        MA_achF(fr,:) = s2 - s1;
    end
    
    
    GR_short{t} = norm3d(MA_ach(frsX(1):fr_mtp,:))./norm3d(MA_ank(frsX(1):fr_mtp,:));
    GR_fix{t} =  norm3d(MA_achF(frsX(1):fr_mtp,:))./norm3d(MA_ankF(frsX(1):fr_mtp,:));
    GR_short_N(t,:) = normaliseNaN(GR_short{t},1,101);
    GR_fix_N(t,:) = normaliseNaN(GR_fix{t},1,101);
    
    GR_fix_frs(t,:) = GR_fix{t}(frs_interest); % all the fixed Mech advantage at key frames
    GR_short_frs(t,:) = GR_short{t}(frs_interest); % all the moving Mech advantage at key frames
    
    GR(t,:) = [GR_fix{t}(end) GR_short{t}(end)];
    
    trialStruct(t).achilles.tib_ach_xrF = tib_ach_pF;
    trialStruct(t).achilles.cal_ach_xrF = cal_ach_dF;
    trialStruct(t).achilles.ach_vecF = ach_vecF;
    trialStruct(t).MA.ankF = MA_ankF;
    trialStruct(t).MA.achF = MA_achF;
    
    
    
    % ------------------velocity of the calc------------------------------
    
    vel_mov{t} = calculateVelocity(cal_ach_d,250);
    vel_fix{t} = calculateVelocity(cal_ach_dF,250);
    
    vel_movP{t} = calculateVelocity(tib_ach_p,250);
    vel_fixP{t} = calculateVelocity(tib_ach_pF,250);
    
    vel_cmp(t,:) = [norm(vel_fix{t}(fr_mtp-3,:)-vel_fixP{t}(fr_mtp-3,:)),norm(vel_mov{t}(fr_mtp-3,:) - vel_movP{t}(fr_mtp-3,:))];
    
    
    figure(600);hold on;
    plot(norm3d(vel_mov{t}(fr_start:fr_mtp,:)),'b'); hold on;
    plot(norm3d(vel_movP{t}(fr_start:fr_mtp,:)),'b:');
    
    plot(norm3d(vel_fix{t}(fr_start:fr_mtp,:)),'k'); hold on;
    plot(norm3d(vel_fixP{t}(fr_start:fr_mtp,:)),'k:');
    legend('moving cal','moving tib','fixed cal','fixed tib')
    title([trialStruct(t).subject ' ' trialStruct(t).trial ' ach vel'])
    
    % title([trialStruct(s).subject ' ' trialStruct(s).trial ' ach vel'])
    
    
    % save the moment arm lengths
    MA_ank_rigid(t,1:2) = norm3d([MA_ankF(fr_start,:);MA_ankF(fr_mtp-1,:)])';
    MA_ank_moving(t,1:2) = norm3d([MA_ank(fr_start,:);MA_ank(fr_mtp-1,:)])';
    MA_ach_rigid(t,1:2) = norm3d([MA_achF(fr_start,:);MA_achF(fr_mtp-1,:)])';
    MA_ach_moving(t,1:2) = norm3d([MA_ach(fr_start,:);MA_ach(fr_mtp-1,:)])';
    
    
    
    p = trialStruct(t).plot;
%     
%     figure(111)
%     subplot(2,1,1)
%     hold on;
%     
%     plot([1,2],[norm3d(MA_ankF(fr_mtp-1,:)),norm3d(MA_ank(fr_mtp-1,:))]-norm3d(MA_ank(fr_mtp-1,:)),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     title('Ankle moment arm fixed (left), moving(right)')
%     
%     subplot(2,1,2)
%     hold on;
%     plot([1,2],[norm3d(MA_achF(fr_mtp-1,:)),norm3d(MA_ach(fr_mtp-1,:))]-norm3d(MA_ach(fr_mtp-1,:)),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     title('Achilles moment arm fixed (left), moving(right)')
%     
%     figure(112)
%     subplot(2,1,1)
%     hold on;
%     
%     plot([1,2],[norm3d(MA_ankF(fr_mtp-1,:)),norm3d(MA_ank(fr_mtp-1,:))],'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     title('Ankle moment arm fixed (left), moving(right)')
%     subplot(2,1,2)
%     hold on;
%     plot([1,2],[norm3d(MA_achF(fr_mtp-1,:)),norm3d(MA_ach(fr_mtp-1,:))],'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     title('Achilles moment arm fixed (left), moving(right)')
%     %       plot(2,norm3d(MA_achF(end-1,:))-norm3d(MA_ach(end-1,:)),'o')
%     
%     %       hold on; plot(norm3d(MA_ach))
%     %       plot(norm3d(MA_ankF),':')
%     %       plot(norm3d(MA_achF),':')
%     %
%     
%     
%     
%     
%     figure;
%     subplot(2,1,1)
%     plot(norm3d(MA_ank))
%     hold on; plot(norm3d(MA_ach))
%     plot(norm3d(MA_ankF),':')
%     plot(norm3d(MA_achF),':')
%     legend('Ankle','Achilles','Ankle Fix','Achilles Fix')
%     title(trialStruct(t).subject)
%     
%     
%     subplot(2,1,2)
%     plot( GR_short{t})
%     hold on; plot( GR_fix{t},':')
%     legend('Shortening GR','Fixed GR')
%     title(trialStruct(t).subject)
%     %
    
end



moving_col = cmap_cb(3,:);%[0,0.5 0.8];
fixed_col = cmap_cb(1,:);
% 
% 
%     figure;
%     hold on;
%     for t = 1:2:ntrials-2
%         p = trialStruct(t).plot;
%         frsX = trialStruct(t).keyFrames.cropFrsXray;
%       
% %         plot(linspace(0,100,diff(frsX)+1 ) , trialStruct(t).kinemat.tib_lean.rot_fix(frsX(1):frsX(2))- trialStruct(t).kinemat.tib_lean.rot_mov(frsX(1):frsX(2)) ,'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%         plot( stancePts{t}, trialStruct(t).kinemat.tib_lean.rot_mov(frsX(1):frsX(2)) ,'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%         plot( stancePts{t}, trialStruct(t).kinemat.tib_lean.rot_fix(frsX(1):frsX(2)) ,'Color',p.col,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%         plot(stancePts{t}(fr_align_s(t)-frsX(1)+1),trialStruct(t).kinemat.tib_lean.rot_fix(fr_align_s(t)),'ko','markerfacecolor','k')
% %         plot(stancePts{t}(fr_align_s(t)-frsX(1)+1),trialStruct(t).kinemat.tib_lean.rot_fix(fr_align_s(t)),'ko','markerfacecolor','k')
%     
%         
%         time_early(t) = (fr_mtp_s(t) - fr_align_s(t))/trialStruct(t).marker_data.FrameRate;
%         perc_early(t) = time_early(t)/(frsX(2)-frsX(1) + 1) * trialStruct(t).marker_data.FrameRate;
%     end
%     ylabel('difference between tibial lean (-) = fixed is more leaned')
%     ylabel('tibial lean (-) = more leaned')
%     xlabel('% stance')
%     makeNicePlotsFunction
%     
%     figure; makeStatsDotPlot(repmat(1:2,1,6), time_early',{'walk','run'})
%     ylabel('Contact time lost with rigid arch [s]')
%     makeNicePlotsFunction
%     figure; makeStatsDotPlot(repmat(1:2,1,6), 100*perc_early',{'walk','run'})
%     ylabel('% of stance lost with rigid arch')
%     makeNicePlotsFunction
%     
% %     close all
%     for t = 2:2:ntrials-2
%          p = trialStruct(t).plot;
%     figure(131);hold on;
%     plot(trialStruct(t).kinemat.ROM_arch,  time_early(t),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.col);
%     
%     
%     xlabel('Arch ROM')
%     ylabel('contact time lost to take off')
% %     
%     figure(132); hold on;
%     
%     plot(trialStruct(t).kinemat.ROM_arch,  trialStruct(t).kinemat.tib_lean.rot_mov(fr_mtp_s(t)),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.col);
%     ylabel('tibial lean at max MTP')
%     xlabel('Arch ROM')
%     
%      figure(133); hold on;
%     
%     plot(trialStruct(t).kinemat.tib_lean.rot_mov(fr_mtp_s(t)),  res.timing.cont_time(t),'Color',p.col,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.col);
%     ylabel('contact time')
%     xlabel('tibial lean at max mtp')
%     end
%     
%     makeNicePlotsFunction
%     
    % the velocity plot over time
% figure; hold on;
% 
% fr_lose = ones(ntrials,1)*3;
% ach_mov_vel = nan(ntrials,101);
% ach_fix_vel = nan(ntrials,101);
% for t = 1:ntrials
%     
%     fr_prop = fr_tib_glob_align(t)-fr_start_s(t);
%     nfr = length(fr_start_s(t):fr_mtp_s(t));
%     npts_interp = round(100*(nfr-fr_lose(t))/nfr)+1;
%     
%     perc_prop(t) = round(fr_prop/nfr*100);
%     perc_end(t) = (nfr-2)/nfr*100;
%     
%     ach_mov = norm3d(vel_mov{t}(fr_start_s(t):fr_mtp_s(t),:)-vel_movP{t}(fr_start_s(t):fr_mtp_s(t),:));
%     ach_fix = norm3d(vel_fix{t}(fr_start_s(t):fr_mtp_s(t),:)-vel_fixP{t}(fr_start_s(t):fr_mtp_s(t),:));
%      
%     
%     p = trialStruct(t).plot;
% %     plot(linspace(0,100,length(ach_mov)),ach_mov,'color',moving_col,'marker','.');hold on;
% %     plot(linspace(0,100,length(ach_fix)),ach_fix,'color',fixed_col,'marker','^');hold on;
%     plot(linspace(0,100,length(ach_fix)),ach_fix-ach_mov,'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);hold on;
%     
%     ach_mov_vel(t,1:npts_interp) = normaliseNaN(adaptiveLowPassButterworth(ach_mov(1:nfr-fr_lose(t))',[14,20],250)',1,npts_interp);%normaliseNaN(ach_mov(1:nfr-fr_lose(s)),1,npts_interp);%
%     ach_fix_vel(t,1:npts_interp) = normaliseNaN(adaptiveLowPassButterworth(ach_fix(1:nfr-fr_lose(t))',[14,20],250)',1,npts_interp);%normaliseNaN(ach_fix(1:nfr-fr_lose(s)),1,npts_interp);%
%     ach_vel_save(t,1:2) = [ ach_mov(nfr-fr_lose(t)),ach_fix(nfr-fr_lose(t))];
% end
% xlabel('% propulsion (between max arch dorsiflexion +max mtp)')
% ylabel('rigid achilles vel - moving achilles vel')
% makeNicePlotsFunction
% grid on;
% 
% figure;
% h = PrettyStdDevGraphs(0:100,mean(ach_mov_vel),std(ach_mov_vel),moving_col,1); hold on;
% 
% h1 = PrettyStdDevGraphs(0:100,mean(ach_fix_vel),std(ach_fix_vel),fixed_col,1);hold on;
% h(2).LineWidth = 5;
% h1(2).LineWidth = 5;
% for t = 1:2:ntrials
%     
%     p = trialStruct(t).plot;
%     h2(t) =  plot(0:100,ach_fix_vel(t,:),'color',fixed_col,'Linewidth',2,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     plot(0:100,ach_mov_vel(t,:),'color',moving_col,'Linewidth',2,'Linewidth',2,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
% end
% legend('moving arch','rigid arch ')
% xlabel('% of Propulsion')
% ylabel('Achilles shortening velocity [mm/s]')
% % legend([h(2) h1(2) h2(7:10 )], {'Moving', 'Rigid',  'RFS - Barefoot','FFS - Barefoot','RFS - Shod','FFS - Shod' })
% makeNicePlotsFunction
% 
% 
% figure;
% h = PrettyStdDevGraphs(0:100,mean(ach_mov_vel),std(ach_mov_vel),moving_col,1); hold on;
% 
% h1 = PrettyStdDevGraphs(0:100,mean(ach_fix_vel),std(ach_fix_vel),fixed_col,1);hold on;
% h(2).LineWidth = 5;
% h1(2).LineWidth = 5;
% for t = 2:2:ntrials
%     
%     p = trialStruct(t).plot;
%     h2(t) =  plot(0:100,ach_fix_vel(t,:),'color',fixed_col,'Linewidth',2,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
%     plot(0:100,ach_mov_vel(t,:),'color',moving_col,'Linewidth',2,'Linewidth',2,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
% end
% legend('moving arch','rigid arch ')
% xlabel('% of Propulsion')
% ylabel('Achilles shortening velocity [mm/s]')
% % legend([h(2) h1(2) h2(7:10 )], {'Moving', 'Rigid',  'RFS - Barefoot','FFS - Barefoot','RFS - Shod','FFS - Shod' })
% makeNicePlotsFunction
% 
% % loook at subject specifci velocity difference
% 
% figure;
% subplot(2,1,1)
% hold on
% 
% for t = 1:2:ntrials
%     
%     p = trialStruct(t).plot;
%     h2(t) =  plot(0:100,ach_fix_vel(t,:)-ach_mov_vel(t,:),'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
% end
% subplot(2,1,2)
% hold on;
% for t = 2:2:ntrials
%     
%     p = trialStruct(t).plot;
%     h2(t) =  plot(0:100,ach_fix_vel(t,:)-ach_mov_vel(t,:),'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle',p.line,'MarkerFaceColor',p.marker_fill);
% end
% 
% xlabel('% of Propulsion')
% ylabel('Achilles shortening velocity difference [mm/s]')
% % legend([h(2) h1(2) h2(7:10 )], {'Moving', 'Rigid',  'RFS - Barefoot','FFS - Barefoot','RFS - Shod','FFS - Shod' })
% makeNicePlotsFunction

% 
% 
% % look at the ACHILLEs tendon velocity at different frames
% ct = 1;
% figure;
% subplot(2,1,1); hold on;title('Ach velocity walk')
% for t = [1:2:ntrials,2:2:ntrials]
%     if ct == 5
%         subplot(2,1,2); hold on;
%         title('Ach velocity run')
%     end
%     [~,Is] = sort(frs_interest_s(t,:));
%     frs_sort = frs_interest_s(t,Is);
%     frsX = trialStruct(t).keyFrames.cropFrsXray;
%     fr_int = frs_sort + frsX(1);
%     vel_fixAch = norm3d(vel_fix{t} - vel_fixP{t});
%     vel_movAch = norm3d(vel_mov{t} - vel_movP{t});
%     
%     p = trialStruct(t).plot;
%  plot(vel_fixAch(fr_int) ,'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle','-','MarkerFaceColor',p.marker_fill);
% plot(vel_movAch(fr_int) ,'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle',':','MarkerFaceColor',p.marker_fill);
%     ct = ct + 1;
%     ylim([0 800])
% end
% 
% % 
% % Look at Mechanical advantage over stance
% % close all
% ct = 1;
% figure;
% subplot(2,1,1); hold on;title('MA diff walk')
% for t = [1:2:ntrials,2:2:ntrials]
%     if ct == ntrials/2+1
%         subplot(2,1,2); hold on;
%         title('MA diff run')
%     end
%     [~,Is] = sort(frs_interest_s(t,:));
%     p = trialStruct(t).plot;
%  plot(stancePts{t}(frs_interest_s(t,Is)),GR_short_frs(t,Is)-GR_fix_frs(t,Is),'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle','-','MarkerFaceColor',p.marker_fill);
% %  plot(frs_interest_s(t,:),GR_short_frs(t,:),'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle',':','MarkerFaceColor',p.marker_fill);
%     ct = ct + 1;
% end
% 
% ct = 1;
% figure;
% subplot(2,1,1); hold on;title('MA @frames of interest walk')
% for t = [1:2:ntrials,2:2:ntrials]
%     if ct == ntrials/2+1
%         subplot(2,1,2); hold on;
%         title('MA @frames of interest run')
%     end
%     [~,Is] = sort(frs_interest_s(t,:));
%     p = trialStruct(t).plot;
%  plot(stancePts{t}(frs_interest_s(t,Is)),GR_short_frs(t,Is),'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle',':','MarkerFaceColor',p.marker_fill);
%  plot(stancePts{t}(frs_interest_s(t,Is)),GR_fix_frs(t,Is),'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle','-','MarkerFaceColor',p.marker_fill);
%     ct = ct + 1;
% end
% 
% ct = 1;
% figure;
% subplot(2,1,1); hold on;title('MA shortening arch for all frames walk')
% for t = [1:2:ntrials,2:2:ntrials]
%     if ct == ntrials/2+1
%         subplot(2,1,2); hold on;
%         title('MA shortening arch for all frames run')
%     end
%     frsX = trialStruct(t).keyFrames.cropFrsXray;
%     
%     p = trialStruct(t).plot;
%  plot(stancePts{t}(1:fr_mtp_s(t)-frsX(1)+1),GR_short{t},'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle','-','MarkerFaceColor',p.marker_fill);
%     ct = ct + 1;
% end

ct = 1;
figure;
subplot(2,1,1); hold on;title('MA difference between moving/rigid for all frames walk')
for t = [1:2:ntrials-2,2:2:ntrials-2]
    if ct >6
        subplot(2,1,2); hold on;
        title('MA difference between moving/rigid for all frames run')
    end
   
    frsX = trialStruct(t).keyFrames.cropFrsXray;
    
    p = trialStruct(t).plot;
 plot(stancePts{t}(1:fr_mtp_s(t)-frsX(1)+1),GR_short{t}-GR_fix{t},'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle',':','MarkerFaceColor',p.marker_fill);
plot([30 100],[0 0],'k')
ct = ct+1;
%     plot(stancePts{t}(1:fr_mtp_s(t)-frsX(1)+1),GR_fix{t},'color',p.col,'Linewidth',2,'marker',p.marker_type,'linestyle','-','MarkerFaceColor',p.marker_fill);
end
makeNicePlotsFunction




%% calculate all the parameters to report:


% look at talar height
makeGroupedStatsDotPlot([tal_ht_fix(:,2:3)-tal_ht_act(:,2:3)],repmat(1:2,1,ntrials/2),{'walk','run'},{'anterior','superior'})

ylabel('Fixed talus position - actual talus position [mm]')
makeNicePlotsFunction
hold on;

plot([0.25,0.5],[tal_ht_fix(1:2:end,2)-tal_ht_act(1:2:end,2),tal_ht_fix(2:2:end,2)-tal_ht_act(2:2:end,2)]','k')
plot([1.25,1.5],[tal_ht_fix(1:2:end,3)-tal_ht_act(1:2:end,3),tal_ht_fix(2:2:end,3)-tal_ht_act(2:2:end,3)]','k')

% 
% makeGroupedStatsDotPlot([tal_ht_fix(:,2:3);tal_ht_act(:,2:3)],[repmat([1,3],1,ntrials/2),repmat([2,4],1,ntrials/2)],{'walk-fixed','walk','run-fixed','run'},{'anterior','superior'})
% ylabel('Distance that the fixed talus is higher and more anterior than the moving arch')
% makeNicePlotsFunction

perc_red_ank = (MA_ank_rigid(:,2)-MA_ank_moving(:,2))./MA_ank_moving(:,2)*100
[mean(perc_red_ank) std(perc_red_ank)]

perc_red_ach = (MA_ach_rigid(:,2)-MA_ach_moving(:,2))./MA_ach_moving(:,2)*100
[mean(perc_red_ach) std(perc_red_ach)]


ank_len = (MA_ank_rigid(:,2)-MA_ank_moving(:,2));
[mean(ank_len) std(ank_len)]
ach_len = (MA_ach_rigid(:,2)-MA_ach_moving(:,2));
[mean(ach_len) std(ach_len)]

perc_gr_red =  (GR(:,1)- GR(:,2))./ GR(:,2)*100;
[mean(perc_gr_red) std(perc_gr_red)]
[mean(GR(:,1)- GR(:,2))  std(GR(:,1)- GR(:,2))]

ttest(MA_rigid(:,2)-MA_moving(:,2))
[p,h,a] = ttest(MA_rigid(:,2)-MA_moving(:,2))

%STATS

[p] = signrank(ach_vel_save(1:7,1),ach_vel_save(1:7,2));
[p] = signrank(GR(1:7,1),GR(1:7,2));
% superior tal
[p] = signrank(tal_ht_act(1:7,3),tal_ht_fix(1:7,3));
[p] = signrank(MA_ank_moving(1:7,2),MA_ank_rigid(1:7,2));
[p] = signrank(MA_ach_moving(1:7,2),MA_ach_rigid(1:7,2));


mean(ach_mov_vel(:,end))
std(ach_mov_vel(:,end))
mean(ach_fix_vel(:,end))
std(ach_fix_vel(:,end))

% achilles velocity at simulated propulsion
mean(ach_mov_vel(:,end))
std(ach_mov_vel(:,end))
for t = 1:nsubj
    vf_1(t) = ach_fix_vel(t,perc_prop(t));
end
mean(vf_1)
std(vf_1)


% [~,p] = ttest(
for t = 1:ntrials
    frsX = trialStruct(t).keyFrames.cropFrsXray;
    cont_time(t) = diff(frsX);%(fr_mtp_s'-fr_start_s');

end
fr_diff_beg = (fr_tib_glob_align'-fr_start_s');
early_t = cont_time-fr_diff_beg;
perc = early_t./cont_time*100;
figure; hold on;
plot(ones(nsubj,1),cont_time(1:2:end)/125,'o')
plot(ones(nsubj,1)*2,cont_time(2:2:end)/250,'x')
xlim([0 3])
ylabel('Contact time')

%% calculate the contact and make a wrist viz for the tibia moved to the same location


for t = 1:ntrials
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    Tm = trialStruct(t).Tm_mm;
    Tfix = trialStruct(t).T_fix;
    nfr = size(Tm.tib,3);
    bone_listM = fields(Tm);
    nBones = length(bone_listM);
    
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    
    fr_tg = trialStruct(t).keyFrames.fr_tg ;
    
    
    dfield = subjStruct(s_ind).bones.tib.dfield;
    pts = subjStruct(s_ind).tal_dome.pts;
    conns = subjStruct(s_ind).tal_dome.cns(:,1:3);
    thresh = surface_dist(s_ind);
    viz = 0;
    
    
    RT_ref =   trialStruct(t).kinemat.T_tib_fix_act_mtp;
    RT_pts = trialStruct(t).T_fix.tal(:,:,fr_mtp);
    
    
    out = JointContact_Dfield(dfield,pts,conns,thresh,1,0,RT_ref,RT_pts);
    trialStruct(t).contact.tal_sim = out;
    
    % get the normal at the contact centroid surface - weight all the
    % normals by their distance
    cc_xr = transformPoints(RT_pts,out.ContactCenter);
    cN_CT = mean(out.Patch.wtdist.* out.Patch.norm);
    cN_xr = unit(transformPoints(RT_pts, cN_CT));
    
    trialStruct(t).contact.tal_sim.cc_xr = cc_xr;
    trialStruct(t).contact.tal_sim.cN_CT= cN_CT;
    trialStruct(t).contact.tal_sim.cN_xr = cN_xr;
    
    
    
    
    anim_dir = fullfile(dir_analy, 'Animation','TalusContact',[trialStruct(t).subject '_' trialStruct(t).trial],filesep);
    if ~exist(anim_dir,'dir')
        mkdir(anim_dir)
    end
    %------ make the simulated tibia/rigid foot iv file -------------------
    
    pat = trialStruct(t).contact.tal_sim.Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSim.iv'))
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSim.iv'),RT_pts(1:3,1:3),RT_pts(1:3,4)')];
    
    for bn = 1:nBones
        
        Tf =  T_fix.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            Tf = RT_ref;
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
    end
    
    fid = fopen(fullfile(anim_dir,'TalusContactSIM_Rigid.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    
    %------ make the actual tibia/rigid foot iv file -------------------
    pat = trialStruct(t).contact.tal(fr_mtp).Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSimNormal.iv'))
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),RT_pts(1:3,1:3),RT_pts(1:3,4)')];
    
    for bn = 1:nBones
        
        Tf =  Tfix.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
    end
    
    fid = fopen(fullfile(anim_dir,'TalusContactRigid.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    
    
    %------ make the actual tibia/rigid foot iv file -------------------
    pat = trialStruct(t).contact.tal(fr_mtp).Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSimNormal.iv'))
    
    T_tal = Tm.tal(:,:,fr_mtp);
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),T_tal(1:3,1:3),T_tal(1:3,4)')];
    
    for bn = 1:nBones
        
        Tf =  Tm.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
    end
    
    fid = fopen(fullfile(anim_dir,'TalusContactMoving.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    
    %------ make the actual tibia/rigid foot iv fileat simulated pushoff -------------------
    %     pat = trialStruct(s).contact.tal(fr_mtp).Patch;
    %     patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSimNormal.iv'))
    
    %     T_tal = Tx.tal(:,:,fr_mtp);
    
    ivstring = createInventorHeader();
    %     ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),T_tal(1:3,1:3),T_tal(1:3,4)')];
    
    for bn = 1:nBones
        
        Tf =  Tm.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        end
        
        Tf =  Tfix.(bone_listM{bn})(:,:,fr_tg); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            Tf =  Tfix.tib_plant(:,:,fr_tg);
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
        
        
        
    end
    
    fid = fopen(fullfile(anim_dir,'Tibia at simulated pushoff.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
end


%% Visualize the contact at any frame - make wrist viz

for t = 1:ntrials
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    Tm = trialStruct(t).Tm_mm;
    Tfix = trialStruct(t).T_fix;
    nfr = size(Tm.tib,3);
    bone_listM = fields(Tm);
    nBones = length(bone_listM);
    
    fr_look = trialStruct(t).keyFrames.max_dors+1; % the frame to look at
    fr_name = 'Max_Ankle_Dors'; % name this frame
%     fr_tg = trialStruct(t).keyFrames.fr_tg ;
    
    anim_dir = fullfile(dir_analy, 'Animation',['TalusContact_' fr_name],[trialStruct(t).subject '_' trialStruct(t).trial],filesep);
    if ~exist(anim_dir,'dir')
        mkdir(anim_dir)
    end

    
    
    %------ make the actual tibia/rigid foot iv file -------------------
    pat = trialStruct(t).contact.tal(fr_look).Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSimNormal.iv'))
    
    T_tal = Tm.tal(:,:,fr_look);
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),T_tal(1:3,1:3),T_tal(1:3,4)')];
    
    for bn = 1:nBones
        
        Tf =  Tm.(bone_listM{bn})(:,:,fr_look); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
            ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
    end
    
    fid = fopen(fullfile(anim_dir,['TalusContactMoving' fr_name '.iv']),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    
    
end

%% Make a visualization of the difference in position for the rigid and non rigid foot + outline graph


for t = 1:ntrials
    
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    
    anim_dir = fullfile(dir_analy, 'Animation','Final_frame_talar_Centroid',trialStruct(t).trial,filesep);
    if ~exist(anim_dir,'dir')
        mkdir(anim_dir)
    end
    Tm = trialStruct(t).Tm_mm;
    Tfix = trialStruct(t).T_fix;
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    
    
    nfr = size(Tm.tib,3);
    bone_listM = fields(Tm);
    nBones = length(bone_listM);
    
    
    ivstring = createInventorHeader();
    
    for bn = 1:nBones % moving arch
        
        Tmm =  Tm.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tmm(1:3,1:3),Tmm(1:3,4)',[0.7 0.7 0.7],0.8)];
    end
    for bn = 1:nBones % when tibia is aligned, the rigid foot
        Tff =  Tfix.(bone_listM{bn})(:,:,fr_tib_glob_align(t)); % T of this frame fixed
        % fixed (faded)
        ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tff(1:3,1:3),Tff(1:3,4)',[0.7 0.7 0.7],0.5)];
        
    end
   for bn = 1:nBones % at MTP max, the rigid foot
        Tff =  Tfix.(bone_listM{bn})(:,:,fr_mtp); % T of this frame fixed
        % fixed (faded)
        ivstring = [ivstring createInventorLink(subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file, Tff(1:3,1:3),Tff(1:3,4)',[0.7 0.7 0.7],0.8)];
        
    end
    
    % add the talar centroid
    
    ivstring = [ivstring createInventorSphere(trialStruct(t).kinemat.tal_cent_short(fr_mtp,:) ,3,[0.25 0.7 0.9],0)];
    ivstring = [ivstring createInventorSphere(trialStruct(t).kinemat.tal_cent_fix(fr_tib_glob_align(t),:) ,3,[0.5 0.7 0.9],0)];
    ivstring = [ivstring createInventorSphere(trialStruct(t).kinemat.tal_cent_fix(fr_mtp,:) ,3,[0 0.7 0.9],0)];
    
    fid = fopen(fullfile(anim_dir,'TalarCentroid.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    fprintf('Wrote file : %s \n',fullfile(anim_dir,'TalarCentroid.iv'))
end

%% make the outlines of the bones - when the fixed and moving tibia are aligned
close all
for t_pick = 2%1:12
% t_pick = 2;
for t = t_pick
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    figure; hold on;
    %     pln_vec =  T_ACS_ph1{s}(1:3,1)';
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    
    Tm = trialStruct(t).Tm_mm;
    Tfix = trialStruct(t).T_fix;
    
    T_ph1 = rotateCoordSys(T_ACS_ph1{t},0,3);
    %     warning('NOT THE ACTUAL PH1 ACS')
    
    bone_listM = fields(Tm);
    ind_fib = contains(bone_listM,'fib');
    bone_listM(ind_fib) = [];
    nBones = length(bone_listM);
    pln_vec = [1,0,0];
    for bn = 1:nBones
        pts = transformPoints(invTranspose(T_ph1)* Tm.(bone_listM{bn})(:,:,fr_mtp),subjStruct(s_ind).bones.(bone_listM{bn}).pts);
        ptsF = transformPoints(invTranspose(T_ph1)* Tfix.(bone_listM{bn})(:,:,fr_tib_glob_align(t)),subjStruct(s_ind).bones.(bone_listM{bn}).pts);
        %
        npts = size(pts,1);
        pts_plane = [];
        for p = 1:npts
            pts_plane(p,:) = closestPointonPlanealongVector(pts(p,:),pln_vec, T_ACS_ph1{t}(1:3,4)',pln_vec);
        end
        k = boundary(pts_plane(:,2),pts_plane(:,3));
        outline.(bone_listM{bn}) = pts_plane(k,:);
        
        
        pts_plane = [];
        for p = 1:npts
            pts_plane(p,:) = closestPointonPlanealongVector(ptsF(p,:),pln_vec, T_ACS_ph1{t}(1:3,4)',pln_vec);
        end
        k = boundary(pts_plane(:,2),pts_plane(:,3));
        
        outlineF.(bone_listM{bn}) = pts_plane(k,:);
        
        h = fill(outline.(bone_listM{bn})(:,2),outline.(bone_listM{bn})(:,3),[.7 .7 .7],'facealpha',0.4,'Edgecolor', [.7 .7 .7],'EdgeAlpha',0.8);%'LineStyle','none');
        fill(outlineF.(bone_listM{bn})(:,2),outlineF.(bone_listM{bn})(:,3),[.7 .7 .7],'facealpha',0.8,'Edgecolor', [.7 .7 .7],'EdgeAlpha',0.8);%,'LineStyle','none');
        %
        %
        %
        %         h = plot(tal_ht_act(s,2),tal_ht_act(s,3),'kx');
        %         h.MarkerSize = 10;
        %         %         h = plot(tal_ht_fix(s,2),tal_ht_fix(s,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        %         h = plot(tal_ht_fix_tg(s,2),tal_ht_fix_tg(s,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        %         h.MarkerSize = 4;
         h = plot(tal_ht_act(t_pick,2),tal_ht_act(t_pick,3) ,'s','MarkerFaceColor','none','MarkerEdgeColor','k');
         h = plot(tal_ht_fix_tg(t_pick,2),tal_ht_fix_tg(t_pick,3) ,'s','MarkerFaceColor','none','MarkerEdgeColor','r');
          
        axis equal
    end
end

for t = [1:t_pick-1,t_pick+1:ntrials-2]
    
        p = trialStruct(t).plot;
    
        %          tal_ht_fix_tg(s,:)
        %
        %       tal_diff_tg(
        %        quiver( tal_ht_act(s,2),tal_ht_act(s,3), abs(tal_diff(su,2)),abs(tal_diff(su,3)),'k')
        if strcmp(trialStruct(t).subject,'SOL001B')
            h = plot(tal_ht_act(t,2)- abs(tal_diff_tg(su,2)),tal_ht_act(t,3)-abs(tal_diff_tg(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','r');% if looking at sim take off
            %              h = plot(tal_ht_act(s,2)+ abs(tal_diff(su,2)),tal_ht_act(s,3)+abs(tal_diff(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','r');
            h.MarkerSize =6;
        else
%             h = plot(tal_ht_act(t,2)- abs(tal_diff_tg(su,2)),tal_ht_act(t,3)-abs(tal_diff_tg(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','k');
%             h = plot(tal_ht_act(t,2),tal_ht_act(t,3),'x','MarkerFaceColor','none','MarkerEdgeColor','k');
%             h = plot(tal_ht_fix_tg(t,2),tal_ht_fix_tg(t,3),'x','MarkerFaceColor','none','MarkerEdgeColor','k');
            h = plot(tal_ht_act(t_pick,2) - tal_diff_tg(t,2),tal_ht_act(t_pick,3) - tal_diff_tg(t,3),p.marker_type,'MarkerFaceColor','none','MarkerEdgeColor','r');
            %        h = plot(tal_ht_act(s,2)+ abs(tal_diff(su,2)),tal_ht_act(s,3)+abs(tal_diff(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','k');
            h.MarkerSize =6;
        end
        
%     end
end

makeNicePlotsFunction
end

%% make the outlines of the bones - at max mtp
close all
for t_pick = 2%1:12
% t_pick = 2;
for t = t_pick
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    figure; hold on;
    %     pln_vec =  T_ACS_ph1{s}(1:3,1)';
    fr_mtp = trialStruct(t).keyFrames.max_mtp;
    
    Tm = trialStruct(t).Tm_mm;
    Tfix = trialStruct(t).T_fix;
    
    T_ph1 = rotateCoordSys(T_ACS_ph1{t},0,3);
    %     warning('NOT THE ACTUAL PH1 ACS')
    
    bone_listM = fields(Tm);
    ind_fib = contains(bone_listM,'fib');
    bone_listM(ind_fib) = [];
    nBones = length(bone_listM);
    pln_vec = [1,0,0];
    for bn = 1:nBones
        pts = transformPoints(invTranspose(T_ph1)* Tm.(bone_listM{bn})(:,:,fr_mtp),subjStruct(s_ind).bones.(bone_listM{bn}).pts);
        ptsF = transformPoints(invTranspose(T_ph1)* Tfix.(bone_listM{bn})(:,:,fr_mtp),subjStruct(s_ind).bones.(bone_listM{bn}).pts);
        %
        npts = size(pts,1);
        pts_plane = [];
        for p = 1:npts
            pts_plane(p,:) = closestPointonPlanealongVector(pts(p,:),pln_vec, T_ACS_ph1{t}(1:3,4)',pln_vec);
        end
        k = boundary(pts_plane(:,2),pts_plane(:,3));
        outline.(bone_listM{bn}) = pts_plane(k,:);
        
        
        pts_plane = [];
        for p = 1:npts
            pts_plane(p,:) = closestPointonPlanealongVector(ptsF(p,:),pln_vec, T_ACS_ph1{t}(1:3,4)',pln_vec);
        end
        k = boundary(pts_plane(:,2),pts_plane(:,3));
        
        outlineF.(bone_listM{bn}) = pts_plane(k,:);
        
        h = fill(outline.(bone_listM{bn})(:,2),outline.(bone_listM{bn})(:,3),[.7 .7 .7],'facealpha',0.4,'Edgecolor', [.7 .7 .7],'EdgeAlpha',0.8);%'LineStyle','none');
        fill(outlineF.(bone_listM{bn})(:,2),outlineF.(bone_listM{bn})(:,3),[.7 .7 .7],'facealpha',0.8,'Edgecolor', [.7 .7 .7],'EdgeAlpha',0.8);%,'LineStyle','none');
        %
        %
        %
        %         h = plot(tal_ht_act(s,2),tal_ht_act(s,3),'kx');
        %         h.MarkerSize = 10;
        %         %         h = plot(tal_ht_fix(s,2),tal_ht_fix(s,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        %         h = plot(tal_ht_fix_tg(s,2),tal_ht_fix_tg(s,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        %         h.MarkerSize = 4;
         h = plot(tal_ht_act(t_pick,2),tal_ht_act(t_pick,3) ,'s','MarkerFaceColor','none','MarkerEdgeColor','k');
         h = plot(tal_ht_fix(t_pick,2),tal_ht_fix(t_pick,3) ,'s','MarkerFaceColor','none','MarkerEdgeColor','b');
          
        axis equal
    end
end

for t = [1:t_pick-1,t_pick+1:ntrials-2]
    
        p = trialStruct(t).plot;
    
        %          tal_ht_fix_tg(s,:)
        %
        %       tal_diff_tg(
        %        quiver( tal_ht_act(s,2),tal_ht_act(s,3), abs(tal_diff(su,2)),abs(tal_diff(su,3)),'k')
        if strcmp(trialStruct(t).subject,'SOL001B')
            h = plot(tal_ht_act(t,2)- abs(tal_diff(su,2)),tal_ht_act(t,3)-abs(tal_diff(su,3)),p.marker_type,'MarkerFaceColor','none','MarkerEdgeColor','b');% if looking at sim take off
            %              h = plot(tal_ht_act(s,2)+ abs(tal_diff(su,2)),tal_ht_act(s,3)+abs(tal_diff(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','r');
            h.MarkerSize =6;
        else
%             h = plot(tal_ht_act(t,2)- abs(tal_diff_tg(su,2)),tal_ht_act(t,3)-abs(tal_diff_tg(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','k');
%             h = plot(tal_ht_act(t,2),tal_ht_act(t,3),'x','MarkerFaceColor','none','MarkerEdgeColor','k');
%             h = plot(tal_ht_fix_tg(t,2),tal_ht_fix_tg(t,3),'x','MarkerFaceColor','none','MarkerEdgeColor','k');
            h = plot(tal_ht_fix(t_pick,2) + tal_diff(t,2),tal_ht_fix(t_pick,3) + tal_diff(t,3),p.marker_type,'MarkerFaceColor','none','MarkerEdgeColor','b');
            %        h = plot(tal_ht_act(s,2)+ abs(tal_diff(su,2)),tal_ht_act(s,3)+abs(tal_diff(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','k');
            h.MarkerSize =6;
        end
        
%     end
end

makeNicePlotsFunction
end
%% Visualize the vectors in Wrist viz


cmap = colormap('parula');

c_ind = round(linspace(1,64,10));
cmap_crop = cmap(c_ind,:);


% clearvars -except trialStruct subjStruct dir_analy cmap_crop cmap nsubj ntrials tal_ht_act tal_ht_fix T_ACS_ph1 plot_cols dir_ref fr_start_s fr_mtp_s demo


for t =  1:ntrials
    
    % get the subject's index in subjectStruct
    s_ind = findInStruct(subjStruct,'subject',trialStruct(t).subject);
    
    fr_start = fr_start_s(t);
    fr_mtp = fr_mtp_s(t);
    
%     fr_start = trialStruct(t).keyFrames.max_dors; % frame of max dorsiflexion, choose as reference frame
    %     max_mtp = trialStruct(s).keyFrames.last_fr_tracked;
%     max_mtp = trialStruct(t).keyFrames.max_mtp;
    
    T = trialStruct(t).Tm_mm;
    T_mt1_fix = T.mt1(:,:,fr_start);
    T_mt1 = T.mt1;
    nfr = size(T_mt1,3);
    bone_listM = fields(T);
    nBones = length(bone_listM);
    
    TanimF = {}; TanimM = {};
    bone_listF = {};
    
    T_fix = trialStruct(t).T_fix;
    
    for bn = 1:nBones
        %         T_bone_fix = T.(bone_listM{bn})(:,:,fr_md);
        ind_fr = 1;
        for fr = 1:nfr
            % give all the fixed transforms with the mt1 relative to the selected max dors position
            %             T_fix.(bone_listM{bn}) = trialStruct(s).T_fix.(bone_listM{bn});% T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * T_bone_fix;
            
            % write this for the animation as well; write ones where there
            % are nans
            if fr >= fr_start && fr <= fr_mtp
                if any(isnan( T_fix.(bone_listM{bn})(:,:,fr)))
                    
                    TanimF{bn}(:,:,ind_fr) = ones(4,4);
                else
                    %                     if strcmp(bone_listM{bn},'tib')
                    %
                    %                          TanimF{bn}(:,:,ind_fr) = T_fix.tib_plant(:,:,fr);
                    %                     else
                    TanimF{bn}(:,:,ind_fr) =  T_fix.(bone_listM{bn})(:,:,fr);
                    %                     end
                end
                
                if any(isnan( T.(bone_listM{bn})(:,:,fr)))
                    
                    TanimM{bn}(:,:,ind_fr) = ones(4,4);
                else
                    TanimM{bn}(:,:,ind_fr) =  T.(bone_listM{bn})(:,:,fr);
                end
                ind_fr = ind_fr+1;
            end
        end
        
        bone_listF{bn} = [bone_listM{bn} 'RIGID'];
    end
    
    % find the average of the achilles plane and project that into x-ray
    % space
    
    mean_ach_xr = trialStruct(t).achilles.cal_ach_xr;

    anim_dir = fullfile(dir_analy, 'Animation','RigidFoot_HelicalAxis',[trialStruct(t).subject '_' trialStruct(t).trial],filesep);
    if ~exist(anim_dir,'dir')
        mkdir(anim_dir)
    end
    
    % write the iv files with the helical axes and the points
    arrow_style = 'arrows_F%i_N%i.iv';
    % create the ini files---------------
    arrowDir = fullfile(anim_dir,'helicalDir');
    if ~exist(arrowDir,'dir')
        mkdir(arrowDir);
    end
    
    ind_fr = 1;
    for fr = fr_start:fr_mtp
        
        if fr > length(trialStruct(t).contact.tal)
            continue
        end
        
        pat = trialStruct(t).contact.tal(fr).Patch;
        patch2iv_vrml2_color(transformPoints(T.tal(:,:,fr),pat.pts),pat.conns,pat.color.dist{2},0.5,fullfile(arrowDir,sprintf(arrow_style,ind_fr,2)))
        
        ivstring = createInventorHeader();
        
        cc_xr = trialStruct(t).contact.tal(fr).cc_xr ;
        cN_CT = trialStruct(t).contact.tal(fr).cN_CT ;
        cN_xr = -trialStruct(t).contact.tal(fr).cN_xr ;
        
        % %
        %
        %          pat = trialStruct(s).contact.tib(fr).Patch;
        %         patch2iv_vrml2_color(transformPoints(T.tib(:,:,fr),pat.pts),pat.conns,pat.color.dist{2},0.5,fullfile(arrowDir,sprintf(arrow_style,ind_fr,2)))
        %
        %          ivstring = createInventorHeader();
        %
        %           cc_xr = trialStruct(s).contact.tib(fr).cc_xr ;
        %           cN_CT = trialStruct(s).contact.tib(fr).cN_CT ;
        %           cN_xr = trialStruct(s).contact.tib(fr).cN_xr ;
        %
        
        
        ivstring = [ivstring createInventorSphere(cc_xr ,1,[1 1 1],0.2)];
        ivstring = [ivstring createInventorArrow(cc_xr,cN_xr,150,2,[1 1 1],0.2)];
        
        
        %          if fr == max_mtp
        %              ivstring = [ivstring createInventorSphere(trialStruct(s).kinemat.tal_cent_short(max_mtp,:) ,3,[0.4 0.4 0.7],0.2)];
        %              ivstring = [ivstring createInventorSphere(trialStruct(s).kinemat.tal_cent_fix(max_mtp,:) ,3,[0.4 0.4 0.7],0.2)];
        %              ivstring = [ivstring createInventorCoordinateSystem( T_ACS_ph1{s},100,3 )];
        %
        %          end
        %           if fr == trialStruct(s).keyFrames.fr_tg
        %               ivFile = trialStruct(s).bones.tib.metadata.orig_file;
        %             ivstring = [ivstring createInventorLink(ivFile,T_fix.tib_plant(1:3,1:3,fr),T_fix.tib_plant(1:3,4,fr),[1 0 0],0.4)];
        %           end
        
        
        % METATARSAL AXIS
        for bp = 2
            if any(isnan([trialStruct(t).helical(bp).s_xrC(fr,:),trialStruct(t).helical(bp).n_xr(fr,:)]))
                continue
            else
                if bp == 2
                    cc = [0 0 1];
                else
                    cc = [0 0 1];
                end
                ivstring = [ivstring createInventorArrow(trialStruct(t).helical(bp).s_xrC(fr,:),trialStruct(t).helical(bp).n_xr(fr,:),200,2,cc,0.2)];
                % centroid of the bone that the HA is trying to be closest to
                ivstring = [ivstring createInventorSphere(trialStruct(t).helical(bp).cent_xr(fr,:),5,cmap(40,:),0.2)];
                
            end
            
        end
        
        
        % MOMENT ARM ARROWS
        MA_achF = trialStruct(t).MA.achF(fr,:);
        MA_ankF = trialStruct(t).MA.ankF(fr,:);
        MA_ach = trialStruct(t).MA.ach(fr,:);
        MA_ank = trialStruct(t).MA.ank(fr,:);
        
        if any(isnan([ MA_achF  MA_ach  MA_ankF  MA_ank]))
            continue
        else
            
            ivstring = [ivstring createInventorArrow(cc_xr,MA_ank,norm(MA_ank),2,cmap(60,:),0.2)];
            %               ivstring = [ivstring createInventorArrow(cc_xr,MA_ankF,norm(MA_ank),2,cmap(60,:),0.2)];
            ivstring = [ivstring createInventorArrow(mean_ach_xr(fr,:),MA_ach,norm(MA_ach),2,cmap(60,:),0.2)];
            %               ivstring = [ivstring createInventorArrow(cc_xr,MA_ank,norm(MA_ank),2,cmap(60,:),0.2)];
            
            
        end
        
        % Achillles arrow
        if any(isnan(mean_ach_xr(fr,:)))
            continue
        else
            %             ivstring = [ivstring createInventorSphere(mean_ach_xr(fr,:),10,cmap(60,:),0.2)];
            ivstring = [ivstring createInventorArrow(mean_ach_xr(fr,:),trialStruct(t).achilles.ach_vec(fr,:),100,2,cmap(60,:),0.2)];
            
            ivstring = [ivstring createInventorSphere(trialStruct(t).achilles.tib_ach_xr(fr,:),5,cmap(10,:),0.2)];
            ivstring = [ivstring createInventorSphere( trialStruct(t).achilles.cal_ach_xr(fr,:),5,cmap(30,:),0.2)];
            
        end
        
        ivstring = [ivstring createInventorGlobalAxes()];
        
        fid = fopen(fullfile(arrowDir,sprintf(arrow_style,ind_fr,1)),'w');
        fprintf(fid, ivstring);
        fclose(fid);
        ind_fr = ind_fr+1;
        
        
    end
    
    arrow_style_ini = strrep(arrow_style,'%i','%d');
    
    create_ini(0,0,1,2,arrowDir,arrow_style_ini, fullfile(anim_dir,'RigidFoot_HelicalAxis.ini'));
    % end of making ini------------------------
    
    % create the RTP files, pos file and linked iv files
    
    write_RTp(bone_listM,TanimM,anim_dir);
    write_RTp(bone_listF,TanimF,anim_dir);
    
    write_pos([bone_listM; bone_listF'],anim_dir,'RigidFoot_HelicalAxis');
    
    rigidiv_dir = fullfile(anim_dir,'rigidiv',filesep);
    if ~exist(rigidiv_dir,'dir')
        mkdir(rigidiv_dir);
    end
    
    for bn = 1:nBones
        % create the normal linked iv files
        ivFile = subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file;
        ivstring = createInventorHeader();
        ivstring = [ivstring createInventorLink(ivFile,eye(3),[0 0 0],[0.7 0.7 0.7],0.8)];
        fid = fopen(fullfile(rigidiv_dir,[bone_listM{bn} '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        % create the rigid foot
        ivFile = subjStruct(s_ind).bones.(bone_listM{bn}).metadata.orig_file;
        ivstring = createInventorHeader();
        ivstring = [ivstring createInventorLink(ivFile,eye(3),[0 0 0],[0.7 0.7 0.7],0.4)];
        fid = fopen(fullfile(rigidiv_dir,[bone_listM{bn} 'RIGID.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        
    end
    
    fprintf('Animation saved as: %s \n',[anim_dir,'RigidFoot_HelicalAxis'])
    
    
    
    
    
    clearvars('T_fix','TanimF','TanimM')
    
end







