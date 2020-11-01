% Arch Levering Foot
close all
clear all
clc
dir_analy = 'E:\ArchSpring_Gearing\';%'C:\Users\Lauren\OneDrive - Queen''s University\Research\Projects\2019 Beaded Mechanisms\Arch Levering\';
bone_list = {'cal','mt1','tal','tib'};
nBones = length(bone_list);

dir_ref = 'E:\ReferenceBones\';% fullfile(dir_analy,'References',filesep);

%% Load data structure
% save([dir_analy 'trialStructWBonesDetails.mat'],'trialStruct')
load([dir_analy 'trialStructWBonesDetails.mat'])
nsubj = length(trialStruct);



%% Re-make the data structure

load([dir_analy 'trialStruct.mat'])
nsubj = length(trialStruct);

for s = 1:nsubj
    
    load(fullfile(trialStruct(s).metadata.subj_dir,'Models','IV','3Aligned Reduced','bonestruct.mat'))
    
        trialStruct(s).bones = bonestruct;

%     dfield = createDfield(trialStruct(s).bones.tal.metadata.orig_file, 0.25, fullfile(dir_analy,'dfields'));
%     trialStruct(s).bones.tal.dfield = dfield;
%     dfield = createDfield(trialStruct(s).bones.tib.metadata.orig_file, 0.25, fullfile(dir_analy,'dfields'));
%     trialStruct(s).bones.tib.dfield = dfield;
    
    filename = dir(fullfile(dir_analy,'dfields',[trialStruct(s).subject '*tal*']));
    load(fullfile(dir_analy,'dfields',filename.name))
    trialStruct(s).bones.tal.dfield = Dfield;
     filename = dir(fullfile(dir_analy,'dfields',[trialStruct(s).subject '*tib*']));
    load(fullfile(dir_analy,'dfields',filename.name))
    trialStruct(s).bones.tib.dfield = Dfield;
end

save([dir_analy 'trialStructWBones.mat'],'trialStruct')

%% set up a colororder based on the subjects 
subj_list_unique = unique({trialStruct(:).subject});
figure;axes();
col_subj = get(gca,'ColorOrder');
rfs_shod = ':';
ffs_shod = '-.';
rfs_bf = '-';
ffs_bf = '--';
hold on;
    plot([1,2],[1 1],rfs_bf,'color',col_subj(1,:))
    plot([1,2],[2 2],ffs_bf,'color',col_subj(1,:))
    plot([1,2],[3 3],rfs_shod,'color',col_subj(1,:))
    plot([1,2],[4 4],ffs_shod,'color',col_subj(1,:))
legend('RFS - Barefoot','FFS - Barefoot','RFS - Shod','FFS - Shod')
ylim([0 5]); xlim([0 3])    

for s = 1:nsubj
    ind_subj = find(strcmp(trialStruct(s).subject,subj_list_unique));
    plot_cols{s,1} = col_subj(ind_subj,:); %'color'
    
    switch trialStruct(s).subject
        case 'PRFC002'% PRFC002 rfs, shod
            plot_cols{s,2} = rfs_shod;
            plot_cols{s,3} = ['RFS - Shod - ' trialStruct(s).subject];
            
        case 'PRFC003'% PRFC003 ffs, shod
            plot_cols{s,2} = ffs_shod;
            plot_cols{s,3} = ['FFS - Shod - ' trialStruct(s).subject];
            
        case 'XBSP00010'% midfoot strike, bf
            plot_cols{s,2} = rfs_bf;
            plot_cols{s,3} = ['RFS - Barefoot - ' trialStruct(s).subject];
            
        case 'XBSP00011' % ffs, bf
            plot_cols{s,2} = ffs_bf;
            plot_cols{s,3} = ['FFS - Barefoot - ' trialStruct(s).subject];
            
        case 'XBSP00012' % ffs, bf
            plot_cols{s,2} = ffs_bf;
            plot_cols{s,3} = ['FFS - Barefoot - ' trialStruct(s).subject];
            
        case 'XBSP00013' % ffs, bf
            plot_cols{s,2} = ffs_bf;
            plot_cols{s,3} = ['FFS - Barefoot - ' trialStruct(s).subject];
            
        case 'SOL001B'
            if contains(trialStruct(s).trial,'rfs')
                if contains(trialStruct(s).trial,'barefoot')
                    plot_cols{s,2} = rfs_bf;
                    plot_cols{s,3} = ['RFS - Barefoot - ' trialStruct(s).subject];
                else
                    
                    plot_cols{s,2} = rfs_shod;
                    plot_cols{s,3} = ['RFS - Shod - ' trialStruct(s).subject];
                end
            elseif contains(trialStruct(s).trial,'ffs')
                if contains(trialStruct(s).trial,'barefoot')
                    plot_cols{s,2} = ffs_bf;
                    plot_cols{s,3} = ['FFS - Barefoot- ' trialStruct(s).subject];
                else
                    
                    plot_cols{s,2} = ffs_shod;
                    plot_cols{s,3} = ['FFS - Shod - ' trialStruct(s).subject];
                end
            end
    end
    %     plot_cols{s,2} = % linestyle based on condition - shod/strike
end


%% Find the frame to reference to for each trial: when the heel comes off the back plate, or the frame of max dorsiflexion of the tibia relative to calc

% also find the frame where the tibia is most upright

for s = 1:nsubj
    nfrs = size(trialStruct(s).Tx_mm.cal,3);
    ang_DP_tibtal = []; align_z = [];ang_DP_mtp = []; ang_DP_calmet = [];
    for fr = 1:nfrs
        [ang_DP_tibtal(fr),ang_SP_tibtal(fr),ang_AD_tibtal(fr)] = eulerYZX(   trialStruct(s).Tx_mm.tib(:,:,fr),trialStruct(s).Tx_mm.tal(:,:,fr) ,...
            trialStruct(s).bones.tib.T_ACS.T_TC,  trialStruct(s).bones.tal.T_ACS.T_TC);
        
        
        [ang_DP_tibcal(fr),ang_SP_tibcal(fr),ang_AD_tibcal(fr)]  = eulerYZX(  trialStruct(s).Tx_mm.tib(:,:,fr),trialStruct(s).Tx_mm.cal(:,:,fr) ,...
            trialStruct(s).bones.tib.T_ACS.T_TC,  trialStruct(s).bones.cal.T_Aligned);  
        
        [ang_DP_calmet(fr),ang_SP_calmet(fr),ang_AD_calmet(fr)]  = eulerYZX(  trialStruct(s).Tx_mm.cal(:,:,fr),trialStruct(s).Tx_mm.mt1(:,:,fr) ,...
            trialStruct(s).bones.cal.T_Aligned,  trialStruct(s).bones.mt1.T_Aligned);
        
        
        [ang_DP_mtp(fr),ang_SP_mtp(fr),ang_AD_mtp(fr)]  = eulerYZX(  trialStruct(s).Tx_mm.mt1(:,:,fr),trialStruct(s).Tx_mm.ph1(:,:,fr) ,...
            trialStruct(s).bones.mt1.T_Aligned,  trialStruct(s).bones.ph1.T_Aligned);
        
        % find when the tibia z axis is most aligned with the global z
        T_acs_tib_x(:,:,fr) =  trialStruct(s).Tx_mm.tib(:,:,fr)* trialStruct(s).bones.tib.T_ACS.T_TC;
        
        align_z(fr) = dot( T_acs_tib_x(1:3,3,fr),[0 0 1]');
    end
    
   
    
    [~,trialStruct(s).keyFrames.max_mtp] = nanmax(ang_DP_mtp);
    [~,trialStruct(s).keyFrames.max_upright] = max(align_z);
    [~,trialStruct(s).keyFrames.max_dors] = max(ang_DP_tibtal);
    
    figure(200);hold on; plot(ang_DP_calmet,'Color', plot_cols{s,1},'Linestyle', plot_cols{s,2});
    trialStruct(s).kinemat.mtp = ang_DP_mtp;
    trialStruct(s).kinemat.mtp_ROM = diff(ang_DP_mtp([trialStruct(s).keyFrames.max_dors,trialStruct(s).keyFrames.max_mtp]));
    
    barvar(s,:) = ang_DP_calmet([trialStruct(s).keyFrames.max_dors,trialStruct(s).keyFrames.max_mtp]);
    ROM_arch(s) = diff(barvar(s,:));
     
    
    if isnan(trialStruct(s).Tx_mm.tib(:,:,trialStruct(s).keyFrames.max_mtp))
        warning(sprintf('Not enough tibia frames tracked %s %s',trialStruct(s).subject,trialStruct(s).trial));
    end
    
    
    figure(201)
    hold on;
    bar(s,ROM_arch(s),'Facecolor', plot_cols{s,1})
    
%     if isfield(trialStruct(s).keyFrames,'fr1moc')
%         s1 = trialStruct(s).keyFrames.fr1moc;
%         yyaxis right
%         hold on
%         plot(trialStruct(s).force_data(1).globForce(3,s1:end))
%         plot(trialStruct(s).force_data(2).globForce(3,s1:end))
%     end
    
    fr_tracked = find(diff(squeeze(trialStruct(s).Tx_mm.tal(1,4,:)))~=0);
    last_fr_tracked = fr_tracked(end);
    trialStruct(s).keyFrames.last_fr_tracked = last_fr_tracked;
    
    figure;
    hold on; plot(ang_DP_mtp,'Color', plot_cols{s,1},'Linestyle', plot_cols{s,2})
    plot(trialStruct(s).keyFrames.max_mtp,ang_DP_mtp(trialStruct(s).keyFrames.max_mtp),'o')
    title([trialStruct(s).subject '_' trialStruct(s).trial])
    
%     plot(align_z)


%     hold on; plot(ang_SP_tibtal)
%     hold on; plot(ang_AD_tibtal)
end




%% Load all the reference meshes for CPD


addpath(genpath('C:\Users\Lauren\Documents\git\cpd_cuda\'))

% choose the reference subject and use that mesh as the reference
% ref_subj = 'XBSP00011';
% ref_ind = find(strcmp(ref_subj,subject_list));

dir_ref_bn{1} = fullfile(dir_ref,'Calcaneus',filesep);
dir_ref_bn{3} = fullfile(dir_ref,'Talus',filesep);
dir_ref_bn{4} = fullfile(dir_ref,'Tibia',filesep);

for bn = [1,3,4]%nBones
    %         T_ACS = trialStruct(ref_ind).bones.(bone_list{bn}).T_Aligned;
    %         ref_pts = trialStruct(ref_ind).bones.(bone_list{bn}).pts;
    %         cns_temp = trialStruct(ref_ind).bones.(bone_list{bn}).cnt;
    bone_file = ls(fullfile(dir_ref_bn{bn},['*' bone_list{bn} '*.iv']));

    [ref_pts, cns_temp] = read_vrml_fast(fullfile(dir_ref_bn{bn},bone_file));
    % load the specifically created ref files
    
    %         trialStruct(s).bones.(bone_list{bn}).pts_AN = transformPoints(T_ACS,ref_pts,-1); % the anatomical points
    cpd_ref.(bone_list{bn}).cns = cns_temp(:,1:3)+1;
    cpd_ref.(bone_list{bn}).pts = ref_pts;
    
end
%% Orient all the bones in the anatomical frame and run CPD

col_map = colormap('summer');
ind = round(linspace(1,64,nsubj));
col_map = col_map(ind,:);

bone_cpd_list = fields(cpd_ref);
for bn = 3:4%1:nBones
    ivstring = createInventorHeader();
    
    for s  = 1:nsubj%6%[1:ref_ind-1,ref_ind+1:nsubj]
        
        
        T_ACS = trialStruct(s).bones.(bone_list{bn}).T_Aligned; % gives the aligned inertial axes
        bone_pts = trialStruct(s).bones.(bone_list{bn}).pts;
        pts_anat_raw = transformPoints(T_ACS,bone_pts,-1);
        if length(pts_anat_raw)>20000
            fprintf('subject %s bone %s\n',trialStruct(s).subject, bone_list{bn})
            continue
        end
                omega = 0.1; beta = 2; lambda = 3; maxIter = 500; tol = 1e-5;
        % CUDA implementation
                [pts_anat,C] = cpd_cuda(pts_anat_raw,cpd_ref.(bone_list{bn}).pts,omega, beta, lambda, maxIter, tol);
        %         Parameter lambda represents the trade off between data fitting and smoothness regularization.
        %       Parameter beta reflects the strength of interaction between points. Small values of produce locally smooth transformation, while large
        % values of beta correspond to nearly pure translation transformation
                pts_CT = transformPoints(T_ACS,pts_anat);
        
        
                save(fullfile(dir_analy,'CPD',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.mat']),'pts_anat','pts_CT','C');
        
        
        load(fullfile(dir_analy,'CPD',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.mat']))
        
        % %
        patch2iv( pts_CT, cpd_ref.(bone_list{bn}).cns ,fullfile(dir_analy,'CPD','Viz',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.iv']),[0.3 0 0.8])
        patch2iv( pts_anat_raw, trialStruct(s).bones.(bone_list{bn}).cnt(:,1:3) ,fullfile(dir_analy,'CPD','Viz',[ trialStruct(s).subject '_' bone_list{bn} '_Anat.iv']))
        patch2iv( pts_anat, cpd_ref.(bone_list{bn}).cns ,fullfile(dir_analy,'CPD','Viz',[ trialStruct(s).subject '_' bone_list{bn} '_CPDAnat.iv']))
        
        ivstring = [ivstring createInventorLink(fullfile(dir_analy,'CPD','Viz',[ trialStruct(s).subject '_' bone_list{bn} '_Anat.iv']),eye(3,3),[0 0 0],col_map(s,:),0.5)];
        
        ivstring = [ivstring createInventorLink(fullfile(dir_analy,'CPD','Viz',[ trialStruct(s).subject '_' bone_list{bn} '_Anat.iv']),eye(3,3),[s*50,0,100],col_map(s,:),0)];
        ivstring = [ivstring createInventorLink(fullfile(dir_analy,'CPD','Viz',[ trialStruct(s).subject '_' bone_list{bn} '_CPDAnat.iv']),eye(3,3),[s*50,0,200],col_map(s,:),0)];
    end
    fid = fopen(fullfile(dir_analy,'CPD','Viz',['ALL_ANAT_' bone_list{bn} '.iv']),'w');
    fprintf(fid,ivstring);
    fclose(fid);
end

%
%
% for s  = [1:ref_ind-1,ref_ind+1:nsubj]
%     for bn = 1:nBones
%
%         load(fullfile(dir_analy,'CPD',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.mat']),'pts_CT','C');
%
%
%         patch2iv( pts_CT, cpd_ref.(bone_list{bn}).cns ,fullfile(dir_analy,'CPD',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.iv']),[0.3 0 0.8])
%
%     end
% end
%% segment the talar surface and tibia surface


[tal_dome.pts,tal_dome.cns] = read_vrml_fast(fullfile(dir_ref,'Talus','ta_dome.iv'));
tal_dome.cns = tal_dome.cns(:,1:3) + 1;
[tib_dome.pts,tib_dome.cns] = read_vrml_fast(fullfile(dir_ref,'Tibia','talocrural_surf.iv'));
tib_dome.cns = tib_dome.cns(:,1:3) + 1;
% calculateSegIVindices(ref_iv,crop_iv)

[~,~,tal_dome.iRef] = intersect(round(tal_dome.pts,2),round(cpd_ref.tal.pts,2),'rows');
[~,~,tib_dome.iRef] = intersect(round(tib_dome.pts,2),round(cpd_ref.tib.pts,2),'rows');

[tal_dome.iRef2,Itd] = sort(tal_dome.iRef);
[tib_dome.iRef2,Itbd] = sort(tib_dome.iRef);
[new_tal_pts,tal_dome.cns] = rewriteCns(cpd_ref.tal.pts,cpd_ref.tal.cns,tal_dome.iRef2);
[~,tib_dome.cns] = rewriteCns(cpd_ref.tib.pts,cpd_ref.tib.cns,tib_dome.iRef2);



bn = 3;% tal
for s = 1:nsubj
    
    % load the raw points
    pts_raw = trialStruct(s).bones.(bone_list{bn}).pts;
    % load the CPD
    load(fullfile(dir_analy,'CPD',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.mat']))
    
    % get the CPD points
    pts_cpd_CT = pts_CT;%pts_raw(C,:);
    
    % segment the achilles surface
    pts_crop = pts_cpd_CT(tal_dome.iRef2,:);
%     pts_crop = pts_crop(Itd,:);
    % save the points of the achilles
    trialStruct(s).tal_dome.pts = pts_crop;
    trialStruct(s).tal_dome.cns = tal_dome.cns;
    
%     figure;
%     patch('faces',tal_dome.cns,'vertices',pts_crop,'facealpha',0.5); hold on  
%     patch('faces',trialStruct(s).bones.(bone_list{bn}).cnt(:,1:3),'vertices',pts_raw,'facealpha',0.1,'facecolor','m'); 

end

bn = 4; % tib
for s = 1:nsubj
    
    % load the raw points
    pts_raw = trialStruct(s).bones.(bone_list{bn}).pts;
    % load the CPD
    load(fullfile(dir_analy,'CPD',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.mat']))
    % get the CPD points
    pts_cpd_CT = pts_CT;%pts_raw(C,:);
    
    % segment the tal-tib surface
    pts_crop = pts_cpd_CT(tib_dome.iRef2,:);
    
    
    % save the points
    trialStruct(s).tib_dome.pts = pts_crop;
    trialStruct(s).tib_dome.cns = tib_dome.cns;
    
    
%     figure;
%     patch('faces',tib_dome.cns,'vertices',pts_crop,'facealpha',0.5); hold on
%     patch('faces',trialStruct(s).bones.(bone_list{bn}).cnt(:,1:3),'vertices',pts_raw,'facealpha',0.1,'facecolor','m'); 
%     figure;
%     pcshow(pointCloud(trialStruct(s).tib_dome.pts))
%     pcshow(pointCloud(trialStruct(s).tal_dome.pts))
end
%% Calculate the plane of the achilles on the calcaneus


cmap = get(gca,'ColorOrder');
[ach_surf.pts,ach_surf.cns] = read_vrml_fast(fullfile(dir_ref,'Calcaneus','Achilles','achilles_surf.iv'));

[~,~,iRef] = intersect(ach_surf.pts,cpd_ref.cal.pts,'rows');

bn = 1; % cal


for s = 1:nsubj
    
    % load the raw points
    pts_raw = trialStruct(s).bones.(bone_list{bn}).pts;
    % load the CPD
    load(fullfile(dir_analy,'CPD',[ trialStruct(s).subject '_' bone_list{bn} '_CPD.mat']))
    % get the CPD points
    pts_cpd_CT = pts_raw(C,:);
    
    % segment the achilles surface
    pts_ach = pts_cpd_CT(iRef,:);
    
    pts_tib_CT = trialStruct(s).bones.tib.pts;
    
    
    % save the points of the achilles
    trialStruct(s).achilles.pts = pts_ach;

     
     
     % move the tibia and calc/achilles insertion to be at the position at
     % max upright
     fr_mu = trialStruct(s).keyFrames.max_upright	;
     
     T_align_tib = trialStruct(s).Tx_mm.tib(:,:,fr_mu);
     T_align_calc = trialStruct(s).Tx_mm.cal(:,:,fr_mu);
     
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
     
     
     
    trialStruct(s).achilles.cal_pt_CT = transformPoints(T_align_calc, mean_ach_pt,-1);
    trialStruct(s).achilles.prox_pt_CT  =  transformPoints(T_align_tib, point_plane(piA,:),-1);
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

for s = 1:nsubj
    fprintf('Working on subject %s...',trialStruct(s).subject)
    fr_md = trialStruct(s).keyFrames.max_dors;
    fr_mtp = trialStruct(s).keyFrames.max_mtp;
    
    dfield = trialStruct(s).bones.tal.dfield;
    pts =  trialStruct(s).tib_dome.pts;%trialStruct(s).bones.tib.pts;
    conns = trialStruct(s).tib_dome.cns;% trialStruct(s).bones.tib.cnt(:,1:3);
    thresh = 5;
    viz = 1;
    
    
    for fr = fr_md:fr_mtp
        RT_ref = trialStruct(s).Tx_mm.tal(:,:,fr);
        RT_pts = trialStruct(s).Tx_mm.tib(:,:,fr);
        
%         if rem(fr,5) == 0
%             out = JointContact_Dfield(dfield,pts,conns,thresh,1,viz,RT_ref,RT_pts);
%         else
            
         out = JointContact_Dfield(dfield,pts,conns,thresh,1,0,RT_ref,RT_pts);
        drawnow
         trialStruct(s).contact.tib(fr).ContactCenter = out.ContactCenter;
         trialStruct(s).contact.tib(fr).ContactCenterDistance = out.ContactCenterDistance;
         trialStruct(s).contact.tib(fr).ContactArea = out.ContactArea;
         trialStruct(s).contact.tib(fr).MinDistance = out.MinDistance;
         trialStruct(s).contact.tib(fr).MinDistanceCentroid = out.MinDistanceCentroid;
         trialStruct(s).contact.tib(fr).AverageDistance = out.AverageDistance;
         trialStruct(s).contact.tib(fr).Penetration = out.Penetration;
         trialStruct(s).contact.tib(fr).Patch = out.Patch;
        
         cc_xr = transformPoints(trialStruct(s).Tx_mm.tib(:,:,fr),out.ContactCenter);
         cN_CT = nanmean(out.Patch.wtdist.* out.Patch.norm);
         cN_xr = unit(transformPoints(trialStruct(s).Tx_mm.tib(:,:,fr), cN_CT));
         
          trialStruct(s).contact.tib(fr).cc_xr = cc_xr;
          trialStruct(s).contact.tib(fr).cN_CT = cN_CT;
          trialStruct(s).contact.tib(fr).cN_xr= cN_xr;
       if out.Penetration == 1
           fprintf(' in frame %i.\n',fr)
       end
    end
   
     dfield = trialStruct(s).bones.tib.dfield;
    pts = trialStruct(s).tal_dome.pts;% trialStruct(s).bones.tal.pts;
    conns =  trialStruct(s).tal_dome.cns;%trialStruct(s).bones.tal.cnt(:,1:3);
    thresh = 4;
    viz = 0;
    
    
    for fr = fr_md:fr_mtp
        RT_ref = trialStruct(s).Tx_mm.tib(:,:,fr);
        RT_pts = trialStruct(s).Tx_mm.tal(:,:,fr);
        
%         figure(1); patch('faces',conns,'vertices',transformPoints(RT_ )
%         if rem(fr,5) == 0
%             out = JointContact_Dfield(dfield,pts,conns,thresh,1,viz,RT_ref,RT_pts);
%         else
            
            out = JointContact_Dfield(dfield,pts,conns,thresh,1,0,RT_ref,RT_pts);
         
         trialStruct(s).contact.tal(fr).ContactCenter = out.ContactCenter;
         trialStruct(s).contact.tal(fr).ContactCenterDistance = out.ContactCenterDistance;
         trialStruct(s).contact.tal(fr).ContactArea = out.ContactArea;
         trialStruct(s).contact.tal(fr).MinDistance = out.MinDistance;
         trialStruct(s).contact.tal(fr).MinDistanceCentroid = out.MinDistanceCentroid;
         trialStruct(s).contact.tal(fr).AverageDistance = out.AverageDistance;
         trialStruct(s).contact.tal(fr).Penetration = out.Penetration;
         trialStruct(s).contact.tal(fr).Patch = out.Patch;
         
         % get the normal at the contact centroid surface - weight all the
         % normals by their distance
         cc_xr = transformPoints(trialStruct(s).Tx_mm.tal(:,:,fr),out.ContactCenter);
         cN_CT = nanmean(out.Patch.wtdist.* out.Patch.norm);
         cN_xr = unit(transformPoints(trialStruct(s).Tx_mm.tal(:,:,fr), cN_CT));
         
          trialStruct(s).contact.tal(fr).cc_xr = cc_xr;
          trialStruct(s).contact.tal(fr).cN_CT= cN_CT;
          trialStruct(s).contact.tal(fr).cN_xr = cN_xr;
    end
  
    
    
    
    fprintf('... Complete!\n')
end

% close all
% 
% for s = 1%:nsubj
%       fr_md = trialStruct(s).keyFrames.max_dors;
%     fr_mtp = trialStruct(s).keyFrames.max_mtp;
%     figure; hold on
%    for fr = fr_md:fr_mtp
%     plot(fr,trialStruct(s).contact.tal(fr).ContactCenter(3),'.')
%    end
%     
% end

%% calculate the talocrural...talocalcaneal? axis
% glb is global, identity 


clearvars -except trialStruct dir_analy cmap_crop cmap nsubj plot_cols dir_ref

load('colorblind_map.mat') % loads the cmap_cb variable

close all
bone_refs = {'tib','tal';'glb','mt1';'cal','mt1'};

for s = 1:nsubj % SOL001
   
    
    
    
    Tx = trialStruct(s).Tx_mm;
    
    
    for br = 1:size(bone_refs,1)
        % total number of frames
         ntot = size(Tx.(bone_refs{br,2}),3);
         
       n_ct = nan(3,ntot); n_xr = nan(ntot,3); 
       s_ct = nan(3,ntot); s_xr = nan(ntot,3); 
        phi = nan(1,ntot); L = nan(1,ntot);
        s_xrC = nan(ntot,3); cent_xr = nan(ntot,3);
        
        fr_md = trialStruct(s).keyFrames.max_dors;
        fr_mtp = trialStruct(s).keyFrames.max_mtp;
        
        fr_md_s(s) = fr_md;
        fr_mtp_s(s) = fr_mtp;
        
% take every  frame between max ankle dorsiflexion and max mtp dorsiflexion

     % !!! added 5 frames to stabilize the helical axis and then removed
     % them after
      T = Tx.(bone_refs{br,2})(:,:,fr_md-5:fr_mtp);
        nfrs = fr_mtp-fr_md+5+1;
%        
        if strcmp(bone_refs{br,1},'glb')
            TR = repmat(eye(4),1,1,nfrs); % reference to global, i.e. identity matrix
        else
            TR = Tx.(bone_refs{br,1})(:,:,fr_md-5:fr_mtp);
        end
        
            
            
        [phi(fr_md-5:fr_mtp),n_ct(:,fr_md-5:fr_mtp),L(fr_md-5:fr_mtp),s_ct(:,fr_md-5:fr_mtp)] = stabilizeHelicalAxis(TR,T,5,9);
              
        % !!! added 5 frames to stabilize the helical axis and then removed
     % them after
        T = Tx.(bone_refs{br,2})(:,:,fr_md:fr_mtp);
        nfrs = fr_mtp-fr_md+1;
%        
        if strcmp(bone_refs{br,1},'glb')
            TR = repmat(eye(4),1,1,nfrs); % reference to global, i.e. identity matrix
        else
            TR = Tx.(bone_refs{br,1})(:,:,fr_md:fr_mtp);
        end
          
        n_ct = n_ct';
        s_ct = s_ct';
        
        n_xr(fr_md:fr_mtp,:) = transformVectors(TR,n_ct(fr_md:fr_mtp,:));
        s_xr(fr_md:fr_mtp,:) = transformPoints(TR,s_ct(fr_md:fr_mtp,:));
        
        cent = trialStruct(s).bones.(bone_refs{br,2}).centroid;
        cent_xr(fr_md:fr_mtp,:) = transformPoints(T,cent);
        %         end
        
        for fr = fr_md:fr_mtp
            s_xrC(fr,:) =  closestPointonVector( cent_xr(fr,:),s_xr(fr,:),n_xr(fr,:),0);
            
        end
        
        trialStruct(s).helical(br).ref_bone = bone_refs{br,1};
        trialStruct(s).helical(br).mv_bone = bone_refs{br,2};
        trialStruct(s).helical(br).phi = phi;
        trialStruct(s).helical(br).phi_sum = nancumtrapz(phi);
        trialStruct(s).helical(br).n_xr = n_xr;
        trialStruct(s).helical(br).s_xr = s_xr;
        trialStruct(s).helical(br).n_ct = n_ct;
        trialStruct(s).helical(br).s_ct = s_ct;
        trialStruct(s).helical(br).cent_xr = cent_xr;
        
        trialStruct(s).helical(br).s_xrC = s_xrC;
        trialStruct(s).helical(br).L = L;
        
        
    end
    ind_nonan = find(~isnan(trialStruct(s).helical(3).phi_sum));
    phi_ROM_arch(s) = trialStruct(s).helical(3).phi_sum(ind_nonan(end));
    
    
    
    
    % calculate the moment arms for the moving and the rigid foot
    
    % MOVING-----------------------------------------------------------
    T_tib = Tx.tib;
    tib_ach_p = transformPoints(T_tib, trialStruct(s).achilles.prox_pt_CT);
    T_cal =Tx.cal;
    cal_ach_d = transformPoints(T_cal, trialStruct(s).achilles.cal_pt_CT	);
    ach_vec = tib_ach_p - cal_ach_d;
    
    MA_ank = nan(ntot,3);
    MA_ach =  nan(ntot,3);
    for fr = fr_md:fr_mtp
        
 % ankle moment arm
        [s1,s2] = closestPointsBtw2Lines( trialStruct(s).contact.tal(fr).cc_xr   , trialStruct(s).helical(2).s_xrC(fr,:),...
             -trialStruct(s).contact.tal(fr).cN_xr, trialStruct(s).helical(2).n_xr(fr,:));
        MA_ank(fr,:) = s2-s1;

        % achilles moment arm
         [s1,s2] = closestPointsBtw2Lines(cal_ach_d(fr,:), trialStruct(s).helical(2).s_xrC(fr,:),...
            ach_vec(fr,:), trialStruct(s).helical(2).n_xr(fr,:));
        MA_ach(fr,:) =s2-s1;
        
      
    end
    
     trialStruct(s).achilles.tib_ach_xr = tib_ach_p;
     trialStruct(s).achilles.cal_ach_xr = cal_ach_d;
     trialStruct(s).achilles.ach_vec = ach_vec;
     trialStruct(s).MA.ank = MA_ank;
      trialStruct(s).MA.ach = MA_ach;
      
    % RIGID------------------------------------------------------------------------------
    
    T_mt1_fix = Tx.mt1(:,:,fr_md);
    T_mt1 = Tx.mt1;
    nfr = size(T_mt1,3);
    
    
    bone_listM = fields(Tx);
    nBones = length(bone_listM);
    for bn = 1:nBones
        T_bone_fix = Tx.(bone_listM{bn})(:,:,fr_md);
        if strcmp(bone_listM{bn},'tib') % give it tibio-talar kinematics
             for fr = 1:ntot
                % give all the fixed transforms with the mt1 relative to the selected max dors position
                T_fix.(bone_listM{bn})(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * Tx.tal(:,:,fr_md) *  invTranspose(Tx.tal(:,:,fr)) * Tx.tib(:,:,fr);
                % lock it in the plantarflexed position
                T_bone_fix = Tx.(bone_listM{bn})(:,:,fr_mtp);
                T_fix.tib_plant(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * Tx.tal(:,:,fr_md) *  invTranspose(Tx.tal(:,:,fr_mtp)) * Tx.tib(:,:,fr_mtp);%T_mt1(:,:,fr) * invTranspose(Tx.mt1(:,:,fr_mtp)) * T_bone_fix;
            end
        elseif  strcmp(bone_listM{bn},'ph1')% give it MTPJ kinematics
            
             for fr = 1:ntot
                % give all the fixed transforms with the mt1 relative to the selected max dors position
                T_fix.(bone_listM{bn})(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * Tx.mt1(:,:,fr_md) *  invTranspose(Tx.mt1(:,:,fr)) * Tx.ph1(:,:,fr);
                
                % THIS IS REDUNDANT BECAUSE IT COULD JUST BE Tx.ph1 -
                % because we're fixing it relative to the bone we're moving
                % everything with
            end
        else
            for fr = 1:ntot
                % give all the fixed transforms with the mt1 relative to the selected max dors position
                T_fix.(bone_listM{bn})(:,:,fr) = T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * T_bone_fix;
            end
        end
    end
    
    trialStruct(s).T_fix    = T_fix;   
    
    
    T_tibF = T_fix.tib;
    tib_ach_pF = transformPoints(T_tibF, trialStruct(s).achilles.prox_pt_CT);
    T_calF = T_fix.cal;
    cal_ach_dF = transformPoints(T_calF, trialStruct(s).achilles.cal_pt_CT	);
    ach_vecF = tib_ach_pF - cal_ach_dF;
    
    
    
    MA_ankF = nan(ntot,3);
    MA_achF =  nan(ntot,3);
    for fr = fr_md:fr_mtp
        % use the talus/ or tibia
        TC_ptF =  transformPoints(T_fix.tal(:,:,fr),trialStruct(s).contact.tal(fr).ContactCenter);
        TC_vecF =   transformVectors(T_fix.tal(:,:,fr),-trialStruct(s).contact.tal(fr).cN_CT);
        %     TC_ptF =  transformPoints(T_fix.tib(:,:,fr),trialStruct(s).contact.tib(fr).ContactCenter);
        %    	TC_vecF =   transformVectors(T_fix.tib(:,:,fr),-trialStruct(s).contact.tib(fr).cN_CT);
        
        
        % ankle moment arm
        [s1,s2] = closestPointsBtw2Lines(TC_ptF, trialStruct(s).helical(2).s_xrC(fr,:),...
           	TC_vecF, trialStruct(s).helical(2).n_xr(fr,:));
        MA_ankF(fr,:) = s1-s2;
        
        % achilles moment arm
         [s1,s2] = closestPointsBtw2Lines(cal_ach_dF(fr,:), trialStruct(s).helical(2).s_xrC(fr,:),...
            ach_vecF(fr,:), trialStruct(s).helical(2).n_xr(fr,:));
        MA_achF(fr,:) = s1-s2;
    end
    
    
      GR_short{s} = norm3d(MA_ach(fr_md:fr_mtp,:))./norm3d(MA_ank(fr_md:fr_mtp,:));
      GR_fix{s} =  norm3d(MA_achF(fr_md:fr_mtp,:))./norm3d(MA_ankF(fr_md:fr_mtp,:));
      GR_short_N(s,:) = normaliseNaN(GR_short{s},1,101);
      GR_fix_N(s,:) = normaliseNaN(GR_fix{s},1,101);
      
      
      GR(s,:) = [GR_fix{s}(end-3) GR_short{s}(end-3)];
      
      trialStruct(s).achilles.tib_ach_xrF = tib_ach_pF;
      trialStruct(s).achilles.cal_ach_xrF = cal_ach_dF;
      trialStruct(s).achilles.ach_vecF = ach_vecF;
      trialStruct(s).MA.ankF = MA_ankF;
      trialStruct(s).MA.achF = MA_achF;
      
      
      %------------- TALAR HEIGHT ANALYSIS---------------------------------
      trialStruct(s).kinemat.tal_cent_fix = transformPoints(T_fix.tal,trialStruct(s).bones.tal.centroid);
      trialStruct(s).kinemat.tal_cent_short = transformPoints(Tx.tal,trialStruct(s).bones.tal.centroid);
      
      T_ACS_ph1{s} = Tx.ph1(:,:,fr_md) * trialStruct(s).bones.ph1.T_Aligned;
      
      tal_ht_init(s,:) = transformPoints(T_ACS_ph1{s} , trialStruct(s).kinemat.tal_cent_fix(fr_md,:),-1);
      tal_ht_act(s,:) = transformPoints(T_ACS_ph1{s} , trialStruct(s).kinemat.tal_cent_short(fr_mtp,:),-1);
      tal_ht_fix(s,:) = transformPoints(T_ACS_ph1{s} , trialStruct(s).kinemat.tal_cent_fix(fr_mtp,:),-1);
      
      
      tal_range(s,:) = tal_ht_act(s,:)- tal_ht_init(s,:);
      tal_diff(s,:) = (tal_ht_fix(s,:) - tal_ht_act(s,:));
      tal_perc(s,:) = tal_diff(s,:)./tal_range(s,:) * 100;
      
      
      % -------------------- push off if tibia is in max plantarflexion----
      T_ACS_tib = trialStruct(s).bones.tib.T_ACS.T_TC;
      tib_pts = trialStruct(s).bones.tib.pts;
      % find the tibia's position
      align_z = [];  
      T_act_plant = Tx.tib(:,:,fr_mtp) * T_ACS_tib;
          for fr = fr_md:fr_mtp
              
           T_fix_plant = T_fix.tib_plant(:,:,fr) * T_ACS_tib;
           
%             pcshow(pointCloud(transformPoints(T_fix.tib_plant(:,:,fr),tib_pts)));hold on;
%             pcshow(pointCloud(transformPoints(Tx.tib(:,:,fr_mtp),tib_pts)));
%             drawnow 
           align_z(fr) = dot(T_fix_plant(1:3,3),T_act_plant(1:3,3));
          end
     
%       figure; plot(align_z)
%       drawnow
      [~,fr_tib_glob_align(s)] = max(align_z);
      mtp_vals =  trialStruct(s).kinemat.mtp([fr_md,fr_tib_glob_align(s),fr_mtp]);
      mtp_perc_tg(s) = (mtp_vals(2)-mtp_vals(1))/ (mtp_vals(3)-mtp_vals(1))*100;
%       tal_ht_act_tg(s,:) = transformPoints(T_ACS_ph1{s} , trialStruct(s).kinemat.tal_cent_short(fr_tib_glob_align(s),:),-1);
      tal_ht_fix_tg(s,:) = transformPoints(T_ACS_ph1{s} , trialStruct(s).kinemat.tal_cent_fix(fr_tib_glob_align(s),:),-1); % look at take off of rigid foot 
      
      tal_diff_tg(s,:) = (tal_ht_act(s,:) - tal_ht_fix_tg(s,:))./(tal_ht_act(s,:) -tal_ht_init(s,:))*100;
      
      red_prop(s) = (fr_mtp - fr_tib_glob_align(s))/(fr_mtp-fr_md);
      prop_time(s) = (fr_mtp - fr_tib_glob_align(s))/250;
      
      trialStruct(s).keyFrames.fr_tg = fr_tib_glob_align(s); % 
      
%       figure; pcshow(pointCloud(transformPoints(T_fix.tib_plant(:,:,fr_md),tib_pts)))
%       hold on; pcshow(pointCloud(transformPoints(Tx.tib(:,:,fr_mtp),tib_pts)))
%       ---------------------------------------------------------------------
       
    % finally, take the tib @ max mtp in the fixed position, and rotate it
    % to align with the global pos of the tib @ max mtp, see where it is on
    % the talus
    T_ACS_tal = trialStruct(s).bones.tal.T_ACS.T_TC;
    
    T_helical = invTranspose(Tx.tib(:,:,fr_mtp)) * T_fix.tib(:,:,fr_mtp);
    hel_params = convertRotation(T_helical,'4x4xn','helical');
    
    trans = 0; 
    q = [0 0 0];
    n = T_ACS_tal(1:3,1)';%hel_params(2:4);
    phi = hel_params(1);
        
    [R,T] = Helical_To_RT(phi, n, trans, q);
    Tt = eye(4);
    Tt(1:3,1:3) = R;
    Tt(1:3,4) = T;
    
    
    T_tib_fix_act_mtp =  T_fix.tib(:,:,fr_mtp)* T_ACS_tal *  Tt * invTranspose(T_ACS_tal);
   
    trialStruct(s).kinemat.T_tib_fix_act_mtp = T_tib_fix_act_mtp;
    
    % ------------------velocity of the calc------------------------------
    
     vel_mov{s} = calculateVelocity(cal_ach_d,250);
     vel_fix{s} = calculateVelocity(cal_ach_dF,250);
     
     vel_movP{s} = calculateVelocity(tib_ach_p,250);
     vel_fixP{s} = calculateVelocity(tib_ach_pF,250);
        
     vel_cmp(s,:) = [norm(vel_fix{s}(fr_mtp-3,:)-vel_fixP{s}(fr_mtp-3,:)),norm(vel_mov{s}(fr_mtp-3,:) - vel_movP{s}(fr_mtp-3,:))];
     
     figure
     plot(norm3d(vel_mov{s}(fr_md:fr_mtp,:)),'b'); hold on;
plot(norm3d(vel_movP{s}(fr_md:fr_mtp,:)),'b:');

     plot(norm3d(vel_fix{s}(fr_md:fr_mtp,:)),'k'); hold on;
plot(norm3d(vel_fixP{s}(fr_md:fr_mtp,:)),'k:');
legend('moving cal','moving tib','fixed cal','fixed tib')
title([trialStruct(s).subject ' ' trialStruct(s).trial ' ach vel'])

% title([trialStruct(s).subject ' ' trialStruct(s).trial ' ach vel'])

        
% save the moment arm lengths 
MA_ank_rigid(s,1:2) = norm3d([MA_ankF(fr_md,:);MA_ankF(fr_mtp-1,:)])';
MA_ank_moving(s,1:2) = norm3d([MA_ank(fr_md,:);MA_ank(fr_mtp-1,:)])';
MA_ach_rigid(s,1:2) = norm3d([MA_achF(fr_md,:);MA_achF(fr_mtp-1,:)])';
MA_ach_moving(s,1:2) = norm3d([MA_ach(fr_md,:);MA_ach(fr_mtp-1,:)])';

      figure(111)
      subplot(2,1,1)
      hold on;
      
      plot([1,2],[norm3d(MA_ankF(fr_mtp-1,:)),norm3d(MA_ank(fr_mtp-1,:))]-norm3d(MA_ank(fr_mtp-1,:)),'o','Color', plot_cols{s,1},'Linestyle', plot_cols{s,2})
      title('Ankle moment arm fixed (left), moving(right)')
      
    subplot(2,1,2)
       hold on;
      plot([1,2],[norm3d(MA_achF(fr_mtp-1,:)),norm3d(MA_ach(fr_mtp-1,:))]-norm3d(MA_ach(fr_mtp-1,:)),'o','Color', plot_cols{s,1},'Linestyle', plot_cols{s,2})
      title('Achilles moment arm fixed (left), moving(right)')   
      
      figure(112)
      subplot(2,1,1)
      hold on;
      
      plot([1,2],[norm3d(MA_ankF(fr_mtp-1,:)),norm3d(MA_ank(fr_mtp-1,:))],'o','Color', plot_cols{s,1},'Linestyle', plot_cols{s,2})
      title('Ankle moment arm fixed (left), moving(right)')
    subplot(2,1,2)
       hold on;
      plot([1,2],[norm3d(MA_achF(fr_mtp-1,:)),norm3d(MA_ach(fr_mtp-1,:))],'o','Color', plot_cols{s,1},'Linestyle', plot_cols{s,2})
      title('Achilles moment arm fixed (left), moving(right)')
%       plot(2,norm3d(MA_achF(end-1,:))-norm3d(MA_ach(end-1,:)),'o')
     
%       hold on; plot(norm3d(MA_ach))
%       plot(norm3d(MA_ankF),':')
%       plot(norm3d(MA_achF),':')
%     




      figure; 
      subplot(2,1,1)
      plot(norm3d(MA_ank))
      hold on; plot(norm3d(MA_ach))
      plot(norm3d(MA_ankF),':')
      plot(norm3d(MA_achF),':')
      legend('Ankle','Achilles','Ankle Fix','Achilles Fix')
      title(trialStruct(s).subject)
      
      
      subplot(2,1,2)
      plot( GR_short{s})
      hold on; plot( GR_fix{s},':')
      legend('Shortening GR','Fixed GR')
      title(trialStruct(s).subject)
%       
      
end
% 
%  figure;
%  h = bar(vel_cmp);
%  h(1).FaceColor = 'k';
%  h(2).FaceColor = 'b';
%  title('shortening velocity')
%  ylabel('achilles vel [mm/s]')
% legend('fixed','moving')
% 
%  figure;
%  h = bar(GR);
%  h(1).FaceColor = 'k';
%  h(2).FaceColor = 'b';
%  title('gear ratio')
%  ylabel('gear ratio')
% legend('fixed','moving')
% 
% 
%     figure;  
% hold on; bar(tal_perc(:,2:3))
%       title('Position difference between fixed and shortening arch')
%       ylabel('% anterior & superior [%]')
%       xlabel('subject')
%       legend('anterior + ','superior + ')

      figure(100)
      hold on; plot( GR_short_N','-')
      hold on; plot(GR_fix_N',':')

      figure; hold on;
      for s = 1:nsubj
      for fr = 1:101
          pv = GR_short_N(s,fr)-GR_short_N(s,1);%-GR_fix_N(s,fr);
%           pv = GR_fix_N(s,fr)-GR_fix_N(s,1);
          if pv > 0
              col = 'g';
          else
              col = 'r';
          end
      hold on; h = plot(fr, (pv)+s,'.','color',col);
      end
      end
      
      for s = 1:nsubj
          for fr = 1:101
%               pv = GR_short_N(s,fr)-GR_short_N(s,1);%-GR_fix_N(s,fr);
                        pv = GR_fix_N(s,fr)-GR_fix_N(s,1);
              if pv > 0
                  col = 'g';
              else
                  col = 'r';
              end
              h1 = plot(fr, (pv)+s,'+','color',col);
          end
      end
    legend([h, h1], {'shortening','fixed'})
    title('gear ratio relative to the first value for each subject')
    
      figure
      hold on;
      for s = 1:nsubj
      h(s) = plot(GR(s,:)','o','Color', plot_cols{s,1},'Linestyle', plot_cols{s,2});
       plot(2,[GR(:,end)]','kx');
      end
%       legend([h(1),h1(1)],{'fixed','shortening'})
      xlim([0,3])
      legend(h,plot_cols{:,3})
      
      figure;
      fr_look = 90;
      hold on; plot([GR_fix_N(:, fr_look),GR_short_N(:, fr_look)]','o','linestyle','-','color',[0.7 0.7 0.7])
      hold on; plot(1,[GR_fix_N(:, fr_look)]','ko')
      hold on; plot(2,[GR_short_N(:, fr_look)]','bo')
   xlim([0,3])
      
      hold on;
      h = plot([]','linestyle','-');

%       cmap = colormap('parula');
     moving_col = cmap_cb(3,:);%[0,0.5 0.8];
     fixed_col = cmap_cb(1,:);
     
      figure
      hold on;
      for s = 1:nsubj
          
         h(s) = plot(GR(s,:)','Color',[0.5 0.5 0.5],'Linestyle', plot_cols{s,2});
         h1(s) = plot(1,GR(s,1),'o','color',fixed_col);
         h2(s) = plot(2,GR(s,end),'o','color',moving_col);
        if s > 6
            h1(s).MarkerFaceColor =  fixed_col;
            h2(s).MarkerFaceColor = moving_col;
        end
      end
      
      xlim([0.5,4])
      ha = gca;
      ha.XTick = [1 2];
      ha.XTickLabel = {'Rigid','Moving'};
      legend([h1(1) h2(1) h1(7) h2(7) h(7:10 )], {'Rigid' 'Moving' 'Rigid - Case study' 'Moving - Case Study','RFS - Barefoot','FFS - Barefoot','RFS - Shod','FFS - Shod' })
      ylabel('Mechanical advantage')
      makeNicePlotsFunction
      
      
      % the velocity plot over time
  
     
    fr_lose = ones(nsubj,1)*3;
%     fr_lose(4) = 3;
%     fr_lose(3) = 2;
%     fr_lose(1) = 1;
    ach_mov_vel = nan(nsubj,101);
    ach_fix_vel = nan(nsubj,101);
     for s = 1:nsubj
        
         fr_prop = fr_tib_glob_align(s)-fr_md_s(s);
         nfr = length(fr_md_s(s):fr_mtp_s(s));
         npts_interp = round(100*(nfr-fr_lose(s))/nfr)+1;
         
         perc_prop(s) = round(fr_prop/nfr*100);
         perc_end(s) = (nfr-2)/nfr*100;
         
         ach_mov = norm3d(vel_mov{s}(fr_md_s(s):fr_mtp_s(s),:)-vel_movP{s}(fr_md_s(s):fr_mtp_s(s),:));
         ach_fix = norm3d(vel_fix{s}(fr_md_s(s):fr_mtp_s(s),:)-vel_fixP{s}(fr_md_s(s):fr_mtp_s(s),:));
         
         
          ach_mov_vel(s,1:npts_interp) = normaliseNaN(adaptiveLowPassButterworth(ach_mov(1:nfr-fr_lose(s))',[14,20],250)',1,npts_interp);%normaliseNaN(ach_mov(1:nfr-fr_lose(s)),1,npts_interp);%
          ach_fix_vel(s,1:npts_interp) = normaliseNaN(adaptiveLowPassButterworth(ach_fix(1:nfr-fr_lose(s))',[14,20],250)',1,npts_interp);%normaliseNaN(ach_fix(1:nfr-fr_lose(s)),1,npts_interp);%
          ach_vel_save(s,1:2) = [ ach_mov(nfr-fr_lose(s)),ach_fix(nfr-fr_lose(s))];
     end
     
     
     figure;
     h = PrettyStdDevGraphs(0:100,mean(ach_mov_vel),std(ach_mov_vel),moving_col,1); hold on;
     
     h1 = PrettyStdDevGraphs(0:100,mean(ach_fix_vel),std(ach_fix_vel),fixed_col,1);hold on;
     h(2).LineWidth = 5;
     h1(2).LineWidth = 5;
     for s = 1:nsubj
         h2(s) =  plot(0:100,ach_fix_vel(s,:),'color',fixed_col,'Linestyle', plot_cols{s,2},'Linewidth',2);
         plot(0:100,ach_mov_vel(s,:),'color',moving_col,'Linestyle', plot_cols{s,2},'Linewidth',2);
     end
     legend('moving arch','rigid arch ')
     xlabel('% of Propulsion')
     ylabel('Achilles shortening velocity [mm/s]')
     legend([h(2) h1(2) h2(7:10 )], {'Moving', 'Rigid',  'RFS - Barefoot','FFS - Barefoot','RFS - Shod','FFS - Shod' })
     makeNicePlotsFunction
%% calculate all the parameters to report:

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
for s = 1:nsubj
vf_1(s) = ach_fix_vel(s,perc_prop(s));
end
mean(vf_1)
std(vf_1)


% [~,p] = ttest(

fr_diff_beg = (fr_tib_glob_align'-fr_md_s')/250;
cont_time = (fr_mtp_s'-fr_md_s')/250;
early_t = cont_time-fr_diff_beg;
perc = early_t./cont_time*100;


%% calculate the contact and make a wrist viz for the tibia moved to the same location


for s = 1:nsubj
    
        
     Tx = trialStruct(s).Tx_mm;
    Tfix = trialStruct(s).T_fix;
    nfr = size(Tx.tib,3);
    bone_listM = fields(Tx);
    nBones = length(bone_listM);
    
    fr_mtp = trialStruct(s).keyFrames.max_mtp;
    
    fr_tg = trialStruct(s).keyFrames.fr_tg ;
    
    
    dfield = trialStruct(s).bones.tib.dfield;
    pts = trialStruct(s).bones.tal.pts;
    conns = trialStruct(s).bones.tal.cnt(:,1:3);
    thresh = 5;
    viz = 0;
    
    
    RT_ref =   trialStruct(s).kinemat.T_tib_fix_act_mtp;
    RT_pts = trialStruct(s).T_fix.tal(:,:,fr_mtp);
    
    
    out = JointContact_Dfield(dfield,pts,conns,thresh,1,0,RT_ref,RT_pts);
    trialStruct(s).contact.tal_sim = out;
    
    % get the normal at the contact centroid surface - weight all the
    % normals by their distance
    cc_xr = transformPoints(RT_pts,out.ContactCenter);
    cN_CT = mean(out.Patch.wtdist.* out.Patch.norm);
    cN_xr = unit(transformPoints(RT_pts, cN_CT));
    
    trialStruct(s).contact.tal_sim.cc_xr = cc_xr;
    trialStruct(s).contact.tal_sim.cN_CT= cN_CT;
    trialStruct(s).contact.tal_sim.cN_xr = cN_xr;
    
    
    
    
    anim_dir = fullfile(dir_analy, 'Animation','TalusContact',[trialStruct(s).subject '_' trialStruct(s).trial],filesep);
    if ~exist(anim_dir,'dir')
        mkdir(anim_dir)
    end
    %------ make the simulated tibia/rigid foot iv file -------------------
    
    pat = trialStruct(s).contact.tal_sim.Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSim.iv'))
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSim.iv'),RT_pts(1:3,1:3),RT_pts(1:3,4)')];
    
     for bn = 1:nBones
       
        Tf =  T_fix.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
            Tf = RT_ref;
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
     end
     
      fid = fopen(fullfile(anim_dir,'TalusContactSIM_Rigid.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    
    %------ make the actual tibia/rigid foot iv file -------------------
    pat = trialStruct(s).contact.tal(fr_mtp).Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSimNormal.iv'))
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),RT_pts(1:3,1:3),RT_pts(1:3,4)')];
    
     for bn = 1:nBones
       
        Tf =  Tfix.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
        
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
     end
     
      fid = fopen(fullfile(anim_dir,'TalusContactRigid.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    
    
    %------ make the actual tibia/rigid foot iv file -------------------
    pat = trialStruct(s).contact.tal(fr_mtp).Patch;
    patch2iv_vrml2_color(pat.pts,pat.conns,pat.color.dist{2},0.5,fullfile(anim_dir,'surfaceSimNormal.iv'))
    
    T_tal = Tx.tal(:,:,fr_mtp);
    
    ivstring = createInventorHeader();
    ivstring = [ivstring createInventorLink(fullfile(anim_dir,'surfaceSimNormal.iv'),T_tal(1:3,1:3),T_tal(1:3,4)')];
    
     for bn = 1:nBones
       
        Tf =  Tx.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
        
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
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
       
        Tf =  Tx.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
        
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        else
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
        end
        
        Tf =  Tfix.(bone_listM{bn})(:,:,fr_tg); % T of this frame
        % write all the bone links
        % moving (faded)
        if strcmp(bone_listM{bn},'tib')
         Tf =  Tfix.tib_plant(:,:,fr_tg);
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        else
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.2)];
        end
        
        
        
        
     end
     
      fid = fopen(fullfile(anim_dir,'Tibia at simulated pushoff.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);

end



%% Make a visualization of the difference in position for the rigid and non rigid foot + outline graph


for s = 1:nsubj
    anim_dir = fullfile(dir_analy, 'Animation','Final_frame_talar_Centroid',trialStruct(s).subject,filesep);
     if ~exist(anim_dir,'dir')
        mkdir(anim_dir) 
     end
    Tx = trialStruct(s).Tx_mm;
    Tfix = trialStruct(s).T_fix;
    fr_mtp = trialStruct(s).keyFrames.max_mtp;
    nfr = size(Tx.tib,3);
    bone_listM = fields(Tx);
    nBones = length(bone_listM);
    
    
    ivstring = createInventorHeader();
    
    for bn = 1:nBones
       
        Tf =  Tx.(bone_listM{bn})(:,:,fr_mtp); % T of this frame
        % write all the bone links
        % moving (faded)
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tf(1:3,1:3),Tf(1:3,4)',[0.7 0.7 0.7],0.8)];
    end
    for bn = 1:nBones
        Tff =  Tfix.(bone_listM{bn})(:,:,fr_tib_glob_align(s)); % T of this frame fixed
        % fixed (faded)
        ivstring = [ivstring createInventorLink(trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file, Tff(1:3,1:3),Tff(1:3,4)',[0.7 0.7 0.7],0.5)];
         
    end
    
    % add the talar centroid
    
    ivstring = [ivstring createInventorSphere(trialStruct(s).kinemat.tal_cent_short(fr_mtp,:) ,3,[0 0.7 0.9],0)];
    ivstring = [ivstring createInventorSphere(trialStruct(s).kinemat.tal_cent_fix(fr_mtp,:) ,3,[0 0.7 0.9],0)];
    
    fid = fopen(fullfile(anim_dir,'TalarCentroid.iv'),'w');
    fprintf(fid,ivstring);
    fclose(fid);

    fprintf('Wrote file : %s \n',fullfile(anim_dir,'TalarCentroid.iv'))
end

%% make the outlines of the bones
close all
s_pick = 6;
for s = s_pick
figure; hold on;
%     pln_vec =  T_ACS_ph1{s}(1:3,1)';
     fr_mtp = trialStruct(s).keyFrames.max_mtp;
     
    Tx = trialStruct(s).Tx_mm;
    Tfix = trialStruct(s).T_fix;
    
    T_ph1 = T_ACS_ph1{s};
    
    bone_listM = fields(Tx);
    ind_fib= contains(bone_listM,'fib');
    bone_listM(ind_fib) = [];
    nBones = length(bone_listM);
   pln_vec = [1,0,0];
    for bn = 1:nBones
        pts = transformPoints(invTranspose(T_ph1)* Tx.(bone_listM{bn})(:,:,fr_mtp),trialStruct(s).bones.(bone_listM{bn}).pts);
        ptsF = transformPoints(invTranspose(T_ph1)* Tfix.(bone_listM{bn})(:,:,fr_tib_glob_align(s)),trialStruct(s).bones.(bone_listM{bn}).pts);
%         
        npts = size(pts,1);
        pts_plane = [];
        for p = 1:npts
        pts_plane(p,:) = closestPointonPlanealongVector(pts(p,:),pln_vec, T_ACS_ph1{s}(1:3,4)',pln_vec);
        end
        k = boundary(pts_plane(:,2),pts_plane(:,3));
         outline.(bone_listM{bn}) = pts_plane(k,:);
         
         
          pts_plane = [];
        for p = 1:npts
        pts_plane(p,:) = closestPointonPlanealongVector(ptsF(p,:),pln_vec, T_ACS_ph1{s}(1:3,4)',pln_vec);
        end
        k = boundary(pts_plane(:,2),pts_plane(:,3));
      
         outlineF.(bone_listM{bn}) = pts_plane(k,:);
         
        h = fill(outline.(bone_listM{bn})(:,2),outline.(bone_listM{bn})(:,3),[.7 .7 .7],'facealpha',0.4,'Edgecolor', [.7 .7 .7],'EdgeAlpha',0.8);%'LineStyle','none');
        fill(outlineF.(bone_listM{bn})(:,2),outlineF.(bone_listM{bn})(:,3),[.7 .7 .7],'facealpha',0.8,'Edgecolor', [.7 .7 .7],'EdgeAlpha',0.8);%,'LineStyle','none');

        h = plot(tal_ht_act(s,2),tal_ht_act(s,3),'kx');
         h.MarkerSize = 10;
%         h = plot(tal_ht_fix(s,2),tal_ht_fix(s,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        h = plot(tal_ht_fix_tg(s,2),tal_ht_fix_tg(s,3),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        h.MarkerSize = 4;
        axis equal
    end
    
    for su = [1:s_pick-1,s_pick+1:nsubj]
        
%          tal_ht_fix_tg(s,:) 
%       
%       tal_diff_tg(
%        quiver( tal_ht_act(s,2),tal_ht_act(s,3), abs(tal_diff(su,2)),abs(tal_diff(su,3)),'k')
       if strcmp(trialStruct(su).subject,'SOL001B')
             h = plot(tal_ht_act(s,2)- abs(tal_diff_tg(su,2)),tal_ht_act(s,3)-abs(tal_diff_tg(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','r');% if looking at sim take off
%              h = plot(tal_ht_act(s,2)+ abs(tal_diff(su,2)),tal_ht_act(s,3)+abs(tal_diff(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','r');
        h.MarkerSize =4;
       else
       h = plot(tal_ht_act(s,2)- abs(tal_diff_tg(su,2)),tal_ht_act(s,3)-abs(tal_diff_tg(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','k');
%        h = plot(tal_ht_act(s,2)+ abs(tal_diff(su,2)),tal_ht_act(s,3)+abs(tal_diff(su,3)),'o','MarkerFaceColor','none','MarkerEdgeColor','k');
        h.MarkerSize =4;
       end
      
    end
end

makeNicePlotsFunction
%% Visualize the vectors 


cmap = colormap('parula');

c_ind = round(linspace(1,64,10));
cmap_crop = cmap(c_ind,:);


clearvars -except trialStruct dir_analy cmap_crop cmap nsubj tal_ht_act tal_ht_fix T_ACS_ph1 plot_cols

for s =  1:nsubj
    
    fr_md = trialStruct(s).keyFrames.max_dors; % frame of max dorsiflexion, choose as reference frame
%     max_mtp = trialStruct(s).keyFrames.last_fr_tracked;
    max_mtp = trialStruct(s).keyFrames.max_mtp;
    
    T = trialStruct(s).Tx_mm;
    T_mt1_fix = T.mt1(:,:,fr_md);
    T_mt1 = T.mt1;
    nfr = size(T_mt1,3);
    bone_listM = fields(T);
    nBones = length(bone_listM);
    
    TanimF = {}; TanimM = {};
    bone_listF = {};
    
    T_fix = trialStruct(s).T_fix;
    
    for bn = 1:nBones
%         T_bone_fix = T.(bone_listM{bn})(:,:,fr_md);
        ind_fr = 1;
        for fr = 1:nfr
            % give all the fixed transforms with the mt1 relative to the selected max dors position
%             T_fix.(bone_listM{bn}) = trialStruct(s).T_fix.(bone_listM{bn});% T_mt1(:,:,fr) * invTranspose(T_mt1_fix) * T_bone_fix;
            
            % write this for the animation as well; write ones where there
            % are nans
            if fr >= fr_md && fr <= max_mtp
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
         ach_pts = trialStruct(s).achilles.pts;
         mean_ach = mean(ach_pts);
         
    for fr = 1:nfr
        % move the achilles average origin 
        mean_ach_xr(fr,:) = transformPoints(T.cal(:,:,fr),mean_ach);
        
    end
    
%     trialStruct(s).T_fix = T_fix;
    
    
    
    anim_dir = fullfile(dir_analy, 'Animation','RigidFoot_HelicalAxis',[trialStruct(s).subject '_' trialStruct(s).trial],filesep);
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
     for fr = fr_md:max_mtp
         
         pat = trialStruct(s).contact.tal(fr).Patch;
        patch2iv_vrml2_color(transformPoints(T.tal(:,:,fr),pat.pts),pat.conns,pat.color.dist{2},0.5,fullfile(arrowDir,sprintf(arrow_style,ind_fr,2)))

         ivstring = createInventorHeader();
         
          cc_xr = trialStruct(s).contact.tal(fr).cc_xr ;
          cN_CT = trialStruct(s).contact.tal(fr).cN_CT ;
          cN_xr = -trialStruct(s).contact.tal(fr).cN_xr ;

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
         for bp = 1:3
             if any(isnan([trialStruct(s).helical(bp).s_xrC(fr,:),trialStruct(s).helical(bp).n_xr(fr,:)]))
                 continue
             else
                 if bp == 2
                     cc = [0 0 1];
                 else
                     cc = [0 0 1];
                 end
                ivstring = [ivstring createInventorArrow(trialStruct(s).helical(bp).s_xrC(fr,:),trialStruct(s).helical(bp).n_xr(fr,:),200,2,cc,0.2)];
             % centroid of the bone that the HA is trying to be closest to
                ivstring = [ivstring createInventorSphere(trialStruct(s).helical(bp).cent_xr(fr,:),5,cmap(40,:),0.2)];
                
             end
             
         end
         
         
         % MOMENT ARM ARROWS
        MA_achF = trialStruct(s).MA.achF(fr,:);
        MA_ankF = trialStruct(s).MA.ankF(fr,:);
        MA_ach = trialStruct(s).MA.ach(fr,:);
        MA_ank = trialStruct(s).MA.ank(fr,:);
         
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
             ivstring = [ivstring createInventorArrow(mean_ach_xr(fr,:),trialStruct(s).achilles.ach_vec(fr,:),100,2,cmap(60,:),0.2)];
             
            ivstring = [ivstring createInventorSphere(trialStruct(s).achilles.tib_ach_xr(fr,:),5,cmap(10,:),0.2)];
            ivstring = [ivstring createInventorSphere( trialStruct(s).achilles.cal_ach_xr(fr,:),5,cmap(30,:),0.2)];

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
        ivFile = trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file;
        ivstring = createInventorHeader();
        ivstring = [ivstring createInventorLink(ivFile,eye(3),[0 0 0],[0.7 0.7 0.7],0.8)];
        fid = fopen(fullfile(rigidiv_dir,[bone_listM{bn} '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        % create the rigid foot
         ivFile = trialStruct(s).bones.(bone_listM{bn}).metadata.orig_file;
        ivstring = createInventorHeader();
        ivstring = [ivstring createInventorLink(ivFile,eye(3),[0 0 0],[0.7 0.7 0.7],0.4)];
        fid = fopen(fullfile(rigidiv_dir,[bone_listM{bn} 'RIGID.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        
     end
     
     fprintf('Animation saved as: %s \n',[anim_dir,'RigidFoot_HelicalAxis'])
     
     
     
     
     
     clearvars('T_fix','TanimF','TanimM')
    
end







