
close all
clear
clc
% code to calculate the midtarsal and subtalar moments, using the various
% methods
subj_name = 'SOL001B';
subj_dir = 'E:\SOL001_VISIT2\';
% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);

list_files = dir([subj_dir '*run*barefoot*']);
% filesnum =  [1,3,4,5];%[7,17,22];op
% list_files = list_files(filesnum);

trialname_list = {list_files(:).name};
nFr = 400;



coreg = csvread(fullfile(subj_dir,'Calibration','Set 1','end_pylon','end_pylon_COREG.csv'));
% load the sync information
cub_marker = load([subj_dir  'Models\IV\Landmarks\cub_marker.stack'])';
nav_marker = load([subj_dir  '\Models\IV\Landmarks\nav_marker.stack'])';
wc = [10,20];
for i = 1:4
    nav_CT{i} = load([subj_dir  '\Models\IV\Landmarks\',sprintf('nav%i.stack',i)])';
    tal_CT{i} = load([subj_dir  '\Models\IV\Landmarks\',sprintf('tal%i.stack',i)])';
end

floor_raw = readmatrix([subj_dir  '\Calibration\Set 1\FloorPosition_fromXMA.csv']);
floor_3xn = reshape(floor_raw,3,6);
floor_mocap = transformPoints(coreg,floor_3xn);

[nfl,V,p] = affine_fit(floor_mocap');
floor_mocap = [floor_mocap'; V(1:3,1)'*400; V(1:3,2)'*400; V(1:3,1)'*400+ V(1:3,2)'*400;-V(1:3,1)'*400; -V(1:3,1)'*400+ V(1:3,2)'*400];
heelpad_33 = readmatrix([subj_dir 'T0033_SOL001_srun_pref_barefoot\XMA_csv\Heelpad_3d.csv']);
heelpad{1} = transformPoints(coreg,heelpad_33(:,1:3));
heelpad{2} = transformPoints(coreg,heelpad_33(:,4:6));
heelpad{3} = transformPoints(coreg,heelpad_33(:,7:9));


load([subj_dir  '\Calibration\Sync\syncFrames.mat'])
ntrials = length(trialname_list);
%% to load saved data

load([subj_dir 'struct_bone.mat']);
load([subj_dir 'struct_data.mat']);

bone_list = {struct_bone(:).bone_name}';
nBones = length(bone_list);
%% recalculate all the raw data
% load the force data
for i = 1:ntrials
    mocap_find_file = ls(fullfile(subj_dir,trialname_list{i},'Mocap','*_tracked*'));
    mocap_file = fullfile(subj_dir,trialname_list{i},'Mocap',mocap_find_file);
    struct_data(i) = MocapDataTransform(mocap_file,'filter','adaptive','forceCutOff',[40 100],'saveProc','off','lab','SOL','resample','force');
    
%     struct_data(i) = MocapDataTransform(mocap_file,'filter','off','saveProc','off','lab','SOL','resample','force');
end

for i = 1:ntrials
    bone_transform_file = ls(fullfile(subj_dir,trialname_list{i},'BoneTransforms','*transformsFILT*'));
    temp_bonesT = load(fullfile(subj_dir,trialname_list{i},'BoneTransforms',bone_transform_file));
    
    struct_data(i).T = temp_bonesT.T; % add the transforms to the data structure
    if i == 1
        bone_list = fields(struct_data(1).T);
        nBones = length(bone_list);
    end
end

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
    
    
    
    
    if bn == 1
        load('E:\SOL001_VISIT2\Models\IV\Reduced\SOL001B_cal_aligned_decimatedDFIELD.mat')
        bone_file = 'E:\SOL001_VISIT2\Models\IV\Reduced\SOL001B_cal_aligned_decimated.iv';
        fprintf('Creating dField. \n')
%         struct_bone(bn).dField = createDfield(bone_file, 0.4);
%         dField = struct_bone(bn).dField ;
%         save('E:\SOL001_VISIT2\Models\IV\Reduced\SOL001B_cal_aligned_decimatedDFIELD.mat','dField')
        struct_bone(bn).dField = dField;
        [pts,cns] = read_vrml_fast( bone_file);
        cns = cns+1;
        struct_bone(bn).ptsRed = pts;% reduced number of points and connections
        struct_bone(bn).cnsRed = cns;
        
    end
    
end

% make a plane manually for cutting the calcaneus mesh

% load the bone transforms + filter the data
for i = 1:ntrials
    
    
    for bn = 1:nBones
            quatT = convertRotation(struct_data(i).T.(bone_list{bn}),'4x4xn','quaternion');
            if all(isnan(quatT),'all')
                continue
            end
                
%                 Tfilt = struct_data(i).T.(bone_list{bn});
%                 quatTfilt = adaptiveLowPassButterworth(quatT',wc,250)';
%                 for fr = 1:nFr
% %                     norm(quatTfilt(fr,2:4))
%                     quatTfilt(fr,2:4) = quatTfilt(fr,2:4) /norm(quatTfilt(fr,2:4));
%                 end
%                 title(bone_list{bn})
%                 drawnow
%                 Tfilt = convertRotation(quatTfilt,'quaternion','4x4xn');
                Tfilt = struct_data(i).T.(bone_list{bn});
        for fr = 1:nFr
            
            
            struct_data(i).Tm.(bone_list{bn})(:,:,fr) = coreg * Tfilt(:,:,fr);
%             struct_data(i).Tm.(bone_list{bn})(:,:,fr) = coreg * struct_data(i).T.(bone_list{bn})(:,:,fr);
            
            
            
            
   
   
    % determine the linear velocity of the centre of mass of the bone
    p_com(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr), struct_bone(bn).cent)/1000;
            
            % move the cuboid and navicular markers conventional in motion
            % capture
            if strcmp('cub',bone_list{bn})
                struct_data(i).landmarks.(bone_list{bn})(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr),cub_marker)/1000;
            elseif strcmp('nav',bone_list{bn})
                for j = 1:4
                    struct_data(i).landmarks.(bone_list{bn}){j}(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr),nav_CT{j})/1000;
                end
            elseif strcmp('tal',bone_list{bn})
                for j = 1:4
                    struct_data(i).landmarks.(bone_list{bn}){j}(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr),tal_CT{j})/1000;
                end
            end
        end
        
         % determine the angular velocity
%         struct_data(i).w.(bone_list{bn})= zeros(3,nfr+1);
        struct_data(i).w.(bone_list{bn}) = calculateRotMatAngularVelocity(struct_data(i).Tm.(bone_list{bn})(1:3,1:3,:),250,'rad');
%         struct_data(i).w.(bone_list{bn})= adaptiveLowPassButterworth(struct_data(i).w.(bone_list{bn}),wc,250);
        % determine the linear velocity of the COM
%         struct_data(i).vcom.(bonesCell{b}) = zeros(3,nfr+1);

%         p_com = adaptiveLowPassButterworth( p_com,wc,250);
        struct_data(i).pcom.(bone_list{bn})  = p_com;
        struct_data(i).vcom.(bone_list{bn})  = calculateVelocity(p_com,250);
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
    
    
    
end


save([subj_dir 'struct_bone.mat'],'struct_bone');
save([subj_dir 'struct_data.mat'],'struct_data');
% 
% figure;
% for bn = [1,5,8,9];
%     patch('faces', struct_bone(bn).cns(1:10:end,1:3),'vertices',struct_bone(bn).pts,'facealpha',0.5);hold on;
%     plotvector3( struct_bone(bn).cent,struct_bone(bn).T_ACS(1:3,1)*50,'r')
%     plotvector3( struct_bone(bn).cent,struct_bone(bn).T_ACS(1:3,2)*50,'g')
%     plotvector3( struct_bone(bn).cent,struct_bone(bn).T_ACS(1:3,3)*50,'b')
% end
% axis equal
bone_list_char = char(bone_list);
%% calculate the cropped calcaneus


calc_plane_pts =  [91.11,56.22,-111.2; 97.44,48.65,-93.67;99.54,61.42,-135.6];
input_points = struct_bone(1).ptsRed;
input_triangles = struct_bone(1).cnsRed(:,1:3);
[nc,~,pc] = affine_fit(calc_plane_pts);

input_plane.Centre = pc;
input_plane.Normal = -nc;
[calc.output_pts,calc.output_cns] = segmentTriangularMesh(input_points,input_triangles,input_plane);
figure;
patch('faces',calc.output_cns,'vertices',calc.output_pts,'facealpha',0.7,'facecolor','b')

%% Compute the helical axes between met-calc, met-tal,cal-tal and tal-tib
close all
bone_pairs = {'cal','mt1'};%'tal','mt1';'tal','cal';'tib','tal'};
nBP = size(bone_pairs,1);
for i = 1:ntrials
    for bp = 1:nBP
        
        T_mR_ct = struct_data(i).Tm.(bone_pairs{bp,1});
        T_mA_ct = struct_data(i).Tm.(bone_pairs{bp,2});
        
        
        [phi_raw,n_raw,L_raw,s_raw] = helicalInstantaneous(T_mR_ct,T_mA_ct);
        %         n = n_raw;
        %         s = s_raw;
        phi = adaptiveLowPassButterworth(phi_raw,[10,30],250);
        L = adaptiveLowPassButterworth(L_raw,[10,30],250);
        %         n = n./norm3d(n);
        %         [phi,n,L,s] = stabilizeHelicalAxis(T_mR_ct,T_mA_ct,4);
        [~,n,~,s] =  helicalInstantaneousGapFill(T_mR_ct,T_mA_ct,6);
        
        
        n_mraw = transformVectors(T_mR_ct,n);
        s_mraw = transformPoints(T_mR_ct,s);
        
        
        n_m = adaptiveLowPassButterworth(n_mraw,[10,30],250);
        s_m = adaptiveLowPassButterworth(s_mraw,[10,30],250);
        n_m = n_m./norm3d(n_m);
        %         for j = 1:3
        %         new_pts = [s_m(j,:)-n_m(j,:); s_m(j,:)+n_m(j,:)];
        %         save_diff(j,:) = diff(new_pts);
        % %         figure; plot(new_pts')
        % %         yyaxis right; plot(diff(new_pts))
        %         end
        %         figure; plot(norm3d(save_diff))
        % %         unstable_pts = norm3d(save_diff) > 1.9;
        %
        ind0 = find(~isnan(phi));
        if i == 1
            ind1 = [126 138 160];
            ind2 = [ 126 147 177];
        elseif i == 2
            
            ind1 = [105 123 153];
            ind2 = [ 113 135 175];
        elseif i == 3
            
            ind1 = [80 91 120];
            ind2 = [82 107 135];
        end
        
        s_m(:,ind2(2):ind1(3)) = NaN;
        figure; plot(s_m'); hold on;
        % %
        %
        % %         for j = 1:3
        %             unstable_pts = [ind2(1):ind1(2),ind2(2):ind1(3)];
        %         ind0 = find(~isnan(norm3d(save_diff)));
        %         save_diff(:,unstable_pts) = NaN;
        %
        %
        %         num_flips = length(ind1);
        %         for j = 1:3
        %
        %             for nf = 1:num_flips
        %
        %                 p1 = ind1(nf) ;
        %                 p2 = ind2(nf);
        %                 if nf == 3%rem(nf,2) == 1
        %                     new_pts = [s_m(j,p1:p2)-n_m(j,p1:p2); s_m(j,p1:p2)+n_m(j,p1:p2)];
        %                 else
        %                     new_pts = [s_m(j,p1:p2)+n_m(j,p1:p2); s_m(j,p1:p2)-n_m(j,p1:p2)];
        %                 end
        %                 save_diff(j,p1:p2) = diff(new_pts);
        %             end
        %         end
        %
        %
        %         figure; plot(save_diff','.')
        
        
        struct_helical(i).trials(bp).refBone = bone_pairs{bp,1};
        struct_helical(i).trials(bp).bone = bone_pairs{bp,2};
        struct_helical(i).trials(bp).phi = phi;
        struct_helical(i).trials(bp).n = n;
        struct_helical(i).trials(bp).L = L;
        struct_helical(i).trials(bp).s = s;
        
        
        %         ind_sm = phi < 0.4;
        %         n(:,ind_sm) = NaN;
        %         s(:,ind_sm) = NaN;
        
        %
        %         n_m(:,ind0) = interpolateNaN(save_diff(:,ind0),ind0,'pchip');
        %         inda = ~isnan(n_m(1,:))
        %         for j = 1:3
        %         interp1(n_m(j,inda),find(inda),ind0,'spline','extrap')
        %         end
        %         n_m = n_m./norm3d(n_m);
        s_m(:,ind0) = interpolateNaN(s_m(:,ind0),ind0,'pchip');
        plot(s_m','--')
        figure; plot(n_m')
        %         figure;
        
        struct_helical(i).trials(bp).n_m = n_m;
        struct_helical(i).trials(bp).s_m = s_m;
        
        bn_ind = findInStruct(struct_bone,'bone_name',bone_pairs{bp,1});
        for fr = 1:nFr-1
            centR_moc = transformPoints(T_mR_ct(:,:,fr),struct_bone(bn_ind).cent);
            s_ref_moc(:,fr) = closestPointonVector(centR_moc,s_m(:,fr),n_m(:,fr));
        end
        struct_helical(i).trials(bp).s_ref_moc = s_ref_moc;
        
        
    end
    
    
end

% for i = 1:3
% figure; plot(struct_helical(i).trials(bp).n_m')
% end

%% animate the force vector in Wrist visualizer

c = [195,48,48;
    0, 153, 153;
    51,51,255;
    153,255,154]/256;
% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
% iv_dir = fullfile(subj_dir,'Models','IV',filesep);
for i = 1:ntrials
    bone_list_write = {};
    anim_dir = fullfile(subj_dir,trialname_list{i},'POS','Filtered_and_GRF',filesep);
    rigidiv_folder = fullfile(anim_dir,'rigidiv',filesep);
    
    % set up the GRF animation folder
    grf_style = 'grf_P%i_F1.iv';
    grf_folder = fullfile(anim_dir,'GRFdir',filesep);
    
    if exist(anim_dir,'dir')==0;     mkdir(anim_dir);     end
    if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end
    if exist(grf_folder,'dir')==0;  mkdir(grf_folder);  end
    
    first_fr = struct_data(i).cropFrsXray(1);
    end_frM   = struct_data(i).cropFrsXray(2);
    
    for bn = 1:nBones
        ivstring = createInventorHeader();
        % make the linked iv file
        ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list_char(bn,:) '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
        
        fid = fopen(fullfile(rigidiv_folder,[bone_list_char(bn,:) '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
    ind = 1;
    for bn = 1:nBones
        if ~isfield(struct_data(i).Tm,(bone_list{bn}))
            continue
        else
            bone_list_write{ind} = bone_list{bn};
            ind = ind+1;
        end
        
        
        for fr = first_fr:end_frM
            Ta.(bone_list{bn})(:,:,fr) = struct_data(i).Tm.(bone_list{bn})(:,:,fr);
        end
        Ta.(bone_list{bn})(isnan(Ta.(bone_list{bn}))) = 1;
        write_RTp(bone_list_char(bn,:) , Ta.(bone_list{bn})(:,:,first_fr:end_frM) , anim_dir)
    end
    
    
    frX = 1;
    for fr = struct_data(i).cropFrsForce(1):struct_data(i).cropFrsForce(2)
        
        grf_filename = sprintf(grf_style,frX);
        frX = frX+1;
        ivstring = createInventorHeader();
        
        for bp = 1:nBP
            if isnan(struct_helical(i).trials(bp).s_m(1,fr-251))
                continue
            end
            ivstring = [ivstring createInventorArrow( struct_helical(i).trials(bp).s_ref_moc(:,fr-251)-50*struct_helical(i).trials(bp).n_m(:,fr-251),...
                struct_helical(i).trials(bp).n_m(:,fr-251),...
                struct_helical(i).trials(bp).phi(fr-251)*100,...
                3,c(bp,:),0)];
        end
        
        
        for fp = 1:2 % each force plate
            COP = struct_data(i).force_data(fp).globCOP(:,fr)*1000;% transformPoints(coreg,struct_data(i).force_data(fp).globCOP(:,fr)*1000,-1);
            force = struct_data(i).force_data(fp).globForce(:,fr);%transformVectors(coreg,struct_data(i).force_data(fp).globForce(:,fr),-1);
            
            %         ivstring = [ivstring createInventorSphere(COP,4,[0.2 0.2 0.8],0)];
            ivstring = [ivstring createInventorArrow(COP,force,norm(force)/5,1.5,[0.2 0.2 0.8],0)];
        end
%         (centroid,width,height,depth,color,transparency)
%         ivstring = [ivstring createInventorCube(p,400,0.1,400,[0.4 0.4 0.4],0.5)
        ivstring = [ivstring createInventorPlanefit(floor_mocap,1,[0.4 0.4 0.4],0.5)];
        fid = fopen(fullfile(grf_folder,grf_filename),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        
    end
    create_ini(0,0,1,1,grf_folder(1:end-1),'grf_P%d_F%d.iv',fullfile(anim_dir,[trialname_list{i} '.ini']))
    
    pos_text = write_pos(bone_list_write,anim_dir,trialname_list{i});
    
    filename = fullfile(anim_dir, [trialname_list{i} '.pos']);
    
    fid = fopen(filename,'w'); % open the file to write
    fprintf(fid,pos_text);
    fclose(fid);
    
    fprintf('Animations are saved in %s.\n',anim_dir)
    
    %
    %     fr = 100;
    %     figure;
    %     hold on;
    %     pts = struct_bone(1).pts;
    %     cns = struct_bone(1).cns;
    %     patch('faces',cns(1:10:end,1:3),'vertices',transformPoints(Tw.(bone_list{1})(:,:,fr),pts) ,'facealpha',0.3)
    %
    %     pts = struct_bone(5).pts;
    %     cns = struct_bone(5).cns;
    %     patch('faces',cns(1:10:end,1:3),'vertices',transformPoints(Tw.(bone_list{5})(:,:,fr),pts) ,'facealpha',0.3)
    %     fr = fr+251
    %     for fp = 1:2
    %     plotvector3(struct_data(i).force_data(fp).globCOP(:,fr)*1000,struct_data(i).force_data(fp).globForce(:,fr))
    %     end
    %
    
    
end


%% Calculate the moments

new_centre = [-0.2 0 0]';

% determine where the force is applied to the foot
fploc1 = struct_data(i).force_data(1).ForcePlateLocation;
fploc2 = struct_data(i).force_data(2).ForcePlateLocation;
[floor_vec1,~,p1] = affine_fit(fploc1);
[floor_vec2,~,p2] = affine_fit(fploc2);

warning('adjust top plate location')

p1 = p';%p1'+[0,0,0.0531]';
p2 = p';%p2'+[0,0,0.0531]';
floor_vec1 = nfl;
floor_vec2 = nfl;
clearvars('s')

for i = 1:ntrials
    
    
    force1 = struct_data(i).force_data(1).globForce;
    force2 = struct_data(i).force_data(2).globForce;
    freemom1 = struct_data(i).force_data(1).globFreeMoment;
    freemom2 = struct_data(i).force_data(2).globFreeMoment;
    COP1 = struct_data(i).force_data(1).globCOP;
    COP2 = struct_data(i).force_data(2).globCOP;
    
    start_frM = struct_data(i).cropFrsForce(1);
    end_frM = struct_data(i).cropFrsForce(2);
    COPf1 = [];
    COPf2 = [];
    for fr = start_frM:end_frM %1:length(force1) % move the COP to the top of the dragonplate
        COPf1(:,fr) = closestPointonPlanealongVector(COP1(:,fr),floor_vec1,p1,force1(:,fr));
        COPf2(:,fr) = closestPointonPlanealongVector(COP2(:,fr),floor_vec2,p2,force2(:,fr));
    end
    % make the processed structure with synced and cropped
    s(i).force{1} = force1(:,start_frM:end_frM);
    s(i).force{2} = force2(:,start_frM:end_frM);
    s(i).COP{1} = COPf1(:,start_frM:end_frM);
    s(i).COP{2} = COPf2(:,start_frM:end_frM);
    s(i).freemom{1} = freemom1(:,start_frM:end_frM);
    s(i).freemom{2} = freemom2(:,start_frM:end_frM);
    
    start_frX = struct_data(i).cropFrsXray(1);
    end_frX = struct_data(i).cropFrsXray(2);
    
    if isfield(struct_data(i).landmarks,'cub')
    s(i).cub = struct_data(i).landmarks.cub(:,start_frX:end_frX);
    end
    %     s(i).nav = struct_data(i).landmarks.nav(:,start_frX:end_frX);
    
    if isfield(struct_data(i).landmarks,'nav')
    for j = 1:8 % number of axes to generate from points
        
        if j <=4
            s(i).marker2{j} = struct_data(i).landmarks.nav{j}(:,start_frX:end_frX);
            
        else
            
            s(i).marker2{j} = struct_data(i).landmarks.tal{j-4}(:,start_frX:end_frX);
        end
    end
    end
    
    for bp = 1:nBP
        s(i).heli(bp).n_m =    struct_helical(i).trials(bp).n_m(:,start_frX:end_frX);
        s(i).heli(bp).s_m =    struct_helical(i).trials(bp).s_m(:,start_frX:end_frX);
        s(i).heli(bp).s_ref_moc =    struct_helical(i).trials(bp).s_ref_moc(:,start_frX:end_frX);
        s(i).heli(bp).phi =    struct_helical(i).trials(bp).phi(start_frX:end_frX);
    end
    for bn = 1:nBones
        if ~isfield(struct_data(i).Tm,(bone_list{bn}))
            continue
        end
        s(i).T.(bone_list{bn}) = struct_data(i).Tm.(bone_list{bn})(:,:, start_frX :end_frX );
        s(i).w.(bone_list{bn}) = struct_data(i).w.(bone_list{bn})(:, start_frX :end_frX );
        s(i).vcom.(bone_list{bn}) = struct_data(i).vcom.(bone_list{bn})(:,start_frX :end_frX );
        s(i).pcom.(bone_list{bn}) = struct_data(i).pcom.(bone_list{bn})(:, start_frX :end_frX );
    end
    
    
    [s(i).ang.archdors,s(i).ang.archabd,s(i).ang.archabd] = eulerYZX( s(i).T.cal,  s(i).T.mt1,struct_bone(1).T_ACS, struct_bone(5).T_ACS);
%     [s(i).ang.toedors,s(i).ang.toeabd,s(i).ang.toeabd] = eulerYZX( s(i).T.mt1,  s(i).T.ph1,struct_bone(5).T_ACS, struct_bone(8).T_ACS);


    % ************** subtalar joint moment arm from GRF*****************
    
    % bone_pairs = {'cal','mt1';'tal','mt1';'tal','cal';'tib','tal'};
    %     for fr = 1:(end_frM - start_frM + 1)
    %
    %
    %         bp = 3;
    %         % find the two closest points (i.e. the points that connect the
    %         % moment arm
    %         [s1,s2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  , s(i).heli(bp).s_m(:,fr)/1000, ...
    %                                     s(i).force{2}(:,fr) , s(i).heli(bp).n_m(:,fr));
    %
    %         r_g_ST = s1-s2;
    %         j_ST = s2;
    %
    %         s(i).r_g_ST(:,fr) = r_g_ST;
    %         s(i).M_ST(:,fr) = cross(r_g_ST,s(i).force{2}(:,fr)) + s(i).freemom{2}(:,fr);
    % %         s(i).M_ST_conv_calc(:,fr) = transformVectors( s(i).T.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_ST(:,fr));
    % %         figure; plotvector3(s(i).COP{1}(:,fr),   s(i).force{1}(:,fr)/norm(  s(i).force{1}(:,fr))*100); hold on; plotvector3(s(i).heli(bp).s_m(:,fr),s(i).heli(bp).n_m(:,fr)/norm(s(i).heli(bp).n_m(:,fr))*100);
    % %         plot3quick_scatter([s1,s2])
    % %         plotvector3(j_ST,r_g_ST); axis equal;
    %
    %
    %
    %         bp = 2;
    %         % find the two closest points (i.e. the points that connect the
    %         % moment arm from the COP under the mt1 to the mt1-talus helical
    %         % axis
    %         [s1,s2] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  , s(i).heli(bp).s_m(:,fr)/1000, ...
    %                                     s(i).force{1}(:,fr) , s(i).heli(bp).n_m(:,fr));
    %
    %         r_g_MtT = s1-s2;
    % %         j_ST = s2;
    %
    %         s(i).r_g_MtT(:,fr) = r_g_MtT;
    %         s(i).M_MtT(:,fr) = cross(r_g_MtT,s(i).force{1}(:,fr)) + s(i).freemom{1}(:,fr);
    %
    %
    % %         s(i).M_MtT_conv_calc(:,fr) =
    % %         transformVectors( s(i).T.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_ST(:,fr));
    %
    %
    %
    %
    %
    %     end
    
    FT = [];
    momT = [];
    r1_com = [];
    r2_com = [];
    M_sum = [];
    cop_new= [];
    %*********** midtarsal joint calculated conventionally + around the axis **********
    for fr = 1:(end_frM - start_frM + 1)
        
        % find the two closest points (i.e. the points that connect the
        % moment arm
%         
%         
%         for j = 1:8
%             midtarsal_JC_a = nanmean([ s(i).cub(:,fr) , s(i).marker2{j}(:,fr) ], 2);
%             
%             axis_vec = diff([ s(i).cub(:,fr) ,  s(i).marker2{j}(:,fr)  ],[],2);
%             axis_vec = axis_vec/norm(axis_vec);
%             [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  ,s(i).marker2{j}(:,fr) , ...
%                 s(i).force{1}(:,fr) ,  axis_vec);
%             
%             [sa2,sb2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  ,s(i).marker2{j}(:,fr) , ...
%                 s(i).force{2}(:,fr) ,  axis_vec);
%             
%             r_g_MT1 = sa1-sb1;
%             r_g_MT2 = sa2-sb2;
%             s(i).s_g_MT1test{j}(:,fr) = sb1;
%             s(i).s_g_MT2test{j}(:,fr) = sb2;
%             s(i).r_g_MT1test{j}(:,fr) = r_g_MT1;
%             s(i).r_g_MT2test{j}(:,fr) = r_g_MT2;
%             s(i).axis_vec_test{j}(:,fr) = axis_vec;
%             s(i).JC_test{j}(:,fr) = midtarsal_JC_a;
%             
%             s(i).M_MT_test{j}(:,fr) = cross(r_g_MT1,s(i).force{1}(:,fr)) + s(i).freemom{1}(:,fr) + cross(r_g_MT2,s(i).force{2}(:,fr)) + s(i).freemom{2}(:,fr);
%             s(i).M_MT_test_calc{j}(:,fr) = transformVectors( s(i).T.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_test{j}(:,fr));
%             
%             
%             
%         end
        
        % calc
            r_com_cop =  s(i).COP{2}(:,fr) - s(i).pcom.cal(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.cal(:,fr) + cross(s(i).w.cal(1:3,fr), r_com_cop);
            PUD(i).cal(:,fr) = dot(s(i).force{2}(:,fr),v_dist) + dot(s(i).freemom{2}(:,fr),s(i).w.cal(:,fr));
            
        % mt1    
            r_com_cop =  s(i).COP{1}(:,fr) - s(i).pcom.mt1(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.mt1(:,fr) + cross(s(i).w.mt1(1:3,fr), r_com_cop);
            PUD(i).mt1(:,fr) = dot(s(i).force{1}(:,fr),v_dist) + dot(s(i).freemom{1}(:,fr),s(i).w.mt1(:,fr));
            
        
        
        % calc
            r_com_cop =  s(i).COP{1}(:,fr) - s(i).pcom.cal(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.cal(:,fr) + cross(s(i).w.cal(1:3,fr), r_com_cop);
            PUD(i).calO(:,fr) = dot(s(i).force{1}(:,fr),v_dist) + dot(s(i).freemom{1}(:,fr),s(i).w.cal(:,fr));
            
        % mt1    
            r_com_cop =  s(i).COP{2}(:,fr) - s(i).pcom.mt1(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.mt1(:,fr) + cross(s(i).w.mt1(1:3,fr), r_com_cop);
            PUD(i).mt1O(:,fr) = dot(s(i).force{2}(:,fr),v_dist) + dot(s(i).freemom{2}(:,fr),s(i).w.mt1(:,fr));
            
            
        
            
        % tal from fp1    
            r_com_cop =  s(i).COP{1}(:,fr) - s(i).pcom.tal(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.tal(:,fr) + cross(s(i).w.tal(1:3,fr), r_com_cop);
            PUD(i).tal1(:,fr) = dot(s(i).force{1}(:,fr),v_dist) + dot(s(i).freemom{1}(:,fr),s(i).w.tal(:,fr));
            
        % tal from fp2 
            r_com_cop =  s(i).COP{2}(:,fr) - s(i).pcom.tal(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.tal(:,fr) + cross(s(i).w.tal(1:3,fr), r_com_cop);
            PUD(i).tal2(:,fr) = dot(s(i).force{2}(:,fr),v_dist) + dot(s(i).freemom{2}(:,fr),s(i).w.tal(:,fr));
%  
        % get the total force
        FT(:,fr) = s(i).force{1}(:,fr) + s(i).force{2}(:,fr);
        momT(:,fr) = s(i).freemom{1}(:,fr) + s(i).freemom{2}(:,fr);
        % calculate the new COP by summing moments around new origin and
        % calculating the new COP
        r1_com(:,fr) = s(i).COP{1}(:,fr) - new_centre;%p_com.(bonesCell{b})(:,f);
        r2_com(:,fr) = s(i).COP{2}(:,fr) - new_centre;%p_com.(bonesCell{b})(:,f);
        
        M_sum(:,fr) = cross(r1_com(:,fr),s(i).force{1}(:,fr)) + cross(r2_com(:,fr),s(i).force{2}(:,fr));
        
        
        cop_new(:,fr) = [-M_sum(2,fr)/FT(3,fr) ;...
            M_sum(1,fr)/FT(3,fr)  ;...
            0];
        
        cop_new(:,fr) = cop_new(:,fr) + new_centre;
        
        % talus
        
            r_com_cop =    cop_new(:,fr) - s(i).pcom.tal(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.tal(:,fr) + cross(s(i).w.tal(1:3,fr), r_com_cop);
            PUD(i).tal(:,fr) = dot(FT(:,fr),v_dist) + dot(momT(:,fr),s(i).w.tal(:,fr));
        
            % tib
        
            r_com_cop =    cop_new(:,fr) - s(i).pcom.tib(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.tib(:,fr) + cross(s(i).w.tib(1:3,fr), r_com_cop);
            PUD(i).tib(:,fr) = dot(FT(:,fr),v_dist) + dot(momT(:,fr),s(i).w.tib(:,fr));
        
                    % distal calc power
        
            r_com_cop =    cop_new(:,fr) - s(i).pcom.cal(:,fr); % second force plate is only one contacting heel
            v_dist = s(i).vcom.cal(:,fr) + cross(s(i).w.cal(1:3,fr), r_com_cop);
            PUD(i).caldist(:,fr) = dot(FT(:,fr),v_dist) + dot(momT(:,fr),s(i).w.cal(:,fr));
        
            %
        %         midtarsal_JC = nanmean([ s(i).cub(:,fr) ,  s(i).nav{2}(:,fr)  ],2);
        %
        %         nav_cub_vec = diff([ s(i).cub(:,fr) ,  s(i).nav(:,fr)  ],[],2);
        %         nav_cub_vec = nav_cub_vec/norm(nav_cub_vec);
        %
        %
        %         r_g1 =  s(i).COP{1}(:,fr) - midtarsal_JC ;
        %         r_g2 =  s(i).COP{2}(:,fr) - midtarsal_JC ;
        % %         sv1(:,fr) = cross(r_g1,s(i).force{1}(:,fr));
        %         s(i).M_MT_conv(:,fr) = cross(r_g1,s(i).force{1}(:,fr)) + s(i).freemom{1}(:,fr) + cross(r_g2,s(i).force{2}(:,fr)) + s(i).freemom{2}(:,fr);
        %
        %         s(i).M_MT_conv_calc(:,fr) = transformVectors( s(i).T.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_conv(:,fr));
        %
        
        
        % using the vector between cuboid and navicular instead of just the point
        %         [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  , s(i).cub(:,fr) , ...
        %             s(i).force{1}(:,fr) , nav_cub_vec);
        %
        %          [sa2,sb2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  , s(i).cub(:,fr) , ...
        %             s(i).force{2}(:,fr) , nav_cub_vec);
        %
        %         r_g_MT1 = sa1-sb1;
        %         r_g_MT2 = sa2-sb2;
        %         s(i).s_g_MT1(:,fr) = sb1;
        %         s(i).s_g_MT2(:,fr) = sb2;
        %
        %
        %         s(i).M_MT_vec(:,fr) = cross(r_g_MT1,s(i).force{1}(:,fr)) + s(i).freemom{1}(:,fr) + cross(r_g_MT2,s(i).force{2}(:,fr)) + s(i).freemom{2}(:,fr);
        %          s(i).M_MT_vec_calc(:,fr) = transformVectors( s(i).T.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_vec(:,fr));
        %
        %
        %         s(i).r_g_MT1(:,fr) = r_g_MT1;
        %         s(i).r_g_MT2(:,fr) = r_g_MT2;
        %
        
        % use the arch axis -
        bp = 1;
        [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  , s(i).heli(bp).s_ref_moc(:,fr)/1000, ...
            s(i).force{1}(:,fr) , s(i).heli(bp).n_m(:,fr));
        r_g_MT1 = sa1-sb1;
        [sa2,sb2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  , s(i).heli(bp).s_ref_moc(:,fr)/1000, ...
            s(i).force{2}(:,fr) , s(i).heli(bp).n_m(:,fr));
        r_g_MT2 = sa2-sb2;
        
        s(i).s_g_MT1ax(:,fr) = sb1;
        s(i).s_g_MT2ax(:,fr) = sb2;
        s(i).r_g_MT1ax(:,fr) = r_g_MT1;
        s(i).r_g_MT2ax(:,fr) = r_g_MT2;
        s(i).M_MT_ax(:,fr) = cross(r_g_MT1,s(i).force{1}(:,fr)) + cross(r_g_MT2,s(i).force{2}(:,fr)) + s(i).freemom{1}(:,fr) + s(i).freemom{2}(:,fr);
        s(i).M_MT_ax_calc(:,fr) = transformVectors( s(i).T.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_ax(:,fr));
        
        
        
        
        %********* calculate the distance of the calcaneus from the floor
 
        if norm(s(i).force{2}(1:3,fr)) > 10 % when the calc force is bigger than 10 N
            n_pts_gr = 50;
            newPts = transformPoints(s(i).T.cal(:,:,fr),calc.output_pts);
            
            limsX = [min(newPts(:,1)) max(newPts(:,1))];
            limsY = [min(newPts(:,2)) max(newPts(:,2))];
            x_grid = linspace(limsX(1),limsX(2),n_pts_gr);
            y_grid = linspace(limsY(1),limsY(2),n_pts_gr);
            [X_grid,Y_grid] = meshgrid(x_grid,y_grid);
            %         limsZ = [min(newPts(:,3)) max(newPts(:,3))];
            findPlaneD = @(n,x,y,z) n(1) * x + n(2) * y + n(3) * z;
            d = findPlaneD(floor_vec2,p2(1),p2(2),p2(3));
            floorPlaneZ = @(n,x,y) (n(1) * x + n(2) * y - d)/(-n(3));
            Z_grid = floorPlaneZ(floor_vec2,X_grid,Y_grid);
            
            for ipx = 1:n_pts_gr
                for ipy = 1:n_pts_gr
                    chkPt = [X_grid(ipx,ipy),Y_grid(ipx,ipy),Z_grid(ipx,ipy)];
                    s(i).dfield_val(ipx,ipy,fr) = (lookUpDfieldPts(struct_bone(1).dField,chkPt,s(i).T.cal(1:3,1:3,fr),s(i).T.cal(1:3,4,fr)));
                end
            end
            
            s(i).X_grid = X_grid;
            s(i).Y_grid = Y_grid;
            s(i).Z_grid = Z_grid;
            df_fr = s(i).dfield_val(:,:,fr);
            df_sort = sort(df_fr(:));
            df_top = df_sort(1:round(n_pts_gr^2/2)); % get the top 50% of compressed points
            s(i).ave_df(fr) = mean(df_top);
            s(i).top_df(fr) = df_top(end);
            
            if i == 3 && fr ==119
%                  lookUpDfieldPts(struct_bone(1).dField, heelpad{2}(119,:),struct_data(i).T.cal(1:3,1:3,119),struct_data(i).T.cal(1:3,4,119))
            end
            
            
            
        end
%**********************************************
        
    end
    fields_list = {'MT','ST','cal','tib','caldist','TC','tal','tal1','tal2','tal12','calO','mt1O','mt1','MTtemp'};
    
    
    PUD(i).TC = PUD(i).tib - PUD(i).tal;
    PUD(i).MT = PUD(i).tal - PUD(i).mt1;
    
    PUD(i).ST = PUD(i).tal - PUD(i).cal;
    PUD(i).tal12 = PUD(i).tal1 + PUD(i).tal2;
    PUD(i).MTtemp = (PUD(i).calO + PUD(i).mt1O -PUD(i).cal - PUD(i).mt1)/2;
    
    for ff = 1:length(fields_list)
        
        nan_ind = ~isnan(PUD(i).(fields_list{ff}));
        npts = (1:length(nan_ind))/250;
        work(i).(fields_list{ff}) = nan(1,length(nan_ind));
        work(i).(fields_list{ff})( nan_ind ) = cumtrapz(npts(nan_ind),PUD(i).(fields_list{ff})( nan_ind ));
        indpos = PUD(i).(fields_list{ff}) >0;
        indneg = PUD(i).(fields_list{ff}) <0;
        nptsP = length(find(indpos));
        nptsN = length(find(indneg));
        poswork(i).(fields_list{ff}) = trapz(npts(indpos),PUD(i).(fields_list{ff})(indpos));
        negwork(i).(fields_list{ff}) = trapz(npts(indneg),PUD(i).(fields_list{ff})(indneg));
    end
    
    
    s(i).FT = FT;
    s(i).momT = momT;
    s(i).cop_new = cop_new;
    
    
    
end


%% the visualization of heel pad deformation
    
    cmap = colormap('jet');
    cmap = flipud(cmap);
for i = 3
figure;
colormap(cmap)
    for fr = -1:size(s(i).dfield_val,3)
        if fr > 0
    surf(s(i).X_grid,s(i).Y_grid,s(i).Z_grid, s(i).dfield_val(:,:,fr))
        else
            
    surf(s(i).X_grid,s(i).Y_grid,s(i).Z_grid, 20*ones(n_pts_gr))
        end
    shading interp
    hold on;
    colorbar
    caxis([10 20])
    
    h = plot3quick(heelpad{1}(120+fr,:)','k','x');hold on;
    h.MarkerSize = 10;
    h = plot3quick(heelpad{2}(120+fr,:)','k','x');   h.MarkerSize = 10;
    h = plot3quick(heelpad{3}(120+fr,:)','k','x');   h.MarkerSize = 10;
    
    newpts = transformPoints(struct_data(i).Tm.cal(:,:,120+fr),struct_bone(1).ptsRed);
    
    h = patch('faces',struct_bone(1).cnsRed(1:5:end,1:3),'vertices',newpts);
    h.FaceAlpha = 0.7;
    h.FaceColor = 'b';
    hold off
    
    axis equal
    view([-81,7])
    drawnow
    pause(0.2)
    end
end

%% plotting
% close all
% bp=1;i=1;
% figure;plot(norm3d(s(i).r_g_MT1ax')); hold on
% yyaxis right; plot(s(i).force{1}')
% yyaxis left; plot(s(i).heli(bp).phi)
%
%
% figure;plot(norm3d(s(i).r_g_MT2ax')); hold on
% yyaxis right; plot(s(i).force{2}')
% yyaxis left; plot(s(i).heli(bp).phi)

% title_str = 'Using navicular-cuboid axis - in talus cosys';
% title_str = 'Using mt1-calc axis-in talus cosys';
%
% figure;
% hold on;
% for i = 1:3
% plot(PUD(i).tal - PUD(i).mt1)
% end
% title('Midtarsal Power')
%
%
% figure;
% hold on;
% for i = 1:3
% plot(PUD(i).tal - PUD(i).cal)
% end
% title('Subtalar Power')
% fields_list = {'MT','ST','cal','tib','caldist','TC'};

    fields_list = {'cal'}%,'MT','ST','tib','caldist','TC','tal','mt1'};%,'tal1','tal2','tal12','calO','mt1O','MTtemp'};
    
fields_names = {'calc'}%,'midtarsal','subtalar'}%,'tib','calc dist','talocrural','tal','mt1','tal1','tal2','tal12','calO','mt1O','MTtemp'}


%     fields_list = {'MT','ST','cal','tib','caldist','TC'};
    
% fields_names = {'midtarsal','subtalar','calc','tib','calc dist','talocrural'}
nfields = length(fields_list);
leg_str = {'slow-decel','slow-const','slow-accel','norm-decel','fast-decel'};
plot_order = [3,4,2]%,1,5];
figure; hold on;
for i = plot_order
    plot(cumtrapz((1/250)*1:length(s(i).FT),s(i).FT(2,:)))
    impulse_save(i) = trapz((1/250)*1:length(s(i).FT),s(i).FT(2,:));
    % plot(s(i).FT(2,:))
end



c = colormap('parula');
liness = {'--','-.','-.','-.','-',};

xsize = 480;
ysize = 360;
pospos = [ 100, 100,xsize,ysize];
% pospos2 = [ 1920,1320+ysize,xsize,ysize];
% pospos3 = [ 1920,1320+ysize*2,xsize,ysize];
for ff = 1:nfields
    
    
    h = figure;
        set(h,'position',pospos)
    hold on;
    for i = plot_order
        dat = normaliseNaN(PUD(i).(fields_list{ff}),2,101);
        
        hold on;plot(dat,'linestyle',liness{2},'linewidth',2)%,'color',c(i*10,:))
    end
    grid on
    ylabel([ fields_names{ff} ' Power (W)'])
    xlabel('% stance')
    legend(leg_str)
    xlim([0 100])
    
    h = figure;
        set(h,'position',pospos)
    bar([[poswork(plot_order).(fields_list{ff})];[negwork(plot_order).(fields_list{ff})];[poswork(plot_order).(fields_list{ff})]+[negwork(plot_order).(fields_list{ff})]])
    title(['pos/neg/net work ' fields_names{ff}])
    legend(leg_str)
    
    h = figure;
        set(h,'position',pospos)
    hold on;
    for i = plot_order
        %     yyaxis left; hold on
        dat = normaliseNaN(work(i).(fields_list{ff}),2,101);
        work_SAVE(i,:) = [nanmin(dat) nanmax(dat)-nanmin(dat)];
        hold on;plot(dat,'linestyle',liness{i},'linewidth',2)%,'color',c(i*10,:))
        
    end
    
    grid on
    ylabel([ fields_names{ff} ' Work (J)'])
    xlabel('% stance')
    legend(leg_str)
    xlim([0 100])
    
    
    
    
end

%% plot the deformation of the heel pad
cmap = get(gca,'Colororder');
% figure;
% plot(s(i).ave_df); hold on;
% plot(s(i).top_df)
% leg_str = {'slow-decel','slow-const','slow-accel'};%,'norm-decel','fast-decel'};
plot_order = [3,4,2,1,5]

figure;hold on;
for i = plot_order
    
npts = length(s(i).ave_df);
d1 = normaliseNaN(s(i).top_df,2,202);
    d2 = normaliseNaN(s(i).force{2}(3,1:npts),2,202);
    energy{i} = cumtrapz(d1,d2)
plot(d1,d2)
end
xlabel('Distance (mm)')
ylabel('Force (N)')
legend(leg_str)

figure;hold on;
for i = plot_order

plot(energy{i}/1000 )
end
xlabel('Frame')
ylabel('Energy')
legend(leg_str)

figure;hold on;
for i =1:5

workE(i,:) = [-((energy{i}(1)-min(energy{i}))/1000 ) ((energy{i}(end)-min(energy{i}))/1000 )];
end
h =bar([workE(plot_order(1),:)',work_SAVE(plot_order(1),:)',workE(plot_order(2),:)',work_SAVE(plot_order(2),:)',workE(plot_order(3),:)',work_SAVE(plot_order(3),:)']);
% h =bar([workE(plot_order(1),:)',work_SAVE(plot_order(1),:)',...
%     workE(plot_order(2),:)',work_SAVE(plot_order(2),:)',...
%     workE(plot_order(3),:)',work_SAVE(plot_order(3),:)',...
%     workE(plot_order(4),:)',work_SAVE(plot_order(4),:)',...
%     workE(plot_order(5),:)',work_SAVE(plot_order(5),:)']);
for i = 1:2
    h(i).FaceColor = cmap(1,:);
end
for i = 3:4
    h(i).FaceColor = cmap(2,:);
end
for i = 5:6
    h(i).FaceColor = cmap(3,:);
end
for i = [2,4,6]
    h(i).LineStyle = '--';
end
for i =1:6
    h(i).LineWidth = 2;
end

ylabel('Energy')
legend(h([1,3,5]),leg_str)

% plot work vs impulse
figure;
plot( impulse_save(plot_order), workE(plot_order,1),'.','markersize',10);hold on
plot( impulse_save(plot_order), workE(plot_order,2),'.','markersize',10)
figure;hold on;
for i = plot_order

npts = length(s(i).top_df);
plot(s(i).top_df,(s(i).force{2}(3,1:npts)))
end
legend(leg_str)
%% plot the velocity: linear/angular
cmap2 = get(gca,'colororder');
xyz_lines= {'-','--','-.'};
bone_fields = {'cal'}%{'tib','cal','tal','mt1'};
for bf =  1:length(bone_fields)
    ind = 0;
    figure;
    for i = plot_order
        ind = ind+1;
        for j = 1:3
            hold on;
            
            dat = normaliseNaN(norm3d(s(i).w.(bone_fields{bf})(:,:)),2,101);
            
            hold on;plot(dat,'linestyle',xyz_lines{j},'linewidth',2,'color',cmap2(ind,:))
        end
    end
    title([bone_fields{bf} ' - Angular Velocity'])
    
    xlim([0 100])
    xlabel('% stance')
    
     ind = 0;
    figure;
    for i = plot_order
        ind = ind+1;
        for j = 1:3
            hold on;
            
            dat = normaliseNaN(norm3d(s(i).vcom.(bone_fields{bf})(:,:)),2,101);
            
            hold on;plot(dat,'linestyle',xyz_lines{j},'linewidth',2,'color',cmap2(ind,:))
        end
    end
    title([bone_fields{bf} ' - COM Velocity'])
    
    
    xlim([0 100])
    
    xlabel('% stance')
    

    
end

       ind = 0;
    figure;
    for i = plot_order
        ind = ind+1;
        for j = 1
            hold on;
            
            dat = normaliseNaN( s(i).ang.archdors(j,:),2,101);
            
            hold on;plot(dat,'linestyle',xyz_lines{j},'linewidth',2,'color',cmap2(ind,:))
        end
    end
    title(['arch dors'])
    
    
    xlim([0 100])
    
    xlabel('% stance')
%%

figure;
hold on;
for i = plot_order
%     yyaxis left; hold on
    dat = normaliseNaN(work(i).ST,2,101);
    dat(end) - min(dat)
    hold on;plot(dat,'linestyle',liness{2},'linewidth',2,'color',c(i*10,:))
%     
%     yyaxis right; hold on;
%     wdat = normaliseNaN( s(i).FT,2,101);
%     for j = 2
%     plot(wdat(j,:),'linestyle',liness{j},'color',c(i*10,:))
%     end
end
grid on
% yyaxis left
ylabel('subtalar Work (J)')
% ylabel('Force (N)')
% yyaxis right
% ylabel('Talus-mt1 angular velocity')
xlabel('% stance')
legend(leg_str)
xlim([0 100])

figure;
hold on;
for i = plot_order
%     yyaxis left; hold on
    dat = normaliseNaN(PUD(i).MT,2,101);
   
    hold on;plot(dat,'linestyle',liness{2},'linewidth',2,'color',c(i*10,:))
%     
%     yyaxis right; hold on;
%     wdat = normaliseNaN( s(i).FT,2,101);
%     for j = 2
%     plot(wdat(j,:),'linestyle',liness{j},'color',c(i*10,:))
%     end
end
grid on
% yyaxis left
ylabel('Midtarsal Power (W)')
% ylabel('Force (N)')
% yyaxis right
% ylabel('Talus-mt1 angular velocity')
xlabel('% stance')
legend(leg_str)
xlim([0 100])

figure;
hold on;
for i = plot_order
    dat = normaliseNaN(PUD(i).ST,2,101);
    
    hold on;plot(dat,'linestyle',liness{i},'linewidth',2)
end
grid on
ylabel('Subtalar Power (W)')
xlabel('% stance')
legend(leg_str)
xlim([0 100])


% 
% 
% figure;
% hold on;
% for i = plot_order
%     yyaxis left; hold on;
%     dat = normaliseNaN(PUD(i).tal,2,101);
%     
%     hold on;plot(dat,'linestyle',liness{i},'linewidth',2,'color',c(i*10,:))
%     
%      yyaxis right; hold on;
%     wdat = normaliseNaN( s(i).w.tal,2,101);
%     for j = 1:3
%     plot(wdat(j,:),'linestyle',liness{j},'color',c(i*10,:))
%     end
%     
%     disp((struct_data(i).cropFrsXray(2)-struct_data(i).cropFrsXray(1))/250);
% end
% grid on
% ylabel('Tal Power (W)')
% xlabel('% stance')
% legend('slow','pref','fast')
% xlim([0 100])

figure;
hold on;
for i = plot_order
%     yyaxis left; hold on;
    dat = normaliseNaN(PUD(i).cal,2,101);
    
    hold on;plot(dat,'linestyle',liness{i},'linewidth',2,'color',c(i*10,:))
%     
%      yyaxis right; hold on;
%     wdat = normaliseNaN( s(i).w.tal,2,101);
%     for j = 1:3
%     plot(wdat(j,:),'linestyle',liness{j},'color',c(i*10,:))
%     end
%     
%     disp((struct_data(i).cropFrsXray(2)-struct_data(i).cropFrsXray(1))/250);
end
grid on
ylabel('Cal Power (W)')
xlabel('% stance')
legend(leg_str)
xlim([0 100])
figure;
hold on;
for i = plot_order
    dat = normaliseNaN(PUD(i).mt1,2,101);
    
    hold on;plot(dat,'linestyle',liness{i},'linewidth',2)
    
    disp((struct_data(i).cropFrsXray(2)-struct_data(i).cropFrsXray(1))/250);
end
grid on
ylabel('mt1 Power (W)')
xlabel('% stance')
legend(leg_str)
xlim([0 100])


figure;
hold on;
for i = plot_order
    dat = normaliseNaN(PUD(i).caldist,2,101);
    
    hold on;plot(dat,'linestyle',liness{i},'linewidth',2)
end
grid on
ylabel('Distal calcaneus Power (W)')
xlabel('% stance')
legend(leg_str)
xlim([0 100])


figure;
hold on;
for i = plot_order
    dat = PUD(i).caldist;
    
    hold on;plot(0:1/250:((length(dat)-1)/250),dat,'linestyle',liness{i},'linewidth',2)
end
grid on
ylabel('Distal calcaneus Power (W)')
xlabel('time (s)')
legend(leg_str)

figure;
hold on;
for i = plot_order
    dat = s(i).w.tib;
    for j = 1:3
    
         hold on;plot(0:1/250:((length(dat(j,:))-1)/250),dat(j,:)','color',c(15*j,:),'linestyle',liness{i},'linewidth',2)
    end
end
grid on
ylabel('tib angular velocity (rad/s)')
xlabel('time (s)')
legend('slow X','slow Y','slow Z','pref X','pref Y','pref Z','fast X','fast Y','fast Z')

figure;
hold on;
for i = plot_order
    dat = s(i).w.cal;
    for j = 1:3
    
         hold on;plot(0:1/250:((length(dat(j,:))-1)/250),dat(j,:)','color',c(15*j,:),'linestyle',liness{i},'linewidth',2)
    end
end
grid on
ylabel('calangular velocity (rad/s)')
xlabel('time (s)')
legend('slow X','slow Y','slow Z','pref X','pref Y','pref Z','fast X','fast Y','fast Z')


figure;
hold on;
for i = plot_order
    dat = normaliseNaN(PUD(i).tib,2,101);
    
    hold on;plot(dat,'linestyle',liness{i},'linewidth',2)
    
    disp((struct_data(i).cropFrsXray(2)-struct_data(i).cropFrsXray(1))/250);
end
grid on
ylabel('tib Power (W)')
xlabel('% stance')
legend(leg_str)
xlim([0 100])


figure;
hold on;
for i = plot_order
    dat = PUD(i).tib;
    
    hold on;plot(0:1/250:((length(dat)-1)/250),dat,'linestyle',liness{i},'linewidth',2)
    
    disp((struct_data(i).cropFrsXray(2)-struct_data(i).cropFrsXray(1))/250);
end
grid on
ylabel('tib Power (W)')
xlabel('time (s)')
legend(leg_str)


    title_str{1} = 'Sagittal axis';
    title_str{8} = 'Most inclined (through talar neck)';
    title_str{1} = 'Using instantaneous helical (skip by 4 frames)';
for jj = 1
    h = figure;
    set(h,'position',[680   100   560   820])
    for i = plot_order
        dat = normaliseNaN(s(i).M_MT_ax_calc,2,101);
        for j = 1:3
            hold on;plot(dat(j,:),'color',c(20*j,:),'linestyle',liness{i},'linewidth',2)
        end
        
        
        
        disp((struct_data(i).cropFrsXray(2)-struct_data(i).cropFrsXray(1))/250);
    end
    xlabel('% Stance')
    ylabel('Midtarsal Moment (Nm)')
    xlim([0 100])
%     ylim([-125 200])
    h=legend('+dors-slow','+int rot-slow','+abduction-slow','+dors-pref','+int rot-pref','+abduction-pref','+dors-fast','+int rot-fast','+abduction-fast');
    set(h,'location','southeast')
    title(title_str{jj})
    grid on
end
print(['C:\Users\Lauren\Documents\School\PhD\Research\AllStudies\Mechanisms\Moment Graphs\' title_str{jj} ],'-dtiffn')



%% animate the force vector in Wrist visualizer

c = [195,48,48;
    0, 153, 153;
    51,51,255;
    153,255,154]/256;
% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
ivdir = fullfile(subj_dir,'Models','IV');
for i = 1:ntrials
    anim_dir = fullfile(subj_dir,trialname_list{i},'POS','ForcesAnimation',filesep);
    rigidiv_folder = fullfile(anim_dir,'rigidiv',filesep);
    
    % set up the GRF animation folder
    grf_style = 'grf_P%i_F1.iv';
    grf_folder = fullfile(anim_dir,'GRFdir',filesep);
    
    if exist(anim_dir,'dir')==0;     mkdir(anim_dir);     end
    if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end
    if exist(grf_folder,'dir')==0;  mkdir(grf_folder);  end
    
    %     first_fr = struct_data(i).cropFrsXray(1);
    %     end_frM   = struct_data(i).cropFrsXray(2);
    
    for bn = 1:nBones
        ivstring = createInventorHeader();
        % make the linked iv file
        ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list_char(bn,:) '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
        
        fid = fopen(fullfile(rigidiv_folder,[bone_list_char(bn,:) '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
    
    for bn = 1:nBones
        %         for fr = first_fr:end_frM
        %         Tw.(bone_list{bn})(:,:,fr) = struct_data(i).Tm.(bone_list{bn})(:,:,fr);
        %         end
        Tw = s(i).T.(bone_list{bn});
        Tw(isnan(Tw)) = 1;
        write_RTp(bone_list_char(bn,:) , Tw, anim_dir)
    end
    
    
    %     frX = 1;
    for fr = 1:length(struct_data(i).cropFrsForce(1):struct_data(i).cropFrsForce(2) )
        
        grf_filename = sprintf(grf_style,fr);
        %         frX = frX+1;
        ivstring = createInventorHeader();
        
        
        %         ************** HELICAL AXES **************************
        for bp = 1%:nBP
            if any(isnan(s(i).heli(bp).s_m(1,fr)),'all')
                continue
            end
            ivstring = [ivstring createInventorArrow( s(i).heli(bp).s_ref_moc(:,fr)-50*s(i).heli(bp).n_m(:,fr),...
                s(i).heli(bp).n_m(:,fr),...
                s(i).heli(bp).phi(fr)*250,...
                2,[1,0, 0.2],0)];
        end
        
        %*********FORCE VECTORS****************
        for fp = 1:2 % each force plate
            COP = s(i).COP{fp}(:,fr) *1000;% transformPoints(coreg,struct_data(i).force_data(fp).globCOP(:,fr)*1000,-1);
            force = s(i).force{fp}(:,fr);%transformVectors(coreg,struct_data(i).force_data(fp).globForce(:,fr),-1);
            
            %         ivstring = [ivstring createInventorSphere(COP,4,[0.2 0.2 0.8],0)];
            ivstring = [ivstring createInventorArrow(COP,force,norm(force)/5,1.5,[1 1 1],0)];
        end
            
            COP = s(i).cop_new(:,fr) *1000;% transformPoints(coreg,struct_data(i).force_data(fp).globCOP(:,fr)*1000,-1);
            force = s(i).FT(:,fr);%transformVectors(coreg,struct_data(i).force_data(fp).globForce(:,fr),-1);
            
            %         ivstring = [ivstring createInventorSphere(COP,4,[0.2 0.2 0.8],0)];
            ivstring = [ivstring createInventorArrow(COP,force,norm(force)/5,2,[.95,.95,.95],0)];

        
        %         *************** MOMENT ARMS **********************
        %
%         for j = 1:8
%             
%             if any(isnan(s(i).axis_vec_test{j}(1,fr)),'all') 
%                 continue
%                
%             end
%             ivstring = [ivstring createInventorArrow( s(i).JC_test{j}(:,fr)*1000,...
%                 s(i).axis_vec_test{j}(:,fr),...
%                 100,...
%                 2,[0,1, 0.5],0)];
%             ivstring = [ivstring createInventorArrow( s(i).JC_test{j}(:,fr)*1000,...
%                 -s(i).axis_vec_test{j}(:,fr),...
%                 100,...
%                 2,[0,1, 0.5],0)];
%             
%             if any(isnan([s(i).s_g_MT1test{j}(:,fr); s(i).r_g_MT1test{j}(:,fr)]),'all') || (s(i).force{1}(3,fr) < 10)
%             else
%                 ivstring = [ivstring createInventorArrow( s(i).s_g_MT1test{j}(:,fr)*1000,...
%                     s(i).r_g_MT1test{j}(:,fr)*1000,...
%                     norm(s(i).r_g_MT1test{j}(:,fr)*1000),...
%                     1,[0 1 0],0)];
%             end
%             
%             if any(isnan([s(i).s_g_MT2test{j}(:,fr); s(i).r_g_MT2test{j}(:,fr)]),'all') || (s(i).force{2}(3,fr) < 10) 
%             else
%                 ivstring = [ivstring createInventorArrow( s(i).s_g_MT2test{j}(:,fr)*1000,...
%                     s(i).r_g_MT2test{j}(:,fr)*1000,...
%                     norm(s(i).r_g_MT2test{j}(:,fr)*1000),...
%                     1,[0 1 0],0)];
%             end
%         end
%         
        
        
        
        %
        %         if any(isnan([s(i).s_g_MT1ax(:,fr); s(i).r_g_MT1ax(:,fr)]),'all')
        %         else
        %             ivstring = [ivstring createInventorArrow( s(i).s_g_MT1ax(:,fr)*1000,...
        %                 s(i).r_g_MT1ax(:,fr)*1000,...
        %                 norm(s(i).r_g_MT1ax(:,fr)*1000),...
        %                 1,[1,0, 0],0)];
        %         end
        %
        %         if any(isnan([s(i).s_g_MT2ax(:,fr); s(i).r_g_MT2ax(:,fr)]),'all')
        %         else
        %             ivstring = [ivstring createInventorArrow( s(i).s_g_MT2ax(:,fr)*1000,...
        %                 s(i).r_g_MT2ax(:,fr)*1000,...
        %                 norm(s(i).r_g_MT2ax(:,fr)*1000),...
        %                 1,[1,0, 0],0)];
        %         end
        %************* WRITE FILE WITH ARROWS *****************
        fid = fopen(fullfile(grf_folder,grf_filename),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        
    end
    create_ini(0,0,1,1,grf_folder(1:end-1),'grf_P%d_F%d.iv',fullfile(anim_dir,[trialname_list{i} '.ini']))
    
    pos_text = write_pos(bone_list,anim_dir,trialname_list{i});
    
    filename = fullfile(anim_dir, [trialname_list{i} '.pos']);
    
    fid = fopen(filename,'w'); % open the file to write
    fprintf(fid,pos_text);
    fclose(fid);
    
    fprintf('Animations are saved in %s.\n',anim_dir)
    
    %
    %     fr = 100;
    %     figure;
    %     hold on;
    %     pts = struct_bone(1).pts;
    %     cns = struct_bone(1).cns;
    %     patch('faces',cns(1:10:end,1:3),'vertices',transformPoints(Tw.(bone_list{1})(:,:,fr),pts) ,'facealpha',0.3)
    %
    %     pts = struct_bone(5).pts;
    %     cns = struct_bone(5).cns;
    %     patch('faces',cns(1:10:end,1:3),'vertices',transformPoints(Tw.(bone_list{5})(:,:,fr),pts) ,'facealpha',0.3)
    %     fr = fr+251
    %     for fp = 1:2
    %     plotvector3(struct_data(i).force_data(fp).globCOP(:,fr)*1000,struct_data(i).force_data(fp).globForce(:,fr))
    %     end
    %
    
    
end

