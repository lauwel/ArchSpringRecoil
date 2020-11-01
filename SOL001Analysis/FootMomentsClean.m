
% code to calculate the midtarsal and subtalar moments, using the various
% methods
close all
clear
clc

subj_name = 'SOL001B';
subj_dir = 'E:\SOL001_VISIT2\';
% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);


load([subj_dir 'struct_bone.mat']);
load([subj_dir 'struct_data.mat']);

bone_list = {struct_bone(:).bone_name}';
nBones = length(bone_list);

% get the trials that we want
trialname_list = { 'T0019_SOL001_nrun_pref_barefoot';...
    'T0024_SOL001_nrun_rfs_barefoot';%...
    'T0033_SOL001_srun_pref_barefoot';...
    'T0034_SOL001_srun_pref_barefoot'; ...
    %                 'T0030_SOL001_srun_pref_barefoot';...
    'T0039_SOL001_frun_pref_barefoot'};
ntrials = length(trialname_list);

for i = 1:ntrials
    ind_file(i) = findInStruct(struct_data,'filename',trialname_list{i});
end

struct_data = struct_data(ind_file);

% total number of frames for force and xray
nFrX = 400;
nFrF = 1600;

% framerates for xray and force
frRateX = 250;
frRateF = 1000;
% filtering frame rate
% wc = [10,30];
fcF = [40 80];

% markers to model the sensitivity of the midtarsal axis *** may have
% changed
cub_marker = load([subj_dir  'Models\IV\Landmarks\cub_marker.stack'])';
nav_marker = load([subj_dir  '\Models\IV\Landmarks\nav_marker.stack'])';

for i = 1:4
    nav_CT{i} = load([subj_dir  '\Models\IV\Landmarks\',sprintf('nav%i.stack',i)])';
    tal_CT{i} = load([subj_dir  '\Models\IV\Landmarks\',sprintf('tal%i.stack',i)])';
end

coreg = csvread(fullfile(subj_dir,'Calibration','Set 1','end_pylon','end_pylon_COREG.csv'));

% measure the floor location

floor_raw = readmatrix([subj_dir  '\Calibration\Set 1\FloorPosition_fromXMA.csv']);
floor_3xn = reshape(floor_raw,3,6);
floor_mocap = transformPoints(coreg,floor_3xn);
% fit a plane to the floor
[nfl,V,p] = affine_fit(floor_mocap');
floor_mocap = [floor_mocap'; V(1:3,1)'*400; V(1:3,2)'*400; V(1:3,1)'*400+ V(1:3,2)'*400;-V(1:3,1)'*400; -V(1:3,1)'*400+ V(1:3,2)'*400];

%% calculate landmark position
% load all the motion capture data


% load the bone transforms + filter the data
for i = 1:ntrials
    
    
    for bn = 1:nBones
        
        for fr = 1:nFrX
            % move the cuboid and navicular markers conventional in motion
            % capture
            if strcmp('cub',bone_list{bn})
                struct_data(i).landmarks.(bone_list{bn})(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr),cub_marker/1000);
            elseif strcmp('nav',bone_list{bn})
                for j = 1:4
                    struct_data(i).landmarks.(bone_list{bn}){j}(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr),nav_CT{j}/1000);
                end
            elseif strcmp('tal',bone_list{bn})
                for j = 1:4
                    struct_data(i).landmarks.(bone_list{bn}){j}(:,fr) = transformPoints(struct_data(i).Tm.(bone_list{bn})(:,:,fr),tal_CT{j}/1000);
                end
            end
        end
        
        
    end
    
    
end





%% Compute the helical axes between met-calc, met-tal,cal-tal and tal-tib
close all
cmap = colormap('parula');
bone_pairs = {'cal','mt1'};%;'tal','mt1';'tal','cal';'tib','tal'};
nBP = size(bone_pairs,1);
for i = 1:ntrials
    for bp = 1:nBP
        
        T_mR_ct = struct_data(i).Tm.(bone_pairs{bp,1});
        T_mA_ct = struct_data(i).Tm.(bone_pairs{bp,2});
        
        spc = 1; % how many frames to skip
        nFrs = size(T_mR_ct,3);
        [phi_raw,n_raw,L_raw,s_raw] = helicalInstantaneous(T_mR_ct(:,:,1:spc:nFrs),T_mA_ct(:,:,1:spc:nFrs));
        %         phi = phi_raw;
        %         n = n_raw;
        %         L = L_raw;
        %         s = s_raw;
        %
        %
        %
        %         n = n./norm3d(n);
        
        
        
        %          [phi_raw,n_raw,L_raw,s_raw] = stabilizeHelicalAxis(T_mR_ct,T_mA_ct,spc);
        %
        %                 fc = [30,40];
        %         n = adaptiveLowPassButterworth(n_raw,fc,250);
        %         s = adaptiveLowPassButterworth(s_raw,fc,250);
        %         phi = adaptiveLowPassButterworth(phi_raw,fc,250);
        %         L = adaptiveLowPassButterworth(L_raw,fc,250);
        
        %         [phi,n,L,s] =  helicalInstantaneousGapFill(T_mR_ct,T_mA_ct,spc);
        
        %         figure;
        %         hold on; plot(n_raw','-.');plot(n')%(:,1:spc:end)')
        %         figure;
        %         hold on; plot(s_raw','-.');plot(s')%(:,1:spc:end)')
        %          figure;
        %         hold on; plot(phi_raw','-.');plot(phi')%(:,1:spc:end)')
        
        n_m = transformVectors(T_mR_ct,n_raw);
        s_m = transformPoints(T_mR_ct,s_raw);
        bn_ind = findInStruct(struct_bone,'bone_name',bone_pairs{bp,1});
        for fr = 1:nFrX-1
            centR_moc = transformPoints(T_mR_ct(:,:,fr),struct_bone(bn_ind).cent/1000);
            s_ref_moc(:,fr) = closestPointonVector(centR_moc,s_m(:,fr),n_m(:,fr));
            
            s_ref_ct(:,fr) = transformPoints(T_mR_ct(:,:,fr), s_ref_moc(:,fr) ,-1);
        end
        
        % find the average helical axis by dividing into groups by the
        % magnitude of phi
        ind_nan = isnan(phi_raw);
        ind_big_phi = phi_raw>0.4;
        ind_small_phi = ~ind_big_phi;
        avg_ind = findGroupsOfLogicals(ind_big_phi);% find the values to average over
        int_ind = findGroupsOfLogicals(ind_small_phi & ~ind_nan);% find the indices that will be interpolated over
        
        n_avg = nan(3,nFrs);
        s_avg = nan(3,nFrs);
        for gi = 1:size(avg_ind,1)
            ii = false(nFrs,1);
            ii(avg_ind(gi,1):avg_ind(gi,2),1) = 1;
            [avgHAM(gi,:),num,variat] = calcAvgMotionLW( [phi_raw(ii)',n_raw(:,ii)',L_raw(ii)',s_ref_ct(:,ii)']);
            for dd = 1:3
                n_avg(dd,ii) = avgHAM(gi,dd+1)';
                s_avg(dd,ii) = avgHAM(gi,dd+5)';
            end
            
            if mean(n_avg(1,ii)) < 0
                n_avg(:,ii) = -n_avg(:,ii);
            end
        end
        figure; plot(n_raw'); hold on; plot(n_avg');
        for gi = 1:size(int_ind,1)
            p1 = int_ind(gi,1);
            p2 = int_ind(gi,2);
            if p2-p1 < 3
                p1 = p1-1;
                p2 = p2+2;
                n_avg(:,p1:p2)= NaN;
                s_avg(:,p1:p2) = NaN;
            end
            
            for dd = 1:3
                n_avg(dd,p1:p2) = interp1(find(ind_big_phi),n_avg(dd,ind_big_phi),p1:p2,'spline'); 
                s_avg(dd,p1:p2) = interp1(find(ind_big_phi),s_avg(dd,ind_big_phi),p1:p2,'spline');
            end
        end
        plot(n_avg')

        
        n_avg_m = transformVectors(T_mR_ct,n_avg);
        s_avg_m = transformPoints(T_mR_ct,s_avg);
        
        
        
        figure
        hold on;
        for dd = 1:3
            plot(n_raw(dd,:)','color',cmap(dd*15,:));  plot(n_avg(dd,:)','--','color',cmap(dd*15,:))
        end
        figure; plot(n_avg_m')
        %
        %         yyaxis right
        %         plot(phi_raw)
        %
        figure
        hold on;
        for dd = 1:3
            plot(s_ref_ct(dd,:)','color',cmap(dd*15,:));  plot(s_avg(dd,:)','--','color',cmap(dd*15,:))
        end
        %
        struct_helical(i).trials(bp).refBone = bone_pairs{bp,1};
        struct_helical(i).trials(bp).bone = bone_pairs{bp,2};
        struct_helical(i).trials(bp).phi = phi_raw;
        struct_helical(i).trials(bp).n = n_raw;
        struct_helical(i).trials(bp).L = L_raw;
        struct_helical(i).trials(bp).s = s_raw;
        
        %
        struct_helical(i).trials(bp).n_m = n_m;
        struct_helical(i).trials(bp).s_m = s_ref_moc;
        %         struct_helical(i).trials(bp).s_ref_moc = s_ref_moc;
        struct_helical(i).trials(bp).s_avg_m = s_avg_m;
        struct_helical(i).trials(bp).n_avg_m = n_avg_m;
        
        
    end
    
    
end


%% animate the force vector and helical axes in Wrist visualizer

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
        ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
        
        fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '.iv']),'w');
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
        write_RTp(bone_list{bn} , Ta.(bone_list{bn})(:,:,first_fr:end_frM) , anim_dir)
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


clearvars('s')

for i = 1:ntrials
    
    new_centre = [-0.2 0 0]';
    
    % determine where the force is applied to the foot
    fploc1 = struct_data(i).force_data(1).ForcePlateLocation;
    fploc2 = struct_data(i).force_data(2).ForcePlateLocation;
    [floor_vec1,~,p1] = affine_fit(fploc1);
    [floor_vec2,~,p2] = affine_fit(fploc2);
    
    warning('adjust top plate location')
    
    p1 = p'/1000;%p1'+[0,0,0.0531]';
    p2 = p'/1000;%p2'+[0,0,0.0531]';
    floor_vec1 = nfl;
    floor_vec2 = nfl;
    
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
        s(i).heli(bp).s_m =    struct_helical(i).trials(bp).s_m(:,start_frX:end_frX);
        s(i).heli(bp).phi =    struct_helical(i).trials(bp).phi(start_frX:end_frX);
        s(i).heli(bp).L =    struct_helical(i).trials(bp).L(start_frX:end_frX);
        s(i).heli(bp).n_avg_m =    struct_helical(i).trials(bp).n_avg_m(:,start_frX:end_frX);
        s(i).heli(bp).s_avg_m =    struct_helical(i).trials(bp).s_avg_m(:,start_frX:end_frX);
    end
    for bn = 1:nBones
        if ~isfield(struct_data(i).Tm,(bone_list{bn}))
            continue
        end
        s(i).Tm.(bone_list{bn}) = struct_data(i).Tm.(bone_list{bn})(:,:, start_frX :end_frX );
        s(i).w.(bone_list{bn}) = struct_data(i).w.(bone_list{bn})(:, start_frX :end_frX );
        s(i).vcom.(bone_list{bn}) = struct_data(i).vcom.(bone_list{bn})(:,start_frX :end_frX );
        s(i).pcom.(bone_list{bn}) = struct_data(i).pcom.(bone_list{bn})(:, start_frX :end_frX );
    end
    
    
    [s(i).ang.archdors,s(i).ang.archabd,s(i).ang.archabd] = eulerYZX( s(i).Tm.cal,  s(i).Tm.mt1,struct_bone(1).T_ACS, struct_bone(5).T_ACS);
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
    
    % frame references (indices)
    nFrs = (end_frM - start_frM + 1);
    frWA = round([nFrs*0.2 nFrs*0.40]); % weight acceptance
    frPO = round([nFrs*0.60 nFrs * 0.8]); % push off
    frMP = round(nFrs/2); % the midpoint
    FT = [];
    momT = [];
    r1_com = [];
    r2_com = [];
    M_sum = [];
    cop_new= [];
    phi_cum = [];
    w_MT = [];
    w_axis = [];
    v_axis =[];
    P_MT_axis_ave = [];
    P_MT_axis_inst = [];
    P_MT_conv = [];
    % arch axes for weight acceptance and push-off - in GLOBAL
    n_calc = transformVectors(s(i).Tm.cal,s(i).heli.n_m,-1); % calcaneal coordinate system
    s_calc = transformPoints(s(i).Tm.cal,s(i).heli.s_m,-1); % calcaneal coordinate system
    
    arch_ax_WA = mean(n_calc(:,frWA(1):frWA(2)),2); % average
    arch_ax_PO = -mean(n_calc(:,frPO(1):frPO(2)),2);
    
    arch_pt_WA = mean(s_calc(:,frWA(1):frWA(2)),2); % average
    arch_pt_PO = mean(s_calc(:,frPO(1):frPO(2)),2);
    
    
    %     arch_pt_sd_WA = distBtwPoints3xN(s_calc(:,frWA(1):frWA(2)),repmat(arch_pt_WA,1,frWA(2)-frWA(1)+1))
    % now drive with calcaneal kinematics
    % n vector
    s(i).arch_ax_ave_calc(:,1:frWA(2)) = transformVectors(s(i).Tm.cal(:,:,1:frWA(2)), arch_ax_WA);
    s(i).arch_ax_ave_calc(:,frPO(1):nFrs) = transformVectors(s(i).Tm.cal(:,:,frPO(1):nFrs), arch_ax_PO);
    arch_interp = [];
    for j = 1:3
        arch_interp(j,:) = interp1([1:frWA(2),frPO(1):nFrs], s(i).arch_ax_ave_calc(j,[1:frWA(2),frPO(1):nFrs]),1:nFrs,'spline');
    end
    %     figure; plot([1:frWA(2),frPO(1):nFrs], s(i).arch_ax_ave_calc(:,[1:frWA(2),frPO(1):nFrs])'); hold on; plot(1:nFrs,arch_interp)
    
    s(i).arch_ax_ave = arch_interp;
    
    % point
    s(i).arch_pt_ave_calc(:,1:frWA(2)) = transformPoints(s(i).Tm.cal(:,:,1:frWA(2)), arch_pt_WA);
    s(i).arch_pt_ave_calc(:,frPO(1):nFrs) = transformPoints(s(i).Tm.cal(:,:,frPO(1):nFrs), arch_pt_PO);
    arch_pt_interp = [];
    for j = 1:3
        arch_pt_interp(j,:) = interp1([1:frWA(2),frPO(1):nFrs], s(i).arch_pt_ave_calc(j,[1:frWA(2),frPO(1):nFrs]),1:nFrs,'spline');
    end
    %     figure; plot([1:frWA(2),frPO(1):nFrs], s(i).arch_pt_ave_calc(:,[1:frWA(2),frPO(1):nFrs])'); hold on; plot(1:nFrs,arch_pt_interp)
    s(i).arch_pt_ave = arch_pt_interp;
    
    
    
    % interpolate the transition point for the instantaneous axis as well
    s(i).arch_ax_inst = [s(i).heli.n_m(:,1:frWA(2)),nan(3,frPO(1)-frWA(2)-1),  -s(i).heli.n_m(:,frPO(1):nFrs)];
    arch_interp = [];
    for j = 1:3
        arch_interp(j,:) = interp1([1:frWA(2),frPO(1):nFrs], s(i).arch_ax_inst(j,[1:frWA(2),frPO(1):nFrs]),1:nFrs,'spline');
    end
    
    s(i).arch_ax_inst(:,frWA(2):frPO(1)) = arch_interp(:,frWA(2):frPO(1)) ;
    %     figure; subplot(2,1,1); plot([1:frWA(2),frPO(1):nFrs],[ s(i).heli.n_m(:,1:frWA(2)),  -s(i).heli.n_m(:,frPO(1):nFrs)],'.')
    %     hold on; plot(1:nFrs, s(i).arch_ax_inst)
    %
    
    
    s(i).arch_pt_inst = [s(i).heli.s_m(:,1:frWA(2)),nan(3,frPO(1)-frWA(2)-1),  s(i).heli.s_m(:,frPO(1):nFrs)];
    arch_pt_interp = [];
    for j = 1:3
        arch_pt_interp(j,:) = interp1([1:frWA(2),frPO(1):nFrs], s(i).arch_pt_inst(j,[1:frWA(2),frPO(1):nFrs]),1:nFrs,'spline');
    end
    
    s(i).arch_pt_inst(:,frWA(2):frPO(1)) = arch_pt_interp(:,frWA(2):frPO(1)) ;
    %      subplot(2,1,2);plot([1:frWA(2),frPO(1):nFrs],[ s(i).heli.s_m(:,1:frWA(2)),  s(i).heli.s_m(:,frPO(1):nFrs)],'.')
    %     hold on on; plot(1:nFrs, s(i).arch_pt_inst)
    
    s(i).arch_pt_conv = (s(i).cub +  s(i).marker2{2} )/2 ;
    v_MT_ave =  calculateVelocity(s(i).arch_pt_ave,frRateX);
    v_MT_conv =  calculateVelocity(s(i).arch_pt_conv,frRateX);
    v_MT_inst=  calculateVelocity(s(i).arch_pt_inst,frRateX);
    %*********** midtarsal joint calculated conventionally + around the axis **********
    for fr = 1:nFrs
        
        %         moments  using fixed cuboid- navicular axis (per Bruening)
        
        
        
        
        % find the two closest points (i.e. the points that connect the
        % moment arm
        
        % Use an array of axes to see how it affects the MT moment
        
        for j = 1:8
            % joint centre
            midtarsal_JC_a = nanmean([ s(i).cub(:,fr) , s(i).marker2{j}(:,fr) ], 2);
            
            axis_vec = diff([ s(i).cub(:,fr) ,  s(i).marker2{j}(:,fr)  ],[],2);
            axis_vec = axis_vec/norm(axis_vec);
            [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  ,s(i).marker2{j}(:,fr) , ...
                s(i).force{1}(:,fr) ,  axis_vec);
            
            [sa2,sb2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  ,s(i).marker2{j}(:,fr) , ...
                s(i).force{2}(:,fr) ,  axis_vec);
            
            r_g_MT1 = sa1-sb1;
            r_g_MT2 = sa2-sb2;
            s(i).s_g_MT1_conv_axis{j}(:,fr) = sb1;
            s(i).s_g_MT2_conv_axis{j}(:,fr) = sb2;
            s(i).r_g_MT1_conv_axis{j}(:,fr) = r_g_MT1;
            s(i).r_g_MT2_conv_axis{j}(:,fr) = r_g_MT2;
            s(i).axis_vec_test{j}(:,fr) = axis_vec;
            s(i).JC_test{j}(:,fr) = midtarsal_JC_a;
            
            
            
            s(i).M_MT_conv_axis{j}(:,fr) = cross(r_g_MT1,s(i).force{1}(:,fr)) + s(i).freemom{1}(:,fr) + cross(r_g_MT2,s(i).force{2}(:,fr)) + s(i).freemom{2}(:,fr);
            s(i).M_MT_conv_axis_calc{j}(:,fr) = transformVectors( s(i).Tm.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_conv_axis{j}(:,fr),-1);
            
        end
        
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
        
        
        
        % Using the joint centre from the marker points (conventional)
        
        midtarsal_JC = nanmean([ s(i).cub(:,fr) ,  s(i).marker2{2}(:,fr)  ],2);
        
        nav_cub_vec = diff([ s(i).cub(:,fr) ,  s(i).marker2{2}(:,fr)  ],[],2);
        nav_cub_vec = nav_cub_vec/norm(nav_cub_vec);
        
        
        r_g1 =  s(i).COP{1}(:,fr) - midtarsal_JC ;
        r_g2 =  s(i).COP{2}(:,fr) - midtarsal_JC ;
        %         sv1(:,fr) = cross(r_g1,s(i).force{1}(:,fr));
        s(i).M_MT_conv_pt(:,fr) = cross(r_g1,s(i).force{1}(:,fr)) + s(i).freemom{1}(:,fr) + cross(r_g2,s(i).force{2}(:,fr)) + s(i).freemom{2}(:,fr);
        
        s(i).M_MT_conv_pt_calc(:,fr) = transformVectors( s(i).Tm.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_conv_pt(:,fr),-1);
        
        
        
        
        % using the average instantaneous arch axis -
        bp = 1;
        
        
        
%         [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  ,s(i).arch_pt_ave(:,fr), ...%s(i).arch_pt_inst(:,fr),...%
%             s(i).force{1}(:,fr) , s(i).arch_ax_ave(:,fr));
                 [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  ,s(i).heli.s_avg_m(:,fr), ...%s(i).arch_pt_inst(:,fr),...%
                    s(i).force{1}(:,fr) , s(i).heli.n_avg_m(:,fr));
        r_g_MT1 = sa1-sb1;
%         [sa2,sb2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  , s(i).arch_pt_ave(:,fr), ...%s(i).arch_pt_inst(:,fr),...%
%             s(i).force{2}(:,fr) , s(i).arch_ax_ave(:,fr));
                [sa2,sb2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  , s(i).heli.s_avg_m(:,fr), ...%s(i).arch_pt_inst(:,fr),...%
                    s(i).force{2}(:,fr) , s(i).heli.n_avg_m(:,fr));
        r_g_MT2 = sa2-sb2;
        
        s(i).s_g_MT1ax_ave(:,fr) = sb1;
        s(i).s_g_MT2ax_ave(:,fr) = sb2;
        s(i).r_g_MT1ax_ave(:,fr) = r_g_MT1;
        s(i).r_g_MT2ax_ave(:,fr) = r_g_MT2;
        s(i).M_MT_ax_ave(:,fr) = cross(r_g_MT1,s(i).force{1}(:,fr)) + cross(r_g_MT2,s(i).force{2}(:,fr)) + s(i).freemom{1}(:,fr) + s(i).freemom{2}(:,fr);
        s(i).M_MT_ax_ave_calc(:,fr) = transformVectors( s(i).Tm.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_ax_ave(:,fr),-1);
        
        
        
        % use the instantaneous axis
        
        %         [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  , s(i).arch_pt_inst(:,fr), ...
        %                                            s(i).force{1}(:,fr) ,   nav_cub_vec (:,fr));
        [sa1,sb1] = closestPointsBtw2Lines(s(i).COP{1}(:,fr)  , s(i).arch_pt_inst(:,fr), ...
            s(i).force{1}(:,fr) ,  s(i).arch_ax_inst(:,fr));
        r_g_MT1 = sa1-sb1;
        [sa2,sb2] = closestPointsBtw2Lines(s(i).COP{2}(:,fr)  , s(i).arch_pt_inst(:,fr), ...
            s(i).force{2}(:,fr) , s(i).arch_ax_inst(:,fr));
        r_g_MT2 = sa2-sb2;
        
        s(i).s_g_MT1ax_inst(:,fr) = sb1;
        s(i).s_g_MT2ax_inst(:,fr) = sb2;
        s(i).r_g_MT1ax_inst(:,fr) = r_g_MT1;
        s(i).r_g_MT2ax_inst(:,fr) = r_g_MT2;
        s(i).M_MT_ax_inst(:,fr) = cross(r_g_MT1,s(i).force{1}(:,fr)) + cross(r_g_MT2,s(i).force{2}(:,fr)) + s(i).freemom{1}(:,fr) + s(i).freemom{2}(:,fr);
        s(i).M_MT_ax_inst_calc(:,fr) = transformVectors( s(i).Tm.tal(:,:,fr) * struct_bone(9).T_ACS, s(i).M_MT_ax_inst(:,fr),-1);
        
        
        w_MT(:,fr) = s(i).w.cal(:,fr)-s(i).w.mt1(:,fr);
        
        
        F_sum = (s(i).force{1}(:,fr) + s(i).force{2}(:,fr));
        w_axis(:,fr) =  s(i).heli.n_m(:,fr).*s(i).heli.phi(:,fr)*pi()/180 *frRateX;
        v_axis(:,fr) = s(i).heli.n_m(:,fr).*s(i).heli.L(:,fr)*frRateX;
        
        P_MT_axis_ave(:,fr) = dot(-s(i).M_MT_ax_ave(:,fr),w_axis(:,fr)) + dot(-F_sum,v_axis(:,fr));
        P_MT_conv(:,fr) = dot(-s(i).M_MT_conv_axis{2}(:,fr),w_axis(:,fr)) + dot(-F_sum,v_axis(:,fr));
        P_MT_axis_inst(:,fr) = dot(-s(i).M_MT_ax_inst(:,fr),w_axis(:,fr)) + dot(-F_sum,v_axis(:,fr));
        
    end
    %     fields_list = {'MT','ST','cal','tib','caldist','TC','tal','tal1','tal2','tal12','calO','mt1O','mt1','MTtemp'};
    
    ind_nonan = ~isnan(s(i).heli.phi);
    phi_cum = s(i).heli.phi;
    phi_cum(ind_nonan) = cumtrapz(s(i).heli.phi(ind_nonan));
    s(i).phi_cum = phi_cum;
    s(i).w_MT = w_MT;
    s(i).w_axis = w_axis;
    s(i).theta_MT = nan(3,length(phi_cum));
    s(i).theta_MT = phi_cum .* s(i).heli.n_m;
    s(i).P.MT_axis_ave = P_MT_axis_ave;
    s(i).P.MT_axis_conv = P_MT_conv;
    s(i).P.MT_axis_inst = P_MT_axis_inst;
    
    s(i).FT = FT;
    s(i).momT = momT;
    s(i).cop_new = cop_new;
    
    
    
end

cmap = colormap('parula');


% Plot the midtarsal power the two different ways
figure;
hold on;
cmap_spec = cmap([10,10,20,20,30],:);
for i = 1:ntrials
     h(1,i)=plot(linspace(0,100,length(s(i).M_MT_conv_pt_calc(1,:))), s(i).P.MT_axis_ave','color',cmap_spec(i,:),'linestyle','-');
     
    pow_min(1,i) = min( s(i).P.MT_axis_ave(10:end));
    pow_max(1,i) = max( s(i).P.MT_axis_ave(10:end));
end
% cmap(i*floor(64/ntrials),:)
hold on;
for i = 1:ntrials
     h(2,i)=plot(linspace(0,100,length(s(i).M_MT_conv_pt_calc(1,:))), s(i).P.MT_axis_conv','color',cmap_spec(i,:),'linestyle','-.');
     
    pow_min(2,i) = min( s(i).P.MT_axis_conv(10:end));
    pow_max(2,i) = max( s(i).P.MT_axis_conv(10:end));
end

for i = 1:ntrials
     h(3,i)=plot(linspace(0,100,length(s(i).M_MT_conv_pt_calc(1,:))), s(i).P.MT_axis_inst','color',cmap_spec(i,:),'linestyle','--');
     
    pow_min(3,i) = min( s(i).P.MT_axis_inst(10:end));
    pow_max(3,i) = max( s(i).P.MT_axis_inst(10:end));
end

legend(h(1:3,1),{'Average axis','conventional joint centre','instantaneous axis'})
xlabel('% stance')
ylabel('Power (W)')
grid on
% legend('Average axis nrun','Average axis srun','Average axis frun','Conventional axis nrun','Conventional axis srun','Conventional axis frun')
makeNicePlotsFunction



contact_times =[.258 .256 .277 .287 .241];% [-13.77 -10.7 -25.92 13.828 -16.49];
[~,i] = sort(contact_times);
in = repmat(i,3,1);
makeStatsDotPlot(repmat([2,1,3],1,5),pow_min(:),{'conventional','Average axis','instantaneous'},cmap([1,10,35],:),'o'); ylabel('minimum power')
makeStatsDotPlot(repmat([2,1,3],1,5),pow_max(:),{'conventional','Average axis','instantaneous'},cmap([1,10,35],:),'o');ylabel('maximum power')
% makeStatsDotPlot(in(:)',pow_max(:),fliplr({'fastest','fast','med','slow','slowest'}),cmap([1,10,20,35,40],:),'o')
% makeStatsDotPlot(in(:)',pow_min(:),fliplr({'fastest','fast','med','slow','slowest'}),cmap([1,10,20,35,40],:),'o')

figure;
hold on;
for i = 1:ntrials
    npts = length(s(i).M_MT_conv_pt_calc(1,:));
    work.ave{1,i} = nancumtrapz((1:npts)/250, s(i).P.MT_axis_ave)';
    first_pk = max(work.ave{1,i}(10:round(npts/2)));
    sec_min = min(work.ave{1,i});
    final_max = max(work.ave{1,i});
    neg_work(1,i) = first_pk - sec_min;
    pos_work(1,i) = final_max - sec_min;
    h(1,i)=plot(linspace(0,100,npts), work.ave{1,i},'color',cmap_spec(i,:),'linestyle','-');
end

hold on;
for i = 1:ntrials
    npts = length(s(i).M_MT_conv_pt_calc(1,:));
     work.ave{2,i} = nancumtrapz((1:npts)/250,s(i).P.MT_axis_conv)';
    first_pk = max(work.ave{2,i}(10:round(npts/2)));
    sec_min = min(work.ave{2,i});
    final_max = max(work.ave{2,i});
    neg_work(2,i) = first_pk - sec_min;
    pos_work(2,i) = final_max - sec_min;
     h(2,i) = plot(linspace(0,100,npts),  work.ave{2,i} ,'color',cmap_spec(i,:),'linestyle','-.');
end

for i = 1:ntrials
    npts = length(s(i).M_MT_conv_pt_calc(1,:)); 
    work.ave{3,i} = nancumtrapz((1:npts)/250,s(i).P.MT_axis_inst)';
    first_pk = max(work.ave{3,i}(10:round(npts/2)));
    sec_min = min(work.ave{3,i});
    final_max = max(work.ave{3,i});
    neg_work(3,i) = first_pk - sec_min;
    pos_work(3,i) = final_max - sec_min;
     h(3,i) = plot(linspace(0,100,npts),  work.ave{3,i} ,'color',cmap_spec(i,:),'linestyle','--');
end
xlabel('% stance')
ylabel('Work [J]')
grid on
legend(h(1:3,1),{'Average axis','conventional joint centre','instantaneous axis'})
makeNicePlotsFunction

% figure;plot(contact_times,pos_work','o');legend('Average axis','conventional','instantaneous');ylabel('positive work')
figure;hold on;
for i = 1:5
plot(1:3,pos_work(:,i),'-o');legend('Average axis','conventional','instantaneous');ylabel('work')
plot(1:3,-neg_work(:,i),'-o');
end
% figure;plot(contact_times,neg_work','o');legend('Average axis','conventional','instantaneous'); ylabel('negative work')

% figure;plot(contact_times,pos_work'+neg_work','o');legend('Average axis','conventional','instantaneous');ylabel('net work')
% 
makeStatsDotPlot(repmat([2,1,3],1,5),pos_work(:),{'conventional','Average axis','instantaneous'},cmap([1,10,35],:),'o')
ylabel('Positive work [J]')
% makeStatsDotPlot([5,5,5,2,2,2,1,1,1,3,3,3,4,4,4]',pos_work(:),fliplr({'fastest','fast','med','slow','slowest'}),cmap([1,10,20,35,40],:),'o')
makeStatsDotPlot(repmat([2,1,3],1,5),-neg_work(:),{'conventional','Average axis','instantaneous'},cmap([1,10,35],:),'o')
ylabel('negative work [J]')
% makeStatsDotPlot([5,5,5,2,2,2,1,1,1,3,3,3,4,4,4]',neg_work(:),fliplr({'fastest','fast','med','slow','slowest'}),cmap([1,10,20,35,40],:),'o')

% Plot the stiffness
% figure;
% hold on;
% for i = 1:ntrials
%     for j = 1:3
%         plot(s(i).theta_MT(j,:)', s(i).M_MT_ax_ave(j,:)','color',cmap(j*15,:))
%     end
% end
%
% legend('dorsiflexion','inversion','adduction')
% figure;
% hold on;
% for i = 1:ntrials
%     for j = 1:3
%         plot(s(i).theta_MT(j,:)','color',cmap(j*15,:))
%     end
% end


figure;
subplot(2,2,1)
hold on;

for i = 1:ntrials
    for j = 1:3
        npts = length(s(i).M_MT_conv_pt_calc(j,:));
        plot(linspace(0,100,npts), s(i).M_MT_conv_pt_calc(j,:)', 'color',cmap(j*15,:),'linestyle','-')
    if j ==1
        [pk_val(1,i), time_pk(1,i)] = max(s(i).M_MT_conv_pt_calc(j,:));
        time_pk(1,i) = time_pk(1,i)/npts;
    end
    end
end
grid on
ylim([-160 160])
legend('dorsiflexion','inversion','adduction')
xlabel('% stance')
ylabel('conventional joint centre')

title('external midtarsal moment in talus co-ordinate system')
subplot(3,1,1)
hold on;
for i = 1:ntrials
    for j = 1:3
        npts = length(s(i).M_MT_conv_axis_calc{2}(j,:));
        plot(linspace(0,100,npts), s(i).M_MT_conv_axis_calc{2}(j,:)', 'color',cmap(j*15,:),'linestyle','-')
    if j ==1
        [pk_val(2,i), time_pk(2,i)] = max(s(i).M_MT_conv_axis_calc{2}(j,:));
        time_pk(2,i) = time_pk(2,i)/npts;
    end
    end
end

grid on
ylim([-160 160])
legend('dorsiflexion','inversion','adduction')
xlabel('% stance')
ylabel(' conventional axis')

subplot(3,1,2)
hold on;

for i = 1:ntrials
    for j = 1:3
        npts = length(s(i).M_MT_ax_ave_calc(j,:));
        plot(linspace(0,100,npts),s(i).M_MT_ax_ave_calc(j,:)', 'color',cmap(j*15,:),'linestyle','-')
    if j ==1
        [pk_val(3,i), time_pk(3,i)] = max(s(i).M_MT_ax_ave_calc(j,:));
        time_pk(3,i) = time_pk(3,i)/npts;
    end
    end
end

grid on
% ylim([-160 160])
legend('dorsiflexion','inversion','adduction')
xlabel('% stance')
ylabel(' average instantaneous axis ')


subplot(3,1,3)
hold on;


for i = 1:ntrials
    for j = 1:3
        npts = length(s(i).M_MT_ax_inst_calc(j,:));
        plot(linspace(0,100,npts),s(i).M_MT_ax_inst_calc(j,:)', 'color',cmap(j*15,:),'linestyle','-')
    if j ==1
        [pk_val(4,i), time_pk(4,i)] = max(s(i).M_MT_ax_inst_calc(j,:));
        time_pk(4,i) = time_pk(4,i)/npts;
    end
    end
end

grid on
% ylim([-160 160])
legend('dorsiflexion','inversion','adduction')
xlabel('% stance')
ylabel('Instantaneous axis')
makeNicePlotsFunction

for i = 1:5
j = 2;
inst_diff{i} = (s(i).M_MT_ax_inst_calc(j,:)'-s(i).M_MT_conv_axis_calc{2}(j,:)')./s(i).M_MT_conv_axis_calc{2}(j,:)';

ave_diff{i} = (s(i).M_MT_ax_ave_calc(j,:)'-s(i).M_MT_conv_axis_calc{2}(j,:)')./s(i).M_MT_conv_axis_calc{2}(j,:)';
end
figure; hold on;
for i = 1:5
    nf = length(inst_diff{i});
    plot(inst_diff{i}(round(0.5*nf):end)*100)
end
figure; hold on;
for i = 1:5
    
    nf = length(inst_diff{i});
    plot(ave_diff{i}(round(0.5*nf):end)*100)
end
[m,st] = makeStatsDotPlot(repmat([2,3,4],1,5),pk_val(:),{'conventional','Average axis','instantaneous'},cmap([1,10,35],:),'o')
ylabel('peak dorsiflexion moment [Nm]')
makeStatsDotPlot(repmat([2,3,4],1,5),time_pk(:),{'conventional','Average axis','instantaneous'},cmap([1,10,35],:),'o')
ylabel('time of peak dorsiflexion moment')
%
for i = 1:3
    subplot(3,1,i)
    xlim([20 90])
    ylim([-150 200])
end

% for i = 1:3
%     for j = 1:3
%     plot(linspace(0,100,length(s(i).M_MT_conv_pt_calc(j,:))), s(i).M_MT_conv_pt_calc(j,:)'- s(i).M_MT_ax_ave_calc(j,:)', 'color',cmap(j*15,:))
%     end
% end
%% animate the force vector +moment arms etc in Wrist visualizer

c = [195,48,48;
    0, 153, 153;
    51,51,255;
    153,255,154]/256;
iha_col = c(1,:);
aha_col = c(2,:);
ca_col = c(3,:);

% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);
for i = 1:ntrials
    anim_dir = fullfile(subj_dir,trialname_list{i},'POS','MomentAnimation',filesep);
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
        ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
        
        fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
    nFrs =length(struct_data(i).cropFrsForce(1):struct_data(i).cropFrsForce(2) );
    for bn = 1:nBones
        %         for fr = first_fr:end_frM
        %         Tw.(bone_list{bn})(:,:,fr) = struct_data(i).Tm.(bone_list{bn})(:,:,fr);
        %         end
        if ~isfield(s(i).Tm,bone_list{bn})
            Tw = nan(4,4,nFrs);
        else
            Tw = s(i).Tm.(bone_list{bn});
            Tw(1:3,4,:) = Tw(1:3,4,:)*1000;
        end
        Tw(isnan(Tw)) = 1;
        write_RTp(bone_list{bn} , Tw, anim_dir)
    end
    
    
    %     frX = 1;
    for fr = 1:nFrs
        
        grf_filename = sprintf(grf_style,fr);
        %         frX = frX+1;
        ivstring = createInventorHeader();
        
        
        %         ************** HELICAL AXES **************************
        for bp = 1%:nBP
            if any(isnan(s(i).heli(bp).s_m(1,fr)),'all')
                continue
            end
            if any(isnan(s(i).heli(bp).s_avg_m(1,fr)),'all')
                continue
            end
            ivstring = [ivstring createInventorArrow( s(i).heli(bp).s_m(:,fr)*1000,...
                -s(i).heli(bp).n_m(:,fr),150, 2,iha_col,0)];
%                 s(i).heli(bp).phi(fr)*250,...
               
            
            ivstring = [ivstring createInventorArrow( s(i).heli(bp).s_avg_m(:,fr)*1000,...
                s(i).heli(bp).n_avg_m(:,fr),150, 2,aha_col,0)];
%                 s(i).heli(bp).phi(fr)*250,...
               
%             
% %             ivstring = [ivstring createInventorArrow( s(i).arch_pt_ave(:,fr)*1000,...
% %                 s(i).arch_ax_ave(:,fr),...
% %                 s(i).heli(bp).phi(fr)*250,...
% %                 2,[0,1, 0.2],0)];
%             
%             ivstring = [ivstring createInventorArrow( s(i).arch_pt_inst(:,fr)*1000,...
%                 s(i).arch_ax_inst(:,fr),...
%                 s(i).heli(bp).phi(fr)*250,...
%                 2,[0,0.2, 0.8],0)];
        end
        
%         LEGEND
%             createInventorText('Prescribed/Interpolated average helical axis',12,[30 15 0],[0.7 0.7 0.7],0),...
%             createInventorArrow([0 30 0],[1 0 0],30,2,[0 0.2 0.8],0),...
text_col = [0 0 0];
        ivstring = [ivstring ,...
            createInventorArrow([0 0 0],[1 0 0],30,2,ca_col,0),...
            createInventorText('conventional axis',12,[30 -5 0],text_col,0),...
            createInventorArrow([0 10 0],[1 0 0],30,2,aha_col,0),...
            createInventorText('average helical axis',12,[30 5 0],text_col,0),...
            createInventorArrow([0 20 0],[1 0 0],30,2,iha_col,0),...
            createInventorText('instantaneous helical axis',12,[30 15 0],text_col,0)];
        %*********FORCE VECTORS****************
        for fp = 1:2 % each force plate
            COP = s(i).COP{fp}(:,fr) *1000;% transformPoints(coreg,struct_data(i).force_data(fp).globCOP(:,fr)*1000,-1);
            force = s(i).force{fp}(:,fr);%transformVectors(coreg,struct_data(i).force_data(fp).globForce(:,fr),-1);
            
            %         ivstring = [ivstring createInventorSphere(COP,4,[0.2 0.2 0.8],0)];
            ivstring = [ivstring createInventorArrow(COP,force,norm(force)/5,1.5,[1 1 1],0)];
        end
        
        
        %         *************** MOMENT ARMS **********************
       
        for j = 2%1:8
            
            if any(isnan(s(i).axis_vec_test{j}(1,fr)),'all')
                continue
                
            end
%             % conventional axis
            ivstring = [ivstring createInventorArrow( s(i).JC_test{j}(:,fr)*1000,...
                s(i).axis_vec_test{j}(:,fr),...
                100,...
                2,ca_col,0)];
            
            ivstring = [ivstring createInventorArrow( s(i).JC_test{j}(:,fr)*1000,...
                -s(i).axis_vec_test{j}(:,fr),...
                100,...
                2,ca_col,0)];
%             % centre
%             ivstring = [ivstring createInventorSphere( s(i).JC_test{j}(:,fr)*1000,...
%                 5,[0.5,0, 0.5],0)];
%             
%             
%             if any(isnan([s(i).s_g_MT1_conv_axis{j}(:,fr); s(i).r_g_MT1_conv_axis{j}(:,fr)]),'all') || (s(i).force{1}(3,fr) < 10)
%             else
%                 ivstring = [ivstring createInventorArrow( s(i).s_g_MT1_conv_axis{j}(:,fr)*1000,...
%                     s(i).r_g_MT1_conv_axis{j}(:,fr)*1000,...
%                     norm(s(i).r_g_MT1_conv_axis{j}(:,fr)*1000),...
%                     1,[0.5,0, 0.5],0)];
%             end
%             
%             if any(isnan([s(i).s_g_MT2_conv_axis{j}(:,fr); s(i).r_g_MT2_conv_axis{j}(:,fr)]),'all') || (s(i).force{2}(3,fr) < 10)
%             else
%                 ivstring = [ivstring createInventorArrow( s(i).s_g_MT2_conv_axis{j}(:,fr)*1000,...
%                     s(i).r_g_MT2_conv_axis{j}(:,fr)*1000,...
%                     norm(s(i).r_g_MT2_conv_axis{j}(:,fr)*1000),...
%                     1,[0.5,0, 0.5],0)];
%             end
%         end
%         
%         % and the instantaneous moment arms
%         %average
%         if any(isnan([s(i).s_g_MT1ax_ave(:,fr); s(i).r_g_MT1ax_ave(:,fr)]),'all') || (s(i).force{1}(3,fr) < 10)
%         else
%             ivstring = [ivstring createInventorArrow(s(i).s_g_MT1ax_ave(:,fr)*1000,...
%                 s(i).r_g_MT1ax_ave(:,fr) *1000,...
%                 norm(s(i).r_g_MT1ax_ave(:,fr) *1000),...
%                 1,[0,1, 0.2],0)];
%             
%         end
%         if any(isnan([s(i).s_g_MT2ax_ave(:,fr); s(i).r_g_MT2ax_ave(:,fr)]),'all') || (s(i).force{2}(3,fr) < 10)
%         else
%             ivstring = [ivstring createInventorArrow(s(i).s_g_MT2ax_ave(:,fr)*1000,...
%                 s(i).r_g_MT2ax_ave(:,fr) *1000,...
%                 norm(s(i).r_g_MT2ax_ave(:,fr) *1000),...
%                 1,[0,1, 0.2],0)];
%             
%         end
%         %instant
%         if any(isnan([s(i).s_g_MT1ax_inst(:,fr); s(i).r_g_MT1ax_inst(:,fr)]),'all') || (s(i).force{1}(3,fr) < 10)
%         else
%             ivstring = [ivstring createInventorArrow(s(i).s_g_MT1ax_inst(:,fr)*1000,...
%                 s(i).r_g_MT1ax_inst(:,fr) *1000,...
%                 norm(s(i).r_g_MT1ax_inst(:,fr) *1000),...
%                 1,[0 0.1 0.8],0)];
%             
%         end
%         if any(isnan([s(i).s_g_MT2ax_inst(:,fr); s(i).r_g_MT2ax_inst(:,fr)]),'all') || (s(i).force{2}(3,fr) < 10)
%         else
%             ivstring = [ivstring createInventorArrow(s(i).s_g_MT2ax_inst(:,fr)*1000,...
%                 s(i).r_g_MT2ax_inst(:,fr) *1000,...
%                 norm(s(i).r_g_MT2ax_inst(:,fr) *1000),...
%                 1,[0 0.1 0.8],0)];
%             
        end
        
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
































