% code to calculate the global axes of the arch 

close all
clear
clc

subj_name = 'SOL001B';
subj_dir = 'E:\SOL001_VISIT2\';
% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);

load([subj_dir 'struct_bone.mat']);
load([subj_dir 'struct_data.mat']);


% get the trials that we want
% 'T0019_SOL001_nrun_pref_barefoot';...
   trialname_list = {                  'T0024_SOL001_nrun_rfs_barefoot';%...
%                  'T0033_SOL001_srun_pref_barefoot';...
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



bone_list = {struct_bone(:).bone_name}';
nBones = length(bone_list);

cmap = colormap('parula');
% make a cropped transform structure & force structure


for i = 1:ntrials
    
    start_frX = struct_data(i).cropFrsXray(1);
    end_frX = struct_data(i).cropFrsXray(2);
    
    for bn = 1:nBones
        if ~isfield(struct_data(i).Tm,bone_list{bn})
             struct_Tm(i).(bone_list{bn})= nan(4,4,end_frX - start_frX + 1);
        else
            
        struct_Tm(i).(bone_list{bn}) = struct_data(i).Tm.(bone_list{bn})(:,:,start_frX:end_frX);
        end        
        
    end
end


% make the force structure cropped

for i = 1:ntrials
     
    start_frF = struct_data(i).cropFrsForce(1);
    end_frF = struct_data(i).cropFrsForce(2);
    
    % crop all the force parameters to sync it up with the helical axes
    for fp = 1:2 % force plates
    struct_force(i).force{fp} = struct_data(i).force_data(fp).globForce(:,start_frF:end_frF);
    struct_force(i).COP{fp} = struct_data(i).force_data(fp).globCOP(:,start_frF:end_frF);
    struct_force(i).freemom{fp} = struct_data(i).force_data(fp).globFreeMoment(:,start_frF:end_frF);
    end
    
    
end
figure; hold on;
plot(struct_force(i).force{1}(3,:));plot(struct_force(i).force{2}(3,:))
%% helical axes of each of the bones (global)
close all
% bone_pairs = {'glob','cal';'glob','tal';'glob','mt1';'glob','cmm';'glob','nav';'glob','cub';'glob','tib'}%;'tal','tib'}%;'mt1','ph1';'mt1','cmm';'mt1','nav';'mt1','tal'};
bone_pairs = {'mt1','cmm';'mt1','nav';'mt1','tal';'mt1','cal'};%{'cal','tal';'tal','nav';'nav','cmm';'cmm','mt1'};%
bone_pairs = {'cal','tal';'tal','nav';'nav','cmm';'cmm','mt1';'tib','tal'};%;'cal','cub'};%;'glob','cub'};
nBP = size(bone_pairs,1);
% [58 65 56] % max mtp frame
%     
% [ 34 36  31] % approx start of propulsion
for bn = 1:nBones
    
    T_ACS.(struct_bone(bn).bone_name) = struct_bone(bn).T_ACS;
end


for i = 1:ntrials
%     figure(i);hold on;
%     figure(10+i);hold on;
                s_ref_moc = [];
                
                s_ct = [];
                
                n_m_mt1 = [];
                phi_n_m_mt1 = [];
                
    start_frX = struct_data(i).cropFrsXray(1);
    end_frX = struct_data(i).cropFrsXray(2);
    
    % project all the angles back into the metatarsal co-ordinate
    % system
    ifr = 0;
    for fr = start_frX:end_frX
        ifr = ifr+1;
        T_ACS_m.mt1(:,:,ifr) = struct_data(i).Tm.mt1(:,:,fr) * T_ACS.mt1;
    end
    
    leg_str = {};
    for bp = 1:nBP
        
%     figure(bp);hold on;
%     figure(10+bp);hold on;

        
%         bone that is being compared
        T_mA_ct = struct_data(i).Tm.(bone_pairs{bp,2})(:,:,start_frX:end_frX);
        
        
       
        if strcmp(bone_pairs{bp,1},'glob')
            T_mR_ct = repmat(eye(4),1,1,end_frX-start_frX+1);
        else    
            T_mR_ct = struct_data(i).Tm.(bone_pairs{bp,1})(:,:,start_frX:end_frX);
        end
        
        [phi_raw,n_raw,L_raw,s_raw] = helicalInstantaneous(T_mR_ct,T_mA_ct);
%         [mean_phi,mean_n,mean_L,mean_s]  = stabilizeHelicalAxis(T_mR_ct,T_mA_ct,5);
        
        
        % find the average helical axis by dividing into groups by the
        % magnitude of phi
        ind_big_phi = phi_raw>0.5;
        ind_small_phi = ~ind_big_phi;
        avg_ind = findGroupsOfLogicals(ind_big_phi);% find the values to average over
        int_ind = findGroupsOfLogicals(ind_big_phi);% find the indices that will be interpolated over
        
        nfrs = length(phi_raw);
        n_avg = nan(3,nfrs);
        s_avg = nan(3,nfrs);
        for gi = 1:size(avg_ind)
            ii = false(nfrs,1);
            ii(avg_ind(gi,1):avg_ind(gi,2),1) = 1;
            [avgHAM(gi,:),num,variat] = calcAvgMotionLW( [phi_raw(ii)',n_raw(:,ii)',L_raw(ii)',s_raw(:,ii)']);
            for dd = 1:3
            n_avg(dd,ii) = avgHAM(gi,dd+1)';
            s_avg(dd,ii) = avgHAM(gi,dd+5)';
            end
       
            if mean(n_avg(1,ii)) < 0
                n_avg(:,ii) = -n_avg(:,ii);
            end
        end
%         figure; plot(n_raw'); hold on; plot(n_avg');
        for gi = 1:size(int_ind,1)
            p1 = int_ind(gi,1);
            p2 = int_ind(gi,2);
            if p2-p1 < 3 && p1>1 && p2 < nfrs-2
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
%         plot(n_avg')
       
        
        % space for further processing
        phi = phi_raw;
        n_m = transformVectors(T_mR_ct,n_raw); % drive with the bone kinematics of bone
        L = L_raw;
        % find the closest point to the centroid
        
        s_m = transformPoints(T_mR_ct,s_raw);
             bn_ind = findInStruct(struct_bone,'bone_name',bone_pairs{bp,2});
            for fr = 1:end_frX-start_frX
                centR_moc = transformPoints(T_mA_ct(:,:,fr),struct_bone(bn_ind).cent/1000);
                s_ref_moc(:,fr) = closestPointonVector(centR_moc,s_m(:,fr),n_m(:,fr));
                
                s_ct(:,fr) = transformPoints(T_mR_ct(:,:,fr),s_ref_moc(:,fr),-1);
                
                n_m_mt1(:,fr) = transformVectors(T_ACS_m.mt1(:,:,fr),n_m(:,fr),-1);
                phi_n_m_mt1(:,fr) = phi(fr) * n_m_mt1(:,fr) ;
                
            end
        
        struct_heli(i).trials(bp).refBone = bone_pairs{bp,1};
        struct_heli(i).trials(bp).bone = bone_pairs{bp,2};
        % all in MOCAP space
        struct_heli(i).trials(bp).phi = phi;
        struct_heli(i).trials(bp).n_m = n_m;
        struct_heli(i).trials(bp).L = L;
        struct_heli(i).trials(bp).s_m = s_ref_moc;
        
        
        struct_heli(i).trials(bp).n_ct = n_raw;
        struct_heli(i).trials(bp).s_ct = s_ct;
        
       
        phi_save{bp}(i,:) = normaliseNaN(phi,2,100);
%         n_m_save{i,bp}(1:3,:) = normaliseNaN(n_m,2,100);
        for dd = 1:3
        phi_n_m_save{bp,dd}(i,:) = normaliseNaN(phi_n_m_mt1(dd,:),2,100);
        end
        L_save{bp}(i,:) = normaliseNaN(L,2,100);
        
%         
%         figure(1);
%         hold on;
%         plot(linspace(0,100,length(phi)),phi.*n_m(1,:),'color',cmap(bp*floor(64/nBP),:))
        
%         title(sprintf('%s relative to %s',bone_pairs{bp,2},bone_pairs{bp,1}))
%         ylim([0 2.5])
%         figure(2);
%         hold on;
%         plot(linspace(0,100,length(phi)),L,'color',cmap(bp*floor(64/nBP),:))
%         ylim([-1 1])
%         
        leg_str = [leg_str sprintf('%s relative to %s',bone_pairs{bp,2},bone_pairs{bp,1})];
        
        
        if ~strcmp(bone_pairs{bp,1},'glob')
            T_ACSA = T_ACS.(bone_pairs{bp,1});
            T_ACSR = T_ACS.(bone_pairs{bp,2});
            [x_ang,y_ang,z_ang] = eulerYZX(T_mR_ct,T_mA_ct, T_ACSR,T_ACSA);
            figure(100);hold on; plot(x_ang,'color',cmap(bp*floor(64/nBP),:))
%             figure(101); hold on; plot(y_ang,'color',cmap(bp*floor(64/nBP),:))
%             figure(102); hold on; plot(z_ang,'color',cmap(bp*floor(64/nBP),:)) 
            title(sprintf('%s relative to % s',bone_pairs{bp,2},bone_pairs{bp,1}));
%             [x_ang,y_ang,z_ang] = eulerYZX(repmat(eye(4,4),1,1,length(T_mA_ct)),T_mA_ct,eye(4,4), TACS.tib);
%             figure(106);hold on; plot(x_ang)
%             figure(107); hold on; plot(y_ang)
%             figure(108); hold on; plot(z_ang)
%             force=norm3d(struct_force(i).force{1}+struct_force(i).force{2});
%             figure(103); hold on; plot( phi,force(1:end-1))
%         elseif bp == 9
%               TACS.mt1 = struct_bone(5).T_ACS;
%             TACS.ph1= struct_bone(8).T_ACS;
%             [x_ang,y_ang,z_ang] = eulerYZX(T_mR_ct,T_mA_ct, TACS.mt1, TACS.ph1);
%             figure(103);hold on; plot(x_ang)
%             figure(104); hold on; plot(y_ang)
%             figure(105); hold on; plot(z_ang)
            
        end
% %         title(sprintf('%s relative to %s',bone_pairs{bp,2},bone_pairs{bp,1}))
    clearvars('s_ref_moc','s_ct')
    end
%     legend(bone_pairs(:,2))
    figure(1);
    legend(leg_str)
%     legend(bone_pairs(:,2))
    

%  
%     for fr = 1:end_frX-start_frX
%         
%         transformVectors(T_ACS_m.mt1(:,:,fr),n_m,-1)
%     end
%     

end


leg_dirs = {'dorsiflexion','inversion','adduction'};
close all
for dd = 1:3
    for bp = 1:nBP
        
        figure(dd)
        hold on;
        PrettyStdDevGraphs(1:100,nanmean(  phi_n_m_save{bp,dd},1),nanstd( phi_n_m_save{bp,dd},1),cmap(bp*floor(64/nBP),:) ,1)
        
        ylabel([leg_dirs{dd} ' frame by frame Rotation between bones'])
        % ylim([-2.5 2])
        legend(leg_str)
        
        
        bp_bar(bp,:) = nanmean(   phi_n_m_save{bp,dd}(:,50:end) ,1);
        
       
    end
        figure(dd+15)
        hold on;
        ind_pos = bp_bar>0;
        ind_neg = bp_bar<0;
        bp_pos = bp_bar; bp_neg = bp_bar;
        bp_pos(ind_neg) = 0;
        bp_neg(ind_pos) = 0;
        hp = bar(50:100,(bp_pos)','stacked');  
        hn = bar(50:100,(bp_neg)','stacked');
        
        for bp = 1:nBP
            hp(bp).FaceColor = cmap(bp*floor(64/nBP),:);
            hn(bp).FaceColor = cmap(bp*floor(64/nBP),:);
        end
        
        ylabel([leg_dirs{dd} ' frame by frame Rotation between bones'])
        hl = legend('subtalar','talonavicular','cuneonavicular','cuneo-metatarsal','talocrural');
        hl.Location = 'northwest';
        
end

for bp = 1:nBP
    
    for dd = 1:3
        figure(dd+4)
        hold on;
        phi_n_sum_m_save{bp,dd}(1:ntrials,1:100) =  zeros(ntrials,100);
        phi_n_sum_m_save{bp,dd}(1:ntrials,50:end) = nancumtrapz( phi_n_m_save{bp,dd}(:,50:end));
        PrettyStdDevGraphs(50:100,nanmean( phi_n_sum_m_save{bp,dd}(:,50:end) ,1),nanstd(  phi_n_sum_m_save{bp,dd}(:,50:end),1),cmap(bp*floor(64/nBP),:) ,1)
%         ylim([-25 25])
        ylabel([leg_dirs{dd} ' summed rotation between bones'])
        legend(leg_str)
        
    end
    
end






bone_pairs{ind_sort,2};
%% Make the animation of the global arch axes in wrist viz


% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);
for i = 1:ntrials
    
    
    anim_dir = fullfile(subj_dir,trialname_list{i},'POS','GlobalArchAxesAnimation_FixedTalus',filesep);
    rigidiv_folder = fullfile(anim_dir,'rigidiv',filesep);
    
    % set up the GRF animation folder
    arch_style = 'arch_P%i_F1.iv';
    arch_folder = fullfile(anim_dir,'Archdir',filesep);
    
    if exist(anim_dir,'dir')==0;     mkdir(anim_dir);     end
    if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end
    if exist(arch_folder,'dir')==0;  mkdir(arch_folder);  end
    
        first_frX = struct_data(i).cropFrsXray(1);
        end_frX   = struct_data(i).cropFrsXray(2);
    
        
        % fix the bones relative to a bone; replace any nan frames with the
        % first or last tracked frames
        T_fix =struct_data(i).Tm.tal(:,:,first_frX:end_frX);%repmat(eye(4,4),1,1,end_frX-first_frX+1);% 
        ind_nan = find(isnan(T_fix(1,1,:)))';
        ind_nonan = find(~isnan(T_fix(1,1,:)))';
        for in = ind_nan
            if in < (end_frX-first_frX+1)/2
                T_fix(:,:,in) = T_fix(:,:,ind_nonan(1));
            else
                T_fix(:,:,in) = T_fix(:,:,ind_nonan(end));
            end
        end
        
        
    for bn = 1:nBones
        ivstring = createInventorHeader();
        % make the linked iv file
        jj = find(strcmp(bone_pairs(:,2),bone_list{bn})); % find the index of the bone pair to color the bone the same as the arrow
        if isempty(jj)
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7] ,0.5)];
        else
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),cmap(jj*floor(64/nBP),:) ,0.5)];
        end
        fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
    for bn = 1:nBones
        if ~isfield(struct_data(i).Tm,bone_list{bn})
            continue
        end
        frWV = 1;
        for fr = first_frX:end_frX
            Ttemp = invTranspose(T_fix(:,:,frWV)) * struct_data(i).Tm.(bone_list{bn})(:,:,fr);
            if any(isnan(Ttemp)) % wrist viz likes 1's where nans are
                Tw.(bone_list{bn})(:,:,frWV) = ones(4,4);
            else
                Ttemp(1:3,4) = Ttemp(1:3,4)*1000;
                Tw.(bone_list{bn})(:,:,frWV) = Ttemp;
            end
            frWV = frWV+1;
        end
        
        write_RTp(bone_list{bn} , Tw.(bone_list{bn}), anim_dir)
    end
    
    
    %     frX = 1;
    for fr = 1:(end_frX-first_frX )
        
        arch_filename = sprintf(arch_style,fr);
        %         frX = frX+1;
        ivstring = createInventorHeader();
        
        
        %         ************** HELICAL AXES **************************
        for bp = 1:nBP
            if any(isnan(struct_heli(i).trials(bp).s_m(:,fr)),'all')
                continue
            end
            
            s_m = transformPoints(T_fix(:,:,fr),struct_heli(i).trials(bp).s_m(:,fr),-1)*1000;
            n_m = transformVectors(T_fix(:,:,fr),struct_heli(i).trials(bp).n_m(:,fr),-1);
            
            ivstring = [ivstring createInventorArrow(s_m,...
                n_m,...
               struct_heli(i).trials(bp).phi(fr)*50,...
                2,cmap(bp*floor(64/nBP),:),0)];
            
        end
        
        %*********FORCE VECTORS****************
        for fp = 1:2 % each force plate
            COP = transformPoints(T_fix(:,:,fr),struct_force(i).COP{fp}(:,fr),-1) *1000;% transformPoints(coreg,struct_data(i).force_data(fp).globCOP(:,fr)*1000,-1);
            force = transformVectors(T_fix(:,:,fr),struct_force(i).force{fp}(:,fr),-1);%transformVectors(coreg,struct_data(i).force_data(fp).globForce(:,fr),-1);
            
            %         ivstring = [ivstring createInventorSphere(COP,4,[0.2 0.2 0.8],0)];
            ivstring = [ivstring createInventorArrow(COP,force,norm(force)/5,1.5,[1 1 1],0)];
        end
        

        %************* WRITE FILE WITH ARROWS *****************
        fid = fopen(fullfile(arch_folder,arch_filename),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        
    end
    create_ini(0,0,1,1,arch_folder(1:end-1),'arch_P%d_F%d.iv',fullfile(anim_dir,[trialname_list{i} '.ini']))
    
    pos_text = write_pos(bone_list,anim_dir,trialname_list{i});
    
    filename = fullfile(anim_dir, [trialname_list{i} '.pos']);
    
    fid = fopen(filename,'w'); % open the file to write
    fprintf(fid,pos_text);
    fclose(fid);
    
    fprintf('Animations are saved in %s.\n',anim_dir)
    clearvars('Tw')
end


%% Fix the arch with the metatarsal and compare ankle moment arms
close all


max_mtp = [58 65 56];
match_arch = [16 18 15];%[20,23,20];%
tal_cent = struct_bone(9).cent/1000;

for i = 1:ntrials
    
    nFrs = length(struct_force(i).force{2}(3,:));
    nFrs_save(i) = nFrs;
    %     find the instant of heel lift
    ind_off = find(struct_force(i).force{2}(3,:) < 10);
    ind_secondhalf = ind_off > 0.30 * nFrs/2; % to avoid the beginning frames
    ind_off = ind_off(ind_secondhalf);
    
    ind_prop(i,:) = [ind_off(1),ind_off(end)]; % propulsion frames
    first_frX = ind_prop(i,1);
    end_frX   = ind_prop(i,2);
    mMTP = max_mtp(i);
    
    
    % put the talocrural axis in the talar co-sys
    ind_ref = findInStruct(struct_heli(i).trials,'refBone','tal');
    ind_comp = findInStruct(struct_heli(i).trials,'bone','tib');
    ind_bp = intersect(ind_ref,ind_comp);
    
    axis_TC_ct = struct_heli(i).trials(ind_bp).n_ct; % axis in global (mocap) space
    pos_TC_ct = struct_heli(i).trials(ind_bp).s_ct;
    
    
    % get the metatarsal axis:
    
    ind_ref = findInStruct(struct_heli(i).trials,'refBone','glob');
    ind_comp = findInStruct(struct_heli(i).trials,'bone','mt1');
    ind_bp = intersect(ind_ref,ind_comp);
    
    axis_mt1_m = struct_heli(i).trials(ind_bp).n_m; % axis in global (mocap) space
    pos_mt1_m = struct_heli(i).trials(ind_bp).s_m;
    
    tal_ind = find(strcmp(bone_list,'tal'));
    
    %initialize all the axes and origins for the talocrural axes
    axis_TC_m{i}= nan(3,ind_prop(i,2));
    pos_TC_m{i}= nan(3,ind_prop(i,2));
    axis_TC_lever{i} = nan(3,ind_prop(i,2));
    pos_TC_lever{i} = nan(3,ind_prop(i,2));
    
    for fr = ind_prop(i,1):ind_prop(i,2)-1
        
        % fixed with the metatarsal, but
        T_fix = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) ) * struct_Tm(i).(bone_list{tal_ind})(:,:,first_frX);
        
        
        % fixed to plantarflexed position @max mtp 
        T_plant = struct_Tm(i).mt1(:,:,fr)* invTranspose( struct_Tm(i).mt1(:,:,first_frX) )  * struct_Tm(i).(bone_list{tal_ind})(:,:,first_frX) ;

        if fr == first_frX
            tal_start(:,i) = transformPoints(T_fix,tal_cent);
            tal_start_act(:,i) =  transformPoints(struct_Tm(i).(bone_list{tal_ind})(:,:,fr),tal_cent);
%             tal_start_act(:,i) =  transformPoints(T_plant,tal_cent);
            
        elseif fr == match_arch(i) + first_frX
%             disp('match arch')
%             disp(fr)
            tal_plantmatch(:,i) = transformPoints(T_plant,tal_cent);
%             disp()
        elseif fr == max_mtp(i)
%             disp('max mtp')
%             disp(fr)
            tal_actual(:,i) = transformPoints(struct_Tm(i).(bone_list{tal_ind})(:,:,fr),tal_cent);
%             disp(transformPoints(struct_Tm(i).(bone_list{tal_ind})(:,:,fr),tal_cent))
            tal_rigid(:,i) = transformPoints(T_fix,tal_cent);
%             disp(transformPoints(T_fix,tal_cent))
             
            
        end
        
        %         axis_TC_tal(:,fr) = transformVectors(struct_bone(tal_ind).T_ACS,axis_TC_ct(:,fr),-1); % put n in talus co-sys for every frame
        %
        %         T_drive = invTranspose(struct_Tm.mt1(:,:,fr)) * struct_Tm.tal(:,:,fr) * struct_bone(tal_ind).T_ACS;
        axis_TC_lever{i}(:,fr) = transformVectors(T_fix,axis_TC_ct(:,fr));
        pos_TC_lever{i}(:,fr) = transformPoints(T_fix,pos_TC_ct(:,fr));
        
        for bp = 1:nBP
            bi = strmatch(bone_pairs{bp,2},bone_list); % index of bone in list
            T.(bone_list{bi}).lev(:,:,fr) = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) ) * struct_Tm(i).(bone_list{bi})(:,:,first_frX);
            T.(bone_list{bi}).act(:,:,fr) = struct_Tm(i).(bone_list{bi})(:,:,fr);
        end
        
        % fix it to the point at max MTP
            T.tib.lev(:,:,fr) = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) ) * struct_Tm(i).(bone_list{bi})(:,:,mMTP);
       
        T_fix = struct_Tm(i).(bone_list{tal_ind})(:,:,fr);
        axis_TC_m{i}(:,fr) = transformVectors(T_fix,axis_TC_ct(:,fr));
        pos_TC_m{i}(:,fr) = transformPoints(T_fix,pos_TC_ct(:,fr));
        
%         from COP
    [s1_lev,s2_lev] = closestPointsBtw2Lines( struct_force(i).COP{1}(:,fr) , pos_TC_lever{i}(:,fr),struct_force(i).force{1}(:,fr) ,axis_TC_lever{i}(:,fr));
    [s1_act,s2_act] = closestPointsBtw2Lines( struct_force(i).COP{1}(:,fr),pos_TC_m{i}(:,fr),struct_force(i).force{1}(:,fr) ,   axis_TC_m{i}(:,fr));
% from met head
    [s1_lev,s2_lev] = closestPointsBtw2Lines( pos_mt1_m(:,fr) , pos_TC_lever{i}(:,fr),axis_mt1_m(:,fr) , axis_TC_lever{i}(:,fr));
    [s1_act,s2_act] = closestPointsBtw2Lines( pos_mt1_m(:,fr) ,pos_TC_m{i}(:,fr),axis_mt1_m(:,fr) ,   axis_TC_m{i}(:,fr));
        mom_arm_lev{i}(:,fr) = s1_lev-s2_lev;
        mom_arm_act{i}(:,fr) = s1_act-s2_act;
        
        
        ank_mom_lev{i}(:,fr) = cross( mom_arm_lev{i}(:,fr),struct_force(i).force{1}(:,fr));
        ank_mom_act{i}(:,fr) = cross( mom_arm_act{i}(:,fr),struct_force(i).force{1}(:,fr));
        
        w{i}(:,fr) =  (struct_heli(i).trials(bp).phi(:,fr) * struct_heli(i).trials(bp).n_m(:,fr)) * frRateX;
        ank_pow_lev{i}(:,fr) = dot( ank_mom_lev{i}(:,fr),w{i}(:,fr));
        ank_pow_act{i}(:,fr) = dot( ank_mom_act{i}(:,fr),w{i}(:,fr));
    end
    
    %          [x_ang,y_ang,z_ang] = eulerYZX(repmat(eye(4,4),1,1,fr),T_tib_lev,eye(4,4),eye(4,4));
    %             figure(106);hold on; plot(x_ang)
    %             figure(107); hold on; plot(y_ang)
    %             figure(108); hold on; plot(z_ang)
    %
    %             [x_ang,y_ang,z_ang] = eulerYZX(repmat(eye(4,4),1,1,fr),T_tib_act,eye(4,4),eye(4,4));
    %             figure(106);hold on; plot(x_ang,'--')
    %             figure(107); hold on; plot(y_ang,'--')
    %             figure(108); hold on; plot(z_ang,'--')
    
    
    
    %             T_fix = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) ) * struct_Tm(i).(bone_list{tal_ind})(:,:,first_frX);
    
    for bp = 1:nBP
        if strcmp(bone_pairs{bp,1},'mt1')
            
            bi = strmatch(bone_pairs{bp,2},bone_list); % index of bone in list
            T_ACSA = T_ACS.(bone_pairs{bp,2});
            T_ACSR = T_ACS.(bone_pairs{bp,1});
            t_vals = linspace(first_frX/nFrs*100,end_frX/nFrs*100,end_frX-first_frX);
            [x_ang,y_ang,z_ang] = eulerYZX( T.(bone_list{bi}).lev(:,:,first_frX:end_frX-1), T.mt1.act(:,:,first_frX:end_frX-1), T_ACSA,  T_ACSR);
            
            dors.lev{bp}(i,:) = nan(1,100);
            dors.lev{bp}(i,round(t_vals(1)):100)= normaliseNaN(x_ang,2,100-round(t_vals(1))+1);
            
            figure(100);hold on; 
            plot(t_vals,x_ang,'color',cmap(bp*floor(64/nBP),:))
            [x_ang,y_ang,z_ang] = eulerYZX(  T.(bone_list{bi}).act(:,:,first_frX:end_frX-1),T.mt1.act(:,:,first_frX:end_frX-1), T_ACSA, T_ACSR);
            plot(t_vals, x_ang,'color',cmap(bp*floor(64/nBP),:),'linestyle','--')
            title(sprintf('%s relative to % s',bone_pairs{bp,2},bone_pairs{bp,1}));
            
            
            dors.act{bp}(i,:) = nan(1,100);
            dors.act{bp}(i,round(t_vals(1)):100)= normaliseNaN(x_ang,2,100-round(t_vals(1))+1);

        end
    end
end


bone_full_name = {'medial cuneiform','navicular','talus','cal'};
 figure;
%  cmap = 
%  nBP = 5
for bp = 1:4
%   subplot(3,1,bp)
   hold on;
   
   PrettyStdDevGraphs(1:100,-1*(nanmean(dors.act{bp})-nanmean(dors.lev{bp})),nanstd(dors.act{bp}),cmap(bp*floor(64/nBP),:),1);
   hold on

    ylim([-5 25])
%     ylabel([bone_full_name{bp}])
    
    xlim([50 90])
end
  h = PrettyStdDevGraphs(1:100,-1*(nanmean(dors.lev{bp})-nanmean(dors.lev{bp})),nanstd(dors.lev{bp}),[0 0 0],1);
   h(2).LineStyle = ':';
    h = legend([bone_full_name,'rigid lever']);
    h.Location = 'northwest';
ylabel('plantarflexion (+) of first metatarsal ')
xlabel(' % stance')
makeNicePlotsFunction

v_plant_match =  tal_plantmatch- tal_start;
v_rigid_actual_mtp = tal_rigid -tal_actual;


v_rigid_mtp_start = tal_rigid -tal_start;
v_actual_mtp = tal_actual- tal_start;



mean((v_rigid_mtp_start-v_actual_mtp) ./ v_actual_mtp,2) * 100
std((v_rigid_mtp_start-v_actual_mtp) ./ v_actual_mtp,[],2)*100

mean((v_plant_match-v_actual_mtp) ./ v_actual_mtp,2) * 100
std((v_plant_match-v_actual_mtp) ./ v_actual_mtp,[],2)*100




mtp_vals = [24.74 ,53.23,64.47;...
        26.91, 54.37 ,69.34;...
        27.48,54.43,65.32];
    
perc_mtp_dors =[ mean((mtp_vals(:,2)-mtp_vals(:,1))./(mtp_vals(:,3)-mtp_vals(:,1))) std((mtp_vals(:,2)-mtp_vals(:,1))./(mtp_vals(:,3)-mtp_vals(:,1)))]
%%


for i = 1:ntrials
figure;
for dd = 1:3
    nfrs = length( mom_arm_lev{i}(dd,:));
    ind_plot = ind_prop(i,1):(ind_prop(i,2)-1);
    time_vals = ind_plot/ind_prop(i,2) * 100;
plot(time_vals, mom_arm_lev{i}(dd,ind_plot)','color',cmap(dd*15,:))
hold on; plot( time_vals, mom_arm_act{i}(dd,ind_plot)','color',cmap(dd*15,:),'linestyle','--')
end

plot(time_vals, norm3d( mom_arm_lev{i}(:,ind_plot)'),'k')
hold on; plot( time_vals,norm3d(mom_arm_act{i}(:,ind_plot)'),'--k')

legend('x - lever','x-actual','y-lever','y-actual','z-lever','z-actual','magnitude-lever','magnitude-- actual')
ylabel('Moment arm [m]')
% xlim([30 60])
end


for i = 1:ntrials
figure;
for dd = 1
    nfrs = length( mom_arm_lev{i}(dd,:));
    ind_plot = ind_prop(i,1):(ind_prop(i,2)-1);
    time_vals = ind_plot/ind_prop(i,2) * 100;
plot(time_vals,ank_pow_lev{i}(dd,ind_plot)','color',cmap(dd*15,:))
hold on; plot( time_vals, ank_pow_act{i}(dd,ind_plot)','color',cmap(dd*15,:),'linestyle','--')

end
end
%% make a wrist viz file with the fixed lever and the mobile arch
bone_pairs = bone_pairs(1:3,:);

% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);
for i = 1:ntrials
    
    
    anim_dir = fullfile(subj_dir,trialname_list{i},'POS','Global_LeverComparison',filesep);
    rigidiv_folder = fullfile(anim_dir,'rigidiv',filesep);
    
    % set up the GRF animation folder
    arch_style = 'arch_P%i_F1.iv';
    arch_folder = fullfile(anim_dir,'Archdir',filesep);
    
    if exist(anim_dir,'dir')==0;     mkdir(anim_dir);     end
    if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end
    if exist(arch_folder,'dir')==0;  mkdir(arch_folder);  end
    
        first_frX = ind_prop(i,1);%+match_arch(i);%struct_data(i).cropFrsXray(1);
        end_frX   = ind_prop(i,2);%struct_data(i).cropFrsXray(2);
    
        
        % fix the bones relative to a bone; replace any nan frames with the
        % first or last tracked frames
        T_fix = repmat(eye(4,4),1,1,end_frX-first_frX+1);%struct_data(i).Tm.nav(:,:,first_frX:end_frX);%
        ind_nan = find(isnan(T_fix(1,1,:)))';
        ind_nonan = find(~isnan(T_fix(1,1,:)))';
        for in = ind_nan
            if in < (end_frX-first_frX+1)/2
                T_fix(:,:,in) = T_fix(:,:,ind_nonan(1));
            else
                T_fix(:,:,in) = T_fix(:,:,ind_nonan(end));
            end
        end
        
        
    for bn = 1:nBones
        ivstring = createInventorHeader();
        % make the linked iv file
        jj = find(strcmp(bone_pairs(:,2),bone_list{bn})); % find the index of the bone pair to color the bone the same as the arrow
        if isempty(jj)
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7] ,0.7)];
        else
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),cmap(jj*floor(64/nBP),:) ,0.7)];
        end
        fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        
        ivstring = createInventorHeader();
        % make the linked iv file with the fixed mt1
        jj = find(strcmp(bone_pairs(:,2),bone_list{bn})); % find the index of the bone pair to color the bone the same as the arrow
        if isempty(jj)
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7] ,0.3)];
        else
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),cmap(jj*floor(64/nBP),:) ,0.3)];
        end
        fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '_fixed.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
    for bn = 1:nBones
        
    frWV = 1; % for the moving bones
        for fr = first_frX:end_frX
            Ttemp = invTranspose(T_fix(:,:,frWV)) * struct_Tm(i).(bone_list{bn})(:,:,fr);
            
            if any(isnan(Ttemp)) % wrist viz likes 1's where nans are
                Tw.(bone_list{bn})(:,:,frWV) = ones(4,4);
            else
                Ttemp(1:3,4) = Ttemp(1:3,4)*1000;
                Tw.(bone_list{bn})(:,:,frWV) = Ttemp;
            end
            frWV = frWV+1;
        end
        
        write_RTp(bone_list{bn} , Tw.(bone_list{bn}), anim_dir)
        
        % for the fixed bones
         frWV = 1;
        for fr = first_frX:end_frX
%             Ttemp = invTranspose(T_fix(:,:,frWV)) * struct_Tm(i).mt1(:,:,fr)*invTranspose(struct_Tm(i).(bone_list{bn})(:,:,first_frX))*struct_Tm(i).(bone_list{bn})(:,:,fr);%* struct_Tm(i).mt1(:,:,fr);
            
            if bn ~= 10
                Ttemp = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) ) * struct_Tm(i).(bone_list{bn})(:,:,first_frX) ;
            else
                Ttemp = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) )  * struct_Tm(i).(bone_list{tal_ind})(:,:,first_frX) * invTranspose( struct_Tm(i).(bone_list{tal_ind})(:,:,fr)) * struct_Tm(i).(bone_list{bn})(:,:,fr);
%                 ind_nonan = find(~isnan(squeeze(struct_Tm(i).(bone_list{bn})(1,1,:))));
%                 [mm,ii] = min(abs(ind_nonan-max_mtp(i)));
%                 new_mtp = ind_nonan(ii);
%                 warning('offset mtp maximum point for tibia transform')
%                 % use this to get the tibia fixed to the most plantarflexed
%                 % position
%                 Ttemp = struct_Tm(i).mt1(:,:,fr)* invTranspose( struct_Tm(i).mt1(:,:,first_frX) )  * struct_Tm(i).(bone_list{tal_ind})(:,:,first_frX) * invTranspose( struct_Tm(i).(bone_list{tal_ind})(:,:,new_mtp)) * struct_Tm(i).(bone_list{bn})(:,:,new_mtp);
%                
            end
            
            if any(isnan(Ttemp)) % wrist viz likes 1's where nans are
                Twf.(bone_list{bn})(:,:,frWV) = ones(4,4);
            else
                Ttemp(1:3,4) = Ttemp(1:3,4)*1000;
                Twf.(bone_list{bn})(:,:,frWV) = Ttemp;
            end
            frWV = frWV+1;
        end
        bone_list_fixed{bn,1} = [bone_list{bn} '_fixed'];
        write_RTp(bone_list_fixed{bn,1} , Twf.(bone_list{bn}), anim_dir)
    end
    
    
        frWV = 1;
    for fr = first_frX:end_frX
        
        arch_filename = sprintf(arch_style,frWV);
        %         frX = frX+1;
        ivstring = createInventorHeader();
%         
%         
%         %         ************** HELICAL AXES **************************
%         for bp = ind_bp%1:nBP
%             if any(isnan(pos_TC_lever{i}(:,fr)),'all')
%                 continue
%             end
%             
%             s =  pos_TC_lever{i}(:,fr)*1000;%transformPoints(T_fix(:,:,fr),struct_heli(i).trials(bp).s(:,fr),-1)*1000;
%             n = axis_TC_lever{i}(:,fr);%transformVectors(T_fix(:,:,fr),struct_heli(i).trials(bp).n(:,fr),-1);
%             
%             ivstring = [ivstring createInventorArrow(s,...
%                 n,...
%                200,...
%                 2,cmap(jj*floor(64/nBP),:),0.8)];
%             
%             
%              s =  pos_TC_m{i}(:,fr)*1000;%transformPoints(T_fix(:,:,fr),struct_heli(i).trials(bp).s(:,fr),-1)*1000;
%             n = axis_TC_m{i}(:,fr);%transformVectors(T_fix(:,:,fr),struct_heli(i).trials(bp).n(:,fr),-1);
%             
%             ivstring = [ivstring createInventorArrow(s,...
%                 n,...
%                200,...
%                 2,cmap(jj*floor(64/nBP),:),0.4)];
%         end
%         
%         %*********FORCE VECTORS****************
%         for fp = 1:2 % each force plate
%             COP = transformPoints(T_fix(:,:,frWV),struct_force(i).COP{fp}(:,fr),-1) *1000;% transformPoints(coreg,struct_data(i).force_data(fp).globCOP(:,fr)*1000,-1);
%             force = transformVectors(T_fix(:,:,frWV),struct_force(i).force{fp}(:,fr),-1);%transformVectors(coreg,struct_data(i).force_data(fp).globForce(:,fr),-1);
%             
%             %         ivstring = [ivstring createInventorSphere(COP,4,[0.2 0.2 0.8],0)];
%             ivstring = [ivstring createInventorArrow(COP,force,norm(force)/5,1.5,[1 1 1],0)];
%         end
%         
%         ivstring = [ivstring createInventorGlobalAxes];
        %************* WRITE FILE WITH ARROWS *****************
        fid = fopen(fullfile(arch_folder,arch_filename),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        frWV = frWV+1;
    end
    create_ini(0,0,1,1,arch_folder(1:end-1),'arch_P%d_F%d.iv',fullfile(anim_dir,[trialname_list{i} '.ini']))
    
    pos_text = write_pos([bone_list;bone_list_fixed],anim_dir,trialname_list{i});
    
    filename = fullfile(anim_dir, [trialname_list{i} '.pos']);
    
    fid = fopen(filename,'w'); % open the file to write
    fprintf(fid,pos_text);
    fclose(fid);
    
    fprintf('Animations are saved in %s.\n',anim_dir)
    clearvars('Tw','Twf')
end

%% Make the abstract figure
bone_pairs = bone_pairs(1:3,:);

% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';
iv_dir = fullfile(subj_dir,'Models','IV',filesep);
for i = 1:ntrials
    
    
    anim_dir = fullfile(subj_dir,trialname_list{i},'POS','Global_LeverComparison_Abstract',filesep);
    rigidiv_folder = fullfile(anim_dir,'rigidiv',filesep);
    
    % set up the GRF animation folder
    arch_style = 'arch_P%i_F1.iv';
    arch_folder = fullfile(anim_dir,'Archdir',filesep);
    
    if exist(anim_dir,'dir')==0;     mkdir(anim_dir);     end
    if exist(rigidiv_folder,'dir')==0;  mkdir(rigidiv_folder);  end
    if exist(arch_folder,'dir')==0;  mkdir(arch_folder);  end
    
        first_frX = ind_prop(i,1);%+match_arch(i);%struct_data(i).cropFrsXray(1);
        end_frX   = ind_prop(i,2);%struct_data(i).cropFrsXray(2);
    
        
        % fix the bones relative to a bone; replace any nan frames with the
        % first or last tracked frames
        T_fix = repmat(eye(4,4),1,1,end_frX-first_frX+1);%struct_data(i).Tm.nav(:,:,first_frX:end_frX);%
        ind_nan = find(isnan(T_fix(1,1,:)))';
        ind_nonan = find(~isnan(T_fix(1,1,:)))';
        for in = ind_nan
            if in < (end_frX-first_frX+1)/2
                T_fix(:,:,in) = T_fix(:,:,ind_nonan(1));
            else
                T_fix(:,:,in) = T_fix(:,:,ind_nonan(end));
            end
        end
        
        
    for bn = [1:3,5:nBones]
        ivstring = createInventorHeader();
        
        % make the linked iv file
        tr = 0.7;% transparency
        
        % moving bones
        jj = find(strcmp(bone_pairs(:,2),bone_list{bn})); % find the index of the bone pair to color the bone the same as the arrow
        if isempty(jj)
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7] ,tr)];
        else
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),cmap(jj*floor(64/nBP),:) ,tr-0.2)];
        end
        fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        
        tr = 0.3;
        ivstring = createInventorHeader();
        % make the linked iv file with the fixed mt1
        jj = find(strcmp(bone_pairs(:,2),bone_list{bn})); % find the index of the bone pair to color the bone the same as the arrow
        if isempty(jj)
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7] ,tr)];
        else
            ivstring = [ivstring createInventorLink([iv_dir subj_name '_'  bone_list{bn} '_aligned.iv'],eye(3,3),zeros(3,1),cmap(jj*floor(64/nBP),:) ,tr)];
        end
        fid = fopen(fullfile(rigidiv_folder,[bone_list{bn} '_fixed.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
    for bn = [1:3,5:nBones]
        
    frWV = 1; % for the moving bones
        for fr = [max_mtp(i),max_mtp(i)]%first_frX:end_frX
            Ttemp = invTranspose(T_fix(:,:,frWV)) * struct_Tm(i).(bone_list{bn})(:,:,fr);
            
            if any(isnan(Ttemp)) % wrist viz likes 1's where nans are
                Tw.(bone_list{bn})(:,:,frWV) = ones(4,4);
            else
                Ttemp(1:3,4) = Ttemp(1:3,4)*1000;
                Tw.(bone_list{bn})(:,:,frWV) = Ttemp;
            end
            frWV = frWV+1;
        end
        
        write_RTp(bone_list{bn} , Tw.(bone_list{bn}), anim_dir)
        
        % for the fixed bones
        frWV = 1;
        for fr = [match_arch(i)+first_frX,max_mtp(i)]%first_frX:end_frX
            %             Ttemp = invTranspose(T_fix(:,:,frWV)) * struct_Tm(i).mt1(:,:,fr)*invTranspose(struct_Tm(i).(bone_list{bn})(:,:,first_frX))*struct_Tm(i).(bone_list{bn})(:,:,fr);%* struct_Tm(i).mt1(:,:,fr);
            if bn == 8 % ph1
                Ttemp = struct_Tm(i).ph1(:,:,fr);
            elseif bn ~= 10 % not the tibia
                Ttemp = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) ) * struct_Tm(i).(bone_list{bn})(:,:,first_frX) ;
                
            else
%                 Ttemp = struct_Tm(i).mt1(:,:,fr)  * invTranspose( struct_Tm(i).mt1(:,:,first_frX) )  * struct_Tm(i).(bone_list{tal_ind})(:,:,first_frX) * invTranspose( struct_Tm(i).(bone_list{tal_ind})(:,:,fr)) * struct_Tm(i).(bone_list{bn})(:,:,fr);
                ind_nonan = find(~isnan(squeeze(struct_Tm(i).(bone_list{bn})(1,1,:))));
                [mm,ii] = min(abs(ind_nonan-max_mtp(i)));
                new_mtp = ind_nonan(ii);
                warning('offset mtp maximum point for tibia transform')
                % use this to get the tibia fixed to the most plantarflexed
                % position
                Ttemp = struct_Tm(i).mt1(:,:,fr)* invTranspose( struct_Tm(i).mt1(:,:,first_frX) )  * struct_Tm(i).(bone_list{tal_ind})(:,:,first_frX) * invTranspose( struct_Tm(i).(bone_list{tal_ind})(:,:,new_mtp)) * struct_Tm(i).(bone_list{bn})(:,:,new_mtp);
               
            end
            
            if any(isnan(Ttemp)) % wrist viz likes 1's where nans are
                Twf.(bone_list{bn})(:,:,frWV) = ones(4,4);
            else
                Ttemp(1:3,4) = Ttemp(1:3,4)*1000;
                Twf.(bone_list{bn})(:,:,frWV) = Ttemp;
            end
            frWV = frWV+1;
        end
        bone_list_fixed{bn,1} = [bone_list{bn} '_fixed'];
        write_RTp(bone_list_fixed{bn,1} , Twf.(bone_list{bn}), anim_dir)
    end
    
    
        frWV = 1;
    for fr = first_frX:end_frX
        
        arch_filename = sprintf(arch_style,frWV);
        %         frX = frX+1;
        ivstring = createInventorHeader();
%         
        %************* WRITE FILE WITH ARROWS *****************
        fid = fopen(fullfile(arch_folder,arch_filename),'w');
        fprintf(fid,ivstring);
        fclose(fid);
        frWV = frWV+1;
    end
    create_ini(0,0,1,1,arch_folder(1:end-1),'arch_P%d_F%d.iv',fullfile(anim_dir,[trialname_list{i} '.ini']))
    
    pos_text = write_pos([bone_list;bone_list_fixed],anim_dir,trialname_list{i});
    
    filename = fullfile(anim_dir, [trialname_list{i} '.pos']);
    
    fid = fopen(filename,'w'); % open the file to write
    fprintf(fid,pos_text);
    fclose(fid);
    
    fprintf('Animations are saved in %s.\n',anim_dir)
    clearvars('Tw','Twf')
end

