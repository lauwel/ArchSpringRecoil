clear
clc

bone_dir = 'E:\ShapeModelling\Models\3Aligned Reduced\';
base_dir = 'E:\ShapeModelling\';




SMbone_list = {'cal','tal','mt1','com'}';

nSM = length(SMbone_list); % number of shape models;

% get the file names and set up the length for each shape model
for sm = 1:nSM-1 % don't need to do the combines one here
    
    
    SMbone = SMbone_list{sm};
    load(['E:\ShapeModelling\Models\3Aligned Reduced\bonestruct' SMbone '.mat'])
    
    nbones.(SMbone) = length(bonestruct(:).(SMbone));
    
    filenames.(SMbone) = {bonestruct(:).(SMbone).filename}';
    
    for bn = 1:nbones.(SMbone)
        path_dec = strsplit(bonestruct.(SMbone)(bn).filename,filesep);
        bone_list.(SMbone){bn,1}= path_dec{end}(1:end-3);
        
         spl_str = strsplit(bone_list.(SMbone){bn},'_');
         subjs.(SMbone){bn,1} = spl_str{1};
    end
   
    ind = round(linspace(1,64,nbones.(SMbone)));
     cmap = colormap('summer');
     cmapSM.(SMbone) =  cmap(ind,:);
    
    
end

[subjs.com,Ia,Ib] = intersect(subjs.tal,subjs.cal);
bone_list.com = [bone_list.tal(Ia), bone_list.cal(Ib)];
nbones.com = length(bone_list.com);
filenames.com = [filenames.tal(Ia), filenames.cal(Ib)];

    % get the reference bone
ref_bone.cal = 9; % SOLH005
ref_bone.tal = 9; % SOLH005
ref_bone.mt1 = 4; % PRFC006
ref_bone.com = 9; % SOLH005

%% create the wristviz file to look at all the bones


for sm = 1:nSM-1 % for each shape model
    SMbone = SMbone_list{sm};
    
    load(['E:\ShapeModelling\Models\3Aligned Reduced\bonestruct' SMbone '.mat'])

    ivstring = createInventorHeader();
    for bn = 1:nbones.(SMbone)
        
        if bn < nbones.(SMbone)/2
            i1 = bn;
            i2 = 1;
        else
            i1 = bn - floor(nbones.(SMbone)/2);
            i2 = 2;
        end
        
        
        T_al = invTranspose(bonestruct.(SMbone)(bn).T_Aligned);
        
        ivstring = [ivstring createInventorLink(bonestruct.(SMbone)(bn).filename,T_al(1:3,1:3),T_al(1:3,4)', cmapSM.(SMbone)(bn,:),0.8)];
        ivstring = [ivstring createInventorLink(bonestruct.(SMbone)(bn).filename,T_al(1:3,1:3),T_al(1:3,4)'+[i1*50,0,100*i2], cmapSM.(SMbone)(bn,:),0.2)];
        
        
    end
    ivstring = [ivstring createInventorGlobalAxes()];
    
    fid = fopen([bone_dir SMbone '_anatAligned.iv'],'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    fprintf('Wrote IV file to : %s\n',[bone_dir SMbone '_anatAligned.iv'])
end



%% Compute the CPD about a reference subject's bone


for sm = 1:nSM-1 % for each shape model
    SMbone = SMbone_list{sm};
    
    load(['E:\ShapeModelling\Models\3Aligned Reduced\bonestruct' SMbone '.mat'])
    
  
    
    T_ACS = bonestruct.(SMbone)(ref_bone.(SMbone)).T_Aligned; % gives the aligned inertial axes
    ref_bone_pts = bonestruct.(SMbone)(ref_bone.(SMbone)).pts;
    ref_pts_anat = transformPoints(T_ACS,ref_bone_pts,-1);
    
    for bn = 1:nbones.(SMbone)
        
        
        if bn == ref_bone.(SMbone)
            pts_CT = ref_bone_pts;
            C = 1:length(ref_bone_pts);
        else
            
            T_ACS = bonestruct.(SMbone)(bn).T_Aligned; % gives the aligned inertial axes
            bone_pts = bonestruct.(SMbone)(bn).pts;
            pts_anat = transformPoints(T_ACS,bone_pts,-1);
            
            %         CUDA implementation
            omega = 0.1; beta = 2; lambda = 3; maxIter = 200; tol = 5e-5;
            %         Parameter lambda represents the trade off between data fitting and smoothness regularization.
            %       Parameter beta reflects the strength of interaction between points. Small values of produce locally smooth transformation, while large
            % values of beta correspond to nearly pure translation transformation
            [pts,C] = cpd_cuda(pts_anat,ref_pts_anat,omega, beta, lambda, maxIter, tol);
            
            pts_CT = transformPoints(T_ACS,pts);
            
        end
        
        
        save(fullfile(bone_dir,'CPD',[ bone_list.(SMbone){bn} '_CPD.mat']),'pts_CT','C');
        
    end
end


%% look at all the cpd bones

for sm = 1:nSM-1 % for each shape model
    SMbone = SMbone_list{sm};
    
    
    load(['E:\ShapeModelling\Models\3Aligned Reduced\bonestruct' SMbone '.mat'])
    
    T_ACS = bonestruct.(SMbone)(ref_bone.(SMbone)).T_Aligned; % gives the aligned inertial axes
    ref_bone_pts = bonestruct.(SMbone)(ref_bone.(SMbone)).pts;
    ref_pts_anat = transformPoints(T_ACS,ref_bone_pts,-1);
    ref_bone_cns = bonestruct.(SMbone)(ref_bone.(SMbone)).cnt;
    
    % create the wristviz file to look at all the bones
    
    corr_pts_all_X = [];
    corr_pts_all_Y = [];
    corr_pts_all_Z = [];
    ivstring = createInventorHeader();
    for bn = 1:nbones.(SMbone)
       
        
        % load the pts_CT and correspondence vectors (pts_CT is deformed
        % mesh)
        load(fullfile(bone_dir,'CPD',[ bone_list.(SMbone){bn} '_CPD.mat']));
        cpd_iv_filename = fullfile(bone_dir,'CPD',[ bone_list.(SMbone){bn} '_CPD.iv']);
        
        % the points of the individual bones
        pts = bonestruct.(SMbone)(bn).pts;
        T_ACS = bonestruct.(SMbone)(bn).T_Aligned; % gives the aligned inertial axes
        
        %corresponding points
        corr_pts = pts(C,:);
%         fprintf('size is : %i x %i\n',size(corr_pts))
        corr_pts_anat = transformPoints(T_ACS,corr_pts,-1); % in the anatomical co sys
        
%         fprintf('anat size is : %i x %i\n',size(corr_pts_anat))
        % write the CPD .iv file in CT space
        
        patch2iv(corr_pts, ref_bone_cns,cpd_iv_filename,[0.3 0 0.3])
        
%         fprintf('Before corr_pts_all_X size is : %i x %i\n',size(corr_pts_all_X))
        corr_pts_all_X(:,bn) = corr_pts_anat(:,1);
        corr_pts_all_Y(:,bn) = corr_pts_anat(:,2);
        corr_pts_all_Z(:,bn) = corr_pts_anat(:,3);
        
        if bn < nbones.(SMbone)/2
            i1 = bn;
            i2 = 1;
        else
            i1 = bn - floor(nbones.(SMbone)/2);
            i2 = 2;
        end
        
        
        T_al = invTranspose(bonestruct.(SMbone)(bn).T_Aligned);
        
        ivstring = [ivstring createInventorLink(cpd_iv_filename,T_al(1:3,1:3),T_al(1:3,4)', cmapSM.(SMbone)(bn,:),0.8)];
        ivstring = [ivstring createInventorLink(cpd_iv_filename,T_al(1:3,1:3),T_al(1:3,4)'+[i1*50,0,100*i2], cmapSM.(SMbone)(bn,:),0.2)];
        
        
    end
    ivstring = [ivstring createInventorGlobalAxes()];
    
    fid = fopen([bone_dir SMbone '_anatAlignedCPD.iv'],'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    
    mean_mesh.(SMbone)(:,1) = mean(corr_pts_all_X,2);
    mean_mesh.(SMbone)(:,2) = mean(corr_pts_all_Y,2);
    mean_mesh.(SMbone)(:,3) = mean(corr_pts_all_Z,2);
    
    
    figure; h = plot3quick_scatter(mean_mesh.(SMbone)','k'); h.Marker = '.';
    axis equal
    axis off
    % write the mean mesh to a file
    patch2iv( mean_mesh.(SMbone), ref_bone_cns,fullfile(bone_dir,'CPD',['meanBone_' SMbone '_CPD.iv']),[0.7 0.1 0.7])
end
%% load the mean mesh (newly renewed in geomagics) and rerun the correspondences
% save over the mean Bone in Geomagics

for sm = 1:nSM-1 % for each shape model
    SMbone = SMbone_list{sm};
    
    
    
    load(['E:\ShapeModelling\Models\3Aligned Reduced\bonestruct' SMbone '.mat'])
    
    
    [ref_bone_pts,ref_cns] = read_vrml_fast(fullfile(bone_dir,'CPD',['meanBone_' SMbone '_CPD_GM.iv']));
    ref_cns = ref_cns(:,1:3)+1;

    % create the ivstring to make the iv file with all the bones in it
    ivstring = createInventorHeader();
    
    for bn = 1:nbones.(SMbone)
        
        T_ACS = bonestruct.(SMbone)(bn).T_Aligned; % gives the aligned inertial axes
        bone_pts = bonestruct.(SMbone)(bn).pts;
        pts_anat = transformPoints(T_ACS,bone_pts,-1);
        
        %         CUDA implementation
        omega = 0.1; beta = 2; lambda = 3; maxIter = 700; tol = 1e-7;
        %         Parameter lambda represents the trade off between data fitting and smoothness regularization.
        %       Parameter beta reflects the strength of interaction between points. Small values of produce locally smooth transformation, while large
        % values of beta correspond to nearly pure translation transformation
        [pts_def,C] = cpd_cuda(pts_anat,ref_bone_pts,omega, beta, lambda, maxIter, tol);
        
        pts_CT = transformPoints(T_ACS,pts_def);
        
        % write the CPD .iv file in CT space
        
        cpd_iv_filename = fullfile(bone_dir,'CPD',[bone_list.(SMbone){bn} '_CPD.iv']);
        patch2iv(pts_def, ref_cns,cpd_iv_filename,[0.3 0 0.3])
%         
%         for fnum = 3:7
%             h = figure(fnum);
%             h.Position = [680 200 1000 800];
%             legend('deformed with cpd','raw points')
%             ha = gca;
%             view([-238,9])
%             savefig(h,['E:\ShapeModelling\InformationAndTesting\CPD parameter effects\' ha.Title.String{1} '.fig'])
%         end
        
        
        %for plotting;
        if bn < nbones.(SMbone)/2
            i1 = bn;
            i2 = 1;
        else
            i1 = bn - floor(nbones.(SMbone)/2);
            i2 = 2;
        end
        
        
        T_al = invTranspose(T_ACS);
        
        ivstring = [ivstring createInventorLink(cpd_iv_filename,T_al(1:3,1:3),T_al(1:3,4)', cmapSM.(SMbone)(bn,:),0.8)]; % the transparent ones at the ACS
        ivstring = [ivstring createInventorLink(cpd_iv_filename,T_al(1:3,1:3),T_al(1:3,4)'+[i1*50,0,100*i2], cmapSM.(SMbone)(bn,:),0.2)];
        
        % save the data
        save(fullfile(bone_dir,'CPD',[bone_list.(SMbone){bn} '_CPD.mat']),'pts_def','pts_CT','C');
        
    end
    
    ivstring = [ivstring createInventorGlobalAxes()];
    ivstring = [ivstring createInventorLink(fullfile(bone_dir,'CPD',['meanBone_' SMbone '_CPD.iv']),eye(3,3),[0 0 0],[0.5 0.5 0.5],.1)];
    fid = fopen(['E:\ShapeModelling\Models\3Aligned Reduced\' SMbone '_anatAlignedCPD.iv'],'w');
    fprintf(fid,ivstring);
    fclose(fid);
    
    clear bonestruct
end

clearvars('ref_cns')
%% Prepare to create the shape model - procrustes



for sm = 1:nSM-1 % for each shape model
    SMbone = SMbone_list{sm};
    
    load(['E:\ShapeModelling\Models\3Aligned Reduced\bonestruct' SMbone '.mat'])
    
    [ref_bone_pts,ref_cns.(SMbone)] = read_vrml_fast(fullfile(bone_dir,'CPD',['meanBone_' SMbone '_CPD_GM.iv']));
    ref_cns.(SMbone) = ref_cns.(SMbone)(:,1:3)+1;
    
    figure(sm);
    
        h = plot3quick_scatter(ref_bone_pts'); h.Marker = '.';h.MarkerFaceColor = 'k'; hold on;
    for bn = 1:nbones.(SMbone)
      
        
        cpd_mat_filename = fullfile(bone_dir,'CPD',[bone_list.(SMbone){bn} '_CPD.mat']);
        load(cpd_mat_filename); % loads the deformed points, CT points and correspondences
        
        npts = size(pts_def,1);
        
        [~,pts_scaled] = procrustes(ref_bone_pts,pts_def,'Scaling',true);
        
        h = plot3quick_scatter(pts_scaled'); h.Marker = '.';hold on;
        axis equal
        
        pca_mat.(SMbone)(bn,:) = reshape(pts_scaled',1,npts*3);
    end
    
    % calculate PCA
    [U.(SMbone),Z.(SMbone),~,~,PctExplained.(SMbone),~] = pca(pca_mat.(SMbone),'Economy',true);
    
    
    % write the points to a file for the mean bone
    npts = length(pca_mat.(SMbone))/3;
    
    
    mean_bone.(SMbone).pts = reshape(mean(pca_mat.(SMbone)),3,npts)';
    mean_bone.(SMbone).cns = ref_cns.(SMbone);
    save(fullfile(base_dir,'Results',['final_mean_' SMbone '.mat']),'mean_bone');
    patch2iv(mean_bone.(SMbone).pts,ref_cns.(SMbone),fullfile(base_dir,'Results',['final_mean_' SMbone '.iv']))
    
end
%%
SMbone = 'com'; % combined talus/calc
ref_cns.com = [ref_cns.tal;(max(ref_cns.tal) +ref_cns.cal)];

npts_cal = size(pca_mat.cal(Ia,:),2); % using the points that were the same from tal and cal
npts_tal = size(pca_mat.tal(Ia,:),2);
% add a negative y and z component to the calcaneus for visualization 
offset_cal = repmat([0 -35 -50],nbones.com,npts_cal/3);
offset_tal = repmat([0 0 0],nbones.com,npts_tal/3);

pca_mat.com = [pca_mat.tal(Ib,:),pca_mat.cal(Ia,:)]+[offset_tal,offset_cal];
[U.com,Z.com,~,~,PctExplained.com,~] = pca(pca_mat.com,'Economy',true);


npts = size(pca_mat.com,2);

mean_bone.(SMbone).pts = reshape(mean(pca_mat.(SMbone)),3,npts/3)';
mean_bone.(SMbone).cns = ref_cns.(SMbone);
save(fullfile(base_dir,'Results',['final_mean_' SMbone '.mat']),'mean_bone');
    
patch2iv(mean_bone.(SMbone).pts,ref_cns.(SMbone),fullfile(base_dir,'Results',['final_mean_' SMbone '.iv']))

save([base_dir 'Results/pca_results.mat'],'pca_mat','U','Z','PctExplained','mean_bone');

%% Animate the shape model

% write to server:
% base_dir = 'P:/Personnel/LaurenWelte/ShapeModelling/';

load('E:/ShapeModelling/Results/pca_results.mat');

spread = -3:0.2:3;
nbv = length(spread); % number of bones visualized across the shape model 

 cmapPC = colormap('summer');
ind = round(linspace(1,64,nbv));
 cmapPC =  cmapPC(ind,:);

for sm = 1:nSM-1 % for each shape model
    SMbone = SMbone_list{sm};
  
    anim_dir = fullfile(base_dir,'Results','Animation',SMbone,filesep);

    Xbar = mean(pca_mat.(SMbone));
    npts = length(Xbar);
    
    bone_style = 'bone_N%i_F%i.iv';
    save_bone_style = [SMbone '_PC%i_spr%0.1f'];
    
    pc = 1;
% for the number of principal components to animate - write until variation
% explained in PC is less than 4%
    while PctExplained.(SMbone)(pc) > 4   
        % prepare the directories
        anim_dirPC = fullfile(anim_dir,sprintf('PC%i',pc),filesep);
        rigidiv_dirPC = fullfile(anim_dirPC,'rigidiv',filesep);
        mesh_dir = fullfile(anim_dirPC,'mesh_dir',filesep);
        if ~exist(rigidiv_dirPC,'dir')
            mkdir(rigidiv_dirPC);
        end
         if ~exist(mesh_dir,'dir')
            mkdir(mesh_dir);
         end
        
        ind = 1;
        % create each bone across the spread of the principal components
        for sc = [1:nbv,nbv:-1:1]
            
            ptsv = Xbar + spread(sc).*std(Z.(SMbone)(:,pc)).*U.(SMbone)(:,pc)';
            
            pts_pc_3xn = reshape(ptsv,3,npts/3)';
            
            % create the ini files
            ini_filename = fullfile(mesh_dir,  sprintf(bone_style,ind,1));
            
            patch2iv(pts_pc_3xn,ref_cns.(SMbone),ini_filename, cmapPC(sc,:));
            
            ind = ind + 1;
            
            if ind <= (nbv + 1)% only write the first time through not the reverse
                save_bone_filename = sprintf(save_bone_style,pc,spread(sc));
                bone_save_dir = fullfile(base_dir,'Results','Bones',SMbone);
                if ~exist(bone_save_dir,'dir')
                    mkdir(bone_save_dir);
                end
                save_bone_filename = fullfile(bone_save_dir,[strrep(save_bone_filename,'.','_') '.iv']);
                
                patch2iv(pts_pc_3xn,ref_cns.(SMbone),save_bone_filename);
            end
        end
        
        % write the ini file
        bone_style_ini = strrep(bone_style,'%i','%d');
        
        create_ini(0,0,1,1,mesh_dir,bone_style_ini, [anim_dirPC 'PC' num2str(pc) '.ini'])
        
        
        % write the mean bone file to a link
        mean_file_name = fullfile(base_dir,'Results',['final_mean_' SMbone '.iv']);
        mean_fname_short = ['mean_' SMbone '.iv'];
        
        ivstring = createInventorHeader();
        ivstring = [ivstring createInventorLink(mean_file_name,eye(3,3),[0 0 0],[0.7 0.7 0.7],0.8)];
        
        fid = fopen(fullfile(rigidiv_dirPC,mean_fname_short),'w');
        fprintf(fid,ivstring);
        fclose(fid)
        
        
        write_pos({mean_fname_short(1:end-3)},anim_dirPC,sprintf('PC%i',pc));
        
        write_RTp({mean_fname_short(1:end-3)},{repmat(eye(4,4),1,1,nbv*2)},anim_dirPC);
        
        pc = pc+1;
        
    end
    
    
end
