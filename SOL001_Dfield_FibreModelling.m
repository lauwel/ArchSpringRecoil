%% Initialize the directories

% clear
close all
clc
% code to pull in Tony's data, create the distance fields and then do fibre
% wrapping for the hopping trials

plot_matlab_flag = 0; %1 = plot, 0 = don't
c = colormap('parula');

addpath('c:\Users\Lauren\Documents/MATLAB/Courses/CISC882')
direct = 'c:\Users\Lauren\Documents/School/PhD/Research/AllStudies/SOL001/Models/IV/';
addpath(genpath(direct));
addpath(genpath('c:\Users\Lauren\Documents/School/PhD/KDrive/Code/TheCollective/'))
addpath(genpath('c:\Users\Lauren\Documents/MATLAB/WorkingFolder/'))
direct_df = 'c:\Users\Lauren\Documents/School/PhD/Research/AllStudies/SOL001/Models/Dfield/';
direct_home = 'c:\Users\Lauren\Documents/School/PhD/Research/AllStudies/SOL001/';
direct_code = 'c:\Users\Lauren\Documents/git/FootFibreModelling/';
addpath(genpath(direct_code));

direct_mvt = 'c:\Users\Lauren\Documents/School/PhD/Research/AllStudies/SOL001/Models/POS/';

% windows directories for making wrist viz folders * NORMAL
% direct_anim_wind = 'E:\\PhDResearch\\Visualization\\';%'P:\\LaurenWelte\\Visualization\\';%
% link_folder_iv_wind = 'E:\PhDResearch\Visualization\IV\';%'P:\LaurenWelte\Visualization\IV\';%
% direct_anim = '/home/lauren/Documents/VirtualDesktopTemp/PhDResearch/Visualization/';%'/home/lauren/Documents/HMRC_NAS/Projects/LaurenWelte/Visualization';%

% Make server files: * FOR MIKE
direct_anim_wind = 'P:\\Personnel\\LaurenWelte\\Visualization\\';
link_folder_iv_wind = 'P:\Personnel\LaurenWelte\Visualization\IV\';
direct_anim = 'P:/Personnel/LaurenWelte/Visualization/';%'/home/lauren/Documents/HMRC_NAS/Projects/LaurenWelte/Visualization';
% %----------------
list_files_iv_raw = dir([direct '*.iv']);
nbones = length(list_files_iv_raw);

addpath(genpath('c:\Users\Lauren/Documents/MATLAB/TalarCoordinateSystemLW/')) % Talar coordinate system code
%% Calculate the distance fields
% for i = 1:nfiles
%     filename{i} = list_files_iv_raw(i).name;
% %     Dfield{i} = createDfield([direct filename{i}],0.25,direct_df);
% end
%
% save([direct_df 'SOL001_Dfield_Nov16_17.mat'],'Dfield')

%% Create the bone_structure
load([direct_df 'SOL001_Dfield_Nov16_17.mat'])
[num,text,raw] = xlsread([direct_home 'TrackedDataAndContactTimes.xlsx']);
% These are used to normalize the frames properly
start_frames = num(:,3);
start_track_frames = num(:,6);
num_track_frames = num(:,10);
end_frames = num(:,5);
end_track_frames = num(:,7);
trial_pct = round(num(:,11)*100);
start_pct = round(num(:,12)*100);
start_diff = start_track_frames-start_frames ;

fr_delay = abs([start_diff(1:3) - max(start_diff(1:3));start_diff(4:6) - max(start_diff(4:6))]);
end_diff = end_track_frames-end_frames ;
fr_end_delay = abs([end_diff(1:3) - max(end_diff(1:3));end_diff(4:6) - max(end_diff(4:6))]);

file_list_mvt = dir([direct_mvt '*norm*']);
n_mvt = length(file_list_mvt);


for i = 1:nbones
    
    [centr,~,~,~,eig_vec,~,~,~,~,~] = mass_properties(list_files_iv_raw(i).name);
    
    bone_struct(i).bone = list_files_iv_raw(i).name(1:3);
    [bone_struct(i).pts cns] = read_vrml_fast([direct list_files_iv_raw(i).name]);
    bone_struct(i).cns = cns+1;
%     bone_struct(i).dfield = Dfield{i};
    bone_struct(i).CT_inert_axes = [eig_vec, centr'; 0 0 0 1];
    
%     
%          [x0n, an, rn, d, sigmah, conv, Vx0n, Van, urn, GNlog, a, R0, R] =lscylinder(bone_struct(i).pts,centr',[0;0;-1],15,1,1);
%         figure
%         hold on; plotvector3(x0n,an*100,'k');
    
    
    % plot the bones
%     figure;
%     h = patch('Faces',bone_struct(i).cns(1:10:end,1:3),'Vertices',bone_struct(i).pts(:,1:3));
%     set(h,'facealpha',0.5,'Facecolor',c(2,:),'Edgecolor',c(2,:))
%     T = bone_struct(i).CT_inert_axes;
    % T(1:3,1:3) = T(1:3,1:3)';
%     plotPointsAndCoordSys1([],T)
    for j = 1:n_mvt % for all the motions
        temp_T = dlmread([direct_mvt,(file_list_mvt(j).name), '/Filtered/', bone_struct(i).bone, '_global.RTp']);
        n_fr = temp_T(1,1); % number of frames of tracked data
        % make the RTP file format into a series of 4x4x nframe matrices to
        % save into the structure
        T_44 = repmat(eye(4,4),1,1,n_fr);
        
        for k = 1:n_fr
            T_44(1:3,1:3,k) = temp_T((4*k-3)+1:(4*k),1:3);
            T_44(1:3,4,k) = temp_T(4*k+1,1:3)';
        end
        bone_struct(i).(file_list_mvt(j).name) = T_44;
        
        for k = 1:n_fr-1 % find all the helical axes of itself
                T_helb(:,:,k) =  T_44(:,:,k+1)*invTranspose(T_44(:,:,k));
                [phi(k),n(:,k),L(k),s(:,k)] = helical(T_helb(:,:,k));
                centr_xr = transformPoints(T_44(:,:,k),centr');
                % make it close to the centroid
                s_cent(:,k) = closestPointonVector(centr_xr,s(:,k),n(:,k));
        end
        
        bone_struct(i).([file_list_mvt(j).name '_n']) = n;
        
        bone_struct(i).([file_list_mvt(j).name '_s']) = s_cent;
        
        
        clearvars('T_helb','phi','n','L','s')
    end
    
    
end
% compute the new talar coordinate system
% bone_struct = makingTalarCoordSys(bone_struct);


%% Create the helical_structure
%determine the helical axes between the bones relative to the calcaneusref = 'cal'; ref_numR = 1;
mh1_ct = load([direct_anim '/IV/MH1_ct.stack'])';
ca_ct = load([direct_anim '/IV/CA_ct.stack'])';
cst_ct = load([direct_anim '/IV/CST_ct.stack'])';
talh_ct = load([direct_anim '/IV/tal_head.stack'])';

%
%         figure;
%     hold on
%     % plot the bones
%     h = patch('Faces',bone_struct(ref_numR).cns(1:10:end,1:3),'Vertices',bone_struct(ref_numR).pts(:,1:3));
%     set(h,'facealpha',0.5,'Facecolor',c(2,:),'Edgecolor',c(2,:))
%     T = bone_struct(ref_numR).CT_inert_axes;
%     plotPointsAndCoordSys1([],T)

% comparing bone A to the reference bone R, in xray space (x) and CT space
% (ct), between time frames i, and i+1 (ip1)


for j = 1:n_mvt
    
    indb = 1;
    
    for bn = 1:nbones % for all the reference bones
        
        ref = bone_struct(bn).bone;
%         ref_numR = bn;
        
        T_ACS = bone_struct(bn).CT_inert_axes;
        T_xR_ct= bone_struct(bn).(file_list_mvt(j).name);
        n_fr = size(T_xR_ct,3);
        
        for i = bn:nbones 
           
            if bn == i % i.e  with reference to itself, then input the global one instead
                  T_ACS = eye(4,4);
                  T_xR_ct = repmat(eye(4,4),1,1,n_fr);
                  ref = 'glb';
            else
                T_ACS = bone_struct(bn).CT_inert_axes;
                T_xR_ct= bone_struct(bn).(file_list_mvt(j).name);
                n_fr = size(T_xR_ct,3);
                ref = bone_struct(bn).bone;
            end
                % initialize all variables
                T_i = zeros(4,4,n_fr-1);
                T_ip1 = zeros(4,4,n_fr-1);
                T_hel = zeros(4,4,n_fr-1);
                T_ACS_xR = zeros(4,4,n_fr-1);
                %              clearvars('phi','n','L','s','T_hel','n_ct','s_ct','s_refx_ct','s_refx_xr','s_xr','n_xr');
                phi = zeros(1,n_fr-1);         n = zeros(3,n_fr-1);        L = zeros(1,n_fr-1);        s = zeros(3,n_fr-1);
                n_ct = zeros(3,n_fr-1);        s_ct = zeros(3,n_fr-1);     s_xr = zeros(3,n_fr-1);     n_xr = zeros(3,n_fr-1);
                s_refx_ct = zeros(3,n_fr-1);   s_refx_xr =zeros(3,n_fr-1); refx_cent_xr = zeros(3,n_fr-1);
                
                
                T_xA_ct = bone_struct(i).(file_list_mvt(j).name);
                
                
                cent1 = bone_struct(bn).CT_inert_axes(1:3,4); %reference bone centroid (CT space)
                cent2 = bone_struct(i).CT_inert_axes(1:3,4); % comp bone centroid (CT space)
                
                if  bn == i % i.e  in global, there is only one bone, so put it as close to the centroid as possible
                    cent1 = cent2;
                end
                
                refx_cent_ct = (cent1 + cent2)/ 2; % find the midpoint of the two and start the axis there
                
                for k = 1:n_fr-1
                    T_ACS_xR(:,:,k) = T_xR_ct(:,:,k) * T_ACS;
                    
                    
                    T_i(:,:,k) = invTranspose(T_xR_ct(:,:,k)) * T_xA_ct(:,:,k); % Register both frames to the reference bone
                    T_ip1(:,:,k) = invTranspose(T_xR_ct(:,:,k+1)) * T_xA_ct(:,:,k+1);
                    T_hel(:,:,k) =  T_ip1(:,:,k)*invTranspose(T_i(:,:,k)); % convert the helical axis matrix
                    
                    [phi(k),n(:,k),L(k),s(:,k)] = helical(T_hel(:,:,k));
                    
                    % transform the points into the xray space
                    
                    n_xr(:,k) = transformVectors(T_xR_ct(:,:,k),n(:,k));
                    s_xr(:,k) = transformPoints(T_xR_ct(:,:,k),s(:,k));
                    
                    % move the reference point into xray space and compute the
                    % closest point
                    refx_cent_xr(:,k) = transformPoints(T_xR_ct(:,:,k),refx_cent_ct);
                    
                    % find the closest point to the reference point that's on
                    % the vector
                    s_refx_xr(:,k) = closestPointonVector(refx_cent_xr(:,k),s_xr(:,k),n_xr(:,k));
                    s_refx_ct(:,k) = transformPoints(T_xR_ct(:,:,k),s_refx_xr(:,k),1); % move it back to ct
                    
                    
                    
                end
                
                mid_fr = floor(n_fr/2);
                firstFr = 1+fr_delay(j):mid_fr-2; % arrays containing the appropriate frames
                lastFr = mid_fr+2:n_fr-fr_end_delay(j)-1;
                % determine the average helical axis, proportional to phi, for
                % each half of the trial
                [avgHam{1},~,variat{1}] = calcAvgMotionLW([phi(firstFr)',n_xr(:,firstFr)',L(firstFr)',s_refx_xr(:,firstFr)']);
                [avgHam{2},~,variat{2}] = calcAvgMotionLW([phi(lastFr)',n_xr(:,lastFr)',L(lastFr)',s_refx_xr(:,lastFr)']);
                
                
                % save everything to the helical structure
                
                helical_struct(j).Name = (file_list_mvt(j).name);
                helical_struct(j).Trial(indb).ref_bone = ref;
                helical_struct(j).Trial(indb).comp_bone = bone_struct(i).bone;
                helical_struct(j).Trial(indb).RT = T_hel;
                helical_struct(j).Trial(indb).phi = phi;
                helical_struct(j).Trial(indb).n = n;
                helical_struct(j).Trial(indb).L = L;
                helical_struct(j).Trial(indb).s = s;
                helical_struct(j).Trial(indb).n_xr = n_xr;
                helical_struct(j).Trial(indb).s_xr = s_xr;
                helical_struct(j).Trial(indb).s_refx_xr = s_refx_xr;
                helical_struct(j).Trial(indb).s_refx_ct = s_refx_ct;
                
                for ax = 1:2 % average helical axis parameters
                    helical_struct(j).Trial(indb).avg_phi{ax} = avgHam{ax}(1);
                    helical_struct(j).Trial(indb).avg_n_xr{ax} = avgHam{ax}(2:4);
                    helical_struct(j).Trial(indb).avg_L{ax} = avgHam{ax}(5);
                    helical_struct(j).Trial(indb).avg_s_xr{ax} = avgHam{ax}(6:8);
                    
                end
                clearvars('avgHAM')
                            helical_struct(j).Trial(indb).HAM_var = variat;
                
            
            indb = indb+1;
           
        end
%         
%         % compute the mla markers in x_ray space then project into
%         % the inertial axes of the calcaneus, then compute MLA angle
%         %     figure;
%         if bn == 1 && i == 5 % mt1 ref to calc
%             T_xA_ct = bone_struct(5).(file_list_mvt(j).name);
%             for k = 1:n_fr-1
%                 mh1_xr(:,k) = transformPoints(T_xA_ct(:,:,k),mh1_ct);
%                 ca_xr(:,k) = transformPoints(T_xR_ct(:,:,k),ca_ct);
%                 cst_xr(:,k) = transformPoints(T_xR_ct(:,:,k),cst_ct);
%                 
%                 ca_cst(:,k) = ca_xr(:,k)-cst_xr(:,k);
%                 mh1_cst(:,k) = mh1_xr(:,k) - cst_xr(:,k);
%                 
%                 sag_plane_norm(:,k) = T_ACS_xR(1:3,3,k); % z component is the sagittal plane normal
%                 
%                 ca_cst_proj(:,k) = vecProjOn2VecPlane(ca_cst(:,k),sag_plane_norm(1:3,k));
%                 ca_cst_proj_norm(:,k) = ca_cst_proj(:,k)/norm(ca_cst_proj(:,k));
%                 mh1_cst_proj(:,k) = vecProjOn2VecPlane(mh1_cst(:,k),sag_plane_norm(1:3,k));
%                 mh1_cst_proj_norm(:,k) = mh1_cst_proj(:,k)/norm(mh1_cst_proj(:,k));
%                 
%                 dot_MLA(k) = dot(ca_cst_proj(:,k),mh1_cst_proj(:,k))/(norm(ca_cst_proj(:,k))*norm(mh1_cst_proj(:,k)));
%                 
%                 MLA(k) = acosd(dot_MLA(k));
%                 %         plot3quick(ca_xr(:,k),'r','o');
%                 % hold on;
%                 % plot3quick(mh1_xr(:,k),'k','o');
%                 % plot3quick(cst_xr(:,k),'b','o');
%                 % plotvector3(cst_xr(:,k),ca_cst(:,k),'k');
%                 % plotvector3(cst_xr(:,k),mh1_cst(:,k),'m');
%                 % plotvector3(ca_xr(:,k),sag_plane_norm(:,k),'c');
%                 % plotvector3(ca_xr(:,k),sag_plane_norm(:,k)*50,'c');
%                 % plotvector3(ca_xr(:,k),sag_plane_norm(:,k)*50,'m');
%                 % axis equal
%                 %
%                 
%                 % hold on
%                 % plotvector3([0 0 0],ca_cst(:,k),'k');
%                 % plotvector3([0 0 0],mh1_cst(:,k),'m');
%                 
%             end
%             
%             helical_struct(j).MLA = MLA;
%             clearvars('MLA')
%         end
    end
end

%% PLOTTING the bones and the co-ordinate system 
% for movement j = , CT SPACE


if plot_matlab_flag == 1
    comp = 4; % which helical axis pair to look at - 4 is cal- mt1
    j = 2;
    
    boneA       = helical_struct(j).Trial(comp).comp_bone;
    ref_numA    = findInStruct(bone_struct,'bone',boneA);
    boneR       = helical_struct(j).Trial(comp).ref_bone;
    ref_numR    = findInStruct(bone_struct,'bone',boneR);
    
    figure(1) % CT space
    hold on;
    h = patch('Faces',bone_struct(ref_numR).cns(1:10:end,1:3),'Vertices',bone_struct(ref_numR).pts);% CT space
    set(h,'facealpha',0.5,'Facecolor',c(2,:),'Edgecolor',c(2,:))
    h = patch('Faces',bone_struct(ref_numA).cns(1:10:end,1:3),'Vertices',bone_struct(ref_numA).pts);
    set(h,'facealpha',0.5,'Facecolor',c(30,:),'Edgecolor',c(30,:))
    hold on
    
    s_ct1 = helical_struct(j).Trial(comp).s(1:3,1)';
    n_ct1 = helical_struct(j).Trial(comp).n(1:3,1)';
    
    
    for i = 1:60
        T_axis = helical_struct(j).Trial(comp).RT(:,:,i);
        % plot the helical axis
        s_ct = helical_struct(j).Trial(comp).s(1:3,i)';
        n_ct = helical_struct(j).Trial(comp).n(1:3,i)'*100;
        
        
        h = plotvector3(s_ct,n_ct,'k');
        
        
        %         h = plotvector3(transformPoints(T_axis,s_ct1),transformVectors( T_axis,n_ct1)*100,'r');
        axis equal
    end
    % PLOT IN XRAY SPACE
    
    T_xR_ct= bone_struct(ref_numR).(file_list_mvt(j).name);
    T_xA_ct = bone_struct(ref_numA).(file_list_mvt(j).name);
    for i = 1%:65
        
        T_axis = helical_struct(j).Trial(comp).RT(:,:,i);
        
        % plot the helical axis - using transforms of ct co-ordinate
        s_xr = helical_struct(j).Trial(comp).s_xr(1:3,1);
        n_xr = helical_struct(j).Trial(comp).n_xr(1:3,1);
        
        s_ct = helical_struct(j).Trial(comp).s(1:3,i);
        n_ct = helical_struct(j).Trial(comp).n(1:3,i);
        
        s_ct1 = helical_struct(j).Trial(comp).s(1:3,1);
        n_ct1 = helical_struct(j).Trial(comp).n(1:3,1);
        
        figure(2) % x-ray space, referenced to first bone
        cla
        T_xiR_x1R =  T_xR_ct(:,:,i) * invTranspose(T_xR_ct(:,:,1));
        
        transBoneRpts = transformPoints( T_xR_ct(:,:,i) , bone_struct(ref_numR).pts',0 )';
        transBoneApts = transformPoints( T_xA_ct(:,:,i) , bone_struct(ref_numA).pts',0 )';
        
        h = patch('Faces',bone_struct(ref_numR).cns(1:10:end,1:3),'Vertices',transBoneRpts);% xray space
        set(h,'facealpha',0.5,'Facecolor',c(2,:),'Edgecolor',c(2,:))
        h = patch('Faces',bone_struct(ref_numA).cns(1:10:end,1:3),'Vertices',transBoneApts);
        set(h,'facealpha',0.5,'Facecolor',c(30,:),'Edgecolor',c(30,:))
        hold on
        %         h = plotvector3(transformPoints( T_axis,s_ct1),transformVectors( T_axis,n_ct1)*100,'r');
        %         set(h,'LineWidth',5)
        %         h = plotvector3(transformPoints(T_xR_ct(:,:,1) ,s_ct),transformVectors(T_xR_ct(:,:,1),n_ct),'r');
        h = plotvector3(transformPoints(T_xR_ct(:,:,i),s_ct),transformVectors(T_xR_ct(:,:,i),n_ct)*50,'b');
        set(h,'LineWidth',5)
        
        % making the plot nice
        hold off
        axis equal
        xlabel('X'); %xlim([-100 100])
        ylabel('Y') ; %ylim([-150 100])
        zlabel('Z'); %zlim([100 270])
        
        view([-115 12])
        title('Bones in xray space')
        drawnow
        pause(0.1)
        
    end
end

%% Wrist Viz 1 - Create the animation files to visualise each trial

% arrow = createInventorArrow(base,orient,length,width,color,transparency)
 
    col = [[57, 142, 31]/256;...
        1 1 1;...
        [142, 46, 178]/256];
            
for j = 1:length(helical_struct) % for every trial
    for i = 5%1:length(helical_struct(j).Trial) % for every bone combination
        trial_name = file_list_mvt(j).name;
        
        ham_style ='helicalaxis_P%i_F1.iv';
        folder_comp = sprintf('%sTO%s',helical_struct(j).Trial(i).comp_bone,helical_struct(j).Trial(i).ref_bone);
        
        trial_base_folder = fullfile(direct_anim,trial_name,folder_comp);
        ham_folder = fullfile(trial_base_folder,'HAdir');
        rigidiv_folder = fullfile(trial_base_folder,'rigidiv');
        
        n_fr = size(helical_struct(j).Trial(i).phi,2);
        mid_fr = floor(n_fr/2);
        % write all the helical axes of motion
        
        s_ref_xr = helical_struct(j).Trial(i).s_refx_xr;
        n_xr = helical_struct(j).Trial(i).n_xr;
        phi_k = helical_struct(j).Trial(i).phi;
        
          c = makeColourMap(col(1:3,:),n_fr);
          
        for k = 1:n_fr % how many frames of data
    
            ham_filename = sprintf(ham_style,k);
            
            % make the helical axis
            ivstring = createInventorHeader();
            % plot in Xray space
            ivstring = [ivstring createInventorArrow(s_ref_xr(:,k),n_xr(:,k),   150* phi_k(k) , 3,c(k,:),0)];
            ivstring = [ivstring createInventorSphere(s_ref_xr(:,k),    4,  c(k,:),0)];
            
            if ismember(k,[1+fr_delay(j):mid_fr-2]) == 1
                % is in the first half of the trial including an average
                % helical axis
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).avg_s_xr{1}',helical_struct(j).Trial(i).avg_n_xr{1}',150, 3,[0.1 0.1 0.5],0)];
            elseif ismember(k,[mid_fr+2:n_fr-fr_end_delay(j)]) == 1
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).avg_s_xr{2}',helical_struct(j).Trial(i).avg_n_xr{2}',150, 3,[0.1 0.1 0.5],0)];
            end
            
            % plot in x-ray but not at the talus intersectionfid
%             ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_xr(:,k),helical_struct(j).Trial(i).n_xr(:,k),150, 3,[1 0 0],0)];
%             ivstring = [ivstring createInventorSphere(helical_struct(j).Trial(i).s_xr(:,k),4,[1 0 0],0)];
           
%           (centroid,radius,color,transparency)
%             ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).comp_centr(:,k),helical_struct(j).Trial(i).n_xr(:,k),150, 3,[1 0 0],0)];
            % Write the extra components to a file
            if exist(ham_folder,'dir')==0
                mkdir(ham_folder);
            end
            fid = fopen(fullfile(ham_folder,ham_filename),'w');
            fprintf(fid,ivstring);
            fclose(fid);
        end
        
        % --- make the RTp files - one for each thing that needs animating
        
        write_RTp({bone_struct(:).bone} , {bone_struct(:).(trial_name)} , trial_base_folder)
        
        % --- Make the pos files
        % if _wind is in the var name, it's a path on windows
       
        local_folder_wind = [direct_anim_wind trial_name '\\' folder_comp '\\'];
        
        pos_text = write_pos({bone_struct(:).bone},local_folder_wind,trial_name);
        
        filename = fullfile(trial_base_folder, [trial_name '.pos']);
        fid = fopen(filename,'w'); % open the file to write
        fprintf(fid, pos_text);
        fclose(fid);
        
        % --------- make the ini file to animate the helical axes
        
        create_ini(0,0,1,1,[local_folder_wind 'HAdir'],'helicalaxis_P%%d_F%%d.iv',fullfile(trial_base_folder,[trial_name '.ini']))
        
        % -----------make the linked IV files
        
        link_trial_iv_wind = [local_folder_wind 'rigidiv'];
        
        for b = 1:length(bone_struct)
            ivstring = createInventorHeader();
            ivstring = [ivstring createInventorLink([link_folder_iv_wind bone_struct(b).bone '.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
            if exist(rigidiv_folder,'dir')==0
                mkdir(rigidiv_folder);
            end
            fid = fopen(fullfile(rigidiv_folder,[bone_struct(b).bone '.iv']),'w');
            fprintf(fid,ivstring);
            fclose(fid);
        end
        disp(['Saved file : ' filename '. '])
    end
end

%% Wrist Viz 2 - Create inertial axes for the bones
       
inert_ax_stat_filename = fullfile(direct_anim,'IV','inert_axes_tib.iv');
% make the inertial axes
b = 8;
ivstring = createInventorHeader();
%             ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_xr(:,k),helical_struct(j).Trial(i).n_xr(:,k),150, 3,[1 0 0],0)];
ivstring = [ivstring createInventorArrow(bone_struct(b).CT_inert_axes(1:3,4),bone_struct(b).CT_inert_axes(1:3,1),80, 3,c(17,:),0)];
ivstring = [ivstring createInventorArrow(bone_struct(b).CT_inert_axes(1:3,4),bone_struct(b).CT_inert_axes(1:3,2),80, 3,c(34,:),0)];
ivstring = [ivstring createInventorArrow(bone_struct(b).CT_inert_axes(1:3,4),bone_struct(b).CT_inert_axes(1:3,3),80, 3,c(51,:),0)];
% Write the extra components to a file

fid = fopen(inert_ax_stat_filename,'w');
fprintf(fid,ivstring);
fclose(fid);

%% Wrist Viz 3 - Butterflies - make helical axes for the CT 
jump_str = {'108','120','132','jog1','jog2','jog4'};
% --------------- write a ct helical axis ----------------
frame_cutoffs = {[4:28], [50:77]; [1:20], [57:76]; [2:21] [40:60]};
% 
%     col = [[29, 94, 26]/256;...
%         1 1 1;...
%         [88, 24, 112]/256];
    
    col = [[57, 142, 31]/256;...
        1 1 1;...
        [142, 46, 178]/256];
for j = 1:3%n_mvt
    

    for i = 5%1:length(helical_struct(j).Trial)
        
        hel_ax_stat_filename = fullfile(direct_anim,'Butterflies',['helical_axes' jump_str{j} '_' helical_struct(j).Trial(i).comp_bone 'TO' helical_struct(j).Trial(i).ref_bone '.iv']);
        
        temp = helical_struct(j).Trial(i);
        % [avgHAM,num] = calcAvgMotion([temp.phi', temp.n',temp.s',temp.L',checkNums)
        % make the helical axis
        ivstring = createInventorHeader();
        transpar = 0.5;
        
          phi_k = helical_struct(j).Trial(i).phi;
        if j < 4
            ind = 1;
            c = makeColourMap(col(1:3,:),frame_cutoffs{j,2}(end) );
            for k = frame_cutoffs{j,1} %6:23
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_refx_ct(:,k),helical_struct(j).Trial(i).n(:,k),100*phi_k(k), 3,c(k,:),transpar)];
                ind = ind +1;
            end
            
%             c = makeColourMap(col(2:3,:),frame_cutoffs{j,2}(end) - frame_cutoffs{j,2}(1)+1);
%             ind = 1;
            for k = frame_cutoffs{j,2} %40:60
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_refx_ct(:,k),helical_struct(j).Trial(i).n(:,k),100*phi_k(k), 3,c(k,:),transpar)];
                ind = ind + 1;
            end
        else
            c = colormap('parula');
            for k = 1:length(helical_struct(j).Trial(i).n)
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_refx_ct(:,k),helical_struct(j).Trial(i).n(:,k),100, 3,c(k,:),transpar)];
            end
        end
        
        for b = 1:length(bone_struct)
            ivstring = [ivstring createInventorLink([link_folder_iv_wind bone_struct(b).bone '.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
        end

        % Write the components to a file
        fid = fopen(hel_ax_stat_filename,'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
end

%% calculate the references



% --relative to the bone previous to it
% bone_line_list = {'cal','tal'; 'tal','nav';'nav','cmm'; 'cmm', 'mt1'};
% along_arch_folder = 'rel_to_each';
% --all relative to cal
% bone_line_list = {'cal','tal'; 'cal','nav';'cal','cmm'; 'cal', 'mt1'};
% along_arch_folder = 'rel_to_cal';
%-- all relative to tib
% bone_line_list = {'tib','cal';'tib','tal'; 'tib','nav';'tib','cmm'; 'tib', 'mt1'};
% along_arch_folder = 'rel_to_tib';
% --all relative to mt1
% bone_line_list = {'mt1','cal';'mt1','tal'; 'mt1','nav';'mt1','cmm'};
% along_arch_folder = 'rel_to_mt1';
%-- looking at the talo-cal-nav joint
bone_line_list = {'tal','cal'; 'tal','nav';'tal','cub'};
along_arch_folder = 'rel_to_tal';
% bone_line_list = {'tal','nav';'cal','cub';'nav','cub'};
% along_arch_folder = 'midtarsal';
%-- looking relative to global
% bone_line_list = {'glb','cal'; 'glb','tal';'glb','nav';'glb','cmm'; 'glb', 'mt1'};
% along_arch_folder = 'rel_to_glb';
% bone_line_list = {'cal', 'mt1'};
% along_arch_folder = 'rel_to_cal';
% % bone_line_list = {'tib','cal';'tib','nav'; 'tib','tal'};
% along_arch_folder = 'rel_to_tib_sarah';
% proj_bone_list = {'cal'; 'tal';'nav'; 'cmm'};
% proj_bone_list = {'cal';'cal';'cal';'cal';'cal'};
% proj_bone_list = {'tal';'tal';'tal';'tal';'tal'};
% proj_bone_list = {'mt1';'mt1';'mt1';'mt1';'mt1'};
% bone_line_list = {'tib','cal'; 'tal','nav';'cal','cub'};
% along_arch_folder = 'random_testing';

leg_str = {};
b_ind = [];
n_bonecomps = size(bone_line_list,1);
for i = 1:n_bonecomps % for the bone combinations in question
    rev_flag(i) = 0;
    % find the indices of the first bone and then use those indices to find
    % the bone of reference
    b1_ind = findInStruct(helical_struct(j).Trial,'ref_bone', bone_line_list{i,1});
    b2_ind = findInStruct(helical_struct(j).Trial(b1_ind),'comp_bone', bone_line_list{i,2});
    if isempty(b2_ind) == 1
        b1_ind = findInStruct(helical_struct(j).Trial,'comp_bone', bone_line_list{i,1});
        b2_ind = findInStruct(helical_struct(j).Trial(b1_ind),'ref_bone', bone_line_list{i,2});
        rev_flag(i) = 1; % if the reference bone and comp_bone are switched
    end
    
    b_ind(i) = b1_ind(b2_ind);
    
    % make a cell array of strings to put on the graphs
    leg_str = [leg_str, [bone_line_list{i,2} ' TO ' bone_line_list{i,1}] ]; 
    
end

%% Wrist Viz 4 - Make the wrist viz files for the arch line of helical axes

c1 = [179, 119, 239;...
    123, 114, 229;...
    229, 114, 170;...
    43, 219, 201]/256;

c1 = [108 182 255;...
    160 66 255;...
    255 128 192]/256;
%  c1 = [160 66 255;...
%     255 128 192]/256;   
% c1 = c([1,12,34,45,54],:);
% showColorMap(c1)
for j = 1:length(helical_struct) % for every trial
    % write to wrist viz
    
    trial_name = file_list_mvt(j).name;
    
    ham_style ='helicalaxis_P%i_F1.iv';
    trial_base_folder = fullfile(direct_anim,'AlongArch',trial_name,along_arch_folder);
    ham_folder = fullfile(trial_base_folder,'HAdir');
    rigidiv_folder = fullfile(trial_base_folder,'rigidiv');
    
    % Double check that the HA directory exists
    if exist(ham_folder,'dir')==0
        mkdir(ham_folder);
    end
    if exist(rigidiv_folder,'dir')==0
        mkdir(rigidiv_folder);
    end
    % write all the helical axes of motion for each pair in the same
    % file
    for k = 1:size(helical_struct(j).Trial(1).phi,2) % how many frames of data
        
        ham_filename = sprintf(ham_style,k);
        
        ivstring = createInventorHeader();
        
        for ii = 1:length(b_ind)
            i = b_ind(ii);
            phi_val = helical_struct(j).Trial(i).phi(k);
            
            if rev_flag(ii) == 0
                
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_refx_xr(:,k),helical_struct(j).Trial(i).n_xr(:,k),100*phi_val, 3,c1(ii,:),0)];
                ivstring = [ivstring createInventorSphere(helical_struct(j).Trial(i).s_refx_xr(:,k),4,c1(ii,:),0)];
%                 [0 0.1+0.2*ii 0.5]
                
            elseif rev_flag(ii) == 1
                
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_refx_xr(:,k),-helical_struct(j).Trial(i).n_xr(:,k),100*phi_val, 3,c1(ii,:),0)];
                ivstring = [ivstring createInventorSphere(helical_struct(j).Trial(i).s_refx_xr(:,k),4,c1(ii,:),0)];
                
            end
        end
        fid = fopen(fullfile(ham_folder,ham_filename),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
        warning('Think through the negative n vector')
    
    % --- make the RTp files - one for each thing that needs animating
    
    write_RTp({bone_struct(:).bone} , {bone_struct(:).(trial_name)} , trial_base_folder)
    
    % --- Make the pos files
    % if _wind is in the var name, it's a path on windows
    local_folder_wind = [direct_anim_wind 'AlongArch\\' trial_name '\\' along_arch_folder '\\' ];
    
    pos_text = write_pos({bone_struct(:).bone},local_folder_wind,trial_name);
    
    filename = fullfile(trial_base_folder, [trial_name '.pos']);
    fid = fopen(filename,'w'); % open the file to write
    fprintf(fid, pos_text);
    fclose(fid);
    
    % --------- make the ini file to animate the helical axes
    
    create_ini(0,0,1,1,[local_folder_wind 'HAdir'],'helicalaxis_P%%d_F%%d.iv',fullfile(trial_base_folder,[trial_name '.ini']))
    
    % -----------make the linked IV files
    link_trial_iv_wind = [local_folder_wind 'rigidiv'];
    
    for b = 1:length(bone_struct)
        ivstring = createInventorHeader();
        
        % mini loop to make the target bone be the same colour as the axis
        for ii = 1:length(b_ind)
            if strcmp(char(bone_line_list{ii,2}),bone_struct(b).bone)
                col = c1(ii,:);%[0 0.1+0.2*ii 0.5];
                break
            else
                col = [0.7 0.7 0.7];
            end
        end
        
        ivstring = [ivstring createInventorLink([link_folder_iv_wind bone_struct(b).bone '.iv'],eye(3,3),zeros(3,1),col,0.5)];
        if exist(rigidiv_folder,'dir')==0
            mkdir(rigidiv_folder);
        end
        fid = fopen(fullfile(rigidiv_folder,[bone_struct(b).bone '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    disp(['Saved file : ' filename '. '])
%     uncomment for phi figures
    figure;
    hold on
    %determine the phis for each relationship
    for ii = 1:length(b_ind)
        i = b_ind(ii);
        phi_tot(j,ii,:) = normalise( cumtrapz(helical_struct(j).Trial(i).phi)' );
        plot(squeeze(phi_tot(j,ii,:)))
    end
    
    title(trial_name)
    legend(leg_str)
    
    % determine the angular velocities and plot
%       figure;
%     hold on
%     %determine the angular phis for each relationship
%     for ii = 1:length(b_ind)
%         i = b_ind(ii);
%         ang_vel_tot(j,ii,:) = normalise(helical_struct(j).Trial(i).phi'*250 );
%         plot(squeeze(ang_vel_tot(j,ii,:)))
%      end
%      
%     title(['Ang vel' trial_name])
%     legend(leg_str)
    
end

%% Wrist Viz 5 - Show the axes of bones in global, not relative to anything else

% look at each bone in global


for j = 1:n_mvt% for every trial
    for i = 1:nbones % for every bone combination
        trial_name = file_list_mvt(j).name;
        
        ham_style ='helicalaxis_P%i_F1.iv';
        folder_comp = sprintf('%s_in_global',bone_struct(i).bone);
        
        trial_base_folder = fullfile(direct_anim,'global',trial_name,folder_comp);
        ham_folder = fullfile(trial_base_folder,'HAdir');
        rigidiv_folder = fullfile(trial_base_folder,'rigidiv');
        
        s_g = bone_struct(i).([file_list_mvt(j).name, '_s']);
        n_g = bone_struct(i).([file_list_mvt(j).name, '_n']);
        c = [101 198 247]/256;
        
        % write all the helical axes of motion
        for k = 1:size(s_g,2) % how many frames of data
            
            ham_filename = sprintf(ham_style,k);
            
            % make the helical axis
            ivstring = createInventorHeader();
            % plot in Xray space
            ivstring = [ivstring createInventorArrow(s_g(:,k),n_g(:,k),150, 3,c,0)];
            ivstring = [ivstring createInventorSphere(s_g(:,k),4,c,0)];
           
            % Write the extra components to a file
            if exist(ham_folder,'dir')==0
                mkdir(ham_folder);
            end
            fid = fopen(fullfile(ham_folder,ham_filename),'w');
            fprintf(fid,ivstring);
            fclose(fid);
        end
        
        % --- make the RTp files - one for each thing that needs animating
        
        write_RTp({bone_struct(:).bone} , {bone_struct(:).(trial_name)} , trial_base_folder)
        
        % --- Make the pos files
        % if _wind is in the var name, it's a path on windows
        local_folder_wind = [direct_anim_wind 'global\\' trial_name '\\' folder_comp '\\'];
        
        pos_text = write_pos({bone_struct(:).bone},local_folder_wind,trial_name);
        
        filename = fullfile(trial_base_folder, [trial_name '.pos']);
        fid = fopen(filename,'w'); % open the file to write
        fprintf(fid, pos_text);
        fclose(fid);
        
        % --------- make the ini file to animate the helical axes
        
        create_ini(0,0,1,1,[local_folder_wind 'HAdir'],'helicalaxis_P%%d_F%%d.iv',fullfile(trial_base_folder,[trial_name '.ini']))
        
        % -----------make the linked IV files
%         link_trial_iv_wind = [local_folder_wind 'rigidiv'];
        
        for b = 1:length(bone_struct)
            ivstring = createInventorHeader();
            ivstring = [ivstring createInventorLink([link_folder_iv_wind bone_struct(b).bone '.iv'],eye(3,3),zeros(3,1),[0.7 0.7 0.7],0.5)];
            if exist(rigidiv_folder,'dir')==0
                mkdir(rigidiv_folder);
            end
            fid = fopen(fullfile(rigidiv_folder,[bone_struct(b).bone '.iv']),'w');
            fprintf(fid,ivstring);
            fclose(fid);
        end
        disp(['Saved file : ' filename '. '])
    end
end

%% Wrist Viz 6 - make along arches but with the average axis 

c1 = [108 182 255;...
    160 66 255;...
    255 128 192]/256;
along_arch_folder_new = [along_arch_folder 'plusAvg'];
trial_dev = [];
phi = [];
L = [];
% showColorMap(c1)
for j = 1:length(helical_struct) % for every trial
    % write to wrist viz
    
    trial_name = file_list_mvt(j).name;
    
    ham_style ='helicalaxis_P%i_F1.iv';
    trial_base_folder = fullfile(direct_anim,'AlongArch',trial_name,along_arch_folder_new);
    ham_folder = fullfile(trial_base_folder,'HAdir');
    rigidiv_folder = fullfile(trial_base_folder,'rigidiv');
    
    % Double check that the HA directory exists
    if exist(ham_folder,'dir')==0
        mkdir(ham_folder);
    end
    if exist(rigidiv_folder,'dir')==0
        mkdir(rigidiv_folder);
    end
    % write all the helical axes of motion for each pair in the same
    % file
    for k = 1:size(helical_struct(j).Trial(1).phi,2) % how many frames of data
        
        ham_filename = sprintf(ham_style,k);
        
        ivstring = createInventorHeader();
        
        for ii = 1:length(b_ind)
            i = b_ind(ii);
            phi_val = helical_struct(j).Trial(i).phi(k);
            
            if rev_flag(ii) == 0
                
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_refx_xr(:,k),helical_struct(j).Trial(i).n_xr(:,k),100*phi_val, 3,c1(ii,:),0)];
                ivstring = [ivstring createInventorSphere(helical_struct(j).Trial(i).s_refx_xr(:,k),4,c1(ii,:),0)];
                ax_temp(:,ii) = helical_struct(j).Trial(i).n_xr(:,k);
                orig_temp(:,ii) = helical_struct(j).Trial(i).s_refx_xr(:,k);
            elseif rev_flag(ii) == 1
                
                ivstring = [ivstring createInventorArrow(helical_struct(j).Trial(i).s_refx_xr(:,k),-helical_struct(j).Trial(i).n_xr(:,k),100*phi_val, 3,c1(ii,:),0)];
                ivstring = [ivstring createInventorSphere(helical_struct(j).Trial(i).s_refx_xr(:,k),4,c1(ii,:),0)];
                ax_temp(:,ii) = -helical_struct(j).Trial(i).n_xr(:,k);
                orig_temp(:,ii) = helical_struct(j).Trial(i).s_refx_xr(:,k);
            end
            
            
            
        phi(ii) = helical_struct(j).Trial(i).phi(k);
        
        L(ii) = helical_struct(j).Trial(i).L(k);
        end
        
        % calculate the average talar axis for each frame
        [avgHam_fr,~,variat_fr] = calcAvgMotionLW([phi', ax_temp, L', orig_temp]);
        %         ax_avg = mean(ax_temp,2);
        if ~isempty(find(avgHam_fr~= 0))
            ax_avg = avgHam_fr(1,2:4);
            orig_avg =  orig_temp;%avgHam_fr(1,6:8);
            %         if phi_val > 0.3
            %             orig_avg = mean(orig_temp,2);
            ivstring = [ivstring createInventorArrow(orig_avg ,ax_avg,80, 3,[0.1 0.1 0.1],0)];
            
            for ii = 1:length(b_ind)
                trial_dev(j,k,ii) = acosd(dot(ax_temp(:,ii),ax_avg)/(norm(ax_temp(:,ii))*norm(ax_avg)));
            end
        end
        
        ivstring = [ivstring createInventorSphere(orig_avg ,4,[0.1 0.1 0.1],0)];
        
        fid = fopen(fullfile(ham_folder,ham_filename),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    
    warning('Think through the negative n vector')
    
    % --- make the RTp files - one for each thing that needs animating
    
    write_RTp({bone_struct(:).bone} , {bone_struct(:).(trial_name)} , trial_base_folder)
    
    % --- Make the pos files
    % if _wind is in the var name, it's a path on windows
    local_folder_wind = [direct_anim_wind 'AlongArch\\' trial_name '\\' along_arch_folder_new '\\' ];
    
    pos_text = write_pos({bone_struct(:).bone},local_folder_wind,trial_name);
    
    filename = fullfile(trial_base_folder, [trial_name '.pos']);
    fid = fopen(filename,'w'); % open the file to write
    fprintf(fid, pos_text);
    fclose(fid);
    
    % --------- make the ini file to animate the helical axes
    
    create_ini(0,0,1,1,[local_folder_wind 'HAdir'],'helicalaxis_P%%d_F%%d.iv',fullfile(trial_base_folder,[trial_name '.ini']))
    
    % -----------make the linked IV files
%     link_trial_iv_wind = [local_folder_wind 'rigidiv'];
    
    for b = 1:length(bone_struct)
        ivstring = createInventorHeader();
        
        % mini loop to make the target bone be the same colour as the axis
        for ii = 1:length(b_ind)
            if strcmp(char(bone_line_list{ii,2}),bone_struct(b).bone)
                col = c1(ii,:);%[0 0.1+0.2*ii 0.5];
                break
            else
                col = [0.7 0.7 0.7];
            end
        end
        
        ivstring = [ivstring createInventorLink([link_folder_iv_wind bone_struct(b).bone '.iv'],eye(3,3),zeros(3,1),col,0.5)];
        if exist(rigidiv_folder,'dir')==0
            mkdir(rigidiv_folder);
        end
        fid = fopen(fullfile(rigidiv_folder,[bone_struct(b).bone '.iv']),'w');
        fprintf(fid,ivstring);
        fclose(fid);
    end
    disp(['Saved file : ' filename '. '])

    
end
%% show the trial deviations from the average helical axis tarsal complex - CSB 2018
figure; hold on; 
for j= 4:6
    for bc = 1:3
        ind = find(trial_dev(j,:,bc)~=0);
          trial_dev_temp = normalise(squeeze(trial_dev(j,ind(1):ind(end),bc))', trial_pct(j));
%           axMeans(j-3,bc) = nanmean(trial_dev(j,ind(1):ind(end),bc));
       
            % add nans to make them the same length
            trial_dev_norm = [nan(start_pct(j),1); trial_dev_temp; nan(100-start_pct(j)-trial_pct(j),1) ];
        axMeans(j-3,bc) = nanmean(trial_dev_norm(15:70));
        plot(trial_dev_norm,'color',c1(bc,:)); 
    end
end
%%
figure;
hold on
for i = 1:3
h(i) = bar(i,mean(axMeans(:,i)),'Facecolor',c1(i,:));
errorbar(i,mean(axMeans(:,i)),std(axMeans(:,i)),'.','color','k')
 set(gca,'xtick',[]);
 set(gca,'xcolor',[1 1 1])
end
makeNicePlotsFunction
%% Helical Contributions - put into the inertial co-ordinate system of the reference bone
close all
clearvars('phi_tot')
% go through each trial and every bone combination in question (b_ind) to compute the
% projection into the reference bone's axes
for j = 4:6%length(helical_struct) % for every trial
  
%     for i = 1:length(helical_struct(j).Trial) % for every bone combination
    for bc = 1:n_bonecomps
         
        ref_numR = findInStruct(bone_struct,'bone',bone_line_list{bc,1});
        b_comp_ind = b_ind(bc);
        
        % project the helical axis for that bone pair into the co-ordinate
        % system of the calcaneus - at this point use the inertial axes
        
%         ref_bone = helical_struct(j).Trial(i).ref_bone;
%         ref_numR = findInStruct(bone_struct,'bone',ref_bone);

        i_hat = bone_struct(ref_numR).CT_inert_axes(1:3,1);
        j_hat = bone_struct(ref_numR).CT_inert_axes(1:3,2);
        k_hat = bone_struct(ref_numR).CT_inert_axes(1:3,3);
        
        for k = 1:size(helical_struct(j).Trial(b_comp_ind).phi,2) % how many frames of data
            
            T_xr_ct = bone_struct(ref_numR).(file_list_mvt(j).name)(:,:,k);
            
            i_hat_xr = transformVectors(T_xr_ct,i_hat);
            j_hat_xr = transformVectors(T_xr_ct,j_hat);
            k_hat_xr = transformVectors(T_xr_ct,k_hat);
            
            if rev_flag(bc) == 1 % reverse the helical axis to compute the correct comp bone to ref bone
                n_xr = -1 * helical_struct(j).Trial(b_comp_ind).n_xr(1:3,k);
            else
                n_xr = helical_struct(j).Trial(b_comp_ind).n_xr(1:3,k);
            end
            
            x_proj(k) = dot(n_xr,i_hat_xr);
            y_proj(k) = dot(n_xr,j_hat_xr);
            z_proj(k) = dot(n_xr,k_hat_xr);
            
        end
        
        helical_struct(j).Trial(b_comp_ind).x_proj = x_proj;
        helical_struct(j).Trial(b_comp_ind).y_proj = y_proj;
        helical_struct(j).Trial(b_comp_ind).z_proj = z_proj;
        
        clearvars('x_proj','y_proj','z_proj')
    end
end



init_nans = start_track_frames-start_frames;
end_nans = end_frames -end_track_frames;
comp_names = {'x_proj','y_proj','z_proj'};

n_components = 3;

cnt_mvt = 1;
count = 1;
legend_str = {};
col_map = colormap('jet');
for j = 4:6% for each trial of interest
    
    
    for bc = 1:n_bonecomps
        
        ref_numR = findInStruct(bone_struct,'bone',bone_line_list{bc,1}); % find the index of the reference bone in each case
        
        % Show the contributions in each axis
        b_comp_ind = b_ind(bc);
        
        if rev_flag(bc) == 1 % reverse the relation comp bone to ref bone
            comp_bone = helical_struct(j).Trial(b_comp_ind).ref_bone;
            ref_bone = helical_struct(j).Trial(b_comp_ind).comp_bone;
        else
            comp_bone = helical_struct(j).Trial(b_comp_ind).comp_bone;
            ref_bone = helical_struct(j).Trial(b_comp_ind).ref_bone;
        end
               
        ref_boneL  = bone_line_list{bc,1};
        comp_boneL = bone_line_list{bc,2};
        
        % put in an error check to make sure I've coded this correctly
        if ~strcmp(comp_bone, comp_boneL) || ~strcmp(ref_bone, ref_boneL)
            error('Not comparing the same bones.')
        end
        fprintf('%s relative to %s - determining contributions for trial %s \n',comp_bone,ref_bone,file_list_mvt(j).name)
        
        phi_temp = helical_struct(j).Trial(b_comp_ind).phi';
        
        for cn = 1:3 % for each component of the reference bone axis
            
            % contribution in a specific axis
            cont = helical_struct(j).Trial(b_comp_ind).(comp_names{cn});
            % compute the angular velocity and total phi
            ang_vel_temp =  cont'.* phi_temp * 250 ;
            phi_tot_temp = cumtrapz(cont'.* phi_temp);
            
            % normalise them to the length of trial they represent (of the
            % tracked data)
            ang_vel_norm_temp = normalise(ang_vel_temp, trial_pct(j));
            phi_tot_norm_temp = normalise(phi_tot_temp, trial_pct(j));
    
            % add nans to make them the same length
            ang_vel{bc}(cnt_mvt,cn,:) = [nan(start_pct(j),1); ang_vel_norm_temp; nan(100-start_pct(j)-trial_pct(j),1) ];
            phi_tot{bc}(cnt_mvt,cn,:) = [nan(start_pct(j),1); phi_tot_norm_temp; nan(100-start_pct(j)-trial_pct(j),1) ];

            figure(bc)
            hold on;
            h = plot(0:100, squeeze(ang_vel{bc}(cnt_mvt,cn,:)));
            set(h,'Color',c(17*cn,:))
            xlabel('% Stance')
            ylabel(sprintf('Angular Vel of %s rel to %s [deg]',comp_bone,ref_bone))

            title(sprintf('Angular Vel of %s rel to %s [deg]',comp_bone,ref_bone))
            %             legend('Inv/Ev','Int ext','Flex')
            %             figure((cnt -1 ) * n_components + j + n_comparisons*n_components);
%             figure(cnt+n_bonecomps*n_components)
            figure(bc+n_bonecomps)
            hold on;
            h = plot(0:100, squeeze(phi_tot{bc}(cnt_mvt,cn,:)));
            set(h,'Color',c(17*cn,:))
            xlabel('% Stance')
            ylabel(sprintf('Phi of %s rel to %s [deg]',comp_bone,ref_bone))
            title(sprintf('Phi of %s rel to %s [deg]',comp_bone,ref_bone))
            %             legend('Inv/Ev','Int ext','Flex')
        end
        
        
        %         get the average axes and the deviations
        
        for ax = 1:2
            ham_var(count,:) = [mean(helical_struct(j).Trial(b_comp_ind).HAM_var{ax}) j bc ax];
            n_vec(count,:) = [helical_struct(j).Trial(b_comp_ind).avg_n_xr{ax} j bc ax];
            %             fprintf('n = %0.3f %0.3f %0.3f \n',n_vec(count,1:3) )
            %             fprintf('theta var = %0.3f \n', ham_var(count,1) )
            if ax == 1
                figure(15)
                hold on
                plotvector3([0 0 0],n_vec(count,:),col_map(count*3,:))
                title('Weight acceptance average axis')
                axis equal
                legend_str = [legend_str sprintf('Trial %i - %s to %s', j, comp_bone,ref_bone)];
                
            else
                figure(16)
                hold on
                plotvector3([0 0 0],-n_vec(count,:),col_map(count*3,:))
                title('pushoff average axis')
                axis equal
            end
            count = count + 1;
            
        end
        
        
    end
    
    cnt_mvt = cnt_mvt+1;
end

figure(15)
legend(legend_str)
figure(16)
legend(legend_str)


%% ang_vel graphs for mt1 vs cal

cnt_mvt = 1:3;
c = colormap('parula');
for bc = 1
    figure;
    hold on
    for cn = 1:3
        meanvel = squeeze(nanmean(ang_vel{bc}(cnt_mvt,cn,:),1))';
        stdvel = nanstd(squeeze(ang_vel{bc}(cnt_mvt,cn,:))); % START HERE
        
        PrettyStdDevGraphs(0:100,meanvel,stdvel,c(17*cn,:),1)
        hold on;
    end
end

%%   ang_vel{bc}(cnt_mvt,cn,:)
diff1 = abs(squeeze(nanmean(ang_vel{3} - ang_vel{1},1)));
diff2 = abs(squeeze(nanmean(ang_vel{2} - ang_vel{1},1)));
for xx = 1:3
figure(xx+10)
plot(diff1(xx,:)')
plot([0 100],[nanmean(diff1(xx,:)') nanmean(diff1(xx,:)')])
hold on
plot(diff2(xx,:)')
plot([0 100],[nanmean(diff2(xx,:)') nanmean(diff2(xx,:)')])
end
%% Make the STJ, NTJ and CTJ graphs of angular velocity
close all
c = colormap('parula');
% cn_names = {'STJ Supination (+)','STJ Dorsiflexion (+)','STJ Internal Rotation(+)'}
cn_names =  {'NTJ Eversion (+)','NTJ Dorsiflexion (+)','NTJ Internal Rotation(+)'};
cn_names =  {'TTJ Eversion (+)','TTJ Dorsiflexion (+)','TTJ Internal Rotation(+)'};
cn_names =  {'Eversion (+)','Dorsiflexion (+)','Internal Rotation(+)'};
cn_names2 = {'STJ','NTJ','CTJ'};
cn_names2 = {'Calcaneus-talus','Navicular-talus','Cuboid-talus'};
figure;
testing = 1:n_bonecomps;

for bc = testing
    ref_boneL  = bone_line_list{bc,1};
    comp_boneL = bone_line_list{bc,2};
%     figure;    
    for cn = 1:3
        figure(cn+10)
%         subplot(3,1,cn)
        hold on
        
        h = PrettyStdDevGraphs( 0:100,squeeze(nanmean(ang_vel{bc}(:,cn,:),1))',squeeze(nanstd(ang_vel{bc}(:,cn,:),1))',c(17*cn,:),1);
%         h = PrettyStdDevGraphs( 0:100, squeeze(nanmean(phi_tot{bc}(:,cn,:),1))',squeeze(nanstd(phi_tot{bc}(:,cn,:),1))',c(17*cn,:),1);
        if bc == 2
            set(h(2),'LineStyle','--')
            perc_diff{1,cn} = (nanmean(ang_vel{bc}(:,cn,:),1)-nanmean(ang_vel{1}(:,cn,:),1))./nanmean(ang_vel{1}(:,cn,:),1)*100;
        elseif bc == 3
            set(h(2),'LineStyle',':')
            
            perc_diff{2,cn} = (nanmean(ang_vel{bc}(:,cn,:),1)-nanmean(ang_vel{1}(:,cn,:),1))./nanmean(ang_vel{1}(:,cn,:),1)*100;
        end
        grid on
        
        xlim([15 70])       
         ylim([-200 200])   
        xlabel('% Stance')
        ylabel(sprintf('%s \n Angular Velocity [deg/s]',cn_names{cn}))
        
        h_line(cn,bc) = h(2);
            
    end

    
end

image_direct = '/home/lauren/Documents/VirtualDesktopTemp/PhDResearch/Abstracts/2018 CSB/Figures/';
for cn = 1:3
    figure(cn+10)
    legend(h_line(cn,testing),cn_names2{testing},'Location','southeast')
    makeNicePlotsFunction
    
%     print(sprintf('%s %s Angular Velocity cub nav.tiff',cn_names{cn}),'-dtiff')
end
% for i = 1:2
%     for cn = 1:3
% figure(cn+4);
% hold on
% plot(squeeze(perc_diff{i,cn}))
% ylim([-300 300])
% end
% end
%% Make the contributions plots (Not this one)


figure;
for j = 1:3
    subplot(3,1,j)
        hold on;
    for i = [4,2,5,6]
        nfr = size(helical_struct(j).Trial(i).phi,2);
        fr_s = 1:nfr;%[4:round(nfr/2)-10,round(nfr/2)+10:nfr-10];
        
        % contribution in a specific axis
        cont = helical_struct(j).Trial(i).y_proj;
        cont_mt1 = helical_struct(j).Trial(4).y_proj;
        
        plot(cont_mt1(fr_s).*helical_struct(j).Trial(4).phi(fr_s),  cont(fr_s).*helical_struct(j).Trial(i).phi(fr_s),'x')

    end
    title(sprintf('Inversion contributions of helical angle (projected in calc axes), %s',file_list_mvt(j).name(6:11)))
    xlabel('Phi of met 1 related to calc [deg]')
    ylabel('Phi of individual bone [deg]') 
    legend({helical_struct(j).Trial([4,2,5,6]).comp_bone})
end
%%
% figure;
% hold on
% for cn = 1:3
%     for j = 1:3
%     perc_vals(j,:,cn) = normalise(abs(squeeze(mean(phi_tot(j,cn,2:end),1)))./total_angles(2:end))';
%     end
%     hold on;
%     PrettyStdDevGraphs( 0:100,mean(perc_vals(:,:,cn),1),std(perc_vals(:,:,cn),1),c(17*cn,:),1)
% %     plot(0:100,perc_vals)
%     disp(max(perc_vals(2:end)))
% % end
% ylim([0 1])
% ylabel('% Contribution to motion')
% xlabel('% Hop')
% legend(cn_names)
% plot(0:100, total_angles)

for cn = 1:3
    disp(cn_names{cn})
    for j = 1:3
        mind(j) = min(ang_vel(j,cn,:));
        maxd(j) = max(ang_vel(j,cn,:));
        
    end
    variat = mean(std(ang_vel(:,cn,:),[],1));
    varivari = std(std(ang_vel(:,cn,:),[],1));
    maxd
    mind
    variat
    varivari
    rangemean = mean(maxd)-mean(mind);
    % perc_er = abs(std(md)/mean(md))+abs(std(maxd)/mean(maxd));
    % abs_er = perc_er*rangemean
    std(mind)+std(maxd);
    
end

%% Look at the variation in the midtarsal axes 
n_xr = zeros(3,2);

    cnt = 1;
for j = 4:6  % for every trial
    
    for k = 1:size(helical_struct(j).Trial(1).phi,2) % how many frames of data
        bc_cnt = 1;
        
        for bc = 2:3 % nav and cub in tarsal
              b_comp_ind = b_ind(bc);
              
            if rev_flag(bc) == 1 % reverse the helical axis to compute the correct comp bone to ref bone
                n_xr(:,bc_cnt) = -1 * helical_struct(j).Trial(b_comp_ind).n_xr(1:3,k);
            else
                n_xr(:,bc_cnt) = helical_struct(j).Trial(b_comp_ind).n_xr(1:3,k);
            end
           
         phi_val(cnt,bc_cnt,k) =  helical_struct(j).Trial(b_comp_ind).phi(k);
          bc_cnt = bc_cnt+1;
        end
        
        ang_dev(cnt,k) = acosd(dot(n_xr(:,1),n_xr(:,2))/(norm(n_xr(:,1))*norm(n_xr(:,2))));
       
    end
        cnt = cnt +1;
   
            
end

for i = 1:cnt-1
    ang_dev_norm(i,:) = normalise(ang_dev(i,:)')';
end

figure
subplot(3,1,1)
plot(ang_dev_norm(:,:)')

subplot(3,1,2)
plot(squeeze(phi_val(:,1,:))')
hold on;
plot(squeeze(phi_val(:,2,:))')

subplot(3,1,3)
% plot(
% PrettyStdDevGraphs(0:52