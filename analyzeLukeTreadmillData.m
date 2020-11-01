close all;
clear;
clc
load('C:\Users\Lauren\Documents\School\PhD\Research\AllStudies\1B_WindlassArchSpringRunning\LukeTreadmill\viz3dexport2.mat')

%% Filter the force and determine gait events (i.e. stance time)

ntrials = 13;
force_threshold = 50;
resamp_fact = 20; % force is sampled 20 times faster than mocap
forceFR = 4000;
fc_FP = 35;
fc_moc = 10;
for i = 1:ntrials
    
    force_total = FP1_F{i} + FP2_F{i};
    
    force_filt{i} = LowPassButterworth(force_total,2,fc_FP,forceFR);

    FP1F_filt{i} = LowPassButterworth(FP1_F{i},2,fc_FP,forceFR);
    FP2F_filt{i} =  LowPassButterworth(FP2_F{i},2,fc_FP,forceFR);
    
    FP1M_filt{i} = LowPassButterworth(FP1_M{i},2,fc_FP,forceFR);
    FP2M_filt{i} =  LowPassButterworth(FP2_M{i},2,fc_FP,forceFR);
    
    cp1 = normaliseNaN(FP1_CP{i},1,length(FP1_CP{i}));
    cp2 = normaliseNaN(FP2_CP{i},1,length(FP2_CP{i}));
    
    FP1CP_filt{i} = LowPassButterworth(cp1,2,fc_FP,forceFR);
    FP2CP_filt{i} =  LowPassButterworth(cp2,2,fc_FP,forceFR);
    
    trapz(1/4000,FP1F_filt{i}(2,:)) 
    
    
    ind_0 = force_filt{i}(:,3) > force_threshold;
    id_start = logical([0 0 1 1]);
    id_end = logical([ 1 1 0 0]);
   
    ind_start = strfind(ind_0',id_start)+2; % find the start points
    ind_end = strfind(ind_0',id_end)+2;
    
    
    
    if ind_start(1)>ind_end(1) % if the start pulse is after the first end, get rid of the first end
        ind_end(1) = [];
    end
    if length(ind_start) > length(ind_end) % if there are not the same number of start and end points
        ind_start(end) = [];
    end
    ind_all = sort([ind_start, ind_end]);
    
    ind_mocap = round(ind_all/resamp_fact);
   
    %contact times
    
    ct{i} = [ind_all(1:2:end)',ind_all(2:2:end)'];
    
    n_contacts = size(ct{i},1);
       delete_inds = [];
    for c = 1:n_contacts
       
        nfrC = length(ct{i}(c,1):ct{i}(c,2)); % number of frames in this "event"
        
        
        if nfrC < 200 % if it's less than 200 frames
            delete_inds = [delete_inds c];
            
        end
        
    end
    tt = [];
    ct{i}(delete_inds,:) = [];
    delete_inds = [];
    n_contacts = size(ct{i},1);
    for c = 1:n_contacts
        
        
        nfrC = length(ct{i}(c,1):ct{i}(c,2)); % number of frames in this "event"
        
        nhalf = round(nfrC/2);
         tt(c) = trapz(force_filt{i}(ct{i}(c,1):ct{i}(c,2)-nhalf,2)); % if the net Y impulse is positive, it's the right foot fall
       
    end
   delete_inds =  find(tt < 0.6*max(tt)) ;% get rid of the foot fall b/c it's left

    ct{i}(delete_inds,:) = [];

    
    % update the # contacts
    n_contacts = size(ct{i},1);
    figure;hold on;
    for c = 1:n_contacts
        nfrC = length(ct{i}(c,1):ct{i}(c,2)); 
        plot(force_filt{i}(ct{i}(c,1):ct{i}(c,2),1:3));
    
    end
title(sprintf('subject %i',i))
drawnow


end
% i = 5;
% plot(force_total); 
% figure;hold on; 
% plot(force_filt{i}); plot(find(force_filt{i}(:,3)>force_threshold),force_filt{i}(force_filt{i}(:,3)>force_threshold,1:3),'.')
%%  determine the MTP angle and the pseudo plantar fascia length
m_fields = {'CA','CA2','CPT','CST','MH1','MH2','MH5','MB1','MB2','MB5','IP1','IP2','IP4','MM','LM','MKNEE','LKNEE','S1','S2','S3','S4','NT'};
i = 2;
figure;
hold on;

plot(FP1_CP{i}(ct{i}(1,1):ct{i}(1,2),1),FP1_CP{i}(ct{i}(1,1):ct{i}(1,2),2),'r.-')
plot(FP2_CP{i}(ct{i}(1,1):ct{i}(1,2),1),FP2_CP{i}(ct{i}(1,1):ct{i}(1,2),2),'.-')

figure;
hold on;

plot(FP1F_filt{i}(ct{i}(1,1):ct{i}(1,2),:),'r.-')
plot(FP2F_filt{i}(ct{i}(1,1):ct{i}(1,2),:),'.-')
plot(norm3d(FP1F_filt{i}(ct{i}(1,1):ct{i}(1,2),:)),'r.-')
plot(norm3d(FP2F_filt{i}(ct{i}(1,1):ct{i}(1,2),:)),'b.-')
figure
hold on;

plot(FP1_CP{i}(ct{i}(1,1):ct{i}(1,2),:),'r.-')
plot(FP2_CP{i}(ct{i}(1,1):ct{i}(1,2),:),'.-')
close all;
n_steps = 1:5; % which steps to take
for i = [1:5,7:ntrials] % number of trials
%     figure;
% plot(1:size(CA2{i},1),CA2{i}(:,3)*1000)
% hold on;
% plot(1:size(CA2{i},1),force_filt{i}(1:20:end,3))

    
    s_ind = 1;
  
        m.CA = CA{i};
        m.CA2 = CA2{i};
        m.CPT = CPT{i};
        m.CST = CST{i};
        m.MH1 = MH1{i};
        m.MH2 = MH2{i};
        m.MH5 = MH5{i};
        m.MB1 = MB1{i};
        m.MB2 = MB2{i};
        m.MB5 = MB5{i};
        m.IP1 = IP1{i};
        m.IP2 = IP2{i};
        m.IP4 = IP4{i};
        m.MM = MM{i};
        m.LM = LM{i};
        m.MKNEE = MKNEE{i};
        m.LKNEE = LKNEE{i};
        m.S1 = S1{i};
        m.S2 = S2{i};
        m.S3 = S3{i};
        m.S4 = S4{i};
        m.NT = NT{i};

      figure; hold on;
     
     nfr = size(CA2{i},1);
     x = 1/200:1/200:(nfr)/200;
%      xq = 1/1000:1/1000:(nfr)/200;
%      xq = 1/200:1/200:(nfr)/200;
     xq = 1/forceFR:1/forceFR:(nfr)/200;
     
     for ff = 1:length(m_fields)
         % filter the data
         m_filt.(m_fields{ff}) = LowPassButterworth(m.(m_fields{ff})(:,1:3),2,fc_moc,200)';
         % interpolate over nans
         temp = normaliseNaN(m_filt.(m_fields{ff}),2,nfr);
         % interpolate over same number of force points
         m_int.(m_fields{ff}) = normaliseNaN(temp,2,nfr*20);
     end
    
  
     
     for s = 1:size(ct{i},1) % number of steps
         
           for ff = 1:length(m_fields)
               m_crop.(m_fields{ff}) = m_int.(m_fields{ff})(1:3,ct{i}(s,1):ct{i}(s,2));
           end

            marker_data.CA_ = m_crop.CA2;
            marker_data.CST = m_crop.CST;
            marker_data.CPT = m_crop.CPT;
            marker_data.NT_ = m_crop.NT;
            marker_data.MH1 =  m_crop.MH1;
            marker_data.MH2 =  m_crop.MH2;
            marker_data.MH5 =  m_crop.MH5;
            marker_data.MB1 = m_crop.MB1;
            marker_data.MB2 = m_crop.MB2;
            marker_data.MB5 = m_crop.MB5;
            marker_data.IP1 = m_crop.IP1;
            marker_data.IP2 = m_crop.IP2;
            marker_data.IP4 = m_crop.IP4;
            marker_data.MM_ = m_crop.MM;
            marker_data.LM_ = m_crop.LM;
            marker_data.ME_ = m_crop.MKNEE;
            marker_data.LE_ = m_crop.LKNEE;
            
            model = createFootModel;
            cal_com = [];
            for f = 1:length(ct{i}(s,1):ct{i}(s,2))

                model.pose.shank(:,:,f) = LeastSquaresPose([m_crop.S1(:,f),m_crop.S2(:,f),m_crop.S3(:,f),m_crop.S4(:,f)]);
                model.pose.cal(:,:,f) = LeastSquaresPose([m_crop.CA(:,f),m_crop.CST(:,f),m_crop.CPT(:,f),m_crop.CA2(:,f)]);

                for di = 1:3
                    cal_com(di,f) = mean([m_crop.CA(di,f),...
                        m_crop.CST(di,f) ,...
                        m_crop.CPT(di,f) ]);
                end
            end
%             figure; plot( squeeze(model.pose.shank(2,3,:)))
%              model.pose.shank
%             
%             subplot(3,1,1)
%              hold on;
%             plot(model.F2Ps)
%             title('toe angle')
%             subplot(3,1,2)
%              hold on;
%             plot(model.elongation)
%             title('elongation')
%               subplot(3,1,3)
%              hold on;
%             plot(model.F2Ps,model.elongation)
%             title('snail')


            com_foot = ((m_crop.MH5+m_crop.MH1)/2 - m_crop.CA)/2 +m_crop.CA;
            ankle = (m_crop.MM+m_crop.LM)/2;
            knee = (m_crop.MKNEE+m_crop.LKNEE)/2;
            com_shank = knee+0.3705*(ankle-knee); % clauser
%             
%             
%             figure; hold on;
%             for ff = 1%:20:nfrC
%                 foot_mark = [m_crop.MH5(:,ff), m_crop.MH1(:,ff), m_crop.CA(:,ff)]; 
%                 plot3quick_scatter(foot_mark,'r');
%                 plot3quick_scatter(com_foot(:,ff),'k')
%             end
%                 plotvector3(copC(:,ff),f_grfC(:,ff)/100)
%                 axis equal
%                 drawnow
%                 pause(0.05)
%             end
            
             f1 = FP1F_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)';
            f2 = FP2F_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)';
            
            cop1 = FP1CP_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)' ;%normaliseNaN(FP1CP_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)',2,nfrC);
            cop2 = FP2CP_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)' ;%normaliseNaN(FP2CP_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)',2,nfrC);
            
            m1 = FP1M_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)';
            m2 = FP2M_filt{i}(ct{i}(s,1):ct{i}(s,2),1:3)';
            
            nfrC = length(ct{i}(s,1):ct{i}(s,2));
            % find when fp1 has contact
            ind_1 = FP1F_filt{i}(ct{i}(s,1):ct{i}(s,2),3) > force_threshold;
            id_start = logical([0 0 1 1]);
            id_end = logical([ 1 1 0 0]);
            
            frms_start1 = strfind(ind_1',id_start)+2; % find the start points
            if isempty(frms_start1)== 1
                continue
            end
            frms_end1 = strfind(ind_1',id_end)+2;
            if isempty(frms_end1)== 1 
                frms_end1 = ct{i}(s,2)-ct{i}(s,1)+1;
            end
            % find when fp2 has contact
            ind_2 = FP2F_filt{i}(ct{i}(s,1):ct{i}(s,2),3) > force_threshold;
          
            id_start = logical([0 0 1 1]);
            id_end = logical([ 1 1 0 0]);
            
            frms_start2 = strfind(ind_2',id_start)+2; % find the start points
            frms_end2 = strfind(ind_2',id_end)+2;
            
            ind_small = frms_start2 < 0.7*nfrC;
           frms_start2(ind_small) = []; % get rid of anything in the beginning of the trial (can't be on second plate to start)
             ind_small = frms_end2 < 0.7*nfrC;
            frms_end2(ind_small) = []; % get rid of anything in the end of the trial
            
            
%             ind_interaction = ind_1 & ind_2;
           
%             ind_interaction(1:round(0.8*nfrC)) = 0;
%             frm_inter = find(ind_interaction == 1);
        ind_inter = logical(zeros(1,nfrC));
        ind_inter(frms_start2:frms_end1) = ones(1,frms_end1-frms_start2+1);
%             
            frm_inter = frms_start2:frms_end1;
                
%            USING V3D THRESHOLD
            dist2COP{1} =  norm3d(com_foot - cop1);
            dist2COP{2} =  norm3d(com_foot - cop2);
            indfrms1 = dist2COP{1} < 0.2; % threshold from v3d - must be this close to the segment
            indfrms2 = dist2COP{2} < 0.2;
            thresh_beg = 0.7*nfrC;
            indbeg = logical(zeros(1,nfrC));
            indbeg(1:round(thresh_beg)) = ones(1,round(thresh_beg));
            
            ind_interp = ~indbeg & ind_inter;% which indices to interpolate over
            

            if norm(diff([cop2(1:3,frms_start2),cop1(1:3,frms_start2)]'))>0.05; % COP issue
                continue
            end
            
            f_grf{i,s} = nan(3,nfrC);
            f_grf{i,s}(1:3,frms_start2:frms_end2) = f2(1:3,frms_start2:frms_end2);
            f_grf{i,s}(1:3,frms_start1:frms_end1) = f1(1:3,frms_start1:frms_end1);
            f_grf{i,s}(1:3,ind_interp) = nan(3,length(find(ind_interp==1)));
            f_grf{i,s} = normaliseNaN(f_grf{i,s},2,nfrC);
            
            
            cop_grf{i,s} = nan(3,nfrC);
            cop_grf{i,s}(1:3,frms_start2:frms_end2) = cop2(1:3,frms_start2:frms_end2)-cop2(1:3,frms_start2)+cop1(1:3,frms_start2);
            cop_grf{i,s}(1:3,frms_start1:frms_end1) = cop1(1:3,frms_start1:frms_end1);
            cop_grf{i,s}(1:3,ind_interp) = nan(3,length(find(ind_interp==1)));
            cop_grf{i,s} = normaliseNaN(cop_grf{i,s},2,nfrC);
             
              
            mom_grf{i,s} = nan(3,nfrC);
            mom_grf{i,s}(1:3,frms_start2:frms_end2) = m2(1:3,frms_start2:frms_end2);
            mom_grf{i,s}(1:3,frms_start1:frms_end1) = m1(1:3,frms_start1:frms_end1);
            mom_grf{i,s}(1:3,ind_interp) = nan(3,length(find(ind_interp==1)));
            mom_grf{i,s} = normaliseNaN(mom_grf{i,s},2,nfrC);
            
            
            

            
%            
% %             --------------------------------
           
%             
%             figure(11); plot(f_grf{i,c}') ; hold on;
%             figure(12); plot(cop_grf{i,c}'); hold on;
%             figure; plot3(cop1(1,ind_start1:ind_end1),cop1(2,ind_start1:ind_end1),cop1(3,ind_start1:ind_end1),'.'); hold on;
%             plot3(cop2(1,ind_start2:ind_end2),cop2(2,ind_start2:ind_end2),cop2(3,ind_start2:ind_end2),'.');
%             plotvector3(cop_grf{i,c}(1:3,1:5:end),f_grf{i,c}(1:3,1:5:end))
%             axis equal
  
            f_grfC =  f_grf{i,s};
            m_grfC = mom_grf{i,s};
            copC = cop_grf{i,s};
%             
%             
%             figure(20+3*i);hold on;
%             plot(f_grfC')
%             
%             figure(21+3*i);hold on;
%             plot(m_grfC')
%             
%             figure(22+3*i);hold on;
%             plot(copC')
%             
            
            
            % determine w of shank and of calc
            w_cal =  calculateRotMatAngularVelocity( model.pose.cal(1:3,1:3,:),forceFR,'rad');
            w_shank = calculateRotMatAngularVelocity( model.pose.shank(1:3,1:3,:),forceFR,'rad');

            r_com_cop = copC - cal_com;
            
            r_com_cop_sh = copC - com_shank ;
             
            
            v_com = diffn(cal_com,2)*forceFR;
            v_com_sh = diffn(com_shank,2)*forceFR;
            v_com_cop = diffn(r_com_cop,2)*forceFR; % this is the global from differentiating the vector
            v_com_copsh = diffn(r_com_cop_sh,2)*forceFR; 
            
            vp_cont = v_com + v_com_cop;
            vp_cont_sh = v_com_sh + v_com_copsh;
            
            vp_UD = v_com + cross(w_cal,r_com_cop);
            vp_UDsh = v_com_sh + cross(w_shank,r_com_cop_sh);
            
            cmap = colormap('parula');
            % calculate as UD with calcaneus
%             
 
P_UD  = dot(f_grfC,vp_UD) + dot(m_grfC,w_cal);
           
%             P_UD  = dot(f_grfC,v_com)  + dot(f_grfC,cross(w_cal,r_com_cop))+ dot(m_grfC,w_cal);
            figure(10+4*i); hold on;
            plot(P_UD,'LineWidth',2,'Color','k')
            plot(dot(f_grfC,v_com) ,'--','Color',cmap(38,:))
            plot(dot(f_grfC,cross(w_cal,r_com_cop)),'--','Color',cmap(52,:))
            plot(dot(m_grfC,w_cal),'--','Color',cmap(22,:))
            legend('Sum','Force- rigid','Force- deform','Moment')
           title('Unified Deformable - calc')
           xlabel('Frame')
           ylabel('Power (W)')
           
%             calculate UD shank
%             P_UDsh = dot(f_grfC,vp_UDsh) + dot(m_grfC,w_shank);
            P_UDsh = dot(f_grfC,vp_UDsh)  +  dot(m_grfC,w_shank);
              figure(11+4*i); hold on;
            plot(P_UDsh,'LineWidth',2,'Color','k')
            plot(dot(f_grfC,v_com_sh),'--','Color',cmap(38,:))
            plot(dot(f_grfC,cross(w_shank,r_com_cop_sh)),'--','Color',cmap(52,:))
            plot(dot(m_grfC,w_shank),'--','Color',cmap(22,:))
            legend('Sum','Force- rigid','Force- deform','Moment')
           title('Shank Unified Deformable')
           xlabel('Frame')
           ylabel('Power (W)')
           
           figure(100);hold on;
           for jj = 1:3
           plot(w_cal(jj,:),'color',cmap(10+10*jj,:))
           end
           drawnow
           
%            calculate as LW with shank
%            
%             P_LWsh = dot(f_grfC,v_com_sh) + dot(f_grfC,v_com_copsh) + dot(m_grfC,w_shank);
%             figure(12+4*i); hold on; 
%             plot(P_LWsh,'LineWidth',2,'Color','k')
%             plot(dot(f_grfC,v_com_sh) ,'--','Color',cmap(38,:))
%             plot(dot(f_grfC,v_com_copsh),'--','Color',cmap(52,:))
%             plot(dot(m_grfC,w_shank),'--','Color',cmap(22,:))
%             legend('Sum','Force- rigid','Force- deform','Moment')
%   
%            title('LW with full COP vel Shank Deformable')
%             xlabel('Frame')
%            ylabel('Power (W)')
%            
%             % calculate LW with the calcaneus
%             
%             P_LW = dot(f_grfC,v_com) +dot(f_grfC,v_com_cop) + dot(m_grfC,w_cal);
%             figure(13+4*i); hold on;
%             plot(P_LW,'LineWidth',2,'Color','k')
%             plot(dot(f_grfC,v_com) ,'--','Color',cmap(38,:))
%             plot(dot(f_grfC,v_com_cop),'--','Color',cmap(52,:))
%             plot(dot(m_grfC,w_cal),'--','Color',cmap(22,:))
%             legend('Sum','Force- rigid','Force- deform','Moment')
%             title('LW- with full COP vel Calc/floor power')
%             xlabel('Frame')
%            ylabel('Power (W)')
%             figure(300);
%             h1 = plot(vp_UD','k');
%             hold on;
%             h2 = plot(vp_cont','b');
%             legend([h1,h2],{'UD velocity','contact pt velocity'})
             s_ind = s_ind+1;
    end
   
    
    
end