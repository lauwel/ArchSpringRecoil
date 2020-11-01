% Processing Toni's data


% close all
clear
clc
trial_name = 'T0024_SOL001_nrun_rfs_barefoot';
subj_dir = 'E:\SOL001_VISIT2\';
calib_dir = [subj_dir 'Calibration\'];
calib_tr_dir = [calib_dir 'Set 1\'];
subj_name = 'SOL001B';

% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\'; % if saving animation to
% server
ivdir = [subj_dir 'Models\IV\'];
trialdir = [subj_dir trial_name '\'];

frRate = 250;
fc = 15;
wc = [20 60];
pretrigFrms = frRate * 1; % 1 second of pre trig
process_mocap = 0;
nfr = 400;
load([calib_dir,'Sync/syncFrames.mat']); % this loads the syncFrs variable so we know where mocap is relative to x-ray
coreg = csvread([calib_tr_dir,'end_pylon\end_pylon_COREG.csv']);
%% process the motion capture file - must have exported a matlab file

if process_mocap == 1
    mocap_file = [trialdir,'Mocap\', trial_name,'.mat'];
    [mocap_data] = MocapDataTransform({mocap_file}, 'filter','on',...
        'forceCutOff',50,...
        'mocapCutOff',50,...
        'lab','SOL',...
        'visualiseForces','on',...
        'resample','force');
    mocap_file_processed = [trialdir,'Mocap\', trial_name,'_processed.mat'];
    save(mocap_file_processed,'mocap_data')
    
else
    load([trialdir,'Mocap\', trial_name,'_processed.mat'])
end
fr1 = pretrigFrms + syncFrs.(trial_name);
fr2 = fr1 + nfr;
for i = 1:2 % each force plate
    force{i} = mocap_data.force_data(i).globForce(1:3,fr1:fr2);
    cop{i} = mocap_data.force_data(i).globCOP(1:3,fr1:fr2);
    moment{i} = mocap_data.force_data(i).globMoment(1:3,fr1:fr2);
    freemom{i} = mocap_data.force_data(i).globFreeMoment(1:3,fr1:fr2);
end

force{1}(:,1:125) = repmat([0;0;0.00001],1,125);
force{2}(:,155:end) = repmat([0;0;0.00001],1,401-154);
freemom{1}(:,1:125) = 0;
freemom{2}(:,155:end) = 0;

%% Load the transforms, transform them to mocap space, calculate angular velocity

bonesCell = {'tib','tal','cal','mt1'};
fr{1} = 121:186; % tracked frames for each bone
fr{2} = 121:186;
fr{3} = 121:186;
fr{4} = 121:186;
clearvars T w
coreg_m = coreg;
coreg_m(1:3,4) = coreg(1:3,4)/1000;
for b = 1:4
    ivfile = [ivdir subj_name '_' bonesCell{b} '_aligned.iv'];
    % import the bone points and connections
    [pts.(bonesCell{b}),cns.(bonesCell{b})] = read_vrml_fast(ivfile);
    cns.(bonesCell{b}) = cns.(bonesCell{b})(:,1:3) + 1;
    pts.(bonesCell{b}) = pts.(bonesCell{b})/1000; % convert to metres
    % calculate the COM
    [com.(bonesCell{b}),~,~,~,axes.(bonesCell{b}),~,~,~,~,~] = mass_properties(pts.(bonesCell{b}),cns.(bonesCell{b}));
    
    autoscoper_dat = dlmread([trialdir,'Autoscoper\Autoscoped\' trial_name '_' bonesCell{b} '_filt_interp.tra']);
    Ttemp = convertRotation(autoscoper_dat,'autoscoper','4x4xn');
    Ttemp(1:3,4,:) = Ttemp(1:3,4,:)/1000;
    
    for f = 1:size(Ttemp,3) %transform every frame into global mocap using the Coregistration
        T.(bonesCell{b})(:,:,f) = (coreg_m) * Ttemp(:,:,f);
        for i = 1:2 % each force plate
            cop_foot{i}(:,f) = closestPointonPlanealongVector(cop{i}(:,f)',[0 0 1],[0 0 0],force{i}(:,f)')'; % because there's a dragon plate on top, but the plane is aligned with the global cosys
        end
        p_com.(bonesCell{b})(:,f) = transformPoints(T.(bonesCell{b})(:,:,f),com.(bonesCell{b}),0);
    end
    w.(bonesCell{b}) = zeros(3,nfr+1);
    % determine the angular velocity
    w.(bonesCell{b})(1:3,fr{b}) = calculateRotMatAngularVelocity(T.(bonesCell{b})(1:3,1:3,fr{b}),250,'rad');
        w.(bonesCell{b})(1:3,fr{b}) = adaptiveLowPassButterworth(w.(bonesCell{b})(1:3,fr{b}),wc,250);
    %     wh.(bonesCell{b}) = calculateHelicalAxisAngularVelocity(repmat(eye(4,4),1,1,401),T.(bonesCell{b}),250,'rad');
    vcom.(bonesCell{b}) = zeros(3,nfr+1);
    % determine the linear velocity of the centre of mass of the bone
    vcom.(bonesCell{b})(1:3,fr{b}) = calculateVelocity(p_com.(bonesCell{b})(:,fr{b}),250);
            vcom.(bonesCell{b})(1:3,fr{b}) = adaptiveLowPassButterworth(vcom.(bonesCell{b})(1:3,fr{b}),wc,250);
end
%
%
% figure;
% cmap = colormap('parula');
% for f =150 %120:185
%     b = 3;
%     plot3quick(mocap_data.force_data(1).ForcePlateLocation([1:4,1],1:3)');hold on;
%     plot3quick(mocap_data.force_data(2).ForcePlateLocation([1:4,1],1:3)');
%     plot3quick_scatter(p_com.(bonesCell{b})(:,f));
%     patch('faces',cns.cal(1:50:end,1:3) ,'vertices',transformPoints(T.cal(:,:,f),pts.cal,0),'facealpha',0.5,'edgecolor',cmap(12,:),'facecolor',cmap(12,:))
%     patch('faces',cns.mt1(1:50:end,1:3) ,'vertices',transformPoints(T.mt1(:,:,f),pts.mt1,0),'facealpha',0.5,'edgecolor',cmap(32,:),'facecolor',cmap(32,:))
%
%     plotvector3(cop_foot{1}(:,f),force{1}(:,f)/1500);
%     plotvector3(cop_foot{2}(:,f),force{2}(:,f)/1500);
%
%
%
%     view(-114,7)
%     axis equal
%     hold off
%     pause(0.5)
% end



%% Calculate the UD power relative to calc, tal and tib


cmap = colormap('parula');
cop_new = zeros(3,nfr+1);

new_centre = [-0.2 0 0]'; % for recalculating the COP
for b = 1:4
    for f = 1:nfr
        
        %    --    individual (two) force plates
        
        for i = 1:2 % which force plate
            
            r_com_cop =  cop_foot{i}(:,f) - p_com.(bonesCell{b})(:,f); % second force plate is only one contacting heel
            v_dist = vcom.(bonesCell{b})(:,f) + cross(w.(bonesCell{b})(1:3,f), r_com_cop);
            PUD(i).(bonesCell{b})(:,f) = dot(force{i}(:,f),v_dist) + dot(freemom{i}(:,f),w.(bonesCell{b})(:,f));
        end
        
        %    --    combined force plate
        
        
        % get the total force
        FT(:,f) = force{1}(:,f) + force{2}(:,f);
        momT(:,f) = freemom{1}(:,f) + freemom{2}(:,f);
        % calculate the new COP by summing moments around new origin and
        % calculating the new COP
        r1_com(:,f) = cop_foot{1}(:,f) - new_centre;%p_com.(bonesCell{b})(:,f);
        r2_com(:,f) = cop_foot{2}(:,f) - new_centre;%p_com.(bonesCell{b})(:,f);
        
        M_sum(:,f) = cross(r1_com(:,f),force{1}(:,f)) + cross(r2_com(:,f),force{2}(:,f));
        
        
        cop_new(:,f) = [-M_sum(2,f)/FT(3,f) ;...
            M_sum(1,f)/FT(3,f)  ;...
            0];
        
        cop_new(:,f) = cop_new(:,f) + new_centre;
        
        r_com_cop =  cop_new(:,f) - p_com.(bonesCell{b})(:,f); % second force plate is only one contacting heel
        v_dist = vcom.(bonesCell{b})(:,f) + cross(w.(bonesCell{b})(1:3,f), r_com_cop);
        PUDsum.(bonesCell{b})(:,f) = dot(FT(:,f),v_dist) + dot(momT(:,f),w.(bonesCell{b})(:,f));
        
    end
    
%     first 20% of stance

    fr_20 = false(1,400);%120 + [1:14]; % 1/5*66 frames
    fr_20(121:134) = true(1,14);
    fr_2050 = false(1,400); 
    fr_2050(135:153) = true(1,19);
    fr_end = false(1,400);
    fr_end(154:186) = true(1,33);
    
    Wpos.(bonesCell{b}) = false(2,400);
    Wneg.(bonesCell{b}) = false(2,400);
    for i = 1:2
         
        pos_ind = PUD(i).(bonesCell{b}) > 0;
        neg_ind = PUD(i).(bonesCell{b}) < 0;
        frp.(bonesCell{b}){1,i}  = length(find(fr_20 & pos_ind));
        frp.(bonesCell{b}){2,i}  = length(find(fr_2050 & pos_ind));
        frp.(bonesCell{b}){3,i}  = length(find(fr_end & pos_ind));
        
        frn.(bonesCell{b}){1,i}  = length(find(fr_20 & neg_ind));
        frn.(bonesCell{b}){2,i}  = length(find(fr_2050 & neg_ind));
        frn.(bonesCell{b}){3,i}  = length(find(fr_end & neg_ind));
%         
%         Wpos.(bonesCell{b})(i,pos_ind) = cumtrapz(PUDsum.(bonesCell{b}));
%         Wneg.(bonesCell{b})(i,neg_ind)  = cumtrapz(PUDsum(i).(bonesCell{b})(neg_ind));
        
        WUD20p(i).(bonesCell{b}) = trapz(PUD(i).(bonesCell{b})(fr_20 & pos_ind));
        WUD50p(i).(bonesCell{b}) = trapz(PUD(i).(bonesCell{b})(fr_2050 & pos_ind));
        WUD100p(i).(bonesCell{b}) = trapz(PUD(i).(bonesCell{b})(fr_end & pos_ind));
        
        WUD20n(i).(bonesCell{b}) = trapz(PUD(i).(bonesCell{b})(fr_20 & neg_ind));
        WUD50n(i).(bonesCell{b}) = trapz(PUD(i).(bonesCell{b})(fr_2050 & neg_ind));
        WUD100n(i).(bonesCell{b}) = trapz(PUD(i).(bonesCell{b})(fr_end & neg_ind));
    end
    
    
    
    figure;
    subplot(3,1,1)
    plot(PUD(1).(bonesCell{b})); hold on;
    plot(PUD(2).(bonesCell{b}));
    plot(PUDsum.(bonesCell{b}));
    plot(PUD(1).(bonesCell{b})+PUD(2).(bonesCell{b}),':')
    
    title(sprintf('UD Power for %s',bonesCell{b}))
    legend('UD Power force anterior','UD Power force plate posterior','UD Power from combined plate', 'Sum of force plate powers')
    xlabel('Frame')
    ylabel('Power [W]')
    xlim([100,200])
    
    subplot(3,1,2)
    plot(w.(bonesCell{b})')
    ylabel('Angular velocity [rad/s]')
    legend('x','y','z')
    xlim([100,200])
    
    subplot(3,1,3)
    plot(vcom.(bonesCell{b})')
    ylabel('Velocity [m/s]')
    legend('x','y','z')
    xlim([100,200])
end


P.MT = PUDsum.tib(:,fr{1})-PUDsum.tal(:,fr{1});
P.ST = PUDsum.tal(:,fr{1})-PUD(2).cal(:,fr{1});
P.TC = PUDsum.tal(:,fr{1})-PUD(1).mt1(:,fr{1});

jointNames = fields(P);

for j = 1:length(jointNames)
    pos_ind = P.(jointNames{j}) > 0;
    neg_ind = P.(jointNames{j}) < 0;
    W.(jointNames{j})(1) = trapz([1:length(find(pos_ind))]/250,P.(jointNames{j})(pos_ind));
    W.(jointNames{j})(2) = trapz([1:length(find(neg_ind))]/250,P.(jointNames{j})(neg_ind));
end



figure;
hold on
h = bar([[W.TC;W.ST;W.MT]';sum([W.TC;W.ST;W.MT],2)']);
set(h(1),'facecolor',cmap(30,:))
set(h(2),'facecolor',cmap(5,:))
set(h(3),'facecolor',cmap(20,:))
legend(h,{'Talocrural','Subtalar','Midtarsal'})
set(gca,'xticklabel',{'+ work','','- work','','net work',''})
ylabel('Work [J]')
makeNicePlotsFunction


stance = linspace(0,100,fr{1}(end)-fr{1}(1)+1);

 % plot all the sums
figure; 
plot(stance,PUDsum.tib(:,fr{1}),'color',cmap(15,:),'linestyle','-');hold on;
plot(stance,PUDsum.tal(:,fr{1}),'color',cmap(30,:),'linestyle','-');hold on;
plot(stance,PUDsum.cal(:,fr{1}),'color',cmap(40,:),'linestyle','-');hold on;
plot(stance,PUDsum.mt1(:,fr{1}),'color',cmap(1,:),'linestyle','-');hold on;
xlabel('% Stance')
ylabel('Power [W]')
% xlim([100,200])
legend('Tibia UD','Talus UD','Calc UD','MT1 UD')
grid on
makeNicePlotsFunction  

% plot the tib and tal sums with mt1 and cal on anterior/posterior plates
figure;
plot(stance,PUDsum.tib(:,fr{1}),'color',cmap(15,:),'linestyle','-');hold on;
plot(stance,PUDsum.tal(:,fr{1}),'color',cmap(30,:),'linestyle','-');
plot(stance,PUD(1).mt1(:,fr{1}),'color',cmap(1,:),'linestyle','-.'); hold on;
plot(stance,PUD(2).cal(:,fr{1}),'color',cmap(40,:),'linestyle','-.');
grid on
% xlim([120 189])
xlabel('% Stance')
ylabel('Power [W]')
legend('UD tib','UD Talus','UD mt1 anterior plate','UD cal posterior plate')%,'sum cal/mt1','midtarsal','diff btw tal and sum')
makeNicePlotsFunction    


% plot the joint powers
figure;
hold on;
plot(stance,PUDsum.tib(:,fr{1})-PUDsum.tal(:,fr{1}),'color',cmap(30,:),'linestyle','-.')
plot(stance,PUDsum.tal(:,fr{1})-PUD(2).cal(:,fr{1}),'color',cmap(5,:),'linestyle','-.')
plot(stance,PUDsum.tal(:,fr{1})-(PUD(1).mt1(:,fr{1})),'color',cmap(20,:),'linestyle','-.')
grid on;

xlabel('% Stance')
ylabel('Power [W]')
legend('Talocrural','Subtalar','Midtarsal')
makeNicePlotsFunction


figure;
plot(cumtrapz(PUDsum.tib(:,fr{1})-PUDsum.tal(:,fr{1})),'color',cmap(30,:),'linestyle','-.'); hold on;
plot(cumtrapz(PUDsum.tal(:,fr{1})-PUD(2).cal(:,fr{1})),'color',cmap(5,:),'linestyle','-.')
plot(cumtrapz(PUDsum.tal(:,fr{1})-(PUD(1).mt1(:,fr{1}))),'color',cmap(20,:),'linestyle','-.')
% figure;
% for f =150 %120:185
%     b = 3;
%     plot3quick(mocap_data.force_data(1).ForcePlateLocation([1:4,1],1:3)');hold on;
%     plot3quick(mocap_data.force_data(2).ForcePlateLocation([1:4,1],1:3)');
%     plot3quick_scatter(p_com.(bonesCell{b})(:,f));
%     patch('faces',cns.cal(1:50:end,1:3) ,'vertices',transformPoints(T.cal(:,:,f),pts.cal,0),'facealpha',0.5,'edgecolor',cmap(12,:),'facecolor',cmap(12,:))
%     patch('faces',cns.mt1(1:50:end,1:3) ,'vertices',transformPoints(T.mt1(:,:,f),pts.mt1,0),'facealpha',0.5,'edgecolor',cmap(32,:),'facecolor',cmap(32,:))
%
%     plotvector3(cop_foot{1}(:,f),force{1}(:,f)/1500);
%     plotvector3(cop_foot{2}(:,f),force{2}(:,f)/1500);
%     plotvector3(cop_new(:,f),FT(:,f)/1500);
%
% %     plotvector3(p_com.(bonesCell{b})(:,f),u);
%
%     view(-114,7)
%     axis equal
%     hold off
%     pause(0.5)
% end
%
%
%
