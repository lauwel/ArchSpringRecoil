
% compute all the transforms for the trials
clc
fc = [12,20];
% walking trials


subjDir = 'E:\SOL001B\';
trialName = 'T0005_SOL001_nwalk_pref_barefoot';

% optimize any bones that need it:
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED.csv']);
beadFile = 'E:\SOL001_VISIT2\Models\bead_positions.txt';
camDir = fullfile(subjDir,'\Calibration\Set 1\Mayacam2\NoQNFilter\');
filterOpts.Type = 'none';
optimizeOpts.Type ='optimize';
optimizeOpts.Bones = {'tib','tal','cal','nav','cmm','cub','mt1','mt5'};

saveOpts.ivDir = fullfile(subjDir,'Models','IV','3Aligned Reduced',filesep);
saveOpts.animDir = fullfile(subjDir,trialName,'POS',filesep);
processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);


% now filter:
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED_optimized.csv']);
filterOpts.Type = '2D';
filterOpts.Fs = 250;
filterOpts.Fc = fc;
optimizeOpts.Type = 'none';
% saveOpts = [];

processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);
%% running
fc = [12,20];

subjDir = 'E:\SOL001B\';

% trialName = 'T0019_SOL001_nrun_pref_barefoot';
% trialName = 'T0024_SOL001_nrun_rfs_barefoot';
% trialName = 'T0033_SOL001_srun_pref_barefoot';
% trialName = 'T0034_SOL001_srun_pref_barefoot';
% trialName = 'T0030_SOL001_srun_pref_barefoot';
% trialName = 'T0039_SOL001_frun_pref_barefoot';
% trialName = 'T0026_SOL001_nrun_ffs_barefoot';
% trialName = 'T0088_SOL001_nrun_rfs_minimal';
trialName = 'T0104_SOL001_nrun_ffs_minimal';

%trialName = 'T0019_jog0001';
% optimize any bones that need it:
% camDir = fullfile(subjDir,'\Calibration\Set 3\Mayacam2\');
camDir = 'E:\SOL001B\Calibration\Set 3\Mayacam2\';

in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED.csv']);
beadFile = 'E:\SOL001B\Models\bead_positions.txt';
filterOpts.Type = 'none';
optimizeOpts.Type ='optimize';
optimizeOpts.Bones = {};%{'tib','tal','cal','nav','cmm','cub','mt1'};
saveOpts.ivDir = fullfile(subjDir,'Models','IV','3Aligned Reduced',filesep);
saveOpts.animDir = fullfile(subjDir,trialName,'POS',filesep);
processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);


% now filter:
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED_optimized.csv']);
filterOpts.Type = '2D';
filterOpts.Fs = 250;
filterOpts.Fc = fc;
optimizeOpts.Type = 'none';
% saveOpts = [];

saveOpts.animDir = fullfile(subjDir,trialName,'POS',filesep);
processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);
% saveOpts.animDir = 'P:\Data\2019-05-02 SOL001_Visit2\';
% saveOpts.ivDir = 'P:\Data\2019-05-02 SOL001_Visit2\Models\IV\';

%% Visit 1 hop


fc = [12,20];
% hop trials


subjDir = 'E:\SOL001_VISIT1\';
trialName = 'T0014_hop108';
trialName = 'T0015_hop120';
trialName = 'T0016_hop132';
% trialName = 'T0005_SOL001_nwalk_pref_barefoot';
in2DFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'UNDISTORTED.csv']);
beadFile = 'E:\SOL001_VISIT1\Models\bead_positions.txt';
camDir = fullfile(subjDir,'Calibration\Calibration 1\Mayacam2\');
filterOpts.Type = '2D';%'none';%'3D';%
filterOpts.Fs = 125;
filterOpts.Fc = fc;
optimizeOpts.Type ='optimize';%'none';% 
optimizeOpts.Bones = {'tal','tib'};
saveOpts = [];
% saveOpts.animDir = 
processBeadData(in2DFile,beadFile,camDir,subjDir,trialName,filterOpts,optimizeOpts,saveOpts);

%% If additional frames have been tracked from autoscoper
bone = 'tib';
auto_dir = fullfile(subjDir,trialName,'Autoscoper','Autoscoped',filesep);
inFile = ls([auto_dir, '*' bone '*interp*.tra']);
autoscoper_dat = dlmread([auto_dir inFile],',');
Ttemp = convertRotation(autoscoper_dat,'autoscoper','4x4xn');

boneTdir = fullfile(subjDir,trialName,'BoneTransforms',filesep);
boneTfile = ls([boneTdir,trialName,'*transforms2DFILT*']);
load([boneTdir boneTfile])


% T.ses = repmat(nan(4,4),1,1,length(T.mt1));
% T.ph1 = repmat(nan(4,4),1,1,length(T.mt1));

dt = squeeze(diff(Ttemp(1,4,:)));
ind_dat = find(dt);
ind_dat = [ind_dat;ind_dat(end)+1];

ind_T = find(~isnan(squeeze(T.(bone)(1,4,:))));

first_fr = min([ind_T;ind_dat]);
end_fr = max([ind_T;ind_dat]);

for i = first_fr:end_fr
    
   if ismember(i,ind_dat) && ~ismember(i,ind_T) % i.e. it is in the autoscoper but not the bone transform file
   
        T.(bone)(:,:,i) = Ttemp(:,:,i);
            
        fprintf('Frame %i was replaced.\n',i)
   end
        
end

save([boneTdir boneTfile],'T')

    fprintf('Bone transform file is rewritten for bone %s: %s \n',bone,[boneTdir boneTfile])
    
%% save the two approaches on the server
% subj_dir = 'P:\Data\2019-05-02 SOL001_Visit2\';

% 
% anim_dir = 'P:\Data\2019-05-02 SOL001_Visit2\T0024_SOL001_nrun_rfs_barefoot\POS\QuaternionFiltered\';
% animateXMATransforms(subj_dir,trial_name,fc,anim_dir)
% 
% 
% anim_dir = 'P:\Data\2019-05-02 SOL001_Visit2\T0024_SOL001_nrun_rfs_barefoot\POS\BeadFiltered\';
% [~] = projectXMA2DPoints(inFile,camFolder,subj_dir,trial_name,fc,'none',[],anim_dir);
%% interpolate a distorted 2D points file
% trial_name = 'T0009_SOL001_nwalk_pref_barefoot';

% trial_name = 'T0024_SOL001_nrun_rfs_barefoot';
inFile = fullfile(subjDir,trialName,'XMA_csv',[trialName 'DISTORTED.csv']);
interpolateXMA2DPoints(inFile,[],'ph1')


