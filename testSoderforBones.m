close all;
clear;
clc;

[file,taldir] = uigetfile('E:/*.iv');
% taldir = 'E:\SOL001_VISIT2\Models\IV\';
% file = 'SOL001B_tal_aligned.iv';
refDir = 'E:\Co-ordinateSystems\TalusRef\';
tal.Ref = fullfile(refDir, 'refTalus.iv');
tal.Dome = fullfile(refDir, 'talarDome.iv');
tal.Calc = fullfile(refDir, 'calcSurf.iv');
tal.Nav = fullfile(refDir, 'navSurf.iv');

surf_names = fields(tal); % get the list of the surfaces we are referencing

[pts.Ref,cns] = read_vrml_fast(tal.Ref);
[pts.Dome,~] = read_vrml_fast(tal.Dome);
[pts.Calc,~] = read_vrml_fast(tal.Calc);
[pts.Nav,~] = read_vrml_fast(tal.Nav);
cns = cns+1;
% get an initial set of inertial axes
[cent,~,~,~,CoM_eigenvectors,~,~,~,~,~] = mass_properties(pts.Ref,cns);

    T_init = eye(4);
    T_init(1:3,1:3) = CoM_eigenvectors;
    T_init(1:3,4) = cent';
    
ds_amt = 1:8:length(pts.Ref); % the amount to downsample
pts.RefDown = pts.Ref(ds_amt,: );
% get the indices of the points in each of the surfaces on the references
for sf = 2:4 % for the dome, calc and nav surfaces
    
    [pts.(surf_names{sf}),~] = read_vrml_fast(tal.(surf_names{sf}));
    [~,~,iRef.(surf_names{sf})] = intersect(pts.(surf_names{sf}),pts.RefDown,'rows');
end

pts.RefDown = transformPoints(T_init,pts.RefDown,-1);

pc.Ref = pointCloud(pts.Ref);
pc.RefDown = pointCloud(pts.RefDown);
pc.Dome = select(pc.RefDown,iRef.Dome);
pc.Calc = select(pc.RefDown,iRef.Calc);
pc.Nav= select(pc.RefDown,iRef.Nav);

tol = 1;
modelNav = pcfitsphere(pc.Nav,tol);
modelCalc = pcfitsphere(pc.Calc,tol);
modelDome = pcfitcylinder(pc.Dome,tol);

% for scaling purposes
xRef = diff(pc.RefDown.XLimits);
yRef = diff(pc.RefDown.YLimits);
zRef = diff(pc.RefDown.ZLimits);

%% Now load in the new talus and segment the new surfaces
pcNew = [];

[ptsNew.Raw,~] = read_vrml_fast(fullfile(taldir,file));

% make the raw point cloud
pcNew.Raw = pointCloud(ptsNew.Raw);

% downsample the point cloud
l1 = length(ptsNew.Raw);
pcNew.RawDown = pointCloud(ptsNew.Raw(1:10:end,:),'color',repmat([.75 .75  .75],length(1:10:l1),1));

% get the initial rotation from the reference bone 
T_affR = eye(4);
T_affR(1:3,1:3) = T_init(1:3,1:3);

% orient the new bone to match the reference one
[T_icp,pcNew.Anat,rmse] = pcregistericp(pcNew.RawDown,pc.RefDown,'Tolerance',[0.001 0.005],'MaxIterations',100,'InitialTransform',affine3d(T_affR));

% determine the average scaling factor
xn = diff(pcNew.Anat.XLimits);
yn = diff(pcNew.Anat.YLimits);
zn = diff(pcNew.Anat.ZLimits);

% scale it to the reference bone
scale_fact = mean([xRef/xn,yRef/yn,zRef/zn]);
scale_mat = eye(4,4);
scale_mat(1:3,1:3) = scale_mat(1:3,1:3) * scale_fact;
T_aff = affine3d(scale_mat);
pcNew.AnatScale = pctransform(pcNew.Anat,T_aff);

% register the points from the new, scaled, anatomically positioned model
% to the 
tic
[TT,pcNew.NewCPD] = pcregistercpd(pc.RefDown,pcNew.AnatScale);
toc

pcNew.Dome = select(pcNew.NewCPD,iRef.Dome);
pcNew.Nav = select(pcNew.NewCPD,iRef.Nav);
pcNew.Calc= select(pcNew.NewCPD,iRef.Calc);

modelNav = pcfitsphere(pcNew.Nav,tol);
modelCalc = pcfitsphere(pcNew.Calc,tol);
modelDome = pcfitcylinder(pcNew.Dome,tol);

modelNav.Radius
modelCalc.Radius
%%
figure;
hold on;
pcshow(pcNew.Dome)
h = plot(modelDome);
h.FaceAlpha = 0.2;

figure; hold on;
pcshow(pcNew.NewCPD)
pcshow(pcNew.Dome,'markersize',50)
pcshow(pcNew.Nav,'markersize',50)
pcshow(pcNew.Calc,'markersize',50)

h = plot(modelNav);
h.FaceAlpha = 0.2;
hold on;
h = plot(modelCalc);
h.FaceAlpha = 0.2;

plot3quick([modelNav.Center;modelCalc.Center]','w')
plot3quick([pcNew.Dome.Location(1,:);pcNew.Nav.Location(1,:)]','m')


STaxis = unit(modelNav.Center-modelCalc.Center);
TCaxis = unit(modelDome.Orientation);

%find an axis that will always be oriented anteriorly
ant_axis = unit(pcNew.Nav.Location(1,:)-pcNew.Dome.Location(1,:));
orient_axis = cross(ant_axis,STaxis);
if dot(orient_axis,TCaxis) < 0 % they are in the opposite direction
    TCaxis = - TCaxis;
end


y = cross(STaxis,TCaxis);
yu = unit(y);

% for TC focused ACS
z = cross(TCaxis,yu);
zu = unit(z);

% for ST focused ACS
x = cross(yu,STaxis);
xu = unit(x);

T_ACS_ST_anat = eye(4);
T_ACS_ST_anat(1:3,1:3) = [xu',yu',STaxis'];

T_ACS_TC_anat = eye(4);
T_ACS_TC_anat(1:3,1:3) = [TCaxis',yu',zu'];

T_icp_R =eye(4);
T_icp_R(1:3,1:3) = T_icp.T(1:3,1:3)';
T_icp_R(1:3,4) =  T_icp.T(4,1:3);

T_ACS_ST_CT = invTranspos9e(T_icp_R) * T_ACS_ST_anat;
T_ACS_TC_CT = invTranspose(T_icp_R) * T_ACS_TC_anat;

figure;
hold on;
pcshow(pcNew.Raw);
plotPointsAndCoordSys1([],T_ACS_ST_CT);
plotPointsAndCoordSys1([],T_ACS_TC_CT);

figure;hold on;
pcshow(pcNew.NewCPD)
plotPointsAndCoordSys1([],T_ACS_TC_anat);
plotPointsAndCoordSys1([],T_ACS_ST_anat);

