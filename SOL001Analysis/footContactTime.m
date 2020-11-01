
filedir = 'C:\Users\Lauren\Desktop\TempMocapSOL001\';

list_files = dir([filedir '*.mat']);

close all

for i = 1:length(list_files)
    mocap_file = [filedir list_files(i).name];
   temp_data = MocapDataTransform(mocap_file,'filter','adaptive','forceCutOff',[40 60],'saveProc','off','lab','SOL','resample','mocap');
   
   
    
    force1 = temp_data.force_data(1).Force;
    force2 = temp_data.force_data(2).Force;
    forceT = force1+force2;
    normForce = norm3d(forceT);
    
    ind = find(normForce(1:8000) > 400) ;%| (force2(3,1:1000) > 100));
    s1 = ind(1) -50; e1 = ind(end) + 50;
    ind = find(normForce(s1:e1) > 20);% | (force2(3,s1:e1) > 20));
    s2 = s1 + ind(1)-1; e2 = s1 + ind(end)+2;
    temp_data.cropFrsForce = [s2 e2];
    
    
    imp = trapz((1/1000)*1:(e2-s2+1),forceT(2,s2:e2));
    
    save_data(i).trial = list_files(i).name;
    save_data(i).cont_time = [e2-s2+1]/1000;
    save_data(i).imp = imp;
    
    
%     
%     figure
%     hold on;
%     plot(normForce')
%     plot([s2 e2],force2(3,[s2 e2]),'o')
%     plot([s2 e2],force1(3,[s2 e2]),'o')
end


