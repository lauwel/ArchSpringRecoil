


subj = {'10','11'};
for ss = 1:2
folder_base = 'P:\VAFootStudy\';
load(['P:\VAFootStudy\Data\XBSP000' subj{ss} '\Models\IV\inertialorientations\bonestructsubject' subj{ss} '.mat'])
Tct(ss).tibACS = bonestruct.tib.inertiaACS;
Tct(ss).talACS = bonestruct.tal.inertiaACS;
tal_temp_file = ['P:\VAFootStudy\Data\XBSP000' subj{ss} '\Trial008\XROMM\Tracking\XBSP000' subj{ss} '_Trial008_tal_interp.tra'];
tib_temp_file = ['P:\VAFootStudy\Data\XBSP000' subj{ss} '\Trial008\XROMM\Tracking\XBSP000' subj{ss} '_Trial008_tib_interp.tra'];
tal_auto = csvread(tal_temp_file);
tib_auto =csvread(tib_temp_file);
Tx(ss).tal = convertRotation(tal_auto,'autoscoper','4x4xn');
Tx(ss).tib = convertRotation(tib_auto,'autoscoper','4x4xn');

[a b c] = eulerYZX(Tx(ss).tib,Tx(ss).tal,Tct(ss).tibACS,Tct(ss).talACS );
figure(1); hold on;
plot(a)
figure(2); hold on;
plot(b)
figure(3); hold on;
plot(c)
end
% 