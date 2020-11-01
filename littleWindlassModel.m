



r = 11.5;
l_arch = 250;
theta = pi()/2 + (0:0.1:pi()/4);

l1 = sqrt(l_arch^2 + r^2);
l2 = r*theta;
pf_length = l1+l2;

figure;

plot((theta-pi()/2)*180/pi(),pf_length)
ylabel('PF length')
xlabel('MTP dorsiflexion')

%% or, how much would the arch change length

pf_length = 270;
l2 = r*theta;
l1 = pf_length - l2;
l_arch = sqrt(l1.^2 - r^2);

figure;

plot((theta-pi()/2)*180/pi(),l_arch)
ylabel('arch length')
xlabel('MTP dorsiflexion')
sens = mean(diff(l_arch)./diff((theta-pi()/2)*180/pi()));
%% for every mm or %strain of plantar fascia, how much MTP dors do we get?
% and how much does that affect arch length?
col_map = colormap('parula');
pf_resting = 265;
figure;
hold on;
ii = 1;
for pf_length = 270:0.1:275
% pf_length = 270;
l2 = r*theta;
l1 = pf_length - l2;
l_arch = sqrt(l1.^2 - r^2);



plot((theta-pi()/2)*180/pi(),l_arch,'color',col_map(ii,:))
ii = ii+1;
end
% lege = num2str([270:0.1:275]'/pf_resting);
caxis([[270,275]'/pf_resting])
colorbar
ylabel('arch length')
xlabel('MTP dorsiflexion')