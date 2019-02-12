pizzi
run Europa_physics.m

delta = 23.56 ; % opening width
opening_rate = 1000/(365*24*3600); % 1km/yr
strain_rate  = 10e-10;
stress_rate = -strain_rate * E; %stress rate
Normal_stress = 40e3; % Stress on Europa
a = 100:10:10e3; b = a/0.5;

F = -0.071-0.535.*(a./b)+0.169.*(a./b).^2-0.09.*(a./b).^3+0.02.*(a./b).^4-1.071.*(b./a).*log(1-a./b);
dF = -0.535.*(1./b)+0.338.*(a./b.^2)-0.27.*(a.^2./b.^3)+0.08.*(a.^3./b.^4)+1.071*(b./a).*(log(1-a./b)./a+1./(b-a));

V = (E*opening_rate/4-stress_rate.*a.*F)./(Normal_stress.*(F+a.*dF));

figure,plot(a,V)