% Solve the HF equations

close all; 
clear all;

% helium
n = 10.;
box = 10.;
dx = box/n;
nelectrons = 2;
% grid3d = zeros(n,n,n);
% flat_grid3d = reshape(grid3d,[1,n^3]);

% electronic hamiltonian
lap = lap3d(n);
ionInt = rIonCoul(n, 2., [box/2.,box/2.,box/2.]/dx); % ion in the centre
H_elec = (-0.5/dx^2)* lap + ionInt;

% coulomb and exchange hamiltonians
guess = ones(n^3,nelectrons/2)./(n.^3);
[HeE, HeState, HeEtot] = HFsolver(n, H_elec, guess, dx, 0.000000005);
disp(['The He total energy:  ', num2str(HeEtot), ' Ha'])
disp(['The He first orbital energy:  ', num2str(HeE(2)), ' Ha'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hydrogen molecule (H2) the coordinates
n = 10;
box = 10.;
dx = box/n;
% grid3d = zeros(n,n,n);
% flat_grid3d = reshape(grid3d,[1,n^3]);

blength = 1.398397232; % bond length in atomic units
% location of ions in the grid
h1 = [(box+blength)/2., box/2., box/2.]./dx;
h2 = [(box-blength)/2., box/2., box/2.]./dx;

lap = lap3d(n);
ionInt = rIonCoul(n, 1., h1) + rIonCoul(n, 1., h2);
H_elec = (-0.5/dx^2)*lap + ionInt;

% coulomb and exchange hamiltonians
guess = ones(n^3,nelectrons/2)./sqrt(n^3);
[H2E, H2State, H2Etot] = HFsolver(n, H_elec, guess, dx, 0.0005);

disp(['The H2 total energy:  ', num2str(H2Etot), ' Ha'])
disp(['The H2 first orbital energy:  ', num2str(H2E(2)), ' Ha'])

% parameters for plott
% 
% for i = 2:4
%     for j = 1:2:n
%         figure;
%         hold on;
%         psi = reshape(HeState(:,i), [n,n,n]);
%         surf((psi(:,:,j)));
%         title(['He atom orbital n = ',num2str(i-1), ' at Z slice ', num2str(j)])
%         view(2)
%         colorbar
%         xlim([1,10]); ylim([1,10]);
%         hold off;
%         saveas(gcf,['He_',num2str(i-1),'slice_',num2str(j),'.png'])
%     end
% end
% 
% for i = 2:4
%     for j = 1:2:n
%         figure;
%         hold on;
%         psi = reshape(H2State(:,i), [n,n,n]);
%         surf(squeeze(psi(:,:,j)));
%         title(['H2 orbital n = ',num2str(i-1), ' at Z slice ', num2str(j)])
%         view(2)
%         xlim([1,10]); ylim([1,10]);
%         colorbar
%         hold off;
%         saveas(gcf,['H2_',num2str(i-1),'slice_',num2str(j),'.png'])
%     end
% end
