%% Load Structure 

close all; clc; clear;

%Loadstructure 
[~,xy,nnod,~,idb,ndof,incid,l,gamma,m,EA,EJ,posit,~,pr]=loadstructure ;

% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

freq = zeros(length(l),1);
fmax = 20; %Hz
omax = fmax*2*pi; %rad/s
i
for i=1:length(l)
    freq(i) = ((pi/l(i))^2)*sqrt((EJ(i)/(m(i))));
    
    if freq(i)/omax <= (2) % safety condition to verify Lmax
 fprintf ('safety factor is less than 2.0 and therefore not ok \n');
    end
end
% %% calcalute the Lm of each element
% scof= 1.2;
% Omax= scof*fmax*2*pi;
% Lm  = sqrt(pi^2/Omax*sqrt(EJ/m))
% 
% %% load input file to obtian the structe data
% [file_i,xy,nnod,sizee,idb,ndof,incidenze,l,gamma,m,EA,EJ,posiz,nbeam,pr] = loadstructure;
% 
% %% Plot undeformed structure
% dis_stru(posiz,l,gamma,xy,pr,idb,ndof);
% xlabel('x [m]'); ylabel('y [m]')
% 
% % Assembling of elements from .inp file 
% [M,K]=assem(incidenze,l,m,EA,EJ,gamma,idb);
% % 
% % check IDB matrix and dof
% disp('n_node  IDB-matix')
% [[1:size(idb,1)]', idb]
% 
% %% contribution of lumped mass to M
% idof_mc = idb(6,2); % lumped mass has ONLY 1 dof-y
% M(idof_mc, idof_mc)=M(idof_mc, idof_mc)+Mc;
% 
% %% contribution of spings to K
% % due to k1
% idof_k1 = idb(3,1);
% K(idof_k1,idof_k1) = K(idof_k1,idof_k1)+k1;
% % due to k2
% Mk2 = [k2 -k2; -k2 k2];
% idof_k2 = [idb(5,2) idb(6,2)];
% K(idof_k2,idof_k2) = K(idof_k2,idof_k2)+Mk2;
% % due to k3
% K_k3_local = [1 0 0 -1 0 0]'*k3*[1 0 0 -1 0 0];% in local reference wrt x
% g = 3/4*pi;% inclination of k3 wrt global frame
% lambda_k3 = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];%3x3
% Lambda_k3 = [lambda_k3 zeros(3,3); zeros(3,3) lambda_k3];
% K_k3_global = Lambda_k3'*K_k3_local*Lambda_k3;%6x6
% idof_k3   = [idb(2,:) idb(4,:)];%1x6
% K(idof_k3,idof_k3) = K(idof_k3,idof_k3)+K_k3_global;
% 
% %% partitioning of matrices 
% MFF=M(1:ndof,1:ndof);
% KFF=K(1:ndof,1:ndof);
% 
% %% natural frequencies and modes of vibration
% [eigenvectors, eigenvalues]=eig(MFF\KFF);% same as eig(inv(MFF)*KFF)
% freq = sqrt(diag(eigenvalues))/2/pi;
% [frqord,ind]=sort(freq);
% 
% % number of modes to be plotted 
% nmodes = 4;
% i_modes= ind(1:nmodes);
% 
% %scaling factor only for visualization (same for all plots, can be adjusted)
% scale_factor = 1;
% 
% for ii=1:nmodes
%     figure
%     diseg2(eigenvectors(:,i_modes(ii)),scale_factor,incidenze,l,gamma,posiz,idb,xy);% plot deformed structure
%     title(['Mode ' num2str(ii) ', Freq f_',num2str(ii),'=',num2str(frqord(ii)),' Hz'])
% end
% 
% %% alfa and beta values for structural damping
% A = [1/(2*2*pi*frqord(1))  2*pi*frqord(1)/2; 
%      1/(2*2*pi*frqord(2))  2*pi*frqord(2)/2; 
%      1/(2*2*pi*frqord(3))  2*pi*frqord(3)/2];
% ab = (A'*A)^-1*A'*[h1;h2;h3];
% 
% alfa = ab(1);
% beta = ab(2);
% C = alfa*M+beta*K;
% %% Contribution due to c2
% Mc2 = [c2 -c2; -c2 c2];
% idof_c2 = [idb(5,2) idb(6,2)];
% C(idof_c2,idof_c2) = C(idof_c2,idof_c2)+Mc2;
% 
% CFF = C(1:ndof,1:ndof);
% %% FRF
% % vector of frequency to be defined
% vett_f = 0:df:fmax; % frequency range of interest
% omega  = vett_f*2*pi; % [Hz]->[rad/s]
% nodA   = 4;% node number of the application point A
% % Forcing vector 
% F0 = zeros(ndof,1);
% idof_Ay = idb(nodA,2); % index of dof at node A(4) in vertical direction
% F0(idof_Ay) = 1;       % unit force
% for k=1:length(vett_f)
%     AA = -omega(k)^2*MFF+sqrt(-1)*omega(k)*CFF+KFF;% inverse of frequency response matrix
%     xx(:,k) = AA\F0;% row-index of dof , column-frequecy range
%     aa(:,k) = -omega(k)^2*xx(:,k);% row-index of dof , column-frequecy range
% end
% 
% %outputs 
% nodA = 4;% node number of the application point A
% nodB = 5;% node number of B
% nodC = 3;% node number of C
% 
% idof_Ay = idb(nodA,2); % index of dof at node A in vertical direction
% idof_By = idb(nodB,2); % index of dof at node B in vertical direction
% idof_Cy = idb(nodC,2); % index of dof at node C in vertical direction
% 
% FRF_dAy = xx(idof_Ay,:);% disp A :driving point(direct) FRF
% FRF_aAy = aa(idof_Ay,:);% acceleration of A
% 
% FRF_dBy = xx(idof_By,:);% disp B: cross FRF
% FRF_aBy = aa(idof_By,:);% acceleration of B
% 
% FRF_dCy = xx(idof_Cy,:);% disp C: cross FRF
% FRF_aCy = aa(idof_Cy,:);% acceleration of C
% 
% magAd = abs(FRF_dAy);
% phsAd = angle(FRF_dAy);
% magAa = abs(FRF_aAy);
% phsAa = angle(FRF_aAy);
% 
% magBd = abs(FRF_dBy);
% phsBd = angle(FRF_dBy);
% magBa = abs(FRF_aBy);
% phsBa = angle(FRF_aBy);
% 
% magCd = abs(FRF_dCy);
% phsCd = angle(FRF_dCy);
% magCa = abs(FRF_aCy);
% phsCa = angle(FRF_aCy);
% %%
% figure
% subplot 311;semilogy(vett_f,magAd,vett_f,magBd,vett_f,magCd);grid;
% hold on
% title('FRF: Vertical Displacements')
% ylabel(['|Y/F| [m/N]'])
% legend('Ay','By','Cy')
% subplot 312;plot(vett_f,phsAd,vett_f,phsBd,vett_f,phsCd);grid
% title('angle')
% hold on
% ylabel(['\Psi [rad]'])
% xlabel('Freq [Hz]')
% subplot 313;
% plot(vett_f,unwrap(phsAd),vett_f,unwrap(phsBd),vett_f,unwrap(phsCd));grid
% ylabel(['\Psi [rad]'])
% xlabel('Freq [Hz]')
% title('unwrap')
% 
% figure
% subplot 211;semilogy(vett_f,magAa,vett_f,magBa,vett_f,magCa);grid;
% title('FRF: Vertical Accelerations')
% ylabel(['|Ypp/F| [m/s^2/N]'])
% legend('Ay','By','Cy')
% subplot 212;plot(vett_f,phsAa,vett_f,phsBa,vett_f,phsCa);grid
% ylabel(['\Psi [rad]'])
% xlabel('Freq [Hz]')
% 
% %% FRF-constraint force
% MCF = M(ndof+1:end,1:ndof);
% KCF = K(ndof+1:end,1:ndof);
% CCF = C(ndof+1:end,1:ndof);
% 
% for k = 1:length(vett_f)
%     A = -omega(k)^2*MFF+sqrt(-1)*omega(k)*CFF+KFF;
%     xx(:,k) = A\F0;% row-index of dof , column-frequecy range
%     rr(:,k) = (-omega(k)^2*MCF+sqrt(-1)*omega(k)*CCF+KCF)*xx(:,k);
% end
% 
% nodC    = 1;
% idof_Dr = idb(nodC,3);% rotation at D
% FRF_DM  = rr(idof_Dr-ndof,:);% rr has the same number of rows as constrained dof 
% mag     = abs(FRF_DM);% constraint moment at D
% phs     = angle(FRF_DM);
% 
% figure
% subplot 211;
% plot(vett_f,mag);grid;
% title('FRF: Clamp Moment (node 1)')
% ylabel(['|M_D/F| [m]'])
% subplot 212;
% plot(vett_f,phs);grid
% ylabel(['\Psi [rad]'])
% xlabel('Freq [Hz]')