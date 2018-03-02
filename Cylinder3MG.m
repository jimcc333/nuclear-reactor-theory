function [] = Cylinder3MG(radius1, thickness2, thickness3, enrichment1, ...
    enrichment2, step_size, H2O_volume_fraction1, H2O_volume_fraction2)
%% 4-GROUP 3-REGION CONCENTRIC CYLINDRICAL INFINE REACTOR
% Author: Cem Bagdatlioglu, 10/2017 
% 
%   This program solves the 4-group reactor eqution for a 3-region infinite
%   cylinder where the two inner zones undergo fission and the outer region
%   does not. The provided material properties for the inner two regions
%   are for uranium (with user-defined enrichment and water fraction) and 
%   the outer zone is a lead reflector. The dimensions of the zones are 
%   also specified by the user. The program returns the flux profile and 
%   multiplication factor. 
%
%   Energy groups:
%     Group 1: > 1 MeV
%     Group 2: 1 keV - 1 MeV
%     Group 3: 1 eV - 1 keV
%     Group 4: < 1 eV

%% Check and assign inputs
if nargin ~= 8      % If all 8 inputs are not provided use default values
    disp('Necessary inputs not specified, using default values')
    R(1) = 20.0;    % outer radius [cm] of Region 1
    R(2) = 35.0;    % outer radius [cm] of Region 2
    R(3) = 60.0;    % outer radius [cm] of Region 3
    x(1) = 0.12;    % enrichment of region 1
    x(2) = 0.10;    % enrichment of region 2
    Delta = 0.1;    % Mesh Spacing [cm]
    h2ovf(1) = 0.70;   % Region 1 v/o occupied by water
    h2ovf(2) = 0.70;   % Region 2 v/o occupied by water
    fprintf('Radii: %.1f, %.1f, %.1f [cm], Delta: %.3f\n', ...
            R(1), R(2), R(3), Delta)
    fprintf('Enrichments: %.2f, %.2f Water volume fractions: %.2f, %.2f \n\n', ...
        x(1), x(2), h2ovf(1), h2ovf(2))
else
    R(1) = radius1;             % outer radius [cm] of Region 1
    R(2) = R(1) + thickness2;   % outer radius [cm] of Region 2
    R(3) = R(2) + thickness3;   % outer radius [cm] of Region 3
    x(1) = enrichment1/100;     % enrichment of region 1
    x(2) = enrichment2/100;     % enrichment of region 2
    Delta = step_size;          % Mesh Spacing [cm]
    h2ovf(1) = H2O_volume_fraction1;    % Region 1 v/o occupied by water
    h2ovf(2) = H2O_volume_fraction2;    % Region 2 v/o occupied by water
    fprintf('Radii: %.1f, %.1f, %.1f [cm], Delta: %.3f\n', ...
            R(1), R(2), R(3), Delta)
    fprintf('Enrichments: %.2f, %.2f Water volume fractions: %.2f, %.2f \n\n', ...
        x(1), x(2), h2ovf(1), h2ovf(2))
end


%% Global Parameters and Microscopic Cross Sections
MAX_ITS = 100;      % Max iterations
TOLER = 0.01;       % fractional tolerance for convergence tests

% Physical properties and group constants
chi = zeros(4,1);     % fission neutron source probability distribution
chi(1) = 0.575;
chi(2) = 0.425;

% Water
N_h2o = 0.030;          % *10^24 atoms/cm^3

sg_h2o(1) = 0;          % (n,gamma) cross section, barns
sg_h2o(2) = 0;          % (n,gamma) cross section, barns
sg_h2o(3) = 0.035;      % (n,gamma) cross section, barns
sg_h2o(4) = 0.57;       % (n,gamma) cross section, barns

str_h2o(1) = 3.08;      % transport cross section,  barns
str_h2o(2) = 10.52;     % transport cross section, barns
str_h2o(3) = 16.55;     % transport cross section, barns
str_h2o(4) = 68.6;      % transport cross section, barns

ss_h2o = zeros(4);      % create 4x4 matrix filled with zeros
ss_h2o(1,2) = 2.81;     % scattering transfer fn, barns
ss_h2o(2,3) = 4.04;   	% ss(i,j) = cross section for scattering
ss_h2o(3,4) = 4.14;     % from group i to group j

sr_h2o(1) = sg_h2o(1);	% removal cross section, barns
sr_h2o(2) = sg_h2o(2);	% removal cross section, barns
sr_h2o(3) = sg_h2o(3);	% removal cross section, barns
sr_h2o(4) = sg_h2o(4);	% removal cross section, barns

for i = 1:4               % add out-of-group scattering to removal cross section
    for j = 1:4           % note: we've defined ss(i,i)=0 for convenience
        sr_h2o(i) = sr_h2o(i) + ss_h2o(i,j);   
    end
end

% Uranium-235
N_235 = 0.048;          % *10^24 atoms/cm^3

sg_235(1) = 0.1;        % (n,gamma) cross section, barns
sg_235(2) = 0.3;        % (n,gamma) cross section, barns
sg_235(3) = 18.0;       % (n,gamma) cross section, barns
sg_235(4) = 97.0;       % (n,gamma) cross section, barns

sf_235(1) = 1.3;        % (n,fission) cross section, barns
sf_235(2) = 1.4;        % (n,fission) cross section, barns
sf_235(3) = 23.0;       % (n,fission) cross section, barns
sf_235(4) = 490.0;      % (n,fission) cross section, barns

nu_235(1) = 2.65;       % neutrons per fission
nu_235(2) = 2.55;       % neutrons per fission
nu_235(3) = 2.5;        % neutrons per fission
nu_235(4) = 2.5;        % neutrons per fission

str_235(1) = 4.7;       % transport cross section, barns
str_235(2) = 7.0;       % transport cross section, barns
str_235(3) = 51.0;      % transport cross section, barns
str_235(4) = 597.0;     % transport cross section, barns

ss_235 = zeros(4);      % create 4x4 matrix filled with zeros
ss_235(1,2) = 1.40;     % scattering transfer fn, barns 
ss_235(2,3) = 0.0;      % ss(i,j) = cross section for scattering
ss_235(3,4) = 0.01;     % from group i to group j

sr_235(1) = sg_235(1) + sf_235(1);  % removal cross section [b]
sr_235(2) = sg_235(2) + sf_235(2);
sr_235(3) = sg_235(3) + sf_235(3);
sr_235(4) = sg_235(4) + sf_235(4);
for i = 1:4             % add out-of-group scattering to removal cross section
    for j = 1:4         % note: defined ss(i,i)=0 for convenience
        sr_235(i) = sr_235(i)+ss_235(i,j);   
    end
end

% Uranium-238
N_238 = 0.048;          % *10^24 atoms/cm^3
sg_238(1) = 0.04;       % (n,gamma) cross section, barns
sg_238(2) = 0.18;       % (n,gamma) cross section, barns
sg_238(3) = 0.80;       % (n,gamma) cross section, barns
sg_238(4) = 2.40;       % (n,gamma) cross section, barns

sf_238(1) = 0.53;       % (n,fission) cross section, barns
sf_238(2) = 0.0;        % (n,fission) cross section, barns
sf_238(3) = 0.0;        % (n,fission) cross section, barns
sf_238(4) = 0.0;        % (n,fission) cross section, barns

nu_238(1) = 2.65;       % neutrons per fission
nu_238(2) = 0.0;        % neutrons per fission
nu_238(3) = 0.0;        % neutrons per fission
nu_238(4) = 0.0;        % neutrons per fission

str_238(1) = 4.7;       % transport cross section, barns
str_238(2) = 7.0;       % transport cross section, barns
str_238(3) = 11.0;      % transport cross section, barns
str_238(4) = 13.0;      % transport cross section, barns

ss_238 = zeros(4);      % create 4x4 matrix filled with zeros
ss_238(1,2) = 2.10;     % scattering transfer fn: barns: 
ss_238(2,3) = 0.0;      % ss(i,j) = cross section for scattering
ss_238(3,4) = 0.01;     % from group i to group j

sr_238(1) = sg_238(1) + sf_238(1);  % removal cross section [b]
sr_238(2) = sg_238(2) + sf_238(2);
sr_238(3) = sg_238(3) + sf_238(3);
sr_238(4) = sg_238(4) + sf_238(4);

for i = 1:4               % add out-of-group scattering to removal cross section
    for j = 1:4           % note: we've defined ss(i,i)=0 for convenience
        sr_238(i) = sr_238(i) + ss_238(i,j);   
    end
end

%% Macroscopic cross sections and other properties by region
% First index, i: region number
% Second index, j: energy group number
for j = 1:4
    for i=1:2
    % --- For the u-235/u-238/water mixtures in regions 1 and 2 ---
    % Macroscopic removal x-s [1/cm]
        Sr(i,j) = (1.-h2ovf(i))*(N_235*x(i)*sr_235(j) + ...
            N_238*(1.-x(i))*sr_238(j)) + N_h2o*h2ovf(i)*sr_h2o(j);
        
    % Macroscopic transport x-s [1/cm]
        Str(i,j) = (1.-h2ovf(i))*(N_235*x(i)*str_235(j) + ...
            N_238*(1.-x(i))*str_238(j)) + N_h2o*h2ovf(i)*str_h2o(j);
        
    % Macroscopic fission x-s [1/cm]
        Sf(i,j) = (1.-h2ovf(i))*(N_235*x(i)*sf_235(j) + ...
            N_238*(1.-x(i))*sf_238(j));
        
    % Fission neutron production rate per unit flux 
    % == neutrons per fission * macroscopic fission cross section [1/cm]
        NuSf(i,j) = (1.-h2ovf(i))*(N_235*x(i)*sf_235(j)*nu_235(j) + ...
            N_238*(1.-x(i))*sf_238(j)*nu_238(j));
        
    % Diffusion coefficient [cm]
        D(i,j) = 1./3./Str(i,j);
        
    % Diffusion area [cm^2]: Use removal cross section
        LSquared(i,j) = D(i,j)/Sr(i,j);
        
    % Group to group scattering cross section: from group j to group e:
        for e=1:4
            Ss(i,j,e) = (1.-h2ovf(i))*(N_235*x(i)*ss_235(j,e) + ...
            N_238*(1.-x(i))*ss_238(j,e)) + N_h2o*h2ovf(i)*ss_h2o(j,e);
        end
    
    end
    
    % --- For the water reflector in region 3 ---
    % Macroscopic removal x-s [1/cm]
    Sr(3,j) = N_h2o*sr_h2o(j);
    
    % Macroscopic transport x-s [1/cm]
    Str(3,j) = N_h2o*str_h2o(j);
    
    % Diffusion coefficient [cm]
    D(3,j) = 1./3./Str(3,j);
    
    % Diffusion area [cm^2]: Use removal cross section
    LSquared(3,j) = D(3,j)/Sr(3,j);    
    
    % Group to group scattering cross section: from group j to group e:
    for e = 1:4
        Ss(3,j,e) = N_h2o*ss_h2o(j,e);
    end
end
    
% Vector of mesh point radii used for plotting
r = linspace(0,R(3),round(R(3)/Delta)+1);

% Number of mesh points in each region
N(1) = round(R(1)/Delta);   
N(2) = round((R(2)-R(1))/Delta);
N(3) = round((R(3)-R(2))/Delta);

% Total number of mesh points:
NTot = N(1)+N(2)+N(3)+1; % include the mesh point at the origin in Reg. 1 

%% Build matrix A
% Make an empty matrix of size NTotal x NTotal:
A = zeros(4*NTot);  

% Since the quantity D/Delta^2 appears frequently assign it a variable
for i = 1:3
    for e = 1:4
        dd2(i,e) = D(i,e)/Delta/Delta; 
    end
end

% Internal mesh points in Region 1. 
for j = 2:N(1)
    for e = 1:4
        % Define ei for indexing the coeffs. according to energy group
        ei = (e-1)*NTot; 
        % Handle terms containing D and Sigma_r
        A(ei+j,ei+j-1) = -dd2(1,e)*(2.*(j-1)-1)/(2.*(j-1));
        A(ei+j,ei+j)   = 2.*dd2(1,e) + Sr(1,e);
        A(ei+j,ei+j+1) = -dd2(1,e)*(2.*(j-1)+1)/(2.*(j-1));
        
        % Handle terms containing Sigma_s(this group -> some other group)
        for e_dest = 1:4
            e_desti = (e_dest-1)*NTot;
            A(e_desti+j,ei+j) = A(e_desti+j,ei+j) - Ss(1,e,e_dest);
        end
    end
end

% Internal mesh points in Region 2. 
for j = N(1)+2:N(1)+N(2)
    for e = 1:4
        % Define ei for indexing the coeffs. according to energy group
        ei = (e-1)*NTot;
        
        % Handle terms containing D and Sigma_r
        A(ei+j,ei+j-1) = -dd2(2,e)*(2.*(j-1)-1)/(2.*(j-1));
        A(ei+j,ei+j)   = 2.*dd2(2,e) + Sr(2,e);
        A(ei+j,ei+j+1) = -dd2(2,e)*(2.*(j-1)+1)/(2.*(j-1));
        
        % Handle terms containing Sigma_s(this group -> some other group)
        for e_dest = 1:4
            e_desti = (e_dest-1)*NTot;
            A(e_desti+j,ei+j) = A(e_desti+j,ei+j) - Ss(2,e,e_dest);
        end 
    end
end

% Internal mesh points in Region 3. 
for j = N(1)+N(2)+2 : N(1)+N(2)+N(3)
    for e = 1:4 
        ei = (e-1)*NTot;
        
        % Handle terms containing D and Sigma_r:
        A(ei+j,ei+j-1) = -dd2(3,e)*(2.*(j-1)-1)/(2.*(j-1));
        A(ei+j,ei+j)   = 2.*dd2(3,e) + Sr(3,e);
        A(ei+j,ei+j+1) = -dd2(3,e)*(2.*(j-1)+1)/(2.*(j-1));
        
        % Handle terms containing Sigma_s(this group -> some other group):
        for e_dest = 1:4
            e_desti = (e_dest-1)*NTot;
            A(e_desti+j,ei+j) = A(e_desti+j,ei+j) - Ss(3,e,e_dest);
        end  
    end
end

% Boundary and interface conditions
% 1.  Symmetry
for e = 1:4 
    ei = (e-1)*NTot;
    A(ei+1,ei+1) = 1;
    A(ei+1,ei+2) = -1;
    
% 2.  Flux Matching - satisfied automatically since the rightmost mesh
% point in region 1 is the same as the leftmost mesh point in region 2

% 3.  Current Matching - Regions 1 and 2
    if N(2) ~= 0
        A(ei+N(1)+1,ei+N(1)) = D(1,e);
        A(ei+N(1)+1,ei+N(1)+1) = -D(1,e)-D(2,e);
        A(ei+N(1)+1,ei+N(1)+2) = D(2,e);
    end
    
% 4.  Flux Matching - satisfied automatically since the rightmost mesh
% point in region 2 is the same as the leftmost mesh point in region 3

% 5.  Current Matching - Regions 2 and 3
    if N(3) ~= 0  
        A(ei+N(1)+N(2)+1,ei+N(1)+N(2)) = D(2,e);
        A(ei+N(1)+N(2)+1,ei+N(1)+N(2)+1) = -D(2,e)-D(3,e);
        A(ei+N(1)+N(2)+1,ei+N(1)+N(2)+2) = D(3,e);
    end
    
% 6.  Vacuum 
    A(ei+N(1)+N(2)+N(3)+1,ei+N(1)+N(2)+N(3)+1) = 1;
end

%% Build fission source term S = F*phi/k
F = zeros((N(1)+N(2)+N(3)+1)*4);

% Fission source coefficients for interior mesh points in region 1 
for j = 2:N(1)
    for e = 1:4   %% calculate contribution to each group:
        ei = (e-1)*NTot;
        for e_source = 1:4  %% add contribution from all groups:
            e_srci = (e_source-1)*NTot;  
            F(ei+j,e_srci+j) = NuSf(1,e_source)*chi(e);
        end 
    end
end

% Fission source coefficients for interior mesh points in region 2 
for j = N(1)+2:N(1)+N(2)
    for e = 1:4
        ei = (e-1)*NTot;
        for e_source = 1:4
            e_srci = (e_source-1)*NTot;  
            F(ei+j,e_srci+j) = NuSf(2,e_source)*chi(e);
        end 
    end
end

%% Initial guesses
% Flux profile initial guess
phi_guess = zeros(4*NTot,1);

for j = 1:4*NTot
    phi_guess(j) = 1;
end

% Initial guess for the eigenvalue k
k_guess = 1.0;

% Source term F*phi corresponding to the initial guess for phi
S_guess = F*phi_guess;

%% Storage matrices for plotting of results later on
phi_stored = zeros(4*NTot,MAX_ITS);
s_stored   = zeros(4*NTot,MAX_ITS);
k_stored   = zeros(MAX_ITS,1);
fiss_rate_stored = zeros(4*NTot,MAX_ITS);

phi_stored(:,1) = phi_guess;
s_stored(:,1)   = S_guess;
k_stored(1)     = k_guess;
for j = 1:N(1)+1
    for e = 1:4
        ei = (e-1)*NTot;
        fiss_rate_stored(ei+j,1) = phi_stored(ei+j,1)*Sf(1,e);
    end
end

for j = N(1)+2:N(1)+N(2)+1
    for e = 1:4
        ei = (e-1)*NTot;
        fiss_rate_stored(ei+j,1)=phi_stored(ei+j,1)*Sf(2,e);
    end
end
fprintf('Starting flux calculation\n\n')

%% Begin iterative solution
for iteration = 2:MAX_ITS
% A*phi = F*phi_guess/k_guess (A*phi = S/k)
phi = (A\S_guess)./k_guess;

% Calculate the fission source term and eigenvalue arising from this flux
S = F*phi;

new_neut_prod_rate = 0.;
old_neut_prod_rate = 0.;

for e = 1:4
    ei = (e-1)*NTot;
    old_neut_prod_rate = old_neut_prod_rate + 1./4.*pi*NuSf(1,e)* ...
        phi_guess(ei+1)*Delta*Delta;
    new_neut_prod_rate = new_neut_prod_rate + 1./4.*pi*NuSf(1,e)* ...
        phi(ei+1)*Delta*Delta;
end

% Integrate over region 1 interior mesh points 
for j = 2:N(1)
    for e = 1:4
        ei = (e-1)*NTot;    
        old_neut_prod_rate = old_neut_prod_rate + 2.*pi*NuSf(1,e)* ...
            phi_guess(ei+j)*(j-1)*Delta*Delta;
        new_neut_prod_rate = new_neut_prod_rate + 2.*pi*NuSf(1,e)* ...
            phi(ei+j)*(j-1)*Delta*Delta;    
    end
end

j = N(1)+1;

% Portion of mesh interval j that belongs to region 1:
for e = 1:4
    ei = (e-1)*NTot;
    old_neut_prod_rate = old_neut_prod_rate + pi*NuSf(1,e)* ...
        phi_guess(ei+j)*Delta*Delta*((j-1)-1./4.);
    new_neut_prod_rate = new_neut_prod_rate + pi*NuSf(1,e)* ...
        phi(ei+j)*Delta*Delta*((j-1)-1./4.);  

% Portion of mesh interval j that belongs to region 2:
    old_neut_prod_rate = old_neut_prod_rate + pi*NuSf(2,e)* ...
        phi_guess(ei+j)*Delta*Delta*((j-1)+1./4.);
    new_neut_prod_rate = new_neut_prod_rate + pi*NuSf(2,e)* ...
        phi(ei+j)*Delta*Delta*((j-1)+1./4.);
end

% Interior region 2 mesh points
for j = N(1)+2:N(1)+N(2)
    for e = 1:4
        ei = (e-1)*NTot;
        old_neut_prod_rate = old_neut_prod_rate + 2.*pi*NuSf(2,e)* ... 
            phi_guess(ei+j)*(j-1)*Delta*Delta;
        new_neut_prod_rate = new_neut_prod_rate + 2.*pi*NuSf(2,e)* ...
            phi(ei+j)*(j-1)*Delta*Delta;    
    end
end

j = N(1)+N(2)+1;

% Portion of mesh interval j that belongs to region 2:
for e = 1:4
    ei = (e-1)*NTot;
    old_neut_prod_rate = old_neut_prod_rate + pi*NuSf(2,e)* ...
        phi_guess(ei+j)*Delta*Delta*((j-1)-1./4.);
    new_neut_prod_rate = new_neut_prod_rate + pi*NuSf(2,e)* ...
        phi(ei+j)*Delta*Delta*((j-1)-1./4.);  
end

% New estimate for k
k = new_neut_prod_rate/old_neut_prod_rate*k_guess;

% Store the new values of k, S, phi and fission rate:
phi_stored(:,iteration) = phi/k;
s_stored(:,iteration)   = S/k;
k_stored(iteration)     = k;
for j = 1:N(1)+1
    for e = 1:4
        ei = (e-1)*NTot;
        fiss_rate_stored(ei+j,iteration)=...
            phi_stored(ei+j,iteration)*Sf(1,e);
    end
end

for j = N(1)+2:N(1)+N(2)+1
    for e = 1:4
        ei = (e-1)*NTot;
        fiss_rate_stored(ei+j,iteration)=...
            phi_stored(ei+j,iteration)*Sf(2,e);
    end
end

% Check for convergence
converged = true;
for j = 1:length(S)
    if S_guess(j) > 0
        if((S(j)/S_guess(j) > 1.+TOLER) || (S(j)/S_guess(j) < 1.-TOLER))
            converged = false;
        end
    end
end
if(k/k_guess > 1.+TOLER || k/k_guess < 1.-TOLER)
    converged = false;
end

fprintf('Iteration %i. k: %.4f \n', iteration, k)

% If converged is still true, break out of the loop
if(converged)
    break;
end
% If not converged, update guess values
k_guess   = k;
phi_guess = phi;
S_guess   = S;

end

disp('Convergence achieved.');
fprintf('k: %.4f \n', k)

%% Plot data
colors=[1, 0, 0; 0.6, 0, 0; .51, .51, 0; 0, 0, 1];

% Plot Fission Rate Density
figure(2) % Figure 2 first so it's behind Figure 1 when they load on screen
clf reset;

fr_e = zeros(NTot, MAX_ITS, 4);
for e=1:4
    ei = (e-1)*NTot;
    fr_e(:,:,e) = fiss_rate_stored(ei+1:ei+NTot,:);
end

hold on
for e=1:4
    plot(r,fr_e(:,iteration,e),'LineWidth',2.0,'Color',colors(e,:));
end

max_y=max(max(fiss_rate_stored));
set(gca,'XLim',[0 R(3)]);
set(gca,'YLim',[0 max_y*1.1]);
xlabel('Radial position r [cm]');
ylabel('Fission Rate Density \Sigma_{f}\phi [#/cm^{3}-s]');
plot([R(1) R(1)],[0,max_y*1.1],'Color','k');
plot([R(2) R(2)],[0,max_y*1.1],'Color','k');
t1=sprintf('Enrichment = %g%%',x(1)*100);
t2=sprintf('Enrichment = %g%%',x(2)*100);
t3=sprintf('H_{2}O Vol Frac = %g%%',h2ovf*100);
text(R(1)/100,max_y/15,'Region 1','Color','k');
text(R(1)*101/100,max_y/15,'Region 2','Color','k');
text(R(2)*101/100,max_y/15,'Region 3','Color','k');
text(R(1)/100,max_y/50,t1,'Color','k');
text(R(1)*101/100,max_y/50,t2,'Color','k');
text(R(2)*101/100,max_y/50,'Reflector','Color','k');
text(R(1)/100,53*max_y/50,t3,'Color','k');
text(R(1)*101/100,53*max_y/50,t3,'Color','k');
text(r(length(r))/1.35,53*max_y/50,'Energy Group 1','Color',colors(1,:));
text(r(length(r))/1.35,50*max_y/50,'Energy Group 2','Color',colors(2,:));
text(r(length(r))/1.35,47*max_y/50,'Energy Group 3','Color',colors(3,:));
text(r(length(r))/1.35,44*max_y/50,'Energy Group 4','Color',colors(4,:));

% Plot Neutron Flux
figure(1)
clf reset;

phi_e = zeros(NTot, MAX_ITS, 4);
for e=1:4
    ei = (e-1)*NTot;
    phi_e(:,:,e) = phi_stored((ei+1):(ei+NTot),:);
end

hold on;

for e=1:4 
    plot(r,phi_e(:,iteration,e),'LineWidth',2.0,'Color',colors(e,:));
end
max_y=max(max(phi_stored));
set(gca,'XLim',[0 R(3)]);
set(gca,'YLim',[0 max_y*1.1]);
xlabel('Radial position r [cm]');
ylabel('Flux \phi [n/cm^{2}-s], arbitrary scaling');
plot([R(1) R(1)],[0,max_y*1.1],'Color','k');
plot([R(2) R(2)],[0,max_y*1.1],'Color','k');
t1=sprintf('Enrichment = %g%%',x(1)*100);
t2=sprintf('Enrichment = %g%%',x(2)*100);
t3=sprintf('H_{2}O Vol Frac = %g%%',h2ovf(1)*100);
t4=sprintf('H_{2}O Vol Frac = %g%%',h2ovf(2)*100);
text(R(1)/100,max_y/15,'Region 1','Color','k');
text(R(1)*101/100,max_y/15,'Region 2','Color','k');
text(R(2)*101/100,max_y/15,'Region 3','Color','k');
text(R(1)/100,max_y/50,t1,'Color','k');
text(R(1)*101/100,max_y/50,t2,'Color','k');
text(R(2)*101/100,max_y/50,'Reflector','Color','k');
text(R(1)/100,53*max_y/50,t3,'Color','k');
text(R(1)*101/100,53*max_y/50,t4,'Color','k');
text(r(length(r))/1.35,53*max_y/50,'Energy Group 1','Color',colors(1,:));
text(r(length(r))/1.35,50*max_y/50,'Energy Group 2','Color',colors(2,:));
text(r(length(r))/1.35,47*max_y/50,'Energy Group 3','Color',colors(3,:));
text(r(length(r))/1.35,44*max_y/50,'Energy Group 4','Color',colors(4,:));


