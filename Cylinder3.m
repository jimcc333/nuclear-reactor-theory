function [] = Cylinder3(radius1, thickness2, thickness3, enrichment1, ...
    enrichment2, step_size)
%% 3-REGION CONCENTRIC CYLINDRICAL INFINE REACTOR
% Author: Cem Bagdatlioglu, 10/2017
% 
%   This program solves the reactor eqution for a 3-region infinite
%   cylinder where the two inner zones undergo fission and the outer region
%   does not. The provided material properties for the inner two regions
%   are for uranium (with user-defined enrichment) and the outer zone is
%   a lead reflector. The dimensions of the zones are also specified by the
%   user. The program returns the flux profile and multiplication factor. 

%% Check and assign inputs
if nargin ~= 6      % If all 6 inputs are not provided use default values
    disp('Necessary inputs not specified, using default values')
    R(1) = 20.0;    % outer radius [cm] of Region 1
    R(2) = 35.0;    % outer radius [cm] of Region 2
    R(3) = 60.0;    % outer radius [cm] of Region 3
    x(1) = 0.12;    % enrichment of region 1
    x(2) = 0.10;    % enrichment of region 2
    Delta = 0.1;    % Mesh Spacing [cm]
    fprintf('Radii: %.1f, %.1f, %.1f [cm], Delta: %.3f\n', ...
            R(1), R(2), R(3), Delta)
    fprintf('Enrichments: %.2f, %.2f \n\n', x(1), x(2))
else
    R(1) = radius1;             % outer radius [cm] of Region 1
    R(2) = R(1) + thickness2;   % outer radius [cm] of Region 2
    R(3) = R(2) + thickness3;   % outer radius [cm] of Region 3
    x(1) = enrichment1/100;     % enrichment of region 1
    x(2) = enrichment2/100;     % enrichment of region 2
    Delta = step_size;          % Mesh Spacing [cm]
end

%% Global parameters
MAX_ITS = 100;      % Max number of iteration steps
TOLER = 0.00005;    % fractional tolerance for flux and k convergence
% Material properties
N_pb = 0.033;                                   % *10^24 atoms/cm^3
sigma_gamma_pb = 0.2;                           % barns
sigma_tr_pb = 3.0;                              % barns

N_u = 0.048;                                    % *10^24 atoms/cm^3

nu_u235 = 2.4;                                  % neutrons/fission
sigma_gamma_u235 = 0.7;                         % barns
sigma_f_u235 = 2.9;                             % barns
sigma_tr_u235 = 2.5;                            % barns
sigma_a_u235 = sigma_f_u235 + sigma_gamma_u235; % barns

nu_u238 = 2.2;                                  % neutrons/fission
sigma_gamma_u238 = 0.3;                         % barns
sigma_f_u238 = 0.05;                            % barns
sigma_tr_u238 = 2.5;                            % barns
sigma_a_u238 = sigma_f_u238 + sigma_gamma_u238; % barns

%% Cross sections
% Fuel regions
for i = 1:2
    % Macroscopic absorption x-s [1/cm]
    Sigma_a(i) = N_u * (x(i)*sigma_a_u235 + (1.-x(i))*sigma_a_u238);
    % Macroscopic transport x-s [1/cm]
    Sigma_tr(i) = N_u * (x(i)*sigma_tr_u235 + (1.-x(i))*sigma_tr_u238);
    % Macroscopic fission x-s [1/cm]
    Sigma_f(i) = N_u * (x(i)*sigma_f_u235 + (1.-x(i))*sigma_f_u238);
    % Fission neutron production rate per unit flux 
    % == neutrons per fission * macroscopic fission cross section [1/cm]
    NuSigma_f(i) = N_u * (x(i)*nu_u235*sigma_f_u235 + ... 
                                           (1.-x(i))*nu_u238*sigma_f_u238);  
    % Diffusion coefficient [cm]
    D(i) = 1./3./Sigma_tr(i);
    % Diffusion area [cm^2]
    LSquared(i) = D(i)/Sigma_a(i);
end
% Outer region
    % Macroscopic absorption x-s [1/cm]
    Sigma_a(3) = N_pb*sigma_gamma_pb;
    % Macroscopic transport x-s [1/cm]
    Sigma_tr(3) = N_pb*sigma_tr_pb;
    % Diffusion coefficient [cm]
    D(3) = 1./3./Sigma_tr(3);
    % Diffusion area [cm^2]
    LSquared(3) = D(3)/Sigma_a(3);

%% Mesh points
N(1) = round(R(1)/Delta);   
N(2) = round((R(2)-R(1))/Delta);
N(3) = round((R(3)-R(2))/Delta);
NTotal = N(1)+N(2)+N(3)+1;  % including the mesh point at the origin

%% Build the inner regions of matrix A
A = zeros(NTotal);  % Start with an empty matrix where everything is zero

% Using 'i' for region index and 'j' for mesh point

% Since the quantity D/Delta^2 appears frequently assign it a variable
for i = 1:3
    dd2(i) = D(i)/Delta/Delta;
end

% Internal mesh points in Region 1
for j = 2:N(1)
    A(j,j-1) = -dd2(1)*(2.*(j-1)-1)/(2.*(j-1));
    A(j,j)   = 2.*dd2(1) + Sigma_a(1);
    A(j,j+1) = -dd2(1)*(2.*(j-1)+1)/(2.*(j-1));
end

% Internal mesh points in Region 2
for j = N(1)+2:N(1)+N(2)
    A(j,j-1) = -dd2(2)*(2.*(j-1)-1)/(2.*(j-1));
    A(j,j)   = 2.*dd2(2) + Sigma_a(2);
    A(j,j+1) = -dd2(2)*(2.*(j-1)+1)/(2.*(j-1));
end

% Internal mesh points in Region 3
for j = N(1)+N(2)+2 : N(1)+N(2)+N(3)
    A(j,j-1) = -dd2(3)*(2.*(j-1)-1)/(2.*(j-1));
    A(j,j)   = 2.*dd2(3) + Sigma_a(3);
    A(j,j+1) = -dd2(3)*(2.*(j-1)+1)/(2.*(j-1));
end

%% Add boundary and interface conditions to matrix A
% 1.  Symmetry
A(1,1) = 1;
A(1,2) = -1;

% 2.  Flux Matching - satisfied automatically since the rightmost mesh
% point in region 1 is the same as the leftmost mesh point in region 2

% 3.  Current Matching - Regions 1 and 2
if N(2) ~= 0    % region 2 could be zero thickness, in which case skip it 
    A(N(1)+1,N(1)) = D(1);
    A(N(1)+1,N(1)+1) = -D(1)-D(2);
    A(N(1)+1,N(1)+2) = D(2);
end

% 4.  Flux Matching - satisfied automatically since the rightmost mesh
% point in region 2 is the same as the leftmost mesh point in region 3

% 5.  Current Matching - Regions 2 and 3
if N(3) ~= 0   % region 3 could be zero thickness, in which case skip it
    A(N(1)+N(2)+1,N(1)+N(2)) = D(2);
    A(N(1)+N(2)+1,N(1)+N(2)+1) = -D(2)-D(3);
    A(N(1)+N(2)+1,N(1)+N(2)+2) = D(3);
end

% 6.  Vacuum 
A(N(1)+N(2)+N(3)+1,N(1)+N(2)+N(3)+1)=1;

%% Build fission source term S = F*phi/k 
F = zeros(N(1)+N(2)+N(3)+1,1);

% Fission source coefficients for interior mesh points in region 1 
for j = 2:N(1)
    F(j) = NuSigma_f(1);
end

% Fission source coefficients for interior mesh points in region 2 
for j = N(1)+2:N(1)+N(2)
    F(j) = NuSigma_f(2);
end

% Leave region 3 terms at zero since there is no fission there

%% Iterative solution
% Flux profile initial guess: This is called phi(0) in the text.
phi_guess = zeros(N(1)+N(2)+N(3)+1,1);
% The final mesh point is left at zero below to satisfy the BCs
for j = 1:N(1)+N(2)+N(3)  
    phi_guess(j) = 1;
end

% Initial guess for the eigenvalue k:
k_guess = 1.0;

% Source term F*phi corresponding to the initial guess for phi:
S_guess = F.*phi_guess;

% Set up storage matrices for plotting of results later on
phi_stored = zeros(N(1)+N(2)+N(3)+1,MAX_ITS);
s_stored   = zeros(N(1)+N(2)+N(3)+1,MAX_ITS);
k_stored   = zeros(MAX_ITS,1);
% Will compute and store the fission rate Sigma_f*phi as well, even
% though it's not directly used in the finite difference solution.
fiss_rate_stored = zeros(N(1)+N(2)+N(3)+1,MAX_ITS);

phi_stored(:,1) = phi_guess;
s_stored(:,1)   = S_guess;
k_stored(1)     = k_guess;
for j = 1:N(1)+1
    fiss_rate_stored(j,1) = phi_stored(j,1)*Sigma_f(1);
end
for j = N(1)+2:N(1)+N(2)+1
    fiss_rate_stored(j,1) = phi_stored(j,1)*Sigma_f(2);
end
fprintf('Iteration: %g (initial guess). k: %f\n',1,k_guess)

%% Begin solution loops
for iteration = 2:MAX_ITS
% given the values phi_guess and k_guess from the previous iteration (or
% the initial guesses if this is the first iteration), solve for the next
% iteration of phi via: A*phi = F*phi_guess/k_guess (A*phi = S/k)
phi = (A\S_guess)./k_guess;

% Now we have a new and improved phi vector. Calculate the fission source
% term and eigenvalue arising from this flux:
S = F.*phi;

% Store the old and new neutron production rates
new_neut_prod_rate = 0.;
old_neut_prod_rate = 0.;

% Numerical integration, with dA = 2*pi*r dr.  Then int(dA NuSigma_f
% phi) over one mesh element (i.e. from (r_j-Delta/2) to (r_j+Delta/2))
% becomes 2*pi*NuSigma_f*phi(j)*j*Delta^2.

% the first mesh element is unique because it's located at r = 0.  So it
% only covers half the radius range of the others, i.e. we integrate from
% 0 to Delta/2:
old_neut_prod_rate = old_neut_prod_rate + 1.0/4.0*pi*NuSigma_f(1)* ...
                     phi_guess(1)*Delta*Delta;
new_neut_prod_rate = new_neut_prod_rate + 1.0/4.0*pi*NuSigma_f(1)* ...
                     phi(1)*Delta*Delta;

% Integrate over region 1 interior mesh points
for j = 2:N(1)
old_neut_prod_rate = old_neut_prod_rate + 2.0*pi*NuSigma_f(1)* ...
                     phi_guess(j)*(j-1)*Delta*Delta;
new_neut_prod_rate = new_neut_prod_rate + 2.0*pi*NuSigma_f(1)* ...
                     phi(j)*(j-1)*Delta*Delta;    
end

% The mesh point sitting on the boundary between region 1 and region 2 is
% associated with a mesh interval that extends a distance Delta/2 into
% region 1 and Delta/2 into region 2.  Account for this:

% Portion of mesh interval j that belongs to region 1:
j = N(1)+1;
old_neut_prod_rate = old_neut_prod_rate + pi*NuSigma_f(1)* ...
                     phi_guess(j)*Delta*Delta*((j-1)-1./4.);
new_neut_prod_rate = new_neut_prod_rate + pi*NuSigma_f(1)* ...
                     phi(j)*Delta*Delta*((j-1)-1./4.);  
% Portion of mesh interval j that belongs to region 2:
old_neut_prod_rate = old_neut_prod_rate + pi*NuSigma_f(2)* ...
                     phi_guess(j)*Delta*Delta*((j-1)+1./4.);
new_neut_prod_rate = new_neut_prod_rate + pi*NuSigma_f(2)* ...
                     phi(j)*Delta*Delta*((j-1)+1./4.);

% Interior region 2 mesh points
for j = N(1)+2:N(1)+N(2)
    old_neut_prod_rate = old_neut_prod_rate + 2.*pi*NuSigma_f(2)* ... 
                         phi_guess(j)*(j-1)*Delta*Delta;
    new_neut_prod_rate = new_neut_prod_rate + 2.*pi*NuSigma_f(2)* ...
                         phi(j)*(j-1)*Delta*Delta;    
end

% Portion of mesh interval j that belongs to region 2:
j = N(1)+N(2)+1;
old_neut_prod_rate = old_neut_prod_rate + pi*NuSigma_f(2)* ...
                     phi_guess(j)*Delta*Delta*((j-1)-1./4.);
new_neut_prod_rate = new_neut_prod_rate + pi*NuSigma_f(2)* ...
                     phi(j)*Delta*Delta*((j-1)-1./4.);

% Region 3 does not produce any neutrons

% By dividing new_neut_prod_rate by the production rate for the previous
% iteration (old_neut_prod_rate/k_guess), obtain a new estimate for k:
k = new_neut_prod_rate/old_neut_prod_rate*k_guess;

% store the new values of k, S, phi and fission rate:
phi_stored(:,iteration) = phi;
s_stored(:,iteration)   = S;
k_stored(iteration)     = k;
for j = 1:N(1)+1
    fiss_rate_stored(j,iteration) = phi_stored(j,iteration)*Sigma_f(1);
end
for j = N(1)+2:N(1)+N(2)+1
    fiss_rate_stored(j,iteration) = phi_stored(j,iteration)*Sigma_f(2);
end

% Check for convergence
converged = true;
for j = 1:length(S)
    if S_guess(j) > 0
        if (S(j)/S_guess(j) > 1.+TOLER) || (S(j)/S_guess(j) < 1.-TOLER)
            converged = false;
        end
    end
end
if k/k_guess > 1.+TOLER || k/k_guess < 1.-TOLER
    converged = false;
end

fprintf('Iteration: %g. k: %f.  k(%g)/k(%g): %f.\n',iteration,k, ...
        iteration,iteration-1,k/k_guess);

% if converged is still equal to true, break out of the loop:
if(converged)
    break;
end

% otherwise, continue. Use computed values of phi, s and k as starting
% points for the next iteration:
k_guess   = k;
phi_guess = phi;
S_guess   = S;
end

disp('Convergence achieved.');

%% Display some metrics about the case
% Show k
fprintf('\nFinal k:                   %f \n\n', k);
% Show flux profile flatness
flatness = ( max(phi_stored(2:N(1)-1,iteration)) - ...
             min(phi_stored(2:N(1)-1,iteration)) ) / ...
             max(phi_stored(2:N(1)-1,iteration));
fprintf('Region 1 flatness:         %.2f percent \n', 100*(1-flatness))

flatness = ( max(phi_stored(N(1)+1:N(1)+N(2)-1,iteration)) - ...
             min(phi_stored(N(1)+1:N(1)+N(2)-1,iteration)) ) / ...
             max(phi_stored(N(1)+1:N(1)+N(2)-1,iteration));
fprintf('Region 2 flatness:         %.2f percent \n', 100*(1-flatness))

flatness = ( max(phi_stored(2:N(1)+N(2)-1,iteration)) - ...
             min(phi_stored(2:N(1)+N(2)-1,iteration)) ) / ...
             max(phi_stored(2:N(1)+N(2)-1,iteration));
fprintf('Fuel flatness:             %.2f percent \n', 100*(1-flatness))

%% Plot results
% vector of mesh point radii (correct units). Used for plotting.
r = linspace(0,R(3),round(R(3)/Delta)+1);
color = linspace(1,0,iteration-1);

% Plot flux, including previous iterations
figure(2)
clf reset;
for i=1:iteration-1
    label = ['Iteration ',num2str(i)];
    plot(r,fiss_rate_stored(:,i),'DisplayName',label,'LineWidth',0.5, ...
         'Color',[0,color(iteration-i),color(i)]);
    hold on;
end
plot(r,fiss_rate_stored(:,iteration),'DisplayName','Final Rate', ...
     'LineWidth',2.0,'Color','r');
legend('show')
max_y=max(max(fiss_rate_stored));
set(gca,'XLim',[0 R(3)]);
set(gca,'YLim',[0 max_y*1.1]);
xlabel('Radial position r [cm]');
ylabel('Fission Rate Density \Sigma_{f}\phi [#/cm^{3}-s]');
h3 = plot([R(1) R(1)],[0,max_y*1.1],'Color','k');
h4 = plot([R(2) R(2)],[0,max_y*1.1],'Color','k');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
t1 = sprintf('Enrichment = %g%%',x(1)*100);
t2 = sprintf('Enrichment = %g%%',x(2)*100);
text(R(1)/100,max_y/15,'Region 1','Color','k');
text(R(1)*101/100,max_y/15,'Region 2','Color','k');
text(R(2)*101/100,max_y/15,'Region 3','Color','k');
text(R(1)/100,max_y/50,t1,'Color','k');
text(R(1)*101/100,max_y/50,t2,'Color','k');
text(R(2)*101/100,max_y/50,'Reflector','Color','k');

figure(1)
clf reset;
for i=1:iteration-1
    label = ['Iteration ',num2str(i)];
    plot(r,phi_stored(:,i),'DisplayName',label,'LineWidth',0.5, ...
         'Color',[0,color(iteration-i),color(i)]);
    hold on;
end
plot(r,phi_stored(:,iteration),'LineWidth',2.0,'Color','r', ...
     'DisplayName','Final Flux');
legend('show')
max_y=max(max(phi_stored));
set(gca,'XLim',[0 R(3)]);
set(gca,'YLim',[0 max_y*1.1]);
xlabel('Radial position r [cm]');
ylabel('Flux \phi [n/cm^{2}-s], arbitrary scaling');
h1 = plot([R(1) R(1)],[0,max_y*1.1],'Color','k');
h2 = plot([R(2) R(2)],[0,max_y*1.1],'Color','k');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
t1=sprintf('Enrichment = %g%%',x(1)*100);
t2=sprintf('Enrichment = %g%%',x(2)*100);
text(R(1)/100,max_y/15,'Region 1','Color','k');
text(R(1)*101/100,max_y/15,'Region 2','Color','k');
text(R(2)*101/100,max_y/15,'Region 3','Color','k');
text(R(1)/100,max_y/50,t1,'Color','k');
text(R(1)*101/100,max_y/50,t2,'Color','k');
text(R(2)*101/100,max_y/50,'Reflector','Color','k');

end