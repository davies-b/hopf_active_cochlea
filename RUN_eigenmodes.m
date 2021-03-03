% RUN_eigenmodes.m
%
% Computes the resonant frequencies and associated eigenmodes for an array
% of N resonators graded in size with size factor s using the multipole
% method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B
%
% Used to create Figure 2 of Ammari & Davies, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

%% Define parameters

L = 35e-3;              % length of cochlea
N = 22;                 % number of resonators
L_end = L;

s = 1.05;
a = 0.0001;

for i = 1:N
    R(i) = a*s^(i-1);
end

%%% Material parameters
rho0 = 1e3;                 % density of water
kappa0 = 2e9;               % bulk modulus of water
v = sqrt(kappa0/rho0);      % speed of sound in water

rho_b = 1.2;                % density of resonators/air
kappa_b = 1e5;              % bulk modulus of resonators/air
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air

% High contrast parameter \delta
delta=rho_b/rho0;

% Define the position of the centres of the resonators
cx = linspace(0,L,N+2);
cx = cx(2:end-1);
cy = zeros(1,N);

% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
% If we use higher order, then accuracy improves. Usually 3 is sufficient.
N_multi = 3;

%%% Plot the geometry
figure, hold on
t = linspace(0,2*pi);
for n = 1:length(R)
    plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
end
daspect([1 1 1])
hold off


%% Compute initial guesses for the resonances
%
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) min(eig((MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy))));

x = linspace(1, 2*pi*22000, 200);
init = [];
y = zeros(1, length(x));
for i = 1:length(x)
    y(i) = abs(f(x(i)));
end
for i = 2:length(x)-1
    if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e0)
        init = [init x(i)];
    end
end

if length(init) < length(R)
    disp('WARNING: fewer than N initial guesses created')
end

init = sort(init);

%% Use Muller's method to find the resonances

distTol = 5e-5; fTol = 1e-5; iterMax = 10;
resonances = [];
n = 1;
for initGuess = init
        
    z0 = initGuess;
    z1 = initGuess - 1i;
    z2 = initGuess - 2i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e0
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi (%.0f Hz) \n'], real(res), imag(res), real(res)/2/pi)
       resonances = [resonances res];
       n = n + 1;
    end
end
% figure
% scatter(real(resonances), imag(resonances), 'x')

%% Computing eigenmodes along x2 axis
figure
Vol = pi*R.^2;
N_res = length(resonances);

% Grid for field
gridN = max(N*10+1,50);         
gridMinX1 = -0.03+R(1);
gridMinX1 = 0;
gridMaxX1 = cx(end)+R(end)+0.03;
g1 = linspace(gridMinX1, gridMaxX1, gridN);
g2 = 0;

[ g1, g2 ] = meshgrid(g1, g2);
gridPoints = [g1(:) g2(:)]';

gridPointsN = length(gridPoints);

u_store = zeros(gridPointsN,N);
ub_store = zeros(gridPointsN,N);
vals_store = zeros(N,N);
source_delta = zeros(N,1);

for m = 1:N_res
omega = resonances(m);
k = omega/v;                          % assumed v = 1
kb = omega/v_b;                         % assumed v_b = 1

A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);
[V, D] = eig(A);
[~,permutation] = sort(diag(D));
D=D(permutation,permutation);V = V(:,permutation);

phi = [];
psi = [];
N_terms = 2*N_multi+1;

for i = 1:N
    phi = [phi, V((i-1)*14+1:i*14-7,1)];
    psi = [psi, V((i-1)*14+8:i*14,1)];
end

% Calculate field
green = @(k,x,y) -1i/4*besselh(0,1,k*sqrt(x.^2+y.^2));

u = zeros(gridPointsN, 1);
u_b = zeros(gridPointsN, 2);

parfor j = 1:gridPointsN
    gridPoint = gridPoints(:, j);
    
    % Determine whether we are inside or outside the domains
    I = (gridPoint(1)*ones(1,length(cx))-cx).^2 + (gridPoint(2)*ones(1,length(cy))).^2  <= R.^2 ;
    if sum( I ) > 0
        S = 0;
        I = find(I);
        fun1 = @(t) green(kb, gridPoint(1)-cx(I)-R(I).*cos(t), gridPoint(2)-cy(I)-R(I).*sin(t)).*...
            (exp(-3i*t)*phi(1,I)+exp(-2i*t)*phi(2,I)+exp(-1i*t)*phi(3,I)+phi(4,I)+exp(1i*t)*phi(5,I)+exp(2i*t)*phi(6,I)+exp(3i*t)*phi(7,I));
        S = integral(fun1, -pi, pi);
        u(j) = S;
        u_b(j,:) = [S, I];     % stores the values that are in the bubbles
    else
        S = 0;
        for i = 1:N
            fun2 = @(t) green(k, gridPoint(1)-cx(i)-R(i).*cos(t), gridPoint(2)-cy(i)-R(i).*sin(t)).*...
                (exp(-3i*t)*psi(1,i)+exp(-2i*t)*psi(2,i)+exp(-1i*t)*psi(3,i)+psi(4,i)+exp(1i*t)*psi(5,i)+exp(2i*t)*psi(6,i)+exp(3i*t)*psi(7,i));
            S = S + integral(fun2, -pi, pi);
        end
        u(j) = S;
        index(j) = 1;
    end
end



subplot(8,3,m+2)
plot(g1,real(u),'k')
hold on
plot([g1(1),g1(end)],[0,0],'k:')
xlim([-0.01,0.045])
ylim(max(abs(ylim)).*[-1 1])
box off
set(gca,'visible','off')
text(0.038,0.5*max(abs(ylim)), [num2str(real(resonances(m)/2/pi),'%.0f'),' Hz'],'interpreter','latex')

end


subplot(10,3,[1.2,2])
scatter(real(resonances)*10^-3/2/pi,imag(resonances)*10^-3/2/pi,'kx')
x = xlabel('Re $\omega/2\pi$ ($\times 10^3$)','interpreter','latex');
ylabel('Im $\omega/2\pi$ ($\times 10^3$)','interpreter','latex')
set(gca, 'XAxisLocation', 'top')
set(gca, 'ticklabelinterpreter','latex')


%% Plot the membrane eigenmodes for comparison
figure

omeg_vals = [1000, 2000, 7000, 18000];
x = linspace(0,L,200);
for n = 1:length(omeg_vals)
    omeg = omeg_vals(n);
    p = membrane(x,omeg);
    subplot(2,ceil(length(omeg_vals)/2),n)
    plot(x,p,'k')
    hold on
    plot([g1(1),g1(end)],[0,0],'k:')
    xlim([0,L])
    ylim(max(abs(ylim)).*[-1 1])
    box off
    set(gca,'visible','off')
    text(0.028,0.5*max(abs(ylim)), [num2str(omeg),' Hz'],'interpreter','latex')
end

