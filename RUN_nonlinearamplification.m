% RUN_nonlinearamplification.m
%
% Computes the resonant frequencies and associated eigenmodes for an array
% of N bubbles graded in size with size factor s using the multipole
% method, then studies the response as a function of the input frequency,
% for different forcing amplitudes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B
%
% Used to create Fig 2 of:
% Modelling the active cochlea as a fully-coupled system of subwavelength
% Hopf resonators (2019) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

%% Define parameters

N = 6;                       % number of domain components / cell bundles
s = 1.05;
for i = 1:N
    R(i) = s^(i-1);
end

source = [-10, 0];             % location of the signal source

% Material parameters
% rho_b, kappa_b : cell bundles
% rho0, kappa0 : background
high=7000;

rho_b=1;
rho0=high*1;
kappa_b=1;
kappa0=high*1;

% High contrast parameter \delta
delta=rho_b/rho0;

cx = 2*R(1)*ones(1,N);
if N > 1
    for i = 2:N
        cx(i) = cx(i-1) + 2*R(i-1) + R(i);
    end
end
cy = zeros(1,N);

% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
% If we use higher order, then accuracy improves. Usually 3 is sufficient.
N_multi = 3;

%% Compute initial guesses for the resonances
%
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) min(eig((MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy))));

x = [linspace(0.0001, 0.005, 200), linspace(0.005, 0.025, 400)];
init = [];
for correction = [0 0.0000001i 0.000001i 0.0001i]
    y = zeros(1, length(x));
    for i = 1:length(x)
        y(i) = abs(f(x(i) - 0.00000001i/x(i) - correction));
    end
    for i = 2:length(x)-1
        if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e-8)
            init = [init x(i)];
        end
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
    z1 = initGuess - 0.00001i;
    z2 = initGuess - 0.00002i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e-7
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi \n'], real(res), imag(res))
       resonances = [resonances res];
       n = n + 1;
    end
end
% clf
% scatter(real(resonances), imag(resonances), 'x')


%% Computing eigenmodes over the plane

Vol = pi*R.^2;
N_res = length(resonances);

% Grid for field
gridN = max(N*10+1,50);         
gridMinX1 = floor(-10+R(1));
gridMaxX1 = ceil(cx(end)+R(end)+10);
gridMinX2 = -5;
gridMaxX2 = -gridMinX2;
g1 = linspace(gridMinX1, gridMaxX1, gridN);
g2 = linspace(gridMinX2, gridMaxX2, gridN);
Dx = (gridMaxX1 - gridMinX1)/(gridN - 1);
Dy = (gridMaxX2 - gridMinX2)/(gridN - 1);
[ g1, g2 ] = meshgrid(g1, g2);
gridPoints = [g1(:) g2(:)]';

gridPointsN = length(gridPoints);

u_store = zeros(gridPointsN,N);
ub_store = zeros(gridPointsN,N);
vals_store = zeros(N,N);
source_delta = zeros(N,1);

for m = 1:N_res
omega = resonances(m);
k = omega;                          % assumed v = 1
kb = omega;                         % assumed v_b = 1

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

% Calculate field at the source point
for i = 1:N
    fun2 = @(t) green(k, source(1)-cx(i)-R(i).*cos(t), source(2)-cy(i)-R(i).*sin(t)).*...
        (exp(-3i*t)*psi(1,i)+exp(-2i*t)*psi(2,i)+exp(-1i*t)*psi(3,i)+psi(4,i)+exp(1i*t)*psi(5,i)+exp(2i*t)*psi(6,i)+exp(3i*t)*psi(7,i));
    source_delta(m) = source_delta(m) + integral(fun2, -pi, pi);
end

% normalisation
avg_val = zeros(N,1);
total_number = 0;
for i = 1:N
    I = (u_b(:,2) == i);
    total_vals = sum(u_b(I,1));
    number = sum(I);
    total_number = total_number + number;
    avg_val(i) = total_vals/number;
end

norm_u = dot(Vol,avg_val.*conj(avg_val));
norm_u = sqrt(norm_u);

u_store(:,m) = u/norm_u;
ub_store(:,m) = u_b(:,1)/norm_u;
vals_store(:,m) = avg_val/norm_u;
end

%% Computing gamma

gamma = zeros(N);
for n = 1:N
    for m = 1:N
        gamma(n,m) = dot(u_store(:,n),u_store(:,m))*Dx*Dy;
    end
end

gamma1 = inv(gamma);


%% Response plot for k-th eigenmode

k = 2;
DA = Dx*Dy;
F_vals = 10.^[-10:2:-6];
alpha = 1e1;

for F = F_vals
    omeg_vals = linspace(0.92*real(resonances(k)),1.08*real(resonances(k)));
    soln = zeros(N,length(omeg_vals));
    for omeg = omeg_vals
        fun = @(x) fN(x,ub_store,gamma1,source_delta,resonances,omeg,F,DA,alpha);
        I = find(omeg_vals == omeg);
        options = optimset('TolFun',10^-15,'TolX',10^-10,'Display','off');
        soln(:,I) = fsolve(fun, ones(N,1), options);
    end
   
    y = log10(abs(soln(k,:)/F));
    plot(omeg_vals/real(resonances(k)),y,'linewidth',2)
    hold on
    xlim([min(omeg_vals)/real(resonances(k)),max(omeg_vals)/real(resonances(k))])
end

numerator = transpose(gamma1)*source_delta;
value = -numerator(k)./(resonances(k)^2-omeg_vals.^2);
plot(omeg_vals/real(resonances(k)), log10(abs(value)),':','linewidth',2)

set(gca,'TickLabelInterpreter','latex')
xlabel(['$\Omega/\Re\omega_',num2str(k),'$'],'interpreter','latex')
ylabel(['$\log_{10} X_{', num2str(k),'}/F$'],'interpreter','latex')
entries = [num2str(log10(F_vals)','$F=10^{%g}$')];
entries = [entries; '$\beta = 0 $'];
legend(entries,'location','northeast','interpreter','latex')

