% RUN_eigenmodes.m
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
% Used to create Figure 3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

%% Define parameters

L = 35e-3;              % length of cochlea
N = 22;                 % number of resonators
L_end = L;

s = 1.05;
a = 0.0001;

% Rad = L_end*(1-s)/(1-s^N)/3;
% for i = 1:N
%     R(i) = Rad*s^(i-1);
% end
for i = 1:N
    R(i) = a*s^(i-1);
end

% source = [-0.01, 0];             % location of the signal source

%%% Material parameters
rho0 = 1e3;                 % density of water
kappa0 = 2e9;               % bulk modulus of water
v = sqrt(kappa0/rho0);      % speed of sound in water

rho_b = 1.2;                % density of resonators/air
kappa_b = 1e5;              % bulk modulus of resonators/air
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air

% High contrast parameter \delta
delta=rho_b/rho0;

% cx = 2*R(1)*ones(1,N);
% if N > 1
%     for i = 2:N
%         cx(i) = cx(i-1) + 2*R(i-1) + R(i);
%     end
% end
cx = linspace(0,L,N+2);
cx = cx(2:end-1);
cy = zeros(1,N);

% Maximum order for multipole expansion (n = -N_multi, ..., -2, -1, 0, +1,
% +2, ..., +N_multi)
% If we use higher order, then accuracy improves. Usually 3 is sufficient.
N_multi = 3;

% %%% Plot the geometry
% figure, hold on
% t = linspace(0,2*pi);
% for n = 1:length(R)
%     plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
% end
% daspect([1 1 1])
% hold off
% 
% return

%% Compute initial guesses for the resonances
%
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) min(eig((MakeA(R,z,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy))));

x = linspace(1, 2*pi*22000, 200);
init = [];
% for correction = 0%[1i 10i 200i 300i]
    y = zeros(1, length(x));
    for i = 1:length(x)
        y(i) = abs(f(x(i)));
    end
    for i = 2:length(x)-1
        if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e0)
            init = [init x(i)];
        end
    end
% end

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
% clf
% scatter(real(resonances), imag(resonances), 'x')

% return
%% Computing eigenmodes along x2 axis

Vol = pi*R.^2;
N_res = length(resonances);

% Grid for field
gridN = max(N*10+1,50);         
gridMinX1 = -0.03+R(1);
gridMaxX1 = cx(end)+R(end)+0.03;
g1 = linspace(gridMinX1, gridMaxX1, gridN);
g2 = 0;
% Dx = (gridMaxX1 - gridMinX1)/(gridN - 1);
% Dy = (gridMaxX2 - gridMinX2)/(gridN - 1);
% DA = Dx*Dy;
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



% % normalisation
% avg_val = zeros(N,1);
% total_number = 0;
% for i = 1:N
%     I = (u_b(:,2) == i);
%     total_vals = sum(u_b(I,1));
%     number = sum(I);
%     total_number = total_number + number;
%     avg_val(i) = total_vals/number;
% end
% 
% norm_u = dot(Vol,avg_val.*conj(avg_val));
% norm_u = sqrt(norm_u);
% 
% u_store(:,m) = u/norm_u;
% ub_store(:,m) = u_b(:,1)/norm_u;
% vals_store(:,m) = avg_val/norm_u;

% plotting
% uTotal = reshape(u, [gridN gridN]);
% hFig = figure(m+1);
% set(hFig, 'Position', [100 100 1200 900]);
% surf(g1, g2, real(uTotal), 'edgecolor', 'none'); xlabel('x_1'); ylabel('x_2'); title(['u_' num2str(m)])
% axis([gridMinX1, gridMaxX1, gridMinX2, gridMaxX2, min(real(uTotal(:))), max(real(uTotal(:))) ]); rotate3d on;

% subplot(5,2,m)
% subplot(4,3,m)
subplot(8,3,m+2)
plot(g1,real(u),'k')
hold on
plot([g1(1),g1(end)],[0,0],'k:')
xlim([-0.01,0.045])
ylim(max(abs(ylim)).*[-1 1])
box off
set(gca,'visible','off')
% if m > 1
    text(0.036,0.8*max(abs(ylim)),[num2str(real(resonances(m)/2/pi),'%.0f'),' Hz'],'interpreter','latex')
% else
%     text(0.036,0.2*max(abs(ylim)),[num2str(real(resonances(m)/2/pi),'%.0f'),' Hz'],'interpreter','latex')
%     plot([0,0.035],[0.8*max(abs(ylim)),0.8*max(abs(ylim))],'k-+')
%     text(0.012,max(abs(ylim)),'3.5 cm','interpreter','latex')
% end
end


subplot(8,3,[1.2,2])
scatter(real(resonances)*10^-3/2/pi,imag(resonances)*10^-3/2/pi,'kx')
x = xlabel('Re $\omega/2\pi$ ($\times 10^3$)','interpreter','latex');
ylabel('Im $\omega/2\pi$ ($\times 10^3$)','interpreter','latex')
set(gca, 'XAxisLocation', 'top')
% set(x, 'position', get(x,'position')+[0,0,0]);
set(gca, 'ticklabelinterpreter','latex')


