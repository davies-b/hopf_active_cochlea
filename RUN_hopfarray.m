% RUN_    .m
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
% Used to create...
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

% High contrast parameter
delta=rho_b/rho0;

cx = linspace(0,L,N+2);
cx = cx(2:end-1);
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
scatter(real(resonances), imag(resonances), 'x')

%% Computing eigenmodes over the plane

N_res = length(resonances);

% Grid for field
gridN1 = max(N*5+1,50);         
gridN2 = 5;
gridMinX1 = -0.03+R(1);
gridMaxX1 = cx(end)+R(end)+0.03;
gridMinX1 = -0.005 - 0.0001;
gridMaxX1 = cx(end) + R(end);
gridMinX2 = -0.05;
gridMaxX2 = -gridMinX2;
g1 = [linspace(gridMinX1, gridMaxX1, gridN1), cx];
g2 = [linspace(gridMinX2, gridMaxX2, gridN2), 0];
Dx = (gridMaxX1 - gridMinX1)/(gridN1 - 1);
Dy = (gridMaxX2 - gridMinX2)/(gridN2 - 1);
DA = Dx*Dy;
[ g1, g2 ] = meshgrid(g1, g2);
gridPoints = [g1(:) g2(:)]';

gridPointsN = length(gridPoints);

us_store = zeros(gridPointsN,N_res);
ub_store = zeros(gridPointsN,N_res);

for m = 1:N_res
omega = resonances(m);
k = omega/v;                          
kb = omega/v_b;

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

u_s = zeros(gridPointsN, 1);
u_b = zeros(gridPointsN, 1);

parfor j = 1:gridPointsN
    gridPoint = gridPoints(:, j);
    
    % Determine whether we are inside or outside the domains
    I = (gridPoint(1)*ones(1,length(cx))-cx).^2 + (gridPoint(2)*ones(1,length(cy))).^2  <= R.^2 ;
    if sum( I ) > 0
        I = find(I);
        fun1 = @(t) green(kb, gridPoint(1)-cx(I)-R(I).*cos(t), gridPoint(2)-cy(I)-R(I).*sin(t)).*...
            (exp(-3i*t)*phi(1,I)+exp(-2i*t)*phi(2,I)+exp(-1i*t)*phi(3,I)+phi(4,I)+exp(1i*t)*phi(5,I)+exp(2i*t)*phi(6,I)+exp(3i*t)*phi(7,I));
        u_b(j) = integral(fun1, -pi, pi);
    else
        S = 0;
        for i = 1:N
            fun2 = @(t) green(k, gridPoint(1)-cx(i)-R(i).*cos(t), gridPoint(2)-cy(i)-R(i).*sin(t)).*...
                (exp(-3i*t)*psi(1,i)+exp(-2i*t)*psi(2,i)+exp(-1i*t)*psi(3,i)+psi(4,i)+exp(1i*t)*psi(5,i)+exp(2i*t)*psi(6,i)+exp(3i*t)*psi(7,i));
            S = S + integral(fun2, -pi, pi);
        end
        u_s(j) = S;
    end
end

% normalisation
u_norm = sqrt(dot(u_b,u_b)*DA);
u_b = u_b/u_norm;
u_s = u_s/u_norm;

us_store(:,m) = u_s;
ub_store(:,m) = u_b;

% % plotting
% uTotal = reshape(u, [gridN gridN]);
% hFig = figure(m+1);
% set(hFig, 'Position', [100 100 1200 900]);
% surf(g1, g2, real(uTotal), 'edgecolor', 'none'); xlabel('x_1'); ylabel('x_2'); title(['u_' num2str(m)])
% axis([gridMinX1, gridMaxX1, gridMinX2, gridMaxX2, min(real(uTotal(:))), max(real(uTotal(:))) ]); rotate3d on;
end


%% Computing gamma and F

gamma = zeros(N);
for n = 1:N
    for m = 1:N
        gamma(n,m) = dot(ub_store(:,n),ub_store(:,m))*Dx*Dy;
    end
end

force = zeros(N_res,1);
for n = 1:N
    force(n) = sum(us_store(:,n))*Dx*Dy + sum(ub_store(:,n))*Dx*Dy;
end




%% Non-linear amplification
beta = 1e5;
k = 11;

%%% Amplitude-frequence response curve
figure

volumes = [20 75 100];           % volumes in dB SPL
F_vals = db2pa(volumes);

res = resonances(k);
mu_c = -2*real(res)*imag(res)*gamma(k,k)/sqrt(real(res)^2-imag(res)^2);
omeg_vals = linspace(0.9999*real(resonances(k)),1.0001*real(resonances(k)),201);
p0 = db2pa(0);

for F = F_vals
soln = zeros(1,length(omeg_vals));
res = real(resonances(k)); tau = imag(resonances(k));
for n = 1:length(omeg_vals)
    omeg = omeg_vals(n);
    
    A = (res^2-tau^2-omeg^2);
    B = 2*res*tau+mu_c*omeg;
    C = beta*omeg^3;
    
    x = roots([C^2, 0, -2*B*C, 0, A^2+B^2, 0, -F^2*abs(force(k))^2]);
    x = x(abs(imag(x))<1e-15 & x>0);
    soln(n) = x;
end
subplot(2,3,[1,1.9,4,4.9])
if F == F_vals(1)
    semilogy(omeg_vals/sqrt(real(res)^2-imag(res)^2),soln/F,'k')
elseif F==F_vals(2)
    semilogy(omeg_vals/sqrt(real(res)^2-imag(res)^2),soln/F,'k-.')
elseif F == F_vals(3)
    semilogy(omeg_vals/sqrt(real(res)^2-imag(res)^2),soln/F,'k--')
end
hold on
xlim([omeg_vals(1)/sqrt(real(res)^2-imag(res)^2),omeg_vals(end)/sqrt(real(res)^2-imag(res)^2)])
end

omeg_vals = linspace(0.9*real(resonances(k)),1.1*real(resonances(k)),100);

delay = zeros(1,length(omeg_vals));
p0 = db2pa(0);

F = F_vals(1);
soln = zeros(1,length(omeg_vals));
res = real(resonances(k)); tau = imag(resonances(k));
for n = 1:length(omeg_vals)
    omeg = omeg_vals(n);
    
    A = (res^2-tau^2-omeg^2);
    B = 2*res*tau+mu_c*omeg;
    C = beta*omeg^3;
    
    x = roots([C^2, 0, -2*B*C, 0, A^2+B^2, 0, -F^2*abs(force(k))^2]);
    x = x(abs(imag(x))<1e-15 & x>0);
    delay(n) = angle(((res+1i*tau)^2-omeg^2)*x + i*omeg*mu_c*x - 1i*beta*omeg^3*x^3) - angle(force(k)) - pi;
end

subplot(2,3,3)
plot(omeg_vals/sqrt(real(res)^2-imag(res)^2),delay/pi/2,'k')
subplot(2,3,6)
group = -gradient(delay/pi/2);
plot(omeg_vals/sqrt(real(res)^2-imag(res)^2),group,'k')


subplot(2,3,[1,1.9,4,4.9])
set(gca,'TickLabelInterpreter','latex')
set(gca, 'OuterPosition', get(gca,'OuterPosition')+[0, 0.001, 0, 0])
xlabel(['$\Omega/\Omega_{',num2str(k),'}^c$'],'FontSize',12,'interpreter','latex')
ylabel(['Amplification $R_{', num2str(k),'}/F$'],'FontSize',12,'interpreter','latex')%,'rotation',0)
entries = num2str(volumes','$%g$ dB SPL');
legend(entries,'location','northeast','interpreter','latex')

subplot(2,3,3)
set(gca,'TickLabelInterpreter','latex')
xlabel(['$\Omega/\Omega_{',num2str(k),'}^c$'],'FontSize',11,'interpreter','latex')
ylabel('Phase delay $\psi/2\pi$','FontSize',11,'interpreter','latex')
ylim([-0.65,0])

subplot(2,3,6)
set(gca,'TickLabelInterpreter','latex')
xlabel(['$\Omega/\Omega_{',num2str(k),'}^c$'],'FontSize',11,'interpreter','latex')
ylabel('Group delay (s)','FontSize',11,'interpreter','latex')



%% Stability of solutions

res = resonances(k);
omeg = -sqrt(real(res)^2-imag(res)^2);
mu_c = 2*real(res)*imag(res)*gamma(k,k)/omeg;

i=1;
figure
for mu = [mu_c-100, mu_c+100]

Rr = linspace(5e-8,1e-6,8);
Ri = linspace(-2.7e-7,2.7e-7,10);

[Rr,Ri] = meshgrid(Rr,Ri);
Rp = zeros(size(Rr));

for n = 1:numel(Rr)
    R = Rr(n)+1i*Ri(n);
    Rp(n) = ((-omeg^2+res^2)*R-1i*mu*omeg*R+beta*1i*omeg^3*R^3)/(2*1i*omeg*-mu+beta*omeg^2*R^2);
end
subplot(1,2,i)
q = quiver(Rr, Ri, real(Rp), imag(Rp),'k');
Rc = sqrt(mu-mu_c)/abs(omeg)/sqrt(beta);

hold on
if mu > mu_c
    scatter(Rc,0,100,'xk')
else
    scatter(0,0,100,'xk')
end
hold off

grid on
set(gca,'TickLabelInterpreter','latex')
xlabel('Re$(R_k)$','FontSize',12,'interpreter','latex')
ylabel('Im$(R_k)$','FontSize',12,'interpreter','latex')
if mu > mu_c
    title('$\mu>\mu_k^c$','interpreter','latex')
else
    title('$\mu<\mu_k^c$','interpreter','latex')
end
i = i+1;
end


%% Backbone curve

omeg_vals = linspace(real(resonances(k)),1.2*real(resonances(k)),400);
F = 0;
soln = zeros(1,length(omeg_vals));
res = real(resonances(k)); tau = imag(resonances(k));
for n = 1:length(omeg_vals)
    omeg = omeg_vals(n);
    
    A = (res^2-tau^2-omeg^2);
    B = 2*res*tau+mu_c*omeg;
    C = beta*omeg^3;
    
    x = roots([C^2, 0, -2*B*C, 0, A^2+B^2]);
    soln(n) = abs(x(1));
end
figure
plot(omeg_vals/sqrt(real(res)^2-imag(res)^2),(soln)*10^6,'k')
xlim([0.8,1.2])
set(gca,'TickLabelInterpreter','latex')
xlabel(['$\Omega/\Omega_{',num2str(k),'}^c$'],'FontSize',12,'interpreter','latex')
ylabel(['$X_{',num2str(k),'}$ ($\mu$Pa)'],'FontSize',12,'interpreter','latex')



%% Fully-coupled modal system
beta = 1e1;
figure
F = 2;
gamma1 = gamma^-1;
numerator = force;

mu_c = zeros(N,1);
for n = 1:N
    mu_c(n) = -2*real(resonances(n))*imag(resonances(n))*gamma(n,n)/sqrt(real(resonances(n))^2-imag(resonances(n))^2);
end

mu = mu_c(11);

I = find(gridPoints(1,:) == cx(3) & gridPoints(2,:) == 0);
I = I(1);
U = ub_store(I,:);

omeg_vals = linspace(0,1.2e5,400);
soln = zeros(N,length(omeg_vals));
phase = zeros(1,length(omeg_vals));
for n = 1:length(omeg_vals)
    omeg = omeg_vals(n);
    fun = @(x) fNN(x,ub_store,gamma1,force,resonances,omeg,F,DA,mu,beta);
    I = find(omeg_vals == omeg);
    options = optimset('TolFun',10^-15,'TolX',10^-10,'Display','iter');
    init = -numerator./(transpose(real(resonances)).^2-omeg^2*ones(N_res,1));
    soln(:,I) = fsolve(fun, init, options);
    phase(n) = angle(U*soln(:,I));
end
   
    for k = 1:N
        y = abs(soln(k,:));
        plot(omeg_vals/2/pi*10^-3,pa2dbref(y),'k')
        hold on
    end

xlabel('Forcing frequency: $\Omega/2\pi$ (kHz)','interpreter','latex')
ylabel('Amplitude (dB)','interpreter','latex')
xlim([min(omeg_vals)/2/pi*10^-3,max(omeg_vals)/2/pi*10^-3])
box off

%% Phase delay curve

omeg_vals = linspace(0,2*pi*3.5*1e3,400);
soln = zeros(N,length(omeg_vals));
phase = zeros(1,length(omeg_vals));
for n = 1:length(omeg_vals)
    omeg = omeg_vals(n);
    fun = @(x) fNN(x,ub_store,gamma1,force,resonances,omeg,F,DA,mu,beta);
    I = find(omeg_vals == omeg);
    options = optimset('TolFun',10^-15,'TolX',10^-10,'Display','iter');
    init = -numerator./(transpose(real(resonances)).^2-omeg^2*ones(N_res,1));
    soln(:,I) = fsolve(fun, init, options);
    phase(n) = angle(U*soln(:,I));
end

cycles = phase;
for n = 2:length(omeg_vals)
    if phase(n-1) - phase(n) > 0.9*2*pi
        cycles(n:end) = cycles(n:end) + 2*pi;
    end
    if phase(n-1) - phase(n) < -0.9*2*pi
        cycles(n:end) = cycles(n:end) - 2*pi;
    end
end

figure
scatter(omeg_vals/2/pi*10^-3,-cycles/2/pi,1,'k')
xlabel('Forcing frequency $\Omega/2\pi$ (kHz)','interpreter','latex')
ylabel('Phase delay (cycles)','interpreter','latex')
xlim([0,max(omeg_vals)/2/pi*10^-3])
