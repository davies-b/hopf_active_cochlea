function out = fNN(x,ub_store,gamma1,force,resonances,omeg,F,DA,mu,alpha)

N = size(ub_store,2);

p = zeros(size(ub_store,1),1);
for n = 1:N
    p = p + x(n)*ub_store(:,n);
end 

nonlinearity = zeros(N,1);
for n = 1:N
    nonlinearity(n) = dot(ub_store(:,n),p.*p.*conj(p))*DA;
end

res = real(resonances);% + 0.01*imag(resonances)*1i;

% out = (ones(N,1)-omeg^2./transpose(resonances).^2).*x*F + ...
%     F*transpose(gamma1)*delta + ...
%     alpha*transpose(gamma1)*nonlinearity*omeg^3*1i.*F^3;

out = (transpose(resonances).^2-omeg^2).*x + ...
    mu*1i*omeg*x - F*force - ...
    alpha*transpose(gamma1)*nonlinearity*omeg^3*1i;

end