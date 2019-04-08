function out = fN(x,ub_store,gamma1,delta,resonances,omeg,F,DA,alpha)

N = size(ub_store,2);

p = zeros(size(ub_store,1),1);
for n = 1:N
    p = p + x(n)*ub_store(:,n);
end 

nonlinearity = zeros(N,1);
for n = 1:N
    nonlinearity(n) = dot(p.*p.*conj(p),ub_store(:,n))*DA;
end

res = real(resonances) + 0.01*imag(resonances)*1i;

% out = (ones(N,1)-omeg^2./transpose(resonances).^2).*x*F + ...
%     F*transpose(gamma1)*delta + ...
%     alpha*transpose(gamma1)*nonlinearity*omeg^3*1i.*F^3;

out = (transpose(resonances).^2-omeg^2).*x + ...
    F*transpose(gamma1)*delta + ...
    alpha*transpose(gamma1)*nonlinearity*omeg^3*1i;


% --------------------

% p = zeros(size(ub_store,1),1);
% for n = 1:N
%     p = p + ub_store(:,n);
% end 
% 
% nonlinearity = zeros(N,1);
% for n = 1:N
%     nonlinearity(n) = dot(p.*p.*conj(p),ub_store(:,n))*DA;
% end
% 
% res = real(resonances) + 0.01*imag(resonances)*1i;
% 
% out = (ones(N,1)-omeg^2./transpose(resonances).^2).*x*F + ...
%     F*transpose(gamma1)*delta + ...
%     transpose(gamma1)*nonlinearity*omeg^3*1i.*abs(x).^2.*x*F^3;
% ------------------------------------
% for n = 1:N
%     vec(n,1) = F*delta(n) + 
%     Gamma(1,1,1,n)*abs(x(1))^2*x(1)...
%         + 2*Gamma(1,2,1,n)*abs(x(1))^2*x(2) + Gamma(2,2,1,n)*conj(x(1))*x(2)^2 ...
%         + Gamma(1,1,2,n)*x(1)^2*conj(x(2)) + 2*Gamma(1,2,2,n)*x(1)*abs(x(2))^2 ...
%         + Gamma(2,2,2,n)*abs(x(2))^2*x(2);
% end
% 
% out = (transpose(resonances).^2-omeg^2).*x + ...
%     transpose(gamma1)*vec;

end