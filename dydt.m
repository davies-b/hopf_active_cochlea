function out = dydt(t,y,gamma1,resonances,ub_store,mu,beta,DA)

res = diag(resonances.^2);
N = size(ub_store,2);

p = zeros(size(ub_store,1),1);
for n = 1:N
    p = p + y(N+n)*ub_store(:,n);
end 

nonlinearity = zeros(N,1);
for n = 1:N
    nonlinearity(n) = dot(ub_store(:,n),p.*p.*conj(p))*DA;
end

nonlinearity = transpose(gamma1)*nonlinearity;

out = zeros(2*N,1);
out(1:N) = y(N+1:2*N); 
out(N+1:2*N) = -res*y(1:N)+mu*y(N+1:2*N)-beta*nonlinearity;

end