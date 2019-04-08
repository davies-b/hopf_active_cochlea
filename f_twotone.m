function out = f_twotone(x,ub_store,gamma,delta,resonances,omeg1,omeg2,F1,F2,DA,alpha)

N = size(ub_store,2);
S10 = 0;
S01 = 0;
S21 = 0;
S12 = 0;
for n = 1:N
    S10 = S10 + omeg1*x(n)*ub_store(:,n);
    S01 = S01 + omeg2*x(N+n)*ub_store(:,n);
    S21 = S21 + (2*omeg1-omeg2)*x(2*N+n)*ub_store(:,n);
    S12 = S12 + (-omeg1+2*omeg2)*x(3*N+n)*ub_store(:,n);
end

C10 = S10.*abs(S10).^2 + 2.*S10.*(abs(S01).^2+abs(S21).^2+abs(S12).^2) ...
    + S01.^2.*conj(S12) + 2.*S01.*S21.*conj(S10) + 2.*S21.*S12.*conj(S01);

C01 = S01.*abs(S01).^2 + 2.*S01.*(abs(S10).^2+abs(S21).^2+abs(S12).^2) ...
    + S10.^2.*conj(S21) + 2.*S10.*S12.*conj(S01) + 2.*S12.*S21.*conj(S10);

C21 = S21.*abs(S21).^2 + 2.*S21.*(abs(S10).^2+abs(S01).^2+abs(S12).^2) ...
    + S10.^2.*conj(S01) + 2.*S10.*S01.*conj(S12);

C12 = S12.*abs(S12).^2 + 2.*S12.*(abs(S10).^2+abs(S01).^2+abs(S21).^2) ...
    + S01.^2.*conj(S10) + 2.*S10.*S01.*conj(S21);


nonlinearity = zeros(N,4);
for n = 1:N
    nonlinearity(n,1) = dot(C10,ub_store(:,n))*DA;
    nonlinearity(n,2) = dot(C01,ub_store(:,n))*DA;
    nonlinearity(n,3) = dot(C21,ub_store(:,n))*DA;
    nonlinearity(n,4) = dot(C12,ub_store(:,n))*DA;
end

out1 = transpose(gamma)*((transpose(resonances).^2-omeg1^2).*x(1:N)) + ...
    F1*delta - alpha*1i*nonlinearity(:,1);
out2 = transpose(gamma)*((transpose(resonances).^2-omeg2^2).*x(N+1:2*N)) + ...
    F2*delta - alpha*1i*nonlinearity(:,2);
out3 = transpose(gamma)*((transpose(resonances).^2-(2*omeg1-omeg2)^2).*x(1:N)) ...
    - alpha*1i*nonlinearity(:,3);
out4 = transpose(gamma)*((transpose(resonances).^2-(-omeg1+2*omeg2)^2).*x(1:N)) ...
    - alpha*1i*nonlinearity(:,4);

out = [out1; out2; out3; out4];

end