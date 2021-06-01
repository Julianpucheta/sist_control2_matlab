function [K,P] = lqrPropio(A,B,Q,R,obs)%obs indica si quiere el controlador para el observador
H = [A -B*inv(R)*B'; 
    -Q         -A']; %Hamiltoniano

[V,D]=eig(H);
M_PM=[];
%armado de [ M
%           PM ] para el calculo de la matiz P
for ii=1:length(eig(H))
    if real(D(ii,ii))<0
        M_PM = [M_PM V(:,ii)];
    end
end
[row,colums] = size(M_PM);

M_ = M_PM(1:row/2,:); % M
PM = M_PM((row/2)+1:end,:);% P*M
P = real(PM * inv(M_));
if ~obs
    K = inv(R)*B'*P;
else
    K = (inv(R)*B'*P)';
end
end