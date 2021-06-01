%ejercicio1
clc;clear all;
syms a t;
A=[0 0 -2; 0 1 0; 1 0 -3];
pol_car= det(A-a*eye(3));
auto_val=solve(pol_car)
V=[1 auto_val(1) auto_val(1)^2; 1 auto_val(2) auto_val(2)^2;1 auto_val(3) auto_val(3)^2]
V_inv= inv(V)
e=[exp(auto_val(1)*t);exp(auto_val(2)*t);exp(auto_val(3)*t)];
fun=V_inv*e
e_At=fun(1)*eye(3)+fun(2)*A+fun(3)*A^2
%% ejercicio 1 con sylvester
clc;clear all; close all;
syms a t I A1 exp_At;
A=[0 0 -2; 0 1 0; 1 0 -3];
auto_val=eig(A)
ex=[exp(auto_val(1)*t);exp(auto_val(2)*t);exp(auto_val(3)*t)]
a=[I A1 A1^2 exp_At]
S=[ones(1,3)',auto_val,auto_val.^2,ex;a]
de=0;
for i=1:4
    S_aux=S;
    ei4=S(i,4);
    S_aux(:,4)=[];
    S_aux(i,:)=[];
    f=((-1)^(i+1));
    de=de+f*ei4*det(S_aux);
end   
de
exp_At=solve(de,exp_At)





