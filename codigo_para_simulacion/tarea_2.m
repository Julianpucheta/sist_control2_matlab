%-----autor alaniz franco
clc;clear all;
syms av;
dim=1000;
t=linspace(0,1,dim);
A=[0 0 -2; 0 1 0; 1 0 -3];
pol_car= det(A-av*eye(3));
val=[-2;-1;1]
V=[1 val(1) val(1)^2; 1 val(2) val(2)^2;1 val(3) val(3)^2]
V_inv= inv(V)
e=[exp(val(1)*t);exp(val(2)*t);exp(val(3)*t)];
fun=V_inv*e;
a0=fun(1,:);
a1=fun(2,:);
a2=fun(3,:);
for i=1:dim
    e_At=fun(1,i)*eye(3)+fun(2,i)*A+fun(3,i)*A^2;
    e_At1(1,i)=e_At(1,1);
    e_At2(1,i)=e_At(1,2);
    e_At3(1,i)=e_At(1,3);
    e_At4(1,i)=e_At(2,1);
    e_At5(1,i)=e_At(2,2);
    e_At6(1,i)=e_At(2,3);
    e_At7(1,i)=e_At(3,1);
    e_At8(1,i)=e_At(3,2);
    e_At9(1,i)=e_At(3,3);
end
figure(2);plot(t,e_At1,t,e_At3,t,e_At5,t,e_At7,t,e_At9)
%% ejercicio1 con sylvester
clc;clear all;close all;
dim=1000;
syms t;
A=[0 0 -2; 0 1 0; 1 0 -3]
I=eye(3)
exp_At=(exp(t)*(A^2 + 3*A + 2*I))/6 - (exp(-t)*(3*A^2 + 3*A - 6*I))/6 - (exp(-2*t)*(- 2*A^2 + 2*I))/6
exp_At';
exp_At=ans(:)
t=linspace(0,1,dim); subs(exp_At,'t',t);exp_At=double(ans);
figure(1);plot(t,exp_At,t,exp_At,t,exp_At,t,exp_At,t,exp_At);


 