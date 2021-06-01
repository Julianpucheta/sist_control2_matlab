%simbolico 
clc;clear all
syms u alfa alfa_p h h_p fi fi_p fi_pp a b c w
alfa_p = a*(fi-alfa);
fi_pp = -w^2*(fi - alfa - b*u);
h_p = c*alfa;

A = [[subs(subs(subs(subs(diff(alfa_p,alfa),'alfa',0),'fi',0),'fi_p',0),'h',0),...
     subs(subs(subs(subs(diff(alfa_p,fi),'alfa',0),'fi',0),'fi_p',0),'h',0),...
     subs(subs(subs(subs(diff(alfa_p,fi_p),'alfa',0),'fi',0),'fi_p',0),'h',0),... 
     subs(subs(subs(subs(diff(alfa_p,h),'alfa',0),'fi',0),'fi_p',0),'h',0)];
     [0 0 1 0];
     [subs(subs(subs(subs(diff(fi_pp,alfa),'alfa',0),'fi',0),'fi_p',0),'h',0),...
     subs(subs(subs(subs(diff(fi_pp,fi),'alfa',0),'fi',0),'fi_p',0),'h',0),...
     subs(subs(subs(subs(diff(fi_pp,fi_p),'alfa',0),'fi',0),'fi_p',0),'h',0),... 
     subs(subs(subs(subs(diff(fi_pp,h),'alfa',0),'fi',0),'fi_p',0),'h',0)];
     [subs(subs(subs(subs(diff(h_p,alfa),'alfa',0),'fi',0),'fi_p',0),'h',0),...
     subs(subs(subs(subs(diff(h_p,fi),'alfa',0),'fi',0),'fi_p',0),'h',0),...
     subs(subs(subs(subs(diff(h_p,fi_p),'alfa',0),'fi',0),'fi_p',0),'h',0),... 
     subs(subs(subs(subs(diff(h_p,h),'alfa',0),'fi',0),'fi_p',0),'h',0)];
    ];

B =[subs(subs(subs(subs(diff(alfa_p,u),'alfa',0),'fi',0),'fi_p',0),'h',0);
    0;
    subs(subs(subs(subs(diff(fi_p,u),'alfa',0),'fi',0),'fi_p',0),'h',0);
    subs(subs(subs(subs(diff(fi_pp,u),'alfa',0),'fi',0),'fi_p',0),'h',0);
    subs(subs(subs(subs(diff(h_p,u),'alfa',0),'fi',0),'fi_p',0),'h',0)];

disp('Matriz A :')
pretty(simplify(A))
disp('Reemplazando')
a=0.05;b=5;c=50;w=2;
A=[-a a 0 0;
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0]
disp('Matriz B :')
pretty(simplify(B))
disp('Reemplazando')
B=[0;
   0;
   b*w^2;
   0]
%% simulacion
%variables
clc;clear all;%close all;
T=150; At=1e-3; Kmax=T/At; t=linspace(0,T,Kmax);
a=0.05;b=5;c=50;w=2;
alfa_p=0;fi_p=0;fi_pp=0;h_p=0;
alfa=zeros(1,int32(Kmax));fi=zeros(1,int32(Kmax));fi_p=zeros(1,int32(Kmax));h=zeros(1,int32(Kmax));
u=linspace(0,0,int32(Kmax));
A=[-a a 0 0;
    0 0 1 0;
    w^2 -w^2 0 0;
    c 0 0 0];
B=[0;
   0;
   b*w^2;
   0];
% C=[1 0 0 0;
%    0 1 0 0;
%    0 0 0 1]
% D=0;
%condiciones iniciales
alfa(1)=0;fi(1)=0;fi_p(1)=0;h(1)=0;u(1)=0;
Xop=[0 0 0 0]';x=[alfa(1) fi(1) fi_p(1) h(1)]';
alfa_l(1)=0;fi_l(1)=0;fi_p_l(1)=0;h_l(1)=0;
ii = 0;
flag=0;
flag1=0;
for i=1:Kmax-1
     ii=ii+At;
     if(ii>=100)
         u(i)=0;
         flag1 = 1;
     end
     if(ii>=90 && ~flag1)
         u(i) = u(i-1) - 2*(1/Kmax);
%          if(u(i)<0)
%              u(i)=0;
%          end
         flag = 1;
     end
     if(ii>=10 && ~flag)
         u(i) = 0;
     else
         u(i+1)=u(i)+ 2*(1/Kmax);
     end
     
     alfa_p    = a*(fi(i) - alfa(i));
     fi_pp     = (-w^2)*(fi(i)-alfa(i)-(b*u(i)));
     h_p       = c*alfa(i);
     alfa(i+1) = alfa(i) + alfa_p*At;
     fi_p(i+1) = fi_p(i) + fi_pp*At;
     fi(i+1)   = fi(i) + fi_p(i)*At;
     h(i+1)    = h(i) + h_p*At;
     %variables de sistema lineal
     xp        =A*(x-Xop)+B*u(i);
     x         =x+xp*At;
     %Y=C*x;
     alfa_l(i+1)=x(1);
     fi_l(i+1)  =x(2);
     fi_p_l(i+1)=x(3);
     h_l(i+1)   =x(4);
end
figure(1)
subplot(3,1,1);hold on;
plot(t,alfa_l,'b--');title('\alpha_t');
subplot(3,1,2);hold on;
plot(t,fi_l,'b--');title('\phi_t');
% subplot(4,1,3);%hold on;
% plot(t,fi_p_l,'b');title('fi_p');
subplot(3,1,3);hold on;
plot(t,h_l,'b');title('altura (h)');
figure(2)
plot(t,u);title('u (timon de profundida)');
