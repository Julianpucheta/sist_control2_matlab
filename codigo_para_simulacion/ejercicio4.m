clc;clear all;
syms M m dpp dp d fipp fip fi L F u g 
%M=0.5;m=.1;F=.1;l=.6;g=9.8;
ang_inicial=0;
disp('Equilibrio iniestable')
ec1=(M+m)*dpp+m*L*fipp*cos(fi)-m*L*fip^2*sin(fi)+F*dp==u;
ec2 = L*fipp - g*sin(fi) + dpp*cos(fi) == 0;
%para el equilibrio inestable fi=0; ==> cos(0)=1, sin(0)=fi
%reemplaznado
ec1_aux= subs(subs(ec1,cos(fi),round(cos(ang_inicial))),sin(fi),fi);
ec2_aux= subs(subs(ec2,cos(fi),round(cos(ang_inicial))),sin(fi),fi);
%encontrar fipp y dpp
fipp_aux=solve(ec2_aux,fipp);
ec1_aux=subs(ec1_aux,'fipp',fipp_aux);
dpp_aux=solve(ec1_aux,dpp);
fipp_aux=subs(fipp_aux,'dpp',dpp_aux);
disp('dpp es igual a :');pretty(simplify(dpp_aux));
disp('fipp es igual a :');pretty(simplify(fipp_aux));
%linealizando para le equilibro inestable
%haciendo taylor para cada ecuacion y evaluando para el punto de operacion Xe=[0 0 0 0]
mat_A=[ [0 1 0 0];
        [subs(subs(subs(subs(diff(dpp_aux,d),'d',0),dp,0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(dpp_aux,dp),'d',0),'dp',0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(dpp_aux,fi),'d',0),'dp',0),'fi',ang_inicial),'fip',0),... 
         subs(subs(subs(subs(diff(dpp_aux,fip),'d',0),'dp',0),'fi',ang_inicial),'fip',0)];
         [0 0 0 1];
         [subs(subs(subs(subs(diff(fipp_aux,d),'d',0),dp,0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,dp),'d',0),'dp',0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,fi),'d',0),'dp',0),'fi',ang_inicial),'fip',0),... 
         subs(subs(subs(subs(diff(fipp_aux,fip),'d',0),'dp',0),'fi',ang_inicial),'fip',0)]];
 mat_B=[ 0;
        subs(subs(subs(subs(diff(dpp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0);
        0;
        subs(subs(subs(subs(diff(fipp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0)];
 disp('Matriz A:')
 pretty(simplify(mat_A))
 disp('Matriz B:')
 pretty(simplify(mat_B))
 %-------------
disp('equilibrio estable')
ang_inicial=pi;
%fi=pi ==> cos(fi)=-1 y sin(fi)=-fi
ec1_aux=subs(subs(ec1,cos(fi),-1),sin(fi),-fi);
ec2_aux=subs(subs(ec2,cos(fi),-1),sin(fi),-fi);
%encontramos dpp y fpp
fipp_aux=solve(ec2_aux,fipp);
ec1_aux=subs(ec1_aux,'fipp',fipp_aux);
dpp_aux=solve(ec1_aux,dpp);
fipp_aux=subs(fipp_aux,'dpp',dpp_aux);
disp('dpp es igual a :');pretty(simplify(dpp_aux))
disp('fipp es igual a :');pretty(simplify(fipp_aux))
%linealizando para le equilibro estable
%haciendo taylor para cada ecuacion y evaluando para el punto de operacion Xe=[0 0 pi 0]
mat_A=[ [0 1 0 0];
        [subs(subs(subs(subs(diff(dpp_aux,d),'d',0),dp,0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(dpp_aux,dp),'d',0),'dp',0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(dpp_aux,fi),'d',0),'dp',0),'fi',ang_inicial),'fip',0),... 
         subs(subs(subs(subs(diff(dpp_aux,fip),'d',0),'dp',0),'fi',ang_inicial),'fip',0)];
         [0 0 0 1];
         [subs(subs(subs(subs(diff(fipp_aux,d),'d',0),dp,0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,dp),'d',0),'dp',0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,fi),'d',0),'dp',0),'fi',ang_inicial),'fip',0),... 
         subs(subs(subs(subs(diff(fipp_aux,fip),'d',0),'dp',0),'fi',ang_inicial),'fip',0)]];
 mat_B=[ 0;
        subs(subs(subs(subs(diff(dpp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0);
        0;
        subs(subs(subs(subs(diff(fipp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0)];
 disp('Matriz A:')
 pretty(simplify(mat_A))
 disp('Matriz B:')
 pretty(simplify(mat_B))
 
 %%
clc;clear all;
%variables
T=10; At=1e-4; Kmax=T/At; t=linspace(0,T,Kmax);
m=0.1;F=0.1;long=1.2;g=9.8;M=0.5;color='k';angulo_inicial=pi;
d=zeros(1,int32(Kmax));dp=zeros(1,int32(Kmax));fi=zeros(1,int32(Kmax));fip=zeros(1,int32(Kmax));
u=linspace(0,0,int32(Kmax));

%equilibrio inestable
A=[ 0 1 0 0;
    0 -F/M -(g*m)/M 0;
    0 0 0 1;
    0 F/(M*long) g*(M+m)/(M*long) 0]
B=[0;
   1/M;
   0;
   -1/(M*long)]

%equilibrio estable
A1=[ 0 1 0 0;
     0 -F/M -(g*m)/M 0;
     0 0 0 1;
     0 -F/(M*long) -g*(M+m)/(M*long) 0]
B1=[ 0;
     1/M;
     0;
     1/(M*long)]
 
%condiciones iniciales
d(1)=0;dp(1)=0;
dp(1)=0;fi(1)=angulo_inicial;fip(1)=0;dpp=0;fipp=0;
u(1)=0;
Xop=[0 0 angulo_inicial 0]';x=[d(1) dp(1) fi(1) fip(1)]';
dl(1)=0;dpl(1)=0;fil(1)=0;fipl(1)=0;
for i=1:Kmax-1
     dpp      = (long*m*sin(fi(i))*fip(i)^2+ u(i) - F * dp(i) - fipp * long * m * cos(fi(i)))/(M+m);
     fipp     = (-dpp * cos(fi(i))+g*sin(fi(i)))*(1/long);
     dp(i+1)  = dp(i)   + dpp *At;
     d(i+1)   = d(i)    + dp(i)*At;
     fip(i+1) = fip(i)  + fipp*At;
     fi(i+1)  = fi(i)   + fip(i)*At;
     %----------variables de sistema lineal------------
      xp        = A1*(x-Xop)+B1*u(i);
      x         = x+xp*At;
      dl(i)     = x(1);
      dpl(i)    = x(2);
      fil(i)    = x(3);
      fipl(i)   = x(4);
     %Y=C*x;
end
dl(i+1)= x(1);dpl(i+1)= x(2);fil(i+1)= x(3);fipl(i+1)= x(4);
figure(1)
subplot(3,2,1);
plot(t,fi,color);hold on;%plot(t,pi*ones(size(t)),'k');hold on;
plot(t,fil,'r--');hold on;legend('No lineal','Lineal');
grid on;title('Angulo, \phi');
subplot(3,2,3);
plot(t,fip,color);hold on;%
plot(t,fipl,'r--');hold on;
grid on;title('Velocida angular, \omega_t');
subplot(3,2,2); 
plot(t,d,color);hold on;%
plot(t,dl,'r--');hold on;
grid on;title('Posición carro, \delta');
subplot(3,2,4);
plot(t,dp,color);hold on;%
plot(t,dpl,'r--');hold on;
grid on;title('Velocidad de carro');
subplot(3,1,3);
plot(t,u,color);
grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;