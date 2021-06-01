clc;clear all;
syms M m dpp dp d fipp fip fi l F u g 
%M=0.5;m=.1;F=.1;l=.6;g=9.8;
ang_inicial=0;
disp('Equilibrio iniestable')
ec1=(M+m)*dpp+m*l*fipp*cos(fi)-m*l*fip^2*sin(fi)+F*dp==u;
ec2 = l*fipp - g*sin(fi) + dpp*cos(fi) == 0;
%para el equilibrio inestable fi=0; ==> cos(0)=1, sin(0)=fi
%reemplaznado
ec1_aux= subs(subs(ec1,cos(fi),1),sin(fi),fi);
ec2_aux= subs(subs(ec2,cos(fi),1),sin(fi),fi);
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
         [0 0 1 0];
         [subs(subs(subs(subs(diff(fipp_aux,d),'d',0),dp,0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,dp),'d',0),'dp',0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,fi),'d',0),'dp',0),'fi',ang_inicial),'fip',0),... 
         subs(subs(subs(subs(diff(fipp_aux,fip),'d',0),'dp',0),'fi',ang_inicial),'fip',0)]];
 mat_B=[ 0;
        subs(subs(subs(subs(diff(dpp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0);
        0;
        subs(subs(subs(subs(diff(fipp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0);
        0];
 pretty(simplify(mat_A))
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
         [0 0 1 0];
         [subs(subs(subs(subs(diff(fipp_aux,d),'d',0),dp,0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,dp),'d',0),'dp',0),'fi',ang_inicial),'fip',0),...
         subs(subs(subs(subs(diff(fipp_aux,fi),'d',0),'dp',0),'fi',ang_inicial),'fip',0),... 
         subs(subs(subs(subs(diff(fipp_aux,fip),'d',0),'dp',0),'fi',ang_inicial),'fip',0)]];
 mat_B=[ 0;
        subs(subs(subs(subs(diff(dpp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0);
        0;
        subs(subs(subs(subs(diff(fipp_aux,u),'d',0),dp,0),'fi',ang_inicial),'fip',0);
        0];
 pretty(simplify(mat_A))
 pretty(simplify(mat_B))
