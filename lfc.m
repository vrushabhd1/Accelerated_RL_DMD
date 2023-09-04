clc; clear all; close all;
warning('off','all')
warning

%Laplacian matrix
L=24;
s = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 17 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 18 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24 24];
t = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 20 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 24 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
b=length(t);
G = digraph(s,t);
Adj = adjacency(G);
Dia = diag(sum(Adj,2));
Lap = Dia - Adj;

%Global continuous-time system dynamics
Tp=20;Tg=0.08;Tt=0.3;Kp=120;R=2.4;Tij=0.015;
A1=[-1/Tp Kp/Tp 0 -Kp/Tp;0 -1/Tt 1/Tt 0;-1/(R*Tg) 0 -1/Tg 0;2*pi*Tij*(L-1) 0 0 0];
A2=[0 0 0 0;0 0 0 0;0 0 0 0;-2*pi*Tij 0 0 0];
% A1=[-1/Tp Kp/Tp 0 -Kp/Tp;0 -1/Tt 1/Tt 0;-1/(R*Tg) 0 -1/Tg 0;0 0 0 0];
% A2=[0 0 0 0;0 0 0 0;0 0 0 0;2*pi*Tij*(L-1) 0 0 0];
A1=kron(eye(L),A1);
A2=kron(Lap,A2);
A2=full(A2);
A=A1+A2;
% A2=kron(eye(L),A2);
% A=A2;
B=[0 0 1/Tg 0]';
B=kron(eye(L),B);
C=[1 zeros(3,1)';zeros(4,1)'];
C=kron(eye(L),C);
D=0;
sysc=ss(A,B,C,D);
%%Global discrete-time system dynamics
sysd = c2d(sysc,1,'zoh');
[A,B,C,D]=ssdata(sysd);
% DS=tf(sysd);
% bode(DS);
% [DS_s,DS_f] = freqsep(DS,4);
% bodeplot(DS_s,DS_f)
% legend('original','slow','fast')

%Graphical-LQR
Q1=[2 0 0 0;0 2 0 0;0 0 2 0;0 0 0 2];
Q1bar=kron(eye(L),Q1);
Q2=[0.01 0 0 0;0 0.01 0 0;0 0 0.01 0;0 0 0 0.01];
Q2bar=kron(Lap,Q2);
Qbar=Q1bar+Q2bar;
Rd=1;
R=kron(eye(L),Rd);
Pstar=dare(A,B,Qbar,R);
Kstar=inv(R+B'*Pstar*B)*B'*Pstar*A;
try chol(Pstar);
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

%initial conditions
x1=[1;-1;1;-1];
x(:,1) = repmat(x1,L,1);
desired=[];
zbar=[];
dd=[];
zbard=[];
K=[1 1 1 1];Kd=K;
K= repmat(K,L,1);
K= repmat(K',L,1);
K=K';
k=0;z=0;iter=15;iterr=800;


% Data collection original x
   tic
for i=1:iterr
    a1=0.001;
    a2=0.001;
    a3=0.001;
    a4=0.001;
    u(:,i)=-Kstar*x(:,i)+a1*sin(9.8*i)+a2*cos(10.2*i)^2+a3*sin(10*i)+a4*cos(10*i);
    x(:,i+1)=A*x(:,i)+B*u(:,i);
    
% Data matrices
    Hxx=kron(x(:,i)',x(:,i)')-kron(x(:,i+1)',x(:,i+1)');
    Hxu=2*kron(x(:,i)',(u(:,i)+K*x(:,i))');
    Huu=-kron((-u(:,i)+K*x(:,i))',(u(:,i)+K*x(:,i))');
   
%     desired1=x(:,i)'*Qbar*x(:,i)+x(:,i)'*K'*R*K*x(:,i);
%     desired=[desired;desired1];
%     J(:,i)=x(:,i)'*Qbar*x(:,i)+x(:,i)'*K'*R*K*x(:,i); 
%     zbar1=[Hxx Hxu Huu]';
%     zbar=[zbar;zbar1'];
     
end
toc

% %DMD
    tic
for j=1:iter
    thresh = 1e-10;
%     thresh = 0.2;
    [Phi,Phio,lambda,zeta, U,S,V,U_r,S_r,V_r,W_r,W,Atilde1,Atilde,r] = DMD(x,thresh);

% Data matrices reduced x
     Hzxx=Hxx*(kron(Phi',Phi'))';
     Hzxu=Hxu*kron(Phi',eye(L))';
     Hzuu=Huu*kron(eye(1),eye(1))';
     Qd=(U_r)'*Qbar*U_r;
     Ad=conj(U_r)'*A*U_r;
     Bd=conj(U_r)'*B;
%      Pstard=dare(Ad,Bd,Qd,Rd);
     Pstard =conj(U_r)'* Pstar*U_r;
%    Kstard=inv(Rd+Bd'*Pstard*Bd)*Bd'*Pstard*Ad;
     Kstard=(pinv(U_r)* Kstar')';
     
     dd1=zeta(:,j)'*Qd*zeta(:,j)+zeta(:,j)'*Kd'*Rd*Kd*zeta(:,j);
    dd=[dd dd1];
    J1(:,j)=zeta(:,j)'*Qd*zeta(:,j)+zeta(:,j)'*Kd'*Rd*Kd*zeta(:,j);   
    zbard1=[Hzxx Hzxu Hzuu]';
    zbard=[zbard;zbard1'];
        
        z=z+1;
        md=zbard'*zbard;qd=zbard'*dd;
        rank(md);
        Hd=(pinv(md))*qd;
          Pd=[Hd(1,1) Hd(5,1) Hd(9,1) Hd(13,1);Hd(2,1) Hd(6,1) Hd(10,1) Hd(14,1);Hd(3,1) Hd(7,1) Hd(11,1) Hd(15,1);Hd(4,1) Hd(8,1) Hd(12,1) Hd(16,1)];
        Kd=inv(Rd+Hd(97,1))*([Hd(17,1) Hd(37,1) Hd(57,1) Hd(77,1)]);  
        Kd= repmat(Kd,L,1);
        ekd(z,:)=norm(Kd-Kstard);
       ekkd(z,:)=norm(Kd);
       epd(z,:)=norm(Pd-Pstard);
%         epd(z,:)=norm(Pd);
        dd=[];
        zbard=[];
        K=Kd*(U_r)';
        X = (U_r)*zeta;
        Y=C*X;
 end
toc
figure(1)
p=plot(G,'Layout','subspace');
p.Marker = '|';
p.NodeColor = 'b';
p.EdgeColor = 'k';
p.MarkerSize = 15;
hold on;

figure(2)
subplot(2,1,1);
y = 1:96;
semilogy(y,diag(S),'x--','Linewidth',1.5);
ylabel('Singular values ($\sigma$)');
 set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
 set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(2,1,2);
yr = 1:4;
semilogy(yr,diag(S_r),'x--','Linewidth',1.5);
ylabel('Singular values ($\sigma$)');xlabel('Number of components');
set(get(gca,'XLabel'),'FontSize',12); set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
 set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');

t1=1:iter;
 figure(3)
 subplot(2,1,1);
 plot(t1,ekd(1:iter,:),'r','Linewidth',2);
 hold on
 plot(t1,ekd(1:iter,:),'o','Linewidth',2);
 ylabel('$||K_r-K_r^*||$');
  set(get(gca,'YLabel'),'FontSize',14);set(get(gca,'legend'),'FontSize',10);
 set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
 subplot(2,1,2);
 plot(t1,epd(1:iter,:),'o','Linewidth',2);
 hold on 
 plot(t1,epd(1:iter,:),'b','Linewidth',2);
 ylabel('$||P_r-P_r^*||$');xlabel('Iteration $j$');
 set(get(gca,'XLabel'),'FontSize',14);set(get(gca,'YLabel'),'FontSize',14);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');

tt=1:15;
figure(4)
subplot(4,1,1);
plot(tt,X(1:4:96,tt),'r','Linewidth',1);
ylabel('$\Delta f$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,2);
plot(tt,X(2:4:96,tt),'b','Linewidth',1);
ylabel('$\Delta P_m$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,3);
plot(tt,X(3:4:96,tt),'k','Linewidth',1);
ylabel('$\Delta P_g$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,4);
plot(tt,X(4:4:96,tt),'Linewidth',1);
ylabel('$\Delta P_{tie}$');xlabel('Iteration $j$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');


figure(5);
subplot(4,1,1); plot(real(U(:, 1:5)),'Linewidth',1.5); 
ylabel('$U$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,2); plot(real(S(:, 1:5)),'Linewidth',1.5);
ylabel('$S$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,3); plot(real(V(:, 1:5)),'Linewidth',1.5);
ylabel('$V$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,4); plot(real(Phio),'Linewidth',1.5);
xlabel('Number of components');ylabel('$\bar{\theta}$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');

figure(6);
subplot(4,1,1); plot(real(U_r(:, 1:4)),'Linewidth',1.5); 
ylabel('$\tilde{U}$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,2); plot(real(S_r(:, 1:4)),'Linewidth',1.5);
ylabel('$\tilde{S}$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,3); plot(real(V_r(:, 1:4)),'Linewidth',1.5);
ylabel('$\tilde{V}$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
subplot(4,1,4); plot(real(Phi),'Linewidth',1.5);
xlabel('Number of components');ylabel('$\theta$');
set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');

figure(12)
plot(tt,X(:,tt),'Linewidth',1.5);
hold on
xlabel('Iteration $j$');legend('$\Delta f$','$\Delta P_m$','$\Delta P_g$','$\Delta P_{tie}$');
set(get(gca,'XLabel'),'FontSize',14);set(get(gca,'YLabel'),'FontSize',14);set(get(gca,'legend'),'FontSize',10);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
% figure(7)
% plot(tt,zeta(:,tt),'Linewidth',1.5);
% hold on
% xlabel('Iteration $j$');ylabel('$\xi$');
% set(get(gca,'XLabel'),'FontSize',14);set(get(gca,'YLabel'),'FontSize',14);set(get(gca,'legend'),'FontSize',10);
% set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
figure(8)
surf(real(Phi*Atilde*Phi'));
figure(9)
surf(U_r*Atilde*U_r');
figure(10)
surf(Atilde1);
figure(11)
plot(diag(S)/sum(diag(S)),'ro'); title('low rank property for DMD (rank = 4, four modes)');

% plot(tt,Y(:,tt),'Linewidth',2);

% tt1=1:200;
% figure(12)
% subplot(2,2,1);
% plot(tt1,x(1:4:96,tt1),'r','Linewidth',1);
% ylabel('$\Delta f$');xlabel('Iteration $j$');
% set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
% set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
% subplot(2,2,2);
% plot(tt1,x(2:4:96,tt1),'b','Linewidth',1);
% ylabel('$\Delta P_m$');xlabel('Iteration $j$');
% set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
% set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
% subplot(2,2,3);
% plot(tt1,x(3:4:96,tt1),'k','Linewidth',1);
% ylabel('$\Delta P_g$');xlabel('Iteration $j$');
% set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
% set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');set(get(gca,'legend'),'Interpreter','latex');
% subplot(2,2,4);
% plot(tt1,x(4:4:96,tt1),'Linewidth',1);
% ylabel('$\Delta P_{tie}$');xlabel('Iteration $j$');
% set(get(gca,'XLabel'),'FontSize',12);set(get(gca,'YLabel'),'FontSize',12);set(get(gca,'legend'),'FontSize',10);
figure(66)
subplot(1,2,1), semilogy(diag(S),'k')
subplot(1,2,2), plot(cumsum(diag(S))/sum(diag(S)),'k')