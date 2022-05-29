clear
close all

%% grid 100
load('LW100.mat');
dLW=d;
uLW=u;
pLW=p;
eLW=e;
load('VL100.mat');
dVL=d;
uVL=u;
pVL=p;
eVL=e;
load('SW100.mat');
dSW=d;
uSW=u;
pSW=p;
eSW=e;
load('SWWENO100.mat');
dSWWENO=d;
uSWWENO=u;
pSWWENO=p;
eSWWENO=e;
load('RP100.mat');
dRP=d;
uRP=u;
pRP=p;
eRP=e;
load('exact.mat')

% density
figure;
plot(s,d,'--k','LineWidth',0.8);hold on
plot(x,dLW,'LineWidth',0.8);hold on
plot(x,dSW,'LineWidth',0.8);hold on
plot(x,dVL,'LineWidth',0.8);hold on
plot(x,dRP,'LineWidth',0.8);hold on
plot(x,dSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$\rho$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off')
grid on

% pressure
figure;
plot(s,p,'--k','LineWidth',0.8);hold on
plot(x,pLW,'LineWidth',0.8);hold on
plot(x,pSW,'LineWidth',0.8);hold on
plot(x,pVL,'LineWidth',0.8);hold on
plot(x,pRP,'LineWidth',0.8);hold on
plot(x,pSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$p$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off')
grid on

% velocity
figure;
plot(s,u,'--k','LineWidth',0.8);hold on
plot(x,uLW,'LineWidth',0.8);hold on
plot(x,uSW,'LineWidth',0.8);hold on
plot(x,uVL,'LineWidth',0.8);hold on
plot(x,uRP,'LineWidth',0.8);hold on
plot(x,uSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off','Location','northwest')
grid on

% internal energy
figure;
plot(s,e,'--k','LineWidth',0.8);hold on
plot(x,eLW,'LineWidth',0.8);hold on
plot(x,eSW,'LineWidth',0.8);hold on
plot(x,eVL,'LineWidth',0.8);hold on
plot(x,eRP,'LineWidth',0.8);hold on
plot(x,eSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$e$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off','Location','northwest')
grid on

close all
%% grid 400
load('LW400.mat');
dLW=d;
uLW=u;
pLW=p;
eLW=e;
load('VL400.mat');
dVL=d;
uVL=u;
pVL=p;
eVL=e;
load('SW400.mat');
dSW=d;
uSW=u;
pSW=p;
eSW=e;
load('SWWENO400.mat');
dSWWENO=d;
uSWWENO=u;
pSWWENO=p;
eSWWENO=e;
load('RP400.mat');
dRP=d;
uRP=u;
pRP=p;
eRP=e;
load('exact.mat')

% density
figure;
plot(s,d,'--k','LineWidth',0.8);hold on
plot(x,dLW,'LineWidth',0.8);hold on
plot(x,dSW,'LineWidth',0.8);hold on
plot(x,dVL,'LineWidth',0.8);hold on
plot(x,dRP,'LineWidth',0.8);hold on
plot(x,dSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$\rho$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off')
grid on

% pressure
figure;
plot(s,p,'--k','LineWidth',0.8);hold on
plot(x,pLW,'LineWidth',0.8);hold on
plot(x,pSW,'LineWidth',0.8);hold on
plot(x,pVL,'LineWidth',0.8);hold on
plot(x,pRP,'LineWidth',0.8);hold on
plot(x,pSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$p$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off')
grid on

% velocity
figure;
plot(s,u,'--k','LineWidth',0.8);hold on
plot(x,uLW,'LineWidth',0.8);hold on
plot(x,uSW,'LineWidth',0.8);hold on
plot(x,uVL,'LineWidth',0.8);hold on
plot(x,uRP,'LineWidth',0.8);hold on
plot(x,uSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off','Location','northwest')
grid on

% internal energy
figure;
plot(s,e,'--k','LineWidth',0.8);hold on
plot(x,eLW,'LineWidth',0.8);hold on
plot(x,eSW,'LineWidth',0.8);hold on
plot(x,eVL,'LineWidth',0.8);hold on
plot(x,eRP,'LineWidth',0.8);hold on
plot(x,eSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$e$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off','Location','northwest')
grid on

close all
%% grid 200
load('LW200.mat');
dLW=d;
uLW=u;
pLW=p;
eLW=e;
load('VL200.mat');
dVL=d;
uVL=u;
pVL=p;
eVL=e;
load('SW200.mat');
dSW=d;
uSW=u;
pSW=p;
eSW=e;
load('SWWENO200.mat');
dSWWENO=d;
uSWWENO=u;
pSWWENO=p;
eSWWENO=e;
load('RP200.mat');
dRP=d;
uRP=u;
pRP=p;
eRP=e;
load('exact.mat')

% density
figure;
plot(s,d,'--k','LineWidth',0.8);hold on
plot(x,dLW,'LineWidth',0.8);hold on
plot(x,dSW,'LineWidth',0.8);hold on
plot(x,dVL,'LineWidth',0.8);hold on
plot(x,dRP,'LineWidth',0.8);hold on
plot(x,dSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$\rho$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off')
grid on

% pressure
figure;
plot(s,p,'--k','LineWidth',0.8);hold on
plot(x,pLW,'LineWidth',0.8);hold on
plot(x,pSW,'LineWidth',0.8);hold on
plot(x,pVL,'LineWidth',0.8);hold on
plot(x,pRP,'LineWidth',0.8);hold on
plot(x,pSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$p$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off')
grid on

% velocity
figure;
plot(s,u,'--k','LineWidth',0.8);hold on
plot(x,uLW,'LineWidth',0.8);hold on
plot(x,uSW,'LineWidth',0.8);hold on
plot(x,uVL,'LineWidth',0.8);hold on
plot(x,uRP,'LineWidth',0.8);hold on
plot(x,uSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$u$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off','Location','northwest')
grid on

% internal energy
figure;
plot(s,e,'--k','LineWidth',0.8);hold on
plot(x,eLW,'LineWidth',0.8);hold on
plot(x,eSW,'LineWidth',0.8);hold on
plot(x,eVL,'LineWidth',0.8);hold on
plot(x,eRP,'LineWidth',0.8);hold on
plot(x,eSWWENO,'LineWidth',0.8);hold on
xlim([-0.5,0.5])
xticks(-0.5:0.1:0.5)
xlabel('$x$','Interpreter','latex','FontSize',16)
ylabel('$e$','Interpreter','latex','FontSize',16)
legend('exact','Lax-Wendroff','Steger-Warming','van Leer','Roe-Pike','WENO5-SW',...
    'FontSize',13,'FontName','Times New Roman','box','off','Location','northwest')
grid on