
clear; close all ; 

%%%Shape of the defect 

%%%%%% bump
% h=@(x) 0.0002*(x>-2).*(x<2).*(x-2).^2.*(x+2).^2+0.1;
% nodes2 = [-15 0; -2 0; 2 0;15  0;15 0.1;2 0.1; -2 0.1; -15 0.1];
% edges2 = {{1,2}, {2,3}, {3,4}, {4,5}, {5,6}, {6,7,'4*t-2','0.0002*(4*t-4).^2.*(4*t).^2+0.1',0,1}, {7,8}, {8,1},{2,7,'lr'},{3,6,'lr'}};

%%%%% linear bump 
% h=@(x) 0.1+(x>=-5).*(x<=0).*(x+5)*0.0025/5+(x>0).*(x<=4).*(-x+4)*0.0025/4;
% hprim=@(x) (x>=-5).*(x<=0).*0.0025/5-(x>0).*(x<=4).*0.0025/4;
% nodes2 = [-15 0; -5 0; 4 0; 15  0;15 0.1; 4 0.1; -5 0.1; -15 0.1];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','0.1+(t>=-5).*(t<=0).*(t+5)*0.0025/5+(t>0).*(t<=4).*(-t+4)*0.0025/4',-5,4}, {7,8},{8,1}};

%%%%% sin bump 
h=@(x) 0.1+0.0025*(x>=-5).*(x<=5).*sin(pi/10*(x+5)); 
hprim=@(x) 0.0025*pi/10*(x>=-5).*(x<=5).*cos(pi/10*(x+5)); 
nodes2 = [-15 0; -5 0; 5 0; 15  0;15 0.1; 5 0.1; -5 0.1; -15 0.1];
edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','0.1+0.0025*sin(pi/10*(t+5))',-5,5}, {7,8},{8,1}};

%%%% hollow 
% h=@(x) 0.1+(x>=-5).*(x<=0).*(-x-5)*0.0025/5+(x>0).*(x<=4).*(x-4)*0.0025/4;
% nodes2 = [-15 0; -5 0; 4 0; 15  0;15 0.1; 4 0.1; -5 0.1; -15 0.1];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','0.1+(t>=-5).*(t<=0).*(-t-5)*0.0025/5+(t>0).*(t<=4).*(t-4)*0.0025/4',-5,4}, {7,8},{8,1}};

%%%%%% mix hollow bump 
% aa=0.0008*sin(4*pi*sqrt(-3.5+5)/sqrt(9));
% bb=-0.0008*sin(4*pi/sqrt(9));
% h=@(x) 0.1+0.0008*sin(4*pi*sqrt(x+5)/sqrt(9)).*(x>=-3.5).*(x<=4)+bb+aa*(x<-3.5);
% nodes2 = [-15 0; -3.5 0; 4 0; 15  0;15 h(10); 4 h(10); -3.5 h(-10); -15 h(-10)];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t',' 0.1+0.0008*sin(4*pi*sqrt(x+5)/sqrt(9))-0.0008*sin(4*pi/sqrt(9))',-3.5,4}, {7,8},{8,1}};


%%%%%% construction of the waveguide 
dx = 0.2; 
domain2 = Domain(nodes2,edges2);
mesh2 = Mesh(domain2,dx);
mesh2=mesh2.submesh;
mesh2=mesh2.submesh;
mesh2=mesh2.submesh;
mesh2=mesh2.submesh;
%%%% for more precision, refine the mesh 
% mesh2=mesh2.submesh;
% mesh2=mesh2.submesh;
% mesh2=mesh2.submesh;


%%%%%% parameters 
%source=dirac source in xs:
xs=0;
sig=0.005;
delta=@(x) 1/sig/sqrt(2*pi)*exp(-(x-xs).^2/2/sig/sig);
%frequency ensemble:
OM=linspace(30.65,31.4,20);
%reconstruction of width h: 
XPG=OM*0; 
XPD=OM*0; 
%mode number:
n=1; 


for i=1:length(OM)
    k=OM(i); 
    %computation of x^star_k:
    xp=fzero(@(x) h(x)^2-pi^2/k^2,0.1); 
    %wavenumber function 
    kn=@(x) sqrt(k^2-pi^2./(h(x).^2));
    %PML parametrization:
    l=10; 
    at=8; 
    alpha=@(x,y) -1i*k./(-1i*k+at*(x-l)).*(x>=l)-1i*k./(-1i*k+at*(-x-l)).*(x<=-l)+1*(x<l).*(x>-l);
    C=mesh2.P0({@(x,y) alpha(x,y), 0, 0 1}); 
    %boundary conditions:
    bc1={[4,8],0,1,0};
    bc2={[1,2,3],1,0,0}; 
    %without boundary sources:
    %     bc3={[5,6,7],1,0,0};
    %with boundary sources:
    bc3={[5,6,7],1,0,@(x,y) delta(x)}; 
    %solution with internal sources:
    %     v = mesh2.solve(C,-k.^2,mesh2.P0(@(x,y) delta(x).*cos(n*pi*y./h(x))),{bc1,bc2,bc3});
    %solution with boundary sources:
    v = mesh2.solve(C,-k.^2,0,{bc1,bc2,bc3});
    %     temp=0*x; 
    XPD(i)=searchxpdroite(v,mesh2,xs,8);
    XPG(i)=searchxpgauche(v,mesh2,xs,-8);
end

%%%%visualization
x=linspace(-6,6);
plot(x,h(x),'k');
hold on 
plot([-6,-5,XPG(end:-1:1),0], [h(-5),h(-5),n*pi./OM(end:-1:1),h(0)])
plot([0,XPD,5,6], [h(0),n*pi./OM,h(5),h(5)])
