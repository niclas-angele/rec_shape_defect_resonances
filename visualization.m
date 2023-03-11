clear; close all ; 

%%%Shape of the defect 

%%%%%%%%%Linear function
% h=@(x) (0.01/30*x+0.1).*(x>=-4).*(x<=4)+(0.1+0.01/30*4)*(x>4)+(0.1-0.01/30*4)*(x<-4);
% nodes2 = [-8 0; -4 0; 4 0; 12  0;12 (0.1+0.01/30*4); 4 (0.1+0.01/30*4); -4 (0.1-0.01/30*4); -8 (0.1-0.01/30*4)];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','0.01/30*t+0.1',-4,4}, {7,8},{8,1}};

%%%%%%%%% x^3 function (1 inflexion point) 
% h=@(x) (0.1+x.^3/4^2*0.01/30).*(x>=-4).*(x<=4)+(0.1+0.01/30*4)*(x>4)+(0.1-0.01/30*4)*(x<-4);
% nodes2 = [-15 0; -4 0; 4 0; 15  0;15 (0.1+0.01/30*4); 4 (0.1+0.01/30*4); -4 (0.1-0.01/30*4); -15 (0.1-0.01/30*4)];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','0.1+t.^3/4^2*0.01/30',-4,4}, {7,8},{8,1}};

%%%%%%%%%% sqrt(x) function (infinite derivates)
h=@(x) (0.1-0.01*4/30+sqrt(x+4)*4*0.01/30/sqrt(2)).*(x>=-4).*(x<=4)+(0.1+0.01/30*4)*(x>4)+(0.1-0.01/30*4)*(x<-4);
nodes2 = [-15 0; -4 0; 4 0; 15  0;15 (0.1+0.01/30*4); 4 (0.1+0.01/30*4); -4 (0.1-0.01/30*4); -15 (0.1-0.01/30*4)];
edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','0.1-0.01*4/30+sqrt(t+4)*4*0.01/30/sqrt(2)',-4,4}, {7,8},{8,1}};


%%%%%%%%%%% regular function with 1 inflexion point
% hprim=@(x) 0.00005*(x.^4-8*x.^3+16*x.^2).*(x>=0).*(x<4)-0.00005*(-x.^4-8*x.^3-16*x.^2).*(x<0).*(x>-4);
% aa=0.00005*(4.^5/5-2*4.^4+16*4.^3/3);
% bb=-0.00005*(4.^5/5-2*4.^4+16*4.^3/3); 
% h=@(x) 0.1+0.00005*(x.^5/5-2*x.^4+16*x.^3/3).*(x>=0).*(x<=4)- 0.00005*(-x.^5/5-2*x.^4-16*x.^3/3).*(x<0).*(x>=-4)+(bb)*(x<-4)+(aa)*(x>4);
% nodes2 = [-15 0; -4 0; 4 0; 15  0;15 0.1+aa; 4 0.1+aa; -4 0.1+bb; -15 0.1+bb];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','0.1+0.00005*(t.^5/5-2*t.^4+16*t.^3/3).*(x>=0)- 0.00005*(-t.^5/5-2*t.^4-16*t.^3/3).*(x<0)',-4,4}, {7,8},{8,1}};


%%%%%%%%%%% regular function without inflexion point
% bb=0.000003*((-4).^5/5-32.*(-4).^3/3+16^2.*(-4));
% aa=0.000003*((4).^5/5-32.*(4).^3/3+16^2.*(4));
% h=@(x) (0.1+0.000003*(x.^5/5-32.*x.^3/3+16^2.*x)).*(x>-4).*(x<4)+(x<-4).*(0.1+bb)+(x>4).*(0.1+aa);
% nodes2 = [-15 0; -4 0; 4 0; 15  0;15 0.1+aa; 4 0.1+aa; -4 0.1+bb; -15 0.1+bb];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7,'t','(0.1+0.000003*(t.^5/5-32.*t.^3/3+16^2.*t))',-4,4}, {7,8},{8,1}};

%%%%%%%%%%%% without defect 
% h=@(x) 0.1+0*x;
% nodes2 = [-15 0; -4 0; 4 0; 15  0;15 0.1; 4 0.1; -4 0.1; -15 0.1];
% edges2 = {{1,2}, {2,3}, {3,4},{4,5}, {5,6},{6,7}, {7,8},{8,1}};


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
xs=6;
sig=0.005;
delta=@(x) 1/sig/sqrt(2*pi)*exp(-(x-xs).^2/2/sig/sig);
%frequency:
k=31.2;
%mode number:
n=1; 
%computation of x^star_k:
xp=fzero(@(x) h(x)^2-pi^2/k^2,0.1); 
%wavenumber function 
kn=@(x) sqrt(k^2-pi^2./(h(x).^2));


%%%% computation of the wavefield 
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

%%%%%%%visualization
x=linspace(-6,6,500); 
y=linspace(0,0.12,200);
[X,Y]=meshgrid(x,y); 
P=P1togrid(mesh2,v,x,y);
surf(X,Y,abs(P),'edgecolor','none')
view(0,90)
xline(xp,'r')
colorbar

