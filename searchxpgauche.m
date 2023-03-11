function xp=searchxpgauche(v,mesh,xs,xmin)
%%% function to use the data and localize an approximation of x^star_k
%%% using data with a source on the right of the waveguide. v is the total
%%% wavefield, mesh its associated mesh, xs the position of the source and
%%% xmin the minimum of the search interval 
    %%%% search interval 
    z=linspace(max(xs-20,xmin),xs,300); 
    %%% measurements at the bottom 
    tempz=P1togrid(mesh,v,z,0); 
    %first approximation of the 3 parameters to minimize 
    [~,b]=max(abs(tempz));
    xapprox=z(b);
    aap=tempz(b)./0.5357;
    %restriction of the interval around x^star_k
    z=linspace(max(xapprox-2,xmin),min(xapprox+2,xs),300); 
    tempz=P1togrid(mesh,v,z,0); 
    %minimization 
    g=@(y) log(norm((y(1)+1i*y(4))*airy(y(2)*(y(3)-z))-tempz));
    a=fminsearch(g,[real(aap),2,xapprox,imag(aap)]); 
    for j=1:10
        a=fminsearch(g,a); 
    end
    xp=a(3);
end