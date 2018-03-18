function PlotSolidCircle(xc, yc, r )
% plot a solid circle with radius r

x=(-r:0.01:r);
y1=sqrt(r^2-x.^2);
y2=-y1;
x = x+xc;
y1 = y1+yc;
y2 = y2+yc;
patch([x x((2*r/0.01+1):-1:1)],[y1 y2((2*r/0.01+1):-1:1)],'k')



end

