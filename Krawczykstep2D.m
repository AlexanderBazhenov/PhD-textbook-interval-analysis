% Krawczykstep2D
% 2023-02-04
L = [2/(3*X1^(1/3)) 2/(3*X2^(1/3)); -1/(X1+1)^2 -1];
Lambda=inv(2*mid(L));
%
LF = [infsup(0.0,0.0) infsup(0.0,0.0);
infsup(0.0,0.0) infsup(0.0,0.0)]
%
d1 = [infsup(0.0,0.0), infsup(0.0,0.0)]
d2 = [infsup(0.0,0.0), infsup(0.0,0.0)]
%
d1(1) = 1-0.691/(X1^(1/3))+0.755/((X1+1)^2)
d1(2) = 0.755-0.691/(X2^(1/3))
d2(1) = 0.175/(X1^(1/3))-0.691/((X1+1)^2)
d2(2) = 0.309+0.175/(X2^(1/3))
%
LF = [ d1; d2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            X0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baumann
##% cut function
##FTI Examples 1D
##midFdiff = mid(Fdiff)
##radFdiff = rad(Fdiff)
##midradFdiff = midFdiff / radFdiff
##p =cut( midradFdiff, infsup(-1,1))
##% points
##clow = mid(x) - p*rad(x)
##cupp = mid(x) + p*rad(x)
midrad11 = mid(d1(1)) / rad(d1(1))
midrad12 = mid(d1(2)) / rad(d1(2))
midrad21 = mid(d2(1)) / rad(d2(1))
midrad22 = mid(d2(2)) / rad(d2(2))
%
p1(1) = cut( midrad11, infsup(-1,1))
p1(2) = cut( midrad12, infsup(-1,1))
p2(1) = cut( midrad21, infsup(-1,1))
p2(2) = cut( midrad22, infsup(-1,1))
% cut
p1 = [1; 0]
c_lo1 (1) = mid (X1) - p1(1)*rad(X1)
c_lo1 (2) = mid (X2) - p1(2)*rad(X2)
c_up1 (1) = mid (X1) + p1(1)*rad(X1)
c_up1 (2) = mid (X2) + p1(2)*rad(X2)
p2 = [0; 1]
c_lo2 (1) = mid (X1) - p2(1)*rad(X1)
c_lo2 (2) = mid (X2) - p2(2)*rad(X2)
c_up2 (1) = mid (X1) + p2(1)*rad(X1)
c_up2(2) = mid (X2) + p2(2)*rad(X2)
%
xlim([0.5 2.0])
ylim([0.5 1.5])
% Points
plot(mid (X1), mid (X2), 'sk')
plot(c_lo1 (1),c_lo1(2), 'sb')
plot(c_up1(1), c_up1(2), 'sb')
plot(c_lo2 (1),c_lo2(2), 'sr')
plot(c_up2(1), c_up2(2), 'sr')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c1_low
FF1_lo = [infsup(0.0,0.0), infsup(0.0,0.0)]
y1=c_lo1 (1); y2=c_lo1 (2);
y = [ y1; y2];
F(1) = y1^(2/3)+y2^(2/3)-infsup(1.8,2.0) ;
F(2) = 1/(y1+1)-y2+0.3;
FF1_lo = y-Lambda*F';
FF1_lo = FF1_lo + d1*(X0-y)
% c1_up
FF1_up = [infsup(0.0,0.0), infsup(0.0,0.0)]
y1=c_up1 (1); y2=c_up1 (2);
y = [ y1; y2];
F(1) = y1^(2/3)+y2^(2/3)-infsup(1.8,2.0) ;
F(2) = 1/(y1+1)-y2+0.3;
FF1_up = y-Lambda*F';
FF1_up = FF1_up + d1*(X0-y)
%
FF1 = intersect(FF1_lo,FF1_up)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c2_low
FF2_lo = [infsup(0.0,0.0), infsup(0.0,0.0)]
y1=c_lo2 (1); y2=c_lo2 (2);
y = [ y1; y2];
F(1) = y1^(2/3)+y2^(2/3)-infsup(1.8,2.0) ;
F(2) = 1/(y1+1)-y2+0.3;
FF2_lo = y-Lambda*F';
FF2_lo = FF2_lo + d2*(X0-y)
% c2_up
FF2_up = [infsup(0.0,0.0), infsup(0.0,0.0)]
y1=c_up2 (1); y2=c_up2 (2);
y = [ y1; y2];
F(1) = y1^(2/3)+y2^(2/3)-infsup(1.8,2.0) ;
F(2) = 1/(y1+1)-y2+0.3;
FF2_up = y-Lambda*F';
FF2_up = FF2_up + d2*(X0-y)
%
FF2 = intersect(FF2_lo,FF2_up)
%
Kbic = [ FF1(1) ; FF2(2) ]

