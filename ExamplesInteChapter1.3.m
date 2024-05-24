% 2022-06-27
% Examples Interval Analysis Chapter 1.3
% 2024-05-17
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExamplePlace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     START    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pkg list
pkg load interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.3.1
% inf sup
disp('Example 1.3.1')
X = infsup(0,2)
X = infsup(1,2)
XD = X^2 - 2*X + 1
F1 = acos((X^2 - 2*X + 1))
F2 = acos(((X - 2)*X + 1))
F3 = acos(((X - 1)^2))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.3.2
% Lipshitz continiousity
disp('Example 1.3.2')
x = infsup(-3,-1)
y = infsup(2,4)

F = x^2+x+2*y

X = [ x; y ]
% Lipshitz continiousity
##    \begin{equation*}
##        f = g + h \quad \rightarrow \quad
##        \begin{aligned}[t]
##            &g = x^2, \quad h = x+y,\\
##            &L_f(\mbf{X}) = L_g + L_h.
##        \end{aligned}
##    \end{equation*}
%

% g = x^2
% h = x + y
Lg = 2 * mag(X(1))^2 * [1 0]
Lh =  [1 0] + [ 0 1 ]
Lf = Lg + Lh

x1 = [-2.8 2.5]
x2 = [-1.8 3.5]
absDFx1x2 =  abs (x1(1)^2+x1(1)+2*x1(2) - (x2(1)^2+x2(1)+2*x2(2)))
Lf * abs(x1-x2)'
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.3.3
% Interval estimation
disp('Example 1.3.3')

x = infsup(0,2)
y = infsup(1,4)
X = [x; y]
% f(x,y) = x^2+sqrt(y)
RANF = x^2+sqrt(y)

% Diff estimate = MV = Mean Value
Fdiff = [2*x 1/(2*sqrt(y)) ]

% points
c1 = [inf(x); inf(y)]
c2 = [inf(x); sup(y)]
c3 = [sup(x); sup(y)]

% MV estimates
Finfc1 = c1(1)^2 + sqrt(c1(2)) + Fdiff*( X - c1 )
Finfc2 = c2(1)^2 + sqrt(c2(2)) + Fdiff*( X - c2 )
Finfc3 = c3(1)^2 + sqrt(c3(2)) + Fdiff*( X - c3 )

% SL estimates
% Fslopex = x +cx
% Fslopey = 1/(sqrt(y)+sqrt(cy))
Fslopec1 = c1(1)^2 + sqrt(c1(2)) +  (x +c1(1))*( X(1) - c1(1) )+1/(sqrt(y)+sqrt(c1(2)))*(y-c1(2))
Fslopec2 = c2(1)^2 + sqrt(c2(2)) +  (x +c2(1))*( X(1) - c2(1) )+1/(sqrt(y)+sqrt(c2(2)))*(y-c2(2))
Fslopec3 = c3(1)^2 + sqrt(c3(2)) +  (x +c3(1))*( X(1) - c3(1) )+1/(sqrt(y)+sqrt(c3(2)))*(y-c3(2))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.3.3
% Bi-centric estimation
disp('Example 1.3.4')

x = infsup(0,2)
y = infsup(1,4)
X = [x; y]
% f(x,y) = x^2+sqrt(y)
RANF = x^2+sqrt(y)

% Diff estimate = MV = Mean Value
Fdiff = [2*x 1/(2*sqrt(y)) ]

% points

midFdiff = mid(Fdiff)
radFdiff = rad(Fdiff)
midradFdiff1 = midFdiff(1) / radFdiff (1)
midradFdiff2 = midFdiff(2) / radFdiff (2)
% cut function
p1 =cut( midradFdiff1, infsup(-1,1))
p2 =cut( midradFdiff2, infsup(-1,1))
% points
clow1 = mid(x) - p1*rad(x)
cup1 = mid(x) + p1*rad(x)
clow2 = mid(y) - p2*rad(y)
cup2 = mid(y) + p2*rad(y)
% MV estimates
Flow = (clow1)^2 + sqrt(clow2) + Fdiff(1) *( x - clow1 ) + Fdiff(2) *( y - clow2 )
Fup = (cup1)^2 + sqrt(cup2) + Fdiff(1) *( x - cup1 ) + Fdiff(2) *( y - cup2 )
% bi
intersect(Flow, Fup)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




