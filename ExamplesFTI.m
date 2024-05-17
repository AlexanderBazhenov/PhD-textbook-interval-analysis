% 2022-06-27
% Examples FTI
% 2024-05-17
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PLACE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirroot = 'D:\Data\ST\2023\T\'
dirki = 'D:\Data\ST\2023\T\kinterval-0.0.1'
dirOld =  'd:\Data\ST\2022\T\'
% HomePC
dirroot = 'D:\ST\2024\T\'
dirki = 'D:\ST\2023\T\kinterval-0.0.1'
dirOld2022 =  'd:\ST\2022\T\'
dirOld2023 =  'd:\ST\2023\T\'
% FTI
dirroot ='E:\Users\Public\Documents\ST\2023\T\'
dirki = 'e:\Users\Public\Documents\ST\2023\T\kinterval-0.0.1'
dirOld =  'e:\Users\Public\Documents\ST\2022\T\'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     START    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(dirroot)
addpath(dirOld2022)
addpath(dirOld2023)
cd(dirroot)
% pkg list
pkg load interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.1.1
a= sqrt(3)
format long
disp(a)

% Example 1.1.2
a = infsup(-2, -1)
b = infsup(3, 4)
a+b
a-b
a*b
a/b
a+(-1)*b

% Example 1.1.3
a = infsup(-1, 1)
a-a
a-a+a-a+a-a


% Example 1.1.4
a = infsup(3, 3)
b = infsup(3.5, 4.5)
format short
b/sqrt(a^2+b^2)
1/sqrt((a/b)^2+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example 1.2.1
x0 = [0, 11]'
phi = midrad(pi/3, 1.6e-3)
R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ]
format short
T = R * R'
format long
T = R * R'

x1 = T * x0

n=100
x=x1
for jj=1:n
  x = T*x;
end

format short
disp(x)

%%%%%%%%%%%%%%% 2023-06-04
phi = midrad(pi/4, 1.6e-3)
R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ]
format short
T = R * R'
format long
T = R * R'
x1 = T * x0

n=100
x=x1
for jj=1:n
  x = T*x;
end

format short
disp(x)
%%%%%%%%%%%%%%% 2023-06-04


% Example 1.2.2
x0 = [ infsup(-0.5, 0.5), infsup(10.5, 11.5) ]'
phi = midrad(pi/3, 1.6e-3)
R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ]
n=6
x=x0
for jj=1:n
  x = R*x;
end
format short
disp(x)

phi = midrad(0, 1.6e-3)
R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ]
x=R * x0
format short
disp(x)


%%%%%%%%%%%%%%% 2023-06-04
phi = midrad(pi/4, 1.6e-3)
x0 = [ infsup(-0.5, 0.5), infsup(10.5, 11.5) ]'
R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ]
n=8
x=x0
for jj=1:n
  x = R*x;
end
format short
disp(x)

phi = midrad(0, 1.6e-3)
R = [ cos(phi) -sin(phi); sin(phi) cos(phi) ]
x=R * x0
format short
disp(x)


%%%%%%%%%%%%%%% 2023-06-04

%%%%%%%%%%%%%%%%%%%%%%%%%% 2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.2.1
% inf sup

a =  [ infsup(-2,1) infsup(5, 7) ]
infa = inf(a)
supa = sup(a)

A = [infsup(0,3) infsup(9, 11) ; infsup(-8,-6) infsup(3, 4)]
infA = inf(A)
supA = sup(A)

% Example 1.2.2
% vertex
a =  [ infsup(-2,1) infsup(5, 7) ]'
vert_a1 = [ inf(a(1)) inf(a(2)) ]'
vert_a2 = [ inf(a(1)) sup(a(2)) ]'
vert_a3 = [ sup(a(1)) inf(a(2)) ]'
vert_a4 = [ sup(a(1)) sup(a(2)) ]'

% Example 1.2.3
% mult A-int b-point
A = [infsup(-1,1) infsup(2, 5) ; infsup(3,8) infsup(-6, -4)]
b =  [ infsup(2,2); infsup(1, 1) ]
Ab = A*b

% Example 1.2.4
% mid(AB) = mid(A)*B, B - point matrix

A = [infsup(1,8) infsup(0, 3) ; infsup(5,6) infsup(-4, -2)]
B = [infsup(1,1) infsup(2, 2) ; infsup(1,1) infsup(4, 4)]
AB = A*B
mid(AB)
mid(A)*B

% Example 1.2.5
% wrappting effect 1

% Example 1.2.6
% wrappting effect 2

% Example 1.2.7
% matrix mult commutativity

A = [ infsup(-2,1) infsup(1,2) ]
B = [ infsup(2,3) ; infsup(5,6) ]
C = [ infsup(-1,1) infsup(4,5) ]

ABC1 = (A*B)*C
ABC2 = A*(B*C)

C = [ infsup(-1,1);  infsup(4,5) ]

A*(B+C)
A*B+A*C


% Example 1.2.8
% vector martix norm

A = [infsup(6,9) infsup(1, 2) ; infsup(-3,3) infsup(5, 7)]
a =  [ infsup(4,9) ; infsup(-2 , -1) ]


norm1_a = mag(a(1))+mag(a(2))
norm2_a = sqrt(mag(a(1))^2+mag(a(2))^2)

norm1_A1 = max ( mag(A(1,1))+mag(A(2,1)), mag(A(1,2))+mag(A(2,2)))
norminf_A = max ( mag(A(1,1))+mag(A(1,2)), mag(A(1,2))+mag(A(2,2)))

AAT = mag(A)*mag(A')
svd(AAT)

Aa=A*a
norm1_AaM = mag(Aa(1))+mag(Aa(2))
norm1_Av = norm1_a * norm1_A1

% Example 1.2.9
% vector martix norm

A1 = [infsup(-3,4) infsup(-1, 1) ; infsup(-1, 1) infsup(-1, 1)]
A2 = [infsup(-4,4) infsup(-1, 1) ; infsup(-1, 1) infsup(-1, 1)]
A3 = [infsup(-2,5) infsup(-1, 1) ; infsup(-1, 1) infsup(-1, 1)]

radA1 = rad(A1)
radA2 = rad(A2)
radA3 = rad(A3)

norminf_A1 = max ( radA1(1,1)+radA1(1,2), radA1(1,2) + radA1(2,2))
norminf_A2 = max ( radA2(1,1)+radA2(1,2), radA2(1,2) + radA2(2,2))
norminf_A3 = max ( radA3(1,1)+radA3(1,2), radA3(1,2) + radA3(2,2))

norminf_A1 = max ( mag(A1(1,1))+mag(A1(1,2)), mag(A1(1,2))+mag(A1(2,2)))
norminf_A2 = max ( mag(A2(1,1))+mag(A2(1,2)), mag(A2(1,2))+mag(A2(2,2)))
norminf_A3 = max ( mag(A3(1,1))+mag(A3(1,2)), mag(A3(1,2))+mag(A3(2,2)))

% Example 1.2.10
% dist vector matrix

a =  [ infsup(-3,2) ; infsup(5, 6) ]
b =  [ infsup(7,10) ; infsup(-2, -1) ]

dist1_ab = max( abs(inf(a(1))-inf(b(1))),  abs(sup(a(1))-sup(b(1)))) + ...
max( abs(inf(a(2))-inf(b(2))),  abs(sup(a(2))-sup(b(2))))


dist2_ab = sqrt((max( abs(inf(a(1))-inf(b(1))),  abs(sup(a(1))-sup(b(1)))))^2 + ...
(max( abs(inf(a(2))-inf(b(2))),  abs(sup(a(2))-sup(b(2)))))^2)

% Example 1.2.11
%  matrix dist properties

A = [infsup(-8,-5) infsup(-8, -1) ; infsup(-6, 8) infsup(-3, 4)]
B = [infsup(-10,6) infsup(-6, 6) ; infsup(-2, 2) infsup(-9, -7)]
C = [infsup(-7,-2) infsup(-4, 2) ; infsup(5, 10) infsup(-9, 1)]
D = [infsup(-5,10) infsup(3, 10) ; infsup(1, 4) infsup(-5, -4)]

AC = A + C
BD = B + D

dist1_AC = max( abs(inf(AC(1,1))-inf(BD(1,1))),  abs(sup(AC(1,1))-sup(BD(1,1)))) + ...
max( abs(inf(AC(1,2))-inf(BD(1,2))),  abs(sup(AC(1,2))-sup(BD(1,2)))) + ...
max( abs(inf(AC(2,1))-inf(BD(2,1))),  abs(sup(AC(2,1))-sup(BD(2,1)))) + ...
max( abs(inf(AC(2,2))-inf(BD(2,2))),  abs(sup(AC(2,2))-sup(BD(2,2))))


dist2_ACBD = sqrt(max( abs(inf(AC(1,1))-inf(BD(1,1))),  abs(sup(AC(1,1))-sup(BD(1,1))))^2 + ...
max( abs(inf(AC(1,2))-inf(BD(1,2))),  abs(sup(AC(1,2))-sup(BD(1,2))))^2 + ...
max( abs(inf(AC(2,1))-inf(BD(2,1))),  abs(sup(AC(2,1))-sup(BD(2,1))))^2 + ...
max( abs(inf(AC(2,2))-inf(BD(2,2))),  abs(sup(AC(2,2))-sup(BD(2,2))))^2)

dist2_AB = sqrt(max( abs(inf(A(1,1))-inf(B(1,1))),  abs(sup(A(1,1))-sup(B(1,1))))^2 + ...
max( abs(inf(A(1,2))-inf(B(1,2))),  abs(sup(A(1,2))-sup(B(1,2))))^2 + ...
max( abs(inf(A(2,1))-inf(B(2,1))),  abs(sup(A(2,1))-sup(B(2,1))))^2 + ...
max( abs(inf(A(2,2))-inf(B(2,2))),  abs(sup(A(2,2))-sup(B(2,2))))^2)

dist2_CD = sqrt(max( abs(inf(C(1,1))-inf(D(1,1))),  abs(sup(C(1,1))-sup(D(1,1))))^2 + ...
max( abs(inf(C(1,2))-inf(D(1,2))),  abs(sup(C(1,2))-sup(D(1,2))))^2 + ...
max( abs(inf(C(2,1))-inf(D(2,1))),  abs(sup(C(2,1))-sup(D(2,1))))^2 + ...
max( abs(inf(C(2,2))-inf(D(2,2))),  abs(sup(C(2,2))-sup(D(2,2))))^2)

dist2_AB + dist2_CD

Dist_AB = [ max( abs(inf(A(1,1))-inf(B(1,1))),  abs(sup(A(1,1))-sup(B(1,1))))  ...
max( abs(inf(A(1,2))-inf(B(1,2))),  abs(sup(A(1,2))-sup(B(1,2)))) ; ...
max( abs(inf(A(2,1))-inf(B(2,1))),  abs(sup(A(2,1))-sup(B(2,1)))) ...
max( abs(inf(A(2,2))-inf(B(2,2))),  abs(sup(A(2,2))-sup(B(2,2)))) ]

Dist_CD = [ max( abs(inf(C(1,1))-inf(D(1,1))),  abs(sup(C(1,1))-sup(D(1,1))))  ...
max( abs(inf(C(1,2))-inf(D(1,2))),  abs(sup(C(1,2))-sup(D(1,2)))) ; ...
max( abs(inf(C(2,1))-inf(D(2,1))),  abs(sup(C(2,1))-sup(D(2,1)))) ...
max( abs(inf(C(2,2))-inf(D(2,2))),  abs(sup(C(2,2))-sup(D(2,2)))) ]

Dist_AB + Dist_CD

Dist_ACBD = [ max( abs(inf(AC(1,1))-inf(BD(1,1))),  abs(sup(AC(1,1))-sup(BD(1,1)))) ...
max( abs(inf(AC(1,2))-inf(BD(1,2))),  abs(sup(AC(1,2))-sup(BD(1,2)))) ; ...
max( abs(inf(AC(2,1))-inf(BD(2,1))),  abs(sup(AC(2,1))-sup(BD(2,1))))  ...
max( abs(inf(AC(2,2))-inf(BD(2,2))),  abs(sup(AC(2,2))-sup(BD(2,2)))) ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.2.4

A = [infsup(1,3) infsup(10, 13) ; infsup(7,9) infsup(15, 17)]

B = A;
B(2,2) = B(2,2)+10

det(A)
det(B)

% Beck
max(eig(abs(inv(mid(A)))*rad(A)))
max(eig(abs(inv(mid(B)))*rad(B)))

% Rump
max(eig(rad(A)))
min(eig(mid(A)))

max(eig(rad(B)))
min(eig(mid(B)))


% Rohn-Rex
min(eig(rad(A)))
max(eig(mid(A)))

min(eig(rad(B)))
max(eig(mid(B)))

% Example 1.2.5

C = A;
C(2,2) = C(2,2)+6

det(C)

det(inv(mid(A))*A)
det(inv(mid(C))*C)
max(eig(abs(inv(mid(C)))*rad(C)))

% Example 1.2.6

A = [infsup(1,3) infsup(10, 13) ; infsup(7,9) infsup(15, 17)]

Ainv = [ A(2,2) -A(1,2) ; -A(2,1) A(1,1) ]
Ainv / det(A)
det(inv(A))
Ainv / det(A)
det(Ainv / det(A))

% Example 1.2.7

A = [infsup(2,2) infsup(-1, -1) ; infsup(-5,-4.95) infsup(3, 3)]
inv(A)

X0 = [infsup(2.2, 3.8) infsup(0.2, 1.8) ; infsup(3.2, 5.8) infsup(0.7, 2.8)]

n=0
X0 = [infsup(2.5, 3.5) infsup(0.5, 1.5) ; infsup(3.5, 5.5) infsup(1, 2.5)]

X1= mid(X0)+X0*(eye(2)-A*mid(X0));
n = n+1
X1
X0 = X1;

A = [infsup(2,2) infsup(-1, -1) ; infsup(-5,-4.995) infsup(3, 3)]
inv(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.3.1

X = infsup(0,2)
X = infsup(1,2)
XD = X^2 - 2*X + 1
F1 = acos((X^2 - 2*X + 1))
F2 = acos(((X - 2)*X + 1))
F3 = acos(((X - 1)^2))

% Example 1.3.2

x = infsup(2,5)
y = infsup(1,3)

F1 = (x^2+2*y)/x
F2 = x + 2*y /x

X = [ x; y ]

X1 = [ infsup(2,2) ; infsup(1, 1) ]
x = X1(1)
y = X1(2)
X2 = [ infsup(5,5) ; infsup(3, 3) ]
x = X2(1)
y = X2(2)

distIR(X1, X2)
X1
X2
FX11 =  (X1(1)^2+2*X1(2))/X1(1)
FX12 =  (X2(1)^2+2*X2(2))/X2(1)
distIR(FX11, FX12)
FX21 =  X1(1)+2*X1(2)/X1(1)
FX22 =  X2(1)+2*X2(2)/X2(1)
distIR(FX21, FX22)

Lf= [0.5 0.85]
Lf*distIR(X1, X2)

norm( 2*rad(X),1)

RANF = infsup(3.0, 6.2)
distIR(F1, RANF )
distIR(F2, RANF )

C1 = distIR(F1, RANF ) / norm( 2*rad(X),1)
C2 = distIR(F2, RANF ) / norm( 2*rad(X),1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.3.3

x = infsup(1,5)

RANF = x^2

% Diff estimate = MV = Mean Value
Fdiff = 2*x
% points
cinf = inf(x)
cmid = mid(x)
csup = sup(x)
% MV estimates
Finf = inf(x)^2 + Fdiff *( x - cinf )
Fmid = mid(x)^2 + Fdiff *( x - cmid )
Fsup = sup(x)^2 + Fdiff *( x - csup )

% Slope estimate
Fslope = x+y
% SL estimates
SLinf = inf(x)^2 + ( x + cinf  ) *( x - cinf )
SLmid = mid(x)^2 + ( x + cmid  ) *( x - cmid )
SLsup = sup(x)^2 + ( x + csup  ) *( x - csup )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023-11-06
% Example 1.3.3 - 2 variables

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.3.4

x = infsup(1,5)

RANF = x^2

Fdiff = 2*x

midFdiff = mid(Fdiff)
radFdiff = rad(Fdiff)
midradFdiff = midFdiff / radFdiff
% cut function
p =cut( midradFdiff, infsup(-1,1))
% points
clow = mid(x) - p*rad(x)
cupp = mid(x) + p*rad(x)
% MV estimates
Flow = inf(x)^2 +  Fdiff *( x - clow )
Fupp = sup(x)^2 +  Fdiff *( x - cupp )
% bi
intersect(Flow, Fupp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-13
% Example  1.5.1

% Example  1.4.1
infA = [ 3 -5; -5 -3 ]
supA= [ 6 2 ; 7 -1 ]
infb = [ -2; -1 ]
supb = [ 2; 1 ]
A=infsup(infA, supA)
b=infsup(infb, supb)
%
midA=mid(A)
Lambda = inv(midA)
Lambda*A
Lambda*b
%
Eta = norm(eye(2)-Lambda*A, inf)
norm(Lambda*b, inf)

infA = [ -6 2; 5 3 ]
supA= [ -5 3 ; 7 10 ]
infb = [ -3; -1 ]
supb = [ 1; 5 ]
A=infsup(infA, supA)
b=infsup(infb, supb)
%
midA=mid(A)
Lambda = inv(midA)
Lambda*A
Lambda*b
%
Eta = norm(eye(2)-Lambda*A, inf)
norm(Lambda*b, inf)
theta = mag(norm(Lambda*b, inf)) / (1 - mag(Eta) )
%
det(inv(mid(A))*A)
%
xhat = inv(mid(A))*mid(b)
% Beck
% \begin{equation}\label{e:BekkEst}
%\Delta=\left(I-\left|(\mid \mbf{A})^{-1}\right| \cdot \rad \mbf{A}\right)^{-1}\left|(\mid \mbf{A})^{-1}\right|(\rad \mbf{A}|\hat{x}|+\rad \mbf{b}).
%\end{equation}
%xhat = [23; 10]*2/203
aa = inv(eye(2)-abs(inv(mid(A)))*rad(A))
bb = abs(inv(mid(A)))
cc = rad(A)*abs(xhat)+rad(b)
delta = aa*bb*cc
delta = inv(eye(2)-abs(inv(mid(A)))*rad(A))*abs(inv(mid(A)))*(rad(A)*abs(xhat)+rad(b))
%
xBeck = infsup(xhat-delta, xhat+delta)

addpath('d:\Data\ST\2022\T\IntLinInc2D\')
% Uni
Aq = ['E' 'E'; 'E' 'E' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
Part_SSW_color(P1,P2,P3,P4, [ 0 1 0])
hold on
%
POS = [ inf(xBeck(1)),inf(xBeck(2)), 2*rad(xBeck(1)),  2*rad(xBeck(2)) ]
rectangle( "Position", POS, "EdgeColor", [0 0 1] )
POS = [ -theta, -theta, 2*theta, 2*theta ]
rectangle( "Position", POS,  "EdgeColor", [1 0 0]  )
xlabel('\it x_1')
ylabel('\it x_2')
set(gca, "Fontsize", 14)

figure_name_out = 'Figure 1.15.png'
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-14
% Example  1.5.2

infA = [ 4 -2; 2 6 ]
supA= [ 5 0-1 ; 3 7 ]
infb = [ 4; 7 ]
supb = [ 5; 8 ]
A=infsup(infA, supA)
b=infsup(infb, supb)

[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
Part_SSW_color(P1,P2,P3,P4, [ 0 1 0])
hold on

X106 = [ infsup(0.848, 1.775); infsup(0.239, 1.051) ]
bb = A*X106
POS = [ inf(X106(1)), inf(X106(2)), 2*rad(X106(1)), 2*rad(X106(2)) ]
rectangle( "Position", POS,  "EdgeColor", [1 0 0]  )

figure_name_out = 'Figure 1.16.png'
print('-dpng', '-r300', figure_name_out), pwd

% Example 1.5.2

infA = [ 4 -2; 2 6]
supA = [5 -1; 3 7]
infb = [4; 7]
supb = [5; 8]
A = infsup(infA, supA)
b= infsup(infb, supb)

%addpath('d:\Data\ST\2022\T\IntLinInc2D\')
% Uni
Aq = ['E' 'E'; 'E' 'E' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
Part_SSW_color(P1,P2,P3,P4, [ 0 1 0])
hold on
%
X0=[infsup(-5, 5) ; infsup(-5 , 5)]
POS = [ inf(X0(1)),inf(X0(2)), 2*rad(X0(1)),  2*rad(X0(2)) ]
rectangle( "Position", POS, "EdgeColor", [0 0 1] )

xlabel('\it x_1')
ylabel('\it x_2')
set(gca, "Fontsize", 14)


% TRY
X0=[infsup(-5, 5) ; infsup(-5 , 5)]
% TRY

X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
%line(Xplot, Yplot,'Color','blue');
thickness=3
line(Xplot, Yplot,'Color','k', 'LineWidth',thickness);
START_str=strcat(' START INI X= ', num2str(inf(X1), '%2.1f'), '.', num2str(sup(X1), '%2.1f'),  ' Y= ', num2str(inf(X2), '%2.1f'), '.', num2str(sup(X2), '%2.1f'))
% /START


% GAUSS-SEIDEL method step 2x2
X1=X0(1); X2=X0(2);

GS=X0

res1=(b(1)-A(1,2)*X0(2))/A(1,1)
res2=(b(2)-A(2,1)*X0(1))/A(2,2)
GZ(1,1)=res1
GZ(2,1)=res2
% /GAUSS-SEIDEL operator 2x2

X1=GZ(1); X2=GZ(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
%line(Xplot, Yplot,'Color','red');
line(Xplot, Yplot,'Color','k');
Xnew=intersect(X0,GZ)
X0=Xnew;
% /GAUSS-SEIDEL  method step

xlim([-5.5, 5.5])
ylim([-5.5, 5.5])
figure_name_out ='Example 1.5.2.png'
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.5.2
% \eqref{e:Karpova3}
infA = [ 4 -2; 2 6]
supA = [5 -1; 3 7]
infb = [4; 7]
supb = [5; 8]
A = infsup(infA, supA)
b= infsup(infb, supb)

addpath('d:\Data\ST\2022\T\IntLinInc2D\')
% Uni
Aq = ['E' 'E'; 'E' 'E' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
Part_SSW_color(P1,P2,P3,P4, [ 0 1 0])
hold on
%
X0=[infsup(-5, 5) ; infsup(-5 , 5)]
POS = [ inf(X0(1)),inf(X0(2)), 2*rad(X0(1)),  2*rad(X0(2)) ]
rectangle( "Position", POS, "EdgeColor", [0 0 1] )

xlabel('\it x_1')
ylabel('\it x_2')
set(gca, "Fontsize", 14)


% TRY
X0=[infsup(-5, 5) ; infsup(-5 , 5)]
% TRY


X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
%line(Xplot, Yplot,'Color','blue');
thickness=3
line(Xplot, Yplot,'Color','k', 'LineWidth',thickness);
START_str=strcat(' START INI X= ', num2str(inf(X1), '%2.1f'), '.', num2str(sup(X1), '%2.1f'),  ' Y= ', num2str(inf(X2), '%2.1f'), '.', num2str(sup(X2), '%2.1f'))
% /START




% GAUSS-SEIDEL method step 2x2
X1=X0(1); X2=X0(2);

GS=X0

res1=(b(1)-A(1,2)*X0(2))/A(1,1)
res2=(b(2)-A(2,1)*X0(1))/A(2,2)
GZ(1,1)=res1
GZ(2,1)=res2
% /GAUSS-SEIDEL operator 2x2

X1=GZ(1); X2=GZ(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
%line(Xplot, Yplot,'Color','red');
line(Xplot, Yplot,'Color','k');
Xnew=intersect(X0,GZ)
X0=Xnew;
% /GAUSS-SEIDEL  method step
xlim([-5.5, 5.5])
ylim([-5.5, 5.5])
figure_name_out ='Example 1.5.2.png'
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.5.3
% Karpova5
AA = [	infsup(2,4) ,infsup(3,5) ; 	infsup(-6,-4) ,infsup(5,10) ]
b = [ infsup(1,2); infsup(2,6)]
mig(AA)

retval = HBR (AA, b)

infA= inf(AA), supA = sup(AA)
infb= inf(b), supb = sup(b)

addpath('d:\Data\ST\2022\T\IntLinInc2D\')
% addpath('e:\Users\Public\Documents\ST\2022\T\IntLinInc2D\')
% Uni
Aq = ['E' 'E'; 'E' 'E' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
Part_SSW_color(P1,P2,P3,P4, [ 0 1 0])
hold on


% Preconditioning
Lambda = inv(mid(AA))
LambdaA = Lambda * AA
Lambdab =Lambda * b
infA= inf(LambdaA), supA = sup(LambdaA)
infb= inf(Lambdab ), supb = sup(Lambdab )
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
Part_SSW_color(P1,P2,P3,P4, [ 0 0 1])

POS = [ inf(retval(1)),inf(retval(2)), 2*rad(retval(1)),  2*rad(retval(2)) ]
rectangle( "Position", POS, "EdgeColor", [0 0 0] )
figure_name_out ='Example 1.5.3.png'
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.5.4
% \eqref{e:Karpova3}
infA = [ 4 -2; 2 6]
supA = [5 -1; 3 7]
infb = [4; 7]
supb = [5; 8]
A = infsup(infA, supA)
b= infsup(infb, supb)

Lambda = inv(mid(A))
retval = HBR (A, b)
Xopt = retval

% HBR + Gauss-Zeidel
Xopt = retval
retval = HBR0 (A, b)
X0 = retval

% GAUSS-SEIDEL method step 2x2
X1=X0(1); X2=X0(2);

GS=X0

res1=(b(1)-A(1,2)*X0(2))/A(1,1)
res2=(b(2)-A(2,1)*X0(1))/A(2,2)
GZ(1,1)=res1
GZ(2,1)=res2
% /GAUSS-SEIDEL operator 2x2

X1=GZ(1); X2=GZ(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
%line(Xplot, Yplot,'Color','red');
line(Xplot, Yplot,'Color','k');
Xnew=intersect(X0,GZ)
X0=Xnew;
% /GAUSS-SEIDEL  method step

Xmulti = diag(X0)

% Uni
Aq = ['E' 'E'; 'E' 'E' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
clr = [ 0 1 0 ]
hold on
      fill(P1(:,1),P1(:,2),clr)
      plot(P1(:,1),P1(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
% Optimal
X0 = Xopt
POS = [ inf(X0(1)),inf(X0(2)), 2*rad(X0(1)),  2*rad(X0(2)) ]
rectangle( "Position", POS, "EdgeColor", [0 0 1] )

xlabel('\it x_1')
ylabel('\it x_2')
set(gca, "Fontsize", 14)
% Multi
X0 = Xmulti
POS = [ inf(X0(1)),inf(X0(2)), 2*rad(X0(1)),  2*rad(X0(2)) ]
rectangle( "Position", POS, "EdgeColor", [1 0 0] )
% start
retval = HBR0 (A, b)
X0 = retval
POS = [ inf(X0(1)),inf(X0(2)), 2*rad(X0(1)),  2*rad(X0(2)) ]
rectangle( "Position", POS, "EdgeColor", [0 0 0] )
figure_name_out ='Example 1.5.4.png'
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-14
% Example  1.5.5

infA = [ 5 0; 0 6 ]
supA= [ 9 0 ; 0 7 ]
infb = [ 1; 2 ]
supb = [ 2; 6 ]
A=infsup(infA, supA)
b=infsup(infb, supb)

% Uni
Aq = ['E' 'E'; 'E' 'E' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
Part_SSW_color(P1,P2,P3,P4, [ 0 1 0])
hold on

Lambda = 0.5 * eye(2)
C = eye(2) - Lambda*A
d=Lambda*b

addpath('d:\Data\ST\2022\T\kinterval-0.0.1\')
% addpath('e:\Users\Public\Documents\ST\2022\T\kinterval-0.0.1\')
x1 = ki(4/5, -1/5)
x2 = ki(5/3, -1/3)
x = [x1 ; x2]

Cki = [ ki(inf( C(1,1) ), sup(C(1,1)) )  ki(inf( C(1,2) ), sup(C(1,2)) ) ; ki(inf( C(2,1) ), sup(C(2,1)) )  ki(inf( C(2,2) ), sup(C(2,2))) ]
dki = [ ki(inf( d(1) ), sup(d(1)) ) ; ki(inf( d(2) ), sup(d(2)) ) ]
Cki*x + dki + opp(x)

Lambda = 0.05 * eye(2)
C = eye(2) - Lambda*A
d=Lambda*b

xx = infsup (min(P1), max(P1) )
xxki = [ ki(inf( xx(1) ), sup(xx(1)) ) ; ki(inf( xx(2) ), sup(xx(2)) ) ]
Cki = [ ki(inf( C(1,1) ), sup(C(1,1)) )  ki(inf( C(1,2) ), sup(C(1,2)) ) ; ki(inf( C(2,1) ), sup(C(2,1)) )  ki(inf( C(2,2) ), sup(C(2,2))) ]
dki = [ ki(inf( d(1) ), sup(d(1)) ) ; ki(inf( d(2) ), sup(d(2)) ) ]
result = Cki*xxki + dki + opp(xxki)
wid(result)
if max(wid(result)) == 0 display('OK') end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example  1.5.6
% Krawczyk linear
% \eqref{e:Karpova5}
A = [	infsup(2,4) ,infsup(3,5) ; 	infsup(-6,-4) ,infsup(5,10) ]
b = [ infsup(1,2); infsup(2,6)]

X0 = HBR (A, b)
X0 = X0'

Y=inv(mid(A));
I = eye(2)
rho = max(eig(mag(I -Y*A)))


figure
X1=X0(1); X=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','b');

% Krawczyk operator
% 2018
% Krawczyk method step
K = Y*b -(I-Y*A)*X0;
X1=K(1); X2=K(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','k');
Xnew=intersect(X0,K);
X0=Xnew;
% /Krawczyk method step




% Esimate 1.5.1
Lambda = Y
Eta = norm(eye(2)-Lambda*A, inf)
norm(Lambda*b, inf)
theta = mag(norm(Lambda*b, inf)) / (1 - mag(Eta) )
%
X0 = theta * [infsup(-1,1) ; infsup(-1,1)   ]
figure
X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','k');
PlotKrawczyk
% Uni
infA= inf(AA), supA = sup(AA)
infb= inf(b), supb = sup(b)
Aq = ['E' 'E'; 'E' 'E' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
clr = [ 0 1 0 ]
hold on
      fill(P1(:,1),P1(:,2),clr)
      plot(P1(:,1),P1(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P2(:,1),P2(:,2),clr)
      plot(P2(:,1),P2(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
% Optimal
retval = HBR (A, b)
Xopt = retval
X0 = Xopt
POS = [ inf(X0(1)),inf(X0(2)), 2*rad(X0(1)),  2*rad(X0(2)) ]
rectangle( "Position", POS, "EdgeColor", [1 0 0] )

figure_name_out ='Example 1.5.6.png'
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1.6.1


% 1.4.8
% 2022-07-11
% Data
% 1.4.8
% 2022-07-08
infA = [ 2 1; -7 6]
supA= [ 5 2 ; -5 7 ]
infb = [ 3; 7]
supb = [ 4; 8 ]
A=infsup(infA, supA)
b=infsup(infb, supb)
[m, n]= size(A)

varargin=1
% ver 2020
[Tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,varargin)


% Mesh
a11 = A(1,1);
a12 = A(1,2);
a21 = A(2,1);
a22 = A(2,2);
b1 = b(1);
b2 = b(2);
xb=-2.8; xe=2.8; yb=-2.8; ye=2.8;
Nx=50; Ny=50;
xs=(xe-xb)/Nx; ys=(ye-yb)/Ny;
[X,Y] = meshgrid(xb:xs:xe, yb:ys:ye);
% [X,Y] = meshgrid(1, -2.8)
% Uni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 TT = min( rad(b1) - mig(mid(b1) - a11*X - a12*Y), ...
             rad(b2) - mig(mid(b2) - a21*X - a22*Y) );
 figure
   h = meshc(X,Y,TT);
   viewTol3d = [   34.666   57.475]
    view(viewTol3d)
      colormap(winter)
    caxis([-40 -7 -2])
    colorbar
    xlim([-2.5 2.5])
    ylim([-2.5 2.5])
    figure_name_out = 'Figureb1.19a.png'
print('-dpng', '-r300', figure_name_out), pwd
 figure
contourf(X,Y,TT, [-40 -7 -2])
   xlabel('\it x_1');
   ylabel('\it x_2');
    figure_name_out = 'Figureb1.19b.png'
    caxis([-40 -0.7])
      colormap(winter)
      axis('square')
      set(gca, 'Fontsize', 14)
hold on
plot(argmax(1), argmax(2),'xk')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-23
% Example 1.6.2

% Data
% Example 1.4.8

% 2022-07-08
infA = [ 2 1; -7 6]
supA= [ 5 2 ; -5 7 ]
infb = [ 3; 7]
supb = [ 4; 8 ]
A=infsup(infA, supA)
b=infsup(infb, supb)
[m, n]= size(A)

varargin=1
% ver 2020
[Tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,varargin)
%
K = 5
bnew = b +K*infsup(-1, 1)
[Tolmaxnew,argmax,envs,ccode] = tolsolvty(infA,supA,inf(bnew),sup(bnew),varargin)
%
a11 = A(1,1);
a12 = A(1,2);
a21 = A(2,1);
a22 = A(2,2);
b1 = b(1);
b2 = b(2);
xb=-2.8; xe=2.8; yb=-2.8; ye=2.8;
Nx=50; Ny=50;
xs=(xe-xb)/Nx; ys=(ye-yb)/Ny;
[X,Y] = meshgrid(xb:xs:xe, yb:ys:ye);
%
v = [ 0.5; 2]
 TTv = min( v(1)*(rad(b1) - mig(mid(b1) - a11*X - a12*Y)), ...
            v(2)*( rad(b2) - mig(mid(b2) - a21*X - a22*Y)) );
 figure
   h = meshc(X,Y,TTv);
   viewTol3d = [   34.666   57.475]
    view(viewTol3d)
      colormap(winter)
    caxis([-40 -7 -2])
    colorbar
    xlim([-2.5 2.5])
    ylim([-2.5 2.5])
    figure_name_out = 'Figureb1.21a.png'
print('-dpng', '-r300', figure_name_out), pwd
 figure
contourf(X,Y,TTv, [-40 -7 -2])
   xlabel('\it x_1');
   ylabel('\it x_2');
    figure_name_out = 'Figureb1.21b.png'
    caxis([-40 -0.7])
      colormap(winter)
      axis('square')
      set(gca, 'Fontsize', 14)
hold on
plot(argmax(1), argmax(2),'xk')
print('-dpng', '-r300', figure_name_out), pwd

bv = b +K*v*infsup(-1, 1)
[Tolmaxv,argmax,envs,ccode] = tolsolvty(infA,supA,inf(bv),sup(bv),varargin)
% D functional (CTL)
D = CTLfun (A, b)
Dnew = CTLfun (A, bnew)
Dv = CTLfun (A, bv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-25
% Example 1.6.3

% Example 1.6.2

% Data
% Example 1.4.8

% 2022-07-08
infA = [ 2 1; -7 6]
supA= [ 5 2 ; -5 7 ]
infb = [ 3; 7]
supb = [ 4; 8 ]
A=infsup(infA, supA)
b=infsup(infb, supb)
[m, n]= size(A)

varargin=1
% ver 2020
[Tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,varargin)

rad(A)
abs(argmax)
abs(Tolmax) / sum(abs(argmax))

e = 0.5
addpath('d:\Data\ST\2022\T\kinterval-0.0.1\')
% addpath('e:\Users\Public\Documents\ST\2022\T\kinterval-0.0.1\')

aa = ki(1,2)
innerminus(aa, aa)

E = e* ones(2,2)*ki(-1,1)
Aki = [ ki(  inf(A(1,1)), sup(A(1,1))  )...
ki(  inf(A(1,2)), sup(A(1,2))  ) ;  ...
 ki(  inf(A(2,1)), sup(A(2,1))  )...
 ki(  inf(A(2,2)), sup(A(2,2))  )]
% 	\label{e:exAmE1}
AmE = innerminus(Aki , E)

Anew =  ikMatrix( AmE )

[Tolmaxnew,argmaxnew,envs,ccode] = tolsolvty(inf(Anew),sup(Anew),infb,supb,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = [ 1.0 ; 0.9 ]
v(1)
E = e* diag(v) *ones(2,2)*ki(-1,1)
AmE2 = innerminus(Aki , E)

[Tolmax2,argmax2,envs,ccode] = tolsolvty(inf(AmE2 ),sup(AmE2 ),infb,supb,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = CTLfun (A, b)
Dnew = CTLfun (AmE, b)
Dv = CTLfun (AmE2, b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------ PLOT -------------------------------------


A = AmE
A = midrad(mid(A), rad(A))
TT1 = TolGraph2 (A, b, X, Y);  figure
contourf(X,Y,TT1, [0 0.22 ])
   xlabel('\it x_1');
   ylabel('\it x_2');
xlim([ 0.25 0.45])
ylim([ 1.4 1.6])
    figure_name_out = 'Figureb1.22a.png'
      set(gca, 'Fontsize', 14)
print('-dpng', '-r300', figure_name_out), pwd

A = AmE2
TT2 = TolGraph2 (A, b, X, Y);
contourf(X,Y,TT2, [0 0.22 ])
   xlabel('\it x_1');
   ylabel('\it x_2');
xlim([ 0.25 0.45])
ylim([ 1.4 1.6])
    figure_name_out = 'Figureb1.22b.png'
      set(gca, 'Fontsize', 14)
print('-dpng', '-r300', figure_name_out), pwd

% ------------------ .PLOT -------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-28
% Example 1.6.2

% Data
% Example 1.4.8

% 2022-07-08
infA = [ 2 1; -7 6]
supA= [ 5 2 ; -5 7 ]
infb = [ -2; 2]
supb = [ 9; 13 ]
A=infsup(infA, supA)
b=infsup(infb, supb)
[m, n]= size(A)
b1 = infb
b2 = supb
a11 = A(1,1);
a12 = A(1,2);
a21 = A(2,1);
a22 = A(2,2);
%
varargin = 1
[Tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,varargin)
%
y = [ -0.2; 0.5]
[ TT, EST ] = Shaidurov(A, b, y)
%
dirIntLin = strcat(dirroot, 'IntLinInc2D\')
%addpath('d:\T\IntLinInc2D\')
addpath(dirIntLin)
% Tol
Aq = ['A' 'A'; 'A' 'A' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
clr = [ 0 1 0 ]
% Fig 1.23a
figure
       plot(y(1), y(2) ,'x','MarkerFaceColor','r','MarkerSize',3)
hold on
      fill(P1(:,1),P1(:,2),clr)
      plot(P1(:,1),P1(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P2(:,1),P2(:,2),clr)
      plot(P2(:,1),P2(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
% EST Shaidurov
out = EST
xx = [ inf(out(1)) sup(out(1))]
yy = [ inf(out(2)) sup(out(2))]
       vertices = [xx(1), xx(1), xx(2), xx(2), xx(1) ; ...
                   yy(1), yy(2), yy(2), yy(1), yy(1)   ]';
       PP = patch (vertices(:,1), vertices(:,2), 0.95*[1 1 0]);
%
plot(argmax(1), argmax(2), 'xr')
%
set(gca, 'Fontsize', 14)
   xlabel('\it x_1');
   ylabel('\it x_2');
xlim([ -1.0 2.0 ] )
ylim([ -1.0 3.0 ] )
axis('equal')
%
figure_name_out = 'Figureb1.23a.png'
print('-dpng', '-r300', figure_name_out), pwd
%
y = argmax
[ TT, EST ] = Shaidurov(A, b, y)
% Tol
Aq = ['A' 'A'; 'A' 'A' ]
bq = [ 'E' ;'E']
relations= [ '=' ; '=' ]
[V,P1,P2,P3,P4] = MixQtr2D(infA, supA, Aq, infb, supb, bq, relations)
clr = [ 0 1 0 ]

% Fig 1.23b
figure
       plot(y(1), y(2) ,'x','MarkerFaceColor','r','MarkerSize',3)
hold on
      fill(P1(:,1),P1(:,2),clr)
      plot(P1(:,1),P1(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P2(:,1),P2(:,2),clr)
      plot(P2(:,1),P2(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
% EST Shaidurov
out = EST
xx = [ inf(out(1)) sup(out(1))]
yy = [ inf(out(2)) sup(out(2))]
       vertices = [xx(1), xx(1), xx(2), xx(2), xx(1) ; ...
                   yy(1), yy(2), yy(2), yy(1), yy(1)   ]';
       PP = patch (vertices(:,1), vertices(:,2), 0.95*[1 1 0]);
%
plot(argmax(1), argmax(2), 'xr')
%
set(gca, 'Fontsize', 14)
   xlabel('\it x_1');
   ylabel('\it x_2');
xlim([ -1.0 2.0 ] )
ylim([ -1.0 3.0 ] )
axis('equal')
%
figure_name_out = 'Figureb1.23b.png'
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-08-06
% Example 1.7.1

% F = 0.5 - (x^3)*(exp(x))
xarray = -3:0.1:1;
for ii = 1:numel(xarray)
  xnow = xarray(ii);
  Fnow = 0.5 - (xnow^(3))*(exp(xnow));
  Fnow1 = - (xnow^(2))*(exp(xnow))*(xnow+3);
  Fnow2 = - (xnow)*(exp(xnow))*(6+xnow*(xnow+6));
Farray(ii) = Fnow;
Farray1(ii) = Fnow1;
Farray2(ii) = Fnow2;
end
max(Farray)
min(Farray)

% Figure 1.24
figure
plot(xarray, Farray, '-r')
hold on
plot(xarray, Farray1, '-g')
plot(xarray, Farray2, '-b')
xxlim =[ -3.0 3.0 ]
xlim(xxlim )
ylim([ -3.0 3.0 ] )
set(gca, 'Fontsize', 14)
   xlabel('\it x_1');
   ylabel('\it x_2');
plot(mid(X0), 0, 'sk')
xx = [ xxlim(1) xxlim(2) ]
yy = [0 0 ]
plot(xx, yy, '--k')
figure_name_out = 'Figureb1.24.png'
print('-dpng', '-r300', figure_name_out), pwd

% ---------------------------------------------
% START
X0= infsup(0.1,1.0)
%
iter = 0
for jj=1:6
NewtonStep1
end
% ---------------------------------------------
% Newton method step

iter = iter + 1
X1=X0;
Xnow=X1;
%
% y=mid(Xnow);
y=inf(Xnow);
% y=sup(Xnow);
%
  Fnow = 0.5 - (y^(3))*(exp(y));
  Fnow1 = - (Xnow^(2))*(exp(Xnow))*(Xnow+3);
%
% Newtonoperator
% N = y - Fnow/Fnow1;
N = y - Fnow/Fnow1
Xnew=intersect(X0,N)
X0=Xnew;
% /Newton method step
% ---------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023-01-05
% Example 1.7.2
##\begin{equation} \label{e:NonLinEx1}
##\left \{ \begin{aligned}
##&0.8x^2 + 1.5y^2 = 1, \\
##&\ln{x} = y,
##\end{aligned} \right.
##\end{equation}

figure
hold on
xx = 0:sqrt(1/0.8)/50: sqrt(1/0.8)
yy = sqrt( (1- 0.8*xx.^2 )/1.5 )
plot(xx, yy,'.k');
plot(xx, yy,'-k');

xx = 0:sqrt(1/0.8)/50: sqrt(1/0.8)
yy = log(xx)
plot(xx, yy,'.k');
plot(xx, yy,'-k');

xlim([0 1.3])
ylim([-0.3 0.3])
% ---------------------------------------------
% START
X0=[infsup(0.8,1.2), infsup(-0.2,0.2)]'

X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','b');
%
##F = X0
##F(1) = 0.8*X1^2 +1.6*X2^2 -1
##F(2) = log(X1)-X2


% Newton method step
X1=X0(1); X2=X0(2);
J=[1.6*X1 3*X2; 1/X1 -1]
Y=inv(J)
y=mid(X0);
y1=y(1); y2=y(2);
F(1) = 0.8*y1^2 +1.5*y2^2 -1;
F(2) = log(y1)-y2;
% Newtonoperator
N = y' - (Y*F)';
X1=N(1); X2=N(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','red');
Xnew=intersect(X0',N)
X0=Xnew';
% /Newton method step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(y1,y2,'sr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023-01-09
% Example 1.7.3
% % Krawczyk non-linear
##\begin{equation} \label{e:NonLinEx2}
##\left \{ \begin{aligned}
##&\frac{1}{x^2-2} = y, \\
##&x = y^2,
##\end{aligned} \right.
##\end{equation}



% ---------------------------------------------
% START
X0=[infsup(1.5,1.9), infsup(1.1,1.5)]'


figure
X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','b');

% Krawczyk operator
% 2018
% Krawczyk method step
X1=X0(1); X2=X0(2);
L = [-2*X1/(X1^2-2)^2 -1; 1 -2*X2];
Lambda=inv(mid(L));
y=mid(X0);
y1=y(1); y2=y(2);
F(1) = 1/(y1^2-2) -y2;
F(2) = y1-y2^2;
I = eye(2);
K = y-Lambda*F' -(I-Lambda*L)*(X0-y);
X1=K(1); X2=K(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','k');
Xnew=intersect(X0,K);
X0=Xnew;
% /Krawczyk method step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2023-01-13
% Example 1.7.5
##????????? ???????????? ????? ???????, ?????? ??????? ????????????? ??????????? ????????? \label{p:NonLin1D}
##\begin{equation} \label{e:NonLinEx3}
##    [2.0,2.2]\sqrt{x} + [0.3,0.4]x - [3.0,3.1] = 0.
##\end{equation}

a2 = infsup(2.0,2.2)
a1 = infsup(0.3,0.4)
a0 = -infsup(3.0,3.1)

% START
X0= infsup(1, 2)
%
% F = a2*sqrt(X0) + a1*X0 +a0
% J = a2/(2*sqrt(X0)) + a1*X0

disp(' mid X')
X0= infsup(1, 2)
% Newton method step
for steps=1:4
steps-1
J=a2/(2*sqrt(X0)) + a1*X0;
Y=inv(J);
y=mid(X0);
F = a2*sqrt(y) + a1*y + a0 ;
% Newtonoperator
N = y - (Y*F)
% N = y - F/J
Xnew=intersect(X0',N)
X0=Xnew';
% /Newton method step
end


disp('inf X')
X0= infsup(1, 2)
% Newton method step
for steps=1:5
  steps
J=a2/(2*sqrt(X0)) + a1*X0;
Y=inv(J);
y=inf(X0);
F = a2*sqrt(y) + a1*y +a0 ;
% Newtonoperator
N = y - (Y*F)
Xnew=intersect(X0',N)
X0=Xnew';
% /Newton method step
end

disp('sup X')
X0= infsup(1, 2)
% Newton method step
for steps=1:7
J=a2/(2*sqrt(X0)) + a1*X0;
Y=inv(J);
y=inf(X0);
F = a2*sqrt(y) + a1*y +a0 ;
% Newtonoperator
N = y - (Y*F)
Xnew=intersect(X0',N)
X0=Xnew';
% /Newton method step
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.7.6
##\begin{equation} \label{e:NonLinEx4}
##    \left \{
##    \begin{aligned}
##    &[-0.75,-0.50]x + [0.50,0.75]y = 0, \\
##    &[1.00,2.00] \arcsin{x} - y = 0.
##    \end{aligned}
##    \right.
##\end{equation}

% ---------------------------------------------
% START
X0=[infsup(0.1,0.7), infsup(0.1,0.7)]'

figure
X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','b');
hold on
%
##F = X0
##F(1) = 0.8*X1^2 +1.6*X2^2 -1
##F(2) = log(X1)-X2

disp(' mid X')
X0=[infsup(0.1,0.7), infsup(0.1,0.7)]'
% Newton method step
X1=X0(1); X2=X0(2);
##J=[1.6*X1 3*X2; 1/X1 -1]
##Y=inv(J)
y=mid(X0);
y1=y(1); y2=y(2);
F(1) = infsup(-0.75, -0.5)*y1 + infsup(0.5, 0.75) -y2;
F(2) = infsup(1, 2)*asin(y1)-y2;
Y = [ -1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) ), 1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) );
-1/( infsup(0.25, 0.75)*sqrt(1-y1^2) + infsup(0.5, 2.0)),  infsup(-0.75, -0.5)/(infsup(0.5, 0.75)+   infsup(1, 2)/sqrt(1-y1^2)) ]
% Newtonoperator
N = y - (Y*F');
X1=N(1); X2=N(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','red');
Xnew=intersect(X0,N)
X0=Xnew;
% /Newton method step

disp(' inf X')
X0=[infsup(0.1,0.7), infsup(0.1,0.7)]'
% Newton method step
X1=X0(1); X2=X0(2);
##J=[1.6*X1 3*X2; 1/X1 -1]
##Y=inv(J)
y=inf(X0);
y1=y(1); y2=y(2);
F(1) = infsup(-0.75, -0.5)*y1 + infsup(0.5, 0.75) -y2;
F(2) = infsup(1, 2)*asin(y1)-y2;
Y = [ -1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) ), 1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) );
-1/( infsup(0.25, 0.75)*sqrt(1-y1^2) + infsup(0.5, 2.0)),  infsup(-0.75, -0.5)/(infsup(0.5, 0.75)+   infsup(1, 2)/sqrt(1-y1^2)) ]
% Newtonoperator
N = y - (Y*F');
X1=N(1); X2=N(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','red');
Xnew=intersect(X0,N)
X0=Xnew;
% /Newton method step


disp(' sup X')
X0=[infsup(0.1,0.7), infsup(0.1,0.7)]'
% Newton method step
X1=X0(1); X2=X0(2);
##J=[1.6*X1 3*X2; 1/X1 -1]
##Y=inv(J)
y=sup(X0);
y1=y(1); y2=y(2);
F(1) = infsup(-0.75, -0.5)*y1 + infsup(0.5, 0.75) -y2;
F(2) = infsup(1, 2)*asin(y1)-y2;
Y = [ -1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) ), 1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) );
-1/( infsup(0.25, 0.75)*sqrt(1-y1^2) + infsup(0.5, 2.0)),  infsup(-0.75, -0.5)/(infsup(0.5, 0.75)+   infsup(1, 2)/sqrt(1-y1^2)) ]
% Newtonoperator
N = y - (Y*F');
X1=N(1); X2=N(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','red');
Xnew=intersect(X0,N)
X0=Xnew;
% /Newton method step

disp(' infsup X')
X0=[infsup(0.1,0.7), infsup(0.1,0.7)]'
% Newton method step
X1=X0(1); X2=X0(2);
##J=[1.6*X1 3*X2; 1/X1 -1]
##Y=inv(J)
y1=inf(X0(1));
y2=sup(X0(2));
y = [ y1; y2];
F(1) = infsup(-0.75, -0.5)*y1 + infsup(0.5, 0.75) -y2;
F(2) = infsup(1, 2)*asin(y1)-y2;
Y = [ -1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) ), 1/(infsup(0.5, 0.75) +  infsup(1, 2)/sqrt(1-y1^2) );
-1/( infsup(0.25, 0.75)*sqrt(1-y1^2) + infsup(0.5, 2.0)),  infsup(-0.75, -0.5)/(infsup(0.5, 0.75)+   infsup(1, 2)/sqrt(1-y1^2)) ]
% Newtonoperator
N = y - (Y*F');
X1=N(1); X2=N(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','red');
Xnew = intersect(X0,N)
X0=Xnew;
% /Newton method step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0=[infsup(0.1,0.7), infsup(0.1,0.7)]'

figure
X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
h =line(Xplot, Yplot,'Color','b');
hold on
XX = [ 0.2964 0.1811 0.4857 0.2835 ]
YY = [ 0.6018 0.3642 0.5072 0.2874 ]
h = plot(XX, YY, 'ob')

% Newton method step
X1=X0(1); X2=X0(2);
y1=inf(X0(1));
y2=inf(X0(2));
y = [ y1; y2];
ExamplesFTI175
N= y - (Y*F');
plotN
set(h, 'color', [0 0 0])
Nii = N;
%

%
y1=inf(X0(1));
y2=sup(X0(2));
y = [ y1; y2];
ExamplesFTI175
N= y - (Y*F');
plotN
set(h, 'color', [0 0 0])
Nis = N;
%set(h, 'linewidth', 2)
N=intersect(Nii,Nis)
plotN
set(h, 'color', [1 0 0])

Xnew = intersect(X0,N)
X0=Xnew;
% /Newton method step

figure_name_out=strcat('Example1.7.5','InfInf cap InfSup', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbaspect([1 1 1])

y1=mid(X0(1));
y2=mid(X0(2));
y = [ y1; y2];
ExamplesFTI175
N= y - (Y*F');
plotN
set(h, 'color', [0 0 0])
Nm = N;
%
N=intersect(Nii,Nm)
plotN
set(h, 'color', [1 0 0])





