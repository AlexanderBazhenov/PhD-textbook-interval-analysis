% 2022-06-27
% Examples Interval Analysis Chapter 1.1
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

% Example 1.2.1
% inf sup
disp('Example 1.2.1')
a =  [ infsup(-2,1) infsup(5, 7) ]'
infa = inf(a)
supa = sup(a)
A =  [ infsup(0,3) infsup(9, 11); infsup(-8,-6) infsup(3, 4) ]
infA = inf(A)
supA = sup(A)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.2
% vertex
disp('Example 1.2.2')
a =  [ infsup(-2,1) infsup(5, 7) ]'
vert_a1 = [ inf(a(1)) inf(a(2)) ]'
vert_a2 = [ inf(a(1)) sup(a(2)) ]'
vert_a3 = [ sup(a(1)) inf(a(2)) ]'
vert_a4 = [ sup(a(1)) sup(a(2)) ]'
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.3
disp('Example 1.2.3')
A = [infsup(-1,1) infsup(2, 5) ; infsup(3,8) infsup(-6, -4)]
B = [infsup(7,9) infsup(1, 3) ; infsup(-3,-2) infsup(-5, -1)]
% sum A B
C = A + B
% mult A B
C = A * B
% mult A-int b-point
b =  [ infsup(2,2); infsup(1, 1) ]
Ab = A*b
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.4
disp('Example 1.2.4')
% mid(AB) = mid(A)*B, B - point matrix
A = [infsup(1,8) infsup(0, 3) ; infsup(5,6) infsup(-4, -2)]
mid(A)
rad(A)
wid(A)
mag(A)
%
B = [infsup(1,1) infsup(2, 2) ; infsup(1,1) infsup(4, 4)]
AB = A*B
mid(AB)
mid(A)*B
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.5
% wrappting effect 1
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.6
% wrappting effect 2
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.7
% matrix mult commutativity
disp('Example 1.2.7')
A = [ infsup(-2,1) infsup(1,2) ]
B = [ infsup(2,3) ; infsup(5,6) ]
C = [ infsup(-1,1) infsup(4,5) ]
ABC1 = (A*B)*C
ABC2 = A*(B*C)
% C - vector
C = [ infsup(-1,1);  infsup(4,5) ]
A*(B+C)
A*B+A*C
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.8
% vector martix norm
disp('Example 1.2.8')
A = [infsup(6,9) infsup(1, 2) ; infsup(-3,3) infsup(5, 7)]
a =  [ infsup(4,9) ; infsup(-2 , -1) ]
%
norm1_a = mag(a(1))+mag(a(2))
norm2_a = sqrt(mag(a(1))^2+mag(a(2))^2)
%
norm1_A1 = max ( mag(A(1,1))+mag(A(2,1)), mag(A(1,2))+mag(A(2,2)))
norminf_A = max ( mag(A(1,1))+mag(A(1,2)), mag(A(1,2))+mag(A(2,2)))
%
AAT = mag(A)*mag(A')
svd(AAT)
%
Aa=A*a
norm1_AaM = mag(Aa(1))+mag(Aa(2))
norm1_Av = norm1_a * norm1_A1
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.9
% vector martix norm
disp('Example 1.2.9')
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
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.10
% dist vector matrix
disp('Example 1.2.10')
a =  [ infsup(-3,2) ; infsup(5, 6) ]
b =  [ infsup(7,10) ; infsup(-2, -1) ]

dist1_ab = max( abs(inf(a(1))-inf(b(1))),  abs(sup(a(1))-sup(b(1)))) + ...
max( abs(inf(a(2))-inf(b(2))),  abs(sup(a(2))-sup(b(2))))


dist2_ab = sqrt((max( abs(inf(a(1))-inf(b(1))),  abs(sup(a(1))-sup(b(1)))))^2 + ...
(max( abs(inf(a(2))-inf(b(2))),  abs(sup(a(2))-sup(b(2)))))^2)
Dist_ab = max( abs(inf(a(1))-inf(b(1))),  abs(sup(a(1))-sup(b(1)))) , ...
max( abs(inf(a(2))-inf(b(2))),  abs(sup(a(2))-sup(b(2))))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% Example 1.2.11
%  matrix dist properties
disp('Example 1.2.11')
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
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.12
disp('Example 1.2.12')
A= [infsup(1,3) infsup(10, 13) ; infsup(7,9) infsup(15, 17)]
B = A;
B(2,2) = B(2,2)+10
%
det(A)
det(B)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.13
disp('Example 1.2.13')
A= [infsup(1,3) infsup(10, 13) ; infsup(7,9) infsup(15, 17)]
B = A;
B(2,2) = B(2,2)+10
%
% Beck
disp('Beck')
max(svd(abs(inv(mid(A)))*rad(A)))
max(svd(abs(inv(mid(B)))*rad(B)))

% Rump
disp('Rump')
max(svd(rad(A)))
min(svd(mid(A)))

max(svd(rad(B)))
min(svd(mid(B)))


% Rohn-Rex
disp('Rohn-Rex')
min(svd(rad(A)))
max(svd(mid(A)))

min(svd(rad(B)))
max(svd(mid(B)))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.14
disp('Example 1.2.14')
A= [infsup(1,3) infsup(10, 13) ; infsup(7,9) infsup(15, 17)]
B = A;
B(2,2) = B(2,2)+10
C = A;
C(2,2) = C(2,2)+6

det(C)

det(inv(mid(A))*A)
det(inv(mid(C))*C)
max(eig(abs(inv(mid(C)))*rad(C)))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.2.15
disp('Example 1.2.15')
A = [infsup(1,3) infsup(10, 13) ; infsup(7,9) infsup(15, 17)]
%
Ainv11 = 1/ ( A(1,1) - A(1,2)*A(2,1) / A(2,2))
Ainv12 = 1/ ( A(2,1) - A(1,1)*A(2,2) / A(1,2))
Ainv21 = 1/ ( A(1,2) - A(1,1)*A(2,2) / A(2,1))
Ainv22 = 1/ ( A(2,2) - A(1,2)*A(2,1) / A(1,1))
%
Ainv = [ Ainv11 Ainv12 ;  Ainv21 Ainv22 ]
det(Ainv)

Ainv / det(A)
det(Ainv / det(A))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% Example 1.2.16
disp('Example 1.2.16')
A = [infsup(2,2) infsup(-1, -1) ; infsup(-5,-4.95) infsup(3, 3)]
inv(A)


X0 = [infsup(2.5, 3.5) infsup(0.5, 1.5) ; infsup(3.5, 5.5) infsup(1, 2.5)]
%
n=0
%
for jj=1:7
X1= mid(X0)+X0*(eye(2)-A*mid(X0));
retval = Dist1AB22 (X1, X0)
X0 = X1;
end
retval = Dist1AB22 (X1, inv(A))
%
A = [infsup(2,2) infsup(-1, -1) ; infsup(-5,-4.995) infsup(3, 3)]
Ainv11 = 1/ ( A(1,1) - A(1,2)*A(2,1) / A(2,2))
Ainv12 = 1/ ( A(2,1) - A(1,1)*A(2,2) / A(1,2))
Ainv21 = 1/ ( A(1,2) - A(1,1)*A(2,2) / A(2,1))
Ainv22 = 1/ ( A(2,2) - A(1,2)*A(2,1) / A(1,1))
%
Ainv = [ Ainv11 Ainv12 ;  Ainv21 Ainv22 ]
retval = Dist1AB22 (X1, Ainv)
% X0 = [infsup(1.5, 3.5) infsup(0.5, 1.5) ; infsup(3.5, 5.5) infsup(1, 2.5)]


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
