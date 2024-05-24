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

% Example 1.1.1
disp('Example 1.1.1')
a= sqrt(3)
format long
disp(a)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.1.2
disp('Example 1.1.2')
a = infsup(-2, -1)
b = infsup(3, 4)
a+b
a-b
a*b
a/b
a+(-1)*b
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.1.3
disp('Example 1.1.3')
a = infsup(-5, 1)
mid(a)
rad(a)
wid(a)
mag(a)
mig(a)
a*a*a
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.1.4
disp('Example 1.1.4')
a = infsup(-12, -8)
chia = sup(a)/inf(a)
mag(a)
mag(a)*(1-chia)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.1.5
disp('Example 1.1.5')
a = infsup(2, 3)
b = infsup(1, 5)
c = infsup(-1, 8)

dist_ab = max(abs(inf(a)-inf(b)), abs((sup(a)-sup(b))))
dist_bc = max(abs(inf(b)-inf(c)), abs((sup(b)-sup(c))))
dist_ac = max(abs(inf(a)-inf(c)), abs((sup(a)-sup(c))))
max(dist_ab, dist_bc)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.1.6
disp('Example 1.1.6')
a = infsup(-1, 1)
a-a
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% Example 1.1.7
disp('Example 1.1.7')
a = infsup(3, 3)
b = infsup(3.5, 4.5)
format short
b/sqrt(a^2+b^2)
1/sqrt((a/b)^2+1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%