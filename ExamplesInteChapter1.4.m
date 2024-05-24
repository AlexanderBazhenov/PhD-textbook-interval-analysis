% 2022-06-27
% Examples Interval Analysis Chapter 1.4
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.4.1
% Newton 1D
disp('Example 1.4.1')

% 2022-08-06
% Example 1.7.1

% F = (x^3)*(exp(x)) -0.5
xarray = -1:0.1:1;
for ii = 1:numel(xarray)
  xnow = xarray(ii);
  Fnow =  (xnow^(3))*(exp(xnow)) - 0.5;
  Fnow1 = (xnow^(2))*(exp(xnow))*(xnow+3);
  Fnow2 = (xnow)*(exp(xnow))*(6+xnow*(xnow+6));
Farray(ii) = Fnow;
Farray1(ii) = Fnow1;
Farray2(ii) = Fnow2;
end
max(Farray)
min(Farray)

% Figure 1.6
figure
plot(xarray, Farray, '-r')
hold on
plot(xarray, Farray1, '-g')
plot(xarray, Farray2, '-b')
xxlim =[ -1.0 1.0 ]
xlim(xxlim )
ylim([ -1.0 2.0 ] )
set(gca, 'Fontsize', 14)
   xlabel('\it x_1');
   ylabel('\it x_2');
plot(mid(X0), 0, 'sk')
xx = [ xxlim(1) xxlim(2) ]
yy = [0 0 ]
plot(xx, yy, '--k')
grid on
figure_name_out = 'Figure1.6.png'
print('-dpng', '-r300', figure_name_out), pwd

% ---------------------------------------------
% START
X0= infsup(0.5,1.0)
%
iter = 0
for jj=1:6
NewtonStepEx141
end
% ---------------------------------------------
% /Newton method step
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.4.2
% Newton 2D
disp('Example 1.4.2')

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
F = [ F(1) F(2)]'
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

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1.4.3
% Krawczyk
disp('Example 1.4.3')
##\begin{equation} \label{e:NonLinEx2}
##\left \{ \begin{aligned}
##&x^{2/3} + y^{2/3} = 2, \\
##&\displaystyle \frac{1}{x+1} + 0.3 = y,
##\end{aligned} \right.
##\end{equation}

% 2023-02-03
% Example 1.7.4
% ---------------------------------------------
% START
X0=[infsup(1.1,1.5), infsup(0.6,0.9)]'
%X0=[infsup(0.7,1.5), infsup(0.6,1.0)]'

figure
X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
h1=line(Xplot, Yplot,'Color','b');
set(h1, 'linewidth', 1.5)
% Krawczyk operator
% 2018
% Krawczyk method step
X1=X0(1); X2=X0(2);
L = [2/(3*X1^(1/3)) 2/(3*X2^(1/3)); -1/(X1+1)^2 -1];
Lambda=inv(mid(L));
y=mid(X0);
y1=y(1); y2=y(2);
F = [infsup(0.0,0.0), infsup(0.0,0.0)]
F(1) = y1^(2/3)+y2^(2/3)-2;
F(2) = 1/(y1+1)-y2+0.3;
I = eye(2);
K = y-Lambda*F' -(I-Lambda*L)*(X0-y);
X1=K(1); X2=K(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','k');
Xnew=intersect(X0,K);
X0=Xnew;
% /Krawczyk method step
figure_name_out=strcat('Example1.7.4','Mid', '.png')
print('-dpng', '-r300', figure_name_out), pwd

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Example 1.4.4
% Krawczyk
disp('Example 1.4.4')
%\cite{Karpova}
\begin{equation} \label{e:NonLinEx2bic}
\left \{ \begin{aligned}
&x^{2/3} + y^{2/3} = [1.8,2.0], \\
&\displaystyle \frac{1}{x+1} + 0.3 = y.
\end{aligned} \right.
\end{equation}

% ---------------------------------------------
% START
%X0=[infsup(1.1,1.5), infsup(0.6,0.9)]'
X0=[infsup(0.7,1.5), infsup(0.6,1.0)]'

figure
X1=X0(1); X2=X0(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
h1=line(Xplot, Yplot,'Color','b');
set(h1, 'linewidth', 1.5)
% Krawczyk operator
% 2018
% Krawczyk method step
X1=X0(1); X2=X0(2);
L = [2/(3*X1^(1/3)) 2/(3*X2^(1/3)); -1/(X1+1)^2 -1];
Lambda=inv(mid(L));
y=mid(X0);
y1=y(1); y2=y(2);
F = [infsup(0.0,0.0), infsup(0.0,0.0)]
F(1) = y1^(2/3)+y2^(2/3)-infsup(1.8,2.0);
F(2) = 1/(y1+1)-y2+0.3;
I = eye(2);
K = y-Lambda*F' -(I-Lambda*L)*(X0-y);
X1=K(1); X2=K(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
line(Xplot, Yplot,'Color','k');
Xnew=intersect(X0,K);
X0=Xnew;
% /Krawczyk method step
figure_name_out=strcat('Example1.7.4','Mid', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0=[infsup(0.7,1.5), infsup(0.6,1.0)]'

figure
[Xplot, Yplot] = RectX1X2 (X0)
h1=line(Xplot, Yplot,'Color','b');
set(h1, 'linewidth', 1.5)
hold on

% Bicentric Krawczyk
%
##\begin{equation*}
##    \Upphi(x) := x - \varLambda F(x),
##\end{equation*}
%
##\cite{Maminov}:
##\begin{equation*}
##    \varLambda = \left(\alpha \ \m \mbf{L} + \beta I\right)^{-1}, \quad \alpha, \beta \in \mathbb{R}, \quad I = \begin{pmatrix}
##      1 & 0\\
##      0 & 1
##    \end{pmatrix}.
##\end{equation*}
##Если производится построение оператора Кравчика \eqref{e:KrawczykOperator},
##то $\alpha = 1$. При вычислении бицентрированной формы оператора Кравчика необходимо выбирать $\alpha = 2$.
##Значение коэффициента $\beta$ зависит от того, является ли интервальная матрица $\mbf{L} = \mbf{L}_{\natural}$ сингулярной:
##если да, то $\beta = 1$; в противном случае $\beta = 0$.
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
figure_name_out=strcat('Example1.7.4','BicX0', '.png')
print('-dpng', '-r300', figure_name_out), pwd


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

