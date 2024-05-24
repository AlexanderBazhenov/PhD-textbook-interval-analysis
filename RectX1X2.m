## Copyright (C) 2023 Папа
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} RectX1X2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Папа <Папа@DESKTOP-8A3H2B3>
## Created: 2023-02-04

function [Xplot, Yplot] = RectX1X2 (X)
X1=X(1); X2=X(2);
Xplot=[inf(X1) inf(X1) sup(X1) sup(X1) inf(X1) ];
Yplot=[inf(X2) sup(X2) sup(X2) inf(X2) inf(X2)];
endfunction
