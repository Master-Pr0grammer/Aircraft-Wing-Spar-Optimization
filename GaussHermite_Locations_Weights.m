function [xi,wi] = GaussHermite_Locations_Weights(n)
%   Compute Gauss Quadrature Hermite location and weights
%   n - Order
% Output:
%   xi -  location points
%   wi -  gaussian weights
format long
syms x;
Hsym = hermiteH (n,x);            % symbolic Hermite polynomial at n (Hn)
Hsym_1 = hermiteH (n-1,x);        % symbolic Hermite polynomial at n-1 (Hn-1)
Hn = sym2poly (Hsym);             % numeric Hermite polynomial at n (Hn)
Hn_1 = sym2poly (Hsym_1);         % numeric Hermite polynomial at n-1 (Hn-1)
xi = roots (Hn);                  % roots of Hn
wi = 2^(n-1) * factorial(n) * ... % weights
    sqrt(pi) ./(n .* polyval(Hn_1,xi)).^2;
end
