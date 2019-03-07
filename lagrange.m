function [y,ydot] = lagrange(X,Y,x,order)
%LAGRANGE Interpolation using Lagrange's method
%
% order 10 for 15 min sp3 gives centimeter accuracy
% order 1 for 30s clocks - few centimeters
%
%[x,xdot] = LAGRANGE(X,Y,x,order) uses Lagrange interpolation to find F(x)
%and F'(x), where X and Y describe the function Y = F(X).  The length of X
%should be equal to the number of rows of Y.  Y can have an arbitrary
%number of columns, such that the size of y (and ydot) is given by
%[length(x) size(Y,2)].
%
%F'(x) is only computed if a second output argument is provided.
%
% Peter F. Cervelli, U.S. Geological Survey, version 1.0, 2014/07/29.
% pcervelli@usgs.gov

    if size(X,2) == 1
        X = X';
    end
    N = length(X);
    nx = length(x);
    T = bsxfun(@minus,X,X') + speye(N);
    J0 = -order:order;
    for i = nx : -1 : 1
        J = find(x(i) <= X,1) + J0;
        if J(1) < 1
            J = 1 : 2*order + 1;
        elseif J(end) > N
            J = N - 2*order : N;
        end
        R = prod(T(J,J));
        delta = x(i) - X(J);
        delta(delta==0) = eps;
        pX = prod(delta)./delta;
        y(i,:) = pX./R*Y(J,:);
        if nargout > 1
            ydot(i,:) = (sum(pX)-pX)./delta./R*Y(J,:);
        end
    end
end
