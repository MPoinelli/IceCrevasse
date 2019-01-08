function [coef,eps] = LSM (A,X,Y)
% Mattia Poinelli
% TU Delft
% October 2015
%
% Least Mean Square coefficients
%
% INPUT
% A: design matrix
% X: data
% Y: data
%
% OUTPUT
% coef: vector of parameters
% eps: residuals

B = ((A'*A)^(-1))*A';
coef = B*Y;
eps = Y-A*coef;

end
