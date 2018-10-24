function alpha = acos2(cos_alpha,beta)
% Mattia Poinelli
% TU Delft
% March 2016
% 
% acos2 function implemented as reported in pg 792 OCMD
%
% INPUT
% 
% cos_alpha: cos(alpha) defined from -1<cos_alpha<1
% beta: angle for relate the hemisphere function  [rad]
%
% OUTPUT
% 
% alpha: evaluated angle mod(2*pi) [rad]


% input check
if cos_alpha < -1 || cos_alpha > 1
    error ('Warning cosine module between [-1,1]')
end

beta = mod (beta,2*pi);

% hemisphere function
H = [];

if beta  <= pi
    H = 1;
elseif beta > pi
    H = -1;
end

alpha = H*acos(cos_alpha);
alpha = mod (alpha,2*pi);

end
