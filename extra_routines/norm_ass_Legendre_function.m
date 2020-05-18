function [P] = norm_ass_Legendre_function (l,t)
% Mattia Poinelli
% TU Delft
% 21 May 2016
% 
% Computation of normalized associated Legendre functions evaluated at a
% certain value t.
%
% INPUT 
% l: degree (0 order is estimated as well)
% t: 'integration' point
% 
% OUTPUT
% P: matrix containing normalized associated Legendre functions, i.e. a
%    sparse matrix with the diagonal the values of P_{l,l} on the right side
%    of the diagonal the other terms for decreasing order m = l-1 ,...0


% check for order 0 or errors
if l < 0
    error ('Degree has to be 0 or higher')
elseif l == 0
    P = 1;
    return
end

% Initialization of P and u (t)
P = zeros(l,2*l); % size [l X 2l]
u = sqrt(1 - t^2);

for i = 1 : l % per degree
    
    % 'diagonal' recursion
    if i == 1
    P(i,i)   = sqrt(3) * u;
    else
    P(i,i) = u * sqrt((2*i+1)/(2*i)) * P(i-1,i-1);
    end
    
    % initialization of parameter m 
    m = i - 1;
    
    % recursion per degree m
    for j = i + 1 : 2 * i 
             
        % evaluation of coefficients a and b
            
        if m == 0 % order 0, i.e. last term
           
            a = 2 / sqrt (2*i*(i+1));
            b = sqrt( ((i+2)*(i-1)) / (2*i*(i+1)) );
    
        else % all other orders
            
            a = 2 * ( m + 1 ) / sqrt( ( i - m ) * ( i + m + 1 ) );
            b = sqrt( ( i + m + 2) * ( i - m - 1 )) / sqrt( ( i - m ) * ( i + m + 1 ) );

        end
        
        if j == i + 1
        
            P(i,j) = a * ( t / u ) * P(i , j - 1);
            
        
        else
             
             P(i,j) = a * ( t / u ) * P(i , j - 1) - b * P (i ,j - 2);
         
        end
        
        clear a b 
        
        % decreasing m until order 0
        m = m - 1;
        
    end
end
