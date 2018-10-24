function [A,B,c] = Side_Angle_Side (a,b)
% Mattia Poinelli
% JPL
% April 2017
%
% Solution of the spherical triangle given two side angle (a,b) and a
% rect 'inner' rotation angles C. Side-Angle-Side problem or Airplane Problem.
% First option is available (short distances).
% Figure is referenced to slide 8, R. Noomen, FullSkyGeometry, Version 4-10,
% Equations taken from Appendix A, Wertz, OCDM pg. 784 
% 
% 
% INPUT
% a: side angle [rad]
% b: side angle  [rad]
% C: 'inner' rotation angle [rad]
% 
% OUTPUT
% A: 'inner' rotation angle [rad] 
% B: 'inner' rotation angle [rad] 
% c: side angle [rad] 

c = acos2 (cos(a)*cos(b),pi/2);
A = acos2 ((cos(a)-cos(b)*cos(c))/(sin(b)*sin(c)),a);
B = acos2 ((cos(b)-cos(a)*cos(c))/(sin(a)*sin(c)),b);

end