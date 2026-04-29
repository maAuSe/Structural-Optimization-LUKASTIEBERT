function [theta] = vtrans_solid(t,vtype)
%VTRANS_SOLID   Transformation matrix for stress and strain components in matrix (Voigt) notation.
%
%   [theta] = vtrans_solid(t)
%   [theta] = vtrans_solid(t,vtype)
%   computes the transformation matrix between the local and the global
%   coordinate system for t stress or strain vector in matrix notation.
%
%   t          Transformation matrix                (3 * 3)
%   vtype      Vector type 'stress' (default) | 'strain'               
%   theta      Stress transformation matrix         (6 * 6)
%
%   See also TRANS_SOLID8, TRANS_SOLID20.

% Miche Jansen
% 2013

if nargin < 2; vtype = 'stress'; end

% transformatie matrix uit Cook p274 en p194 Zienkiewicz deel1 (6th ed.)
switch vtype
    case 'stress' % Tsigma
theta =[t(:,:).^2 2*t(:,1).*t(:,2) 2*t(:,2).*t(:,3) 2*t(:,3).*t(:,1);
        t(1,:).*t(2,:) t(1,1)*t(2,2)+t(1,2)*t(2,1) t(1,2)*t(2,3)+t(1,3)*t(2,2) t(1,3)*t(2,1)+t(1,1)*t(2,3);
        t(2,:).*t(3,:) t(2,1)*t(3,2)+t(2,2)*t(3,1) t(2,2)*t(3,3)+t(2,3)*t(3,2) t(2,3)*t(3,1)+t(2,1)*t(3,3);
        t(3,:).*t(1,:) t(3,1)*t(1,2)+t(3,2)*t(1,1) t(3,2)*t(1,3)+t(3,3)*t(1,2) t(3,3)*t(1,1)+t(3,1)*t(1,3)];
    
    case 'strain' % Tepsilon
theta = [t(:,:).^2 t(:,1).*t(:,2) t(:,2).*t(:,3) t(:,3).*t(:,1);
          2*t(1,:).*t(2,:) t(1,1)*t(2,2)+t(1,2)*t(2,1) t(1,2)*t(2,3)+t(1,3)*t(2,2) t(1,3)*t(2,1)+t(1,1)*t(2,3);
          2*t(2,:).*t(3,:) t(2,1)*t(3,2)+t(2,2)*t(3,1) t(2,2)*t(3,3)+t(2,3)*t(3,2) t(2,3)*t(3,1)+t(2,1)*t(3,3);
          2*t(3,:).*t(1,:) t(3,1)*t(1,2)+t(3,2)*t(1,1) t(3,2)*t(1,3)+t(3,3)*t(1,2) t(3,3)*t(1,1)+t(3,1)*t(1,3)];
end