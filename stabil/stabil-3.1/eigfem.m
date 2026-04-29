function [phi,omega]=eigfem(K,M,nMode)

%EIGFEM  Compute the eigenmodes and eigenfrequencies of the finite element model.
%
%   [phi,omega]=eigfem(K,M,nMode)
%   [phi,omega]=eigfem(K,M)
%   computes the eigenmodes and eigenfrequencies of the finite element model.
%
%   K          Stiffness matrix (nDOF * nDOF)
%   M          Mass matrix (nDOF * nDOF)
%   nMode      Number of eigenmodes and eigenfrequencies (default: all)
%   phi        Eigenmodes (in columns) (nDOF * nMode)
%   omega      Eigenfrequencies [rad/s] (nMode * 1)

% David Dooms
% March 2008

% Compute eigenfrequencies and mode shapes
if (nargin < 3) || (nMode==length(K))
    [phi,lambda]=eig(full(K),full(M));
    omega=sqrt(diag(lambda));
    [omega,index]=sort(omega);
    phi=phi(:,index);
else
    try                                                    
        [phi,lambda]=eigs(K,M,nMode,'sm');
        [omega,index]=sort(sqrt(diag(lambda)));
        phi=phi(:,index);
    catch                                                   
        [phi,lambda]=eig(full(K),full(M));
        omega=sqrt(diag(lambda));
        [omega,index]=sort(omega);
        phi=phi(:,index);
        if nargin>2
            omega=omega(1:nMode);
            phi=phi(:,1:nMode); 
        end
    end 
end
% Normalize mode shapes
n=sqrt(diag(phi.'*M*phi));                          
phi=phi./n.';   
end
