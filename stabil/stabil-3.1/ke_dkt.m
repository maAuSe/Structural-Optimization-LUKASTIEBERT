function [Ke,Me] = ke_dkt(Node,h,E,nu,rho)

%KE_DKT   DKT plate element stiffness and mass matrix.
%
%   [Ke,Me] = ke_dkt(Node,h,E,nu,rho)
%    Ke     = ke_dkt(Node,h,E,nu) 
%   returns the element stiffness and mass matrix in the global coordinate system 
%   for a three node plate element (isotropic material) in the xy-plane.
%
%   Node       Node definitions                              [x y] (3 * 2)
%   h          Plate thickness (uniform or defined in nodes) [h]/[h1 h2 h3]
%   E          Young's modulus
%   nu         Poisson coefficient
%   rho        Mass density          
%   Ke         Element stiffness matrix (24 * 24)
%   Me         Element mass matrix (24 * 24)
%
%   This element is a Discrete Kirchoff Triangle plate element. This code 
%   is a MATLAB code based on the E-2 element FORTRAN coding, described in:
%   Construction of new efficient three-node triangular thin plate bending
%   elements, C. Jeyachandrabose and J. Kirkhope, Computers & Structures
%   Vol. 23, No. 5, pp. 587-603, 1986.
%
%   See also Q_DKT, SH_T, KE_DKT4.




D = E/(12*(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

IT =[420 140 140 70 35 70;
    420 210 105 126 42 42;
    420 105 210 42 42 126;
    210 126 42 84 21 14;
    210 84 84 42 28 42;
    210 42 126 14 21 84;
    42 28 7 20 4 2;
    42 21 14 12 6 6;
    42 14 21 6 6 12;
    42 7 28 2 4 20].';


% variabele dikte?

b = [Node(2,2)-Node(3,2);
     Node(3,2)-Node(1,2);
     Node(1,2)-Node(2,2)];
c = [Node(3,1)-Node(2,1);
     Node(1,1)-Node(3,1);
     Node(2,1)-Node(1,1)];

deter = b(2)*c(3)-b(3)*c(2);


if length(h) > 1
h21 = h(2)-h(1);
h31 = h(3)-h(1);
else
h21 = 0;
h31 = 0;
end
gamma = [h(1)^3
         h(1)*h(1)*h21
         h(1)*h(1)*h31
         h(1)*h21*h21
         h(1)*h21*h31
         h(1)*h31*h31
         h21^3
         h21*h21*h31
         h21*h31*h31
         h31^3];

     
% E matrix
jj =1;
if length(h) > 1
    jj=10;
if h(1) == h(2) && h(1) == h(3)
    jj=1;
end
end

XX = zeros(3,3);
E = zeros(9,9); 

k=1;
for i=1:3
    for j=i:3
        for l=1:jj
            XX(i,j)=XX(i,j)+IT(k,l)*gamma(l);
        end
            k=k+1;
            XX(j,i)=XX(i,j);
    end
end

for i=1:3
    for j=i:3
        for k1=1:3
            ii = (i-1)*3+k1;
            for k2=1:3
                jj=(j-1)*3+k2;
                E(ii,jj)=D(i,j)*XX(k1,k2)/deter/840;
                E(jj,ii)=E(ii,jj);
            end
        end
    end
end
% Q matrix

Q = q_dkt(b,c,deter);

Ke = Q.'*E*Q;
%Ke = T.'*Q.'*E*Q*T;

if nargout > 1
Me = zeros(9,9);
% gausspunten uit ansys theory
x = 1/6*[4 1 1;1 4 1;1 1 4];
H = [1 1 1]/3;
for iGauss=1:3
    Ni = sh_t(x(:,iGauss),b,c);
    Me = Me + Ni.'*rho*Ni*H(iGauss)*h(1)*deter/2;
end
end
end


        
        




    
    
    

