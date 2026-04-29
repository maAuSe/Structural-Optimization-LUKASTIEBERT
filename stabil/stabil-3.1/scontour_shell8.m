function Scontour = scontour_shell8(Node,Stress,Svalues)

%SCONTOUR_SHELL8   Matrix to plot contours in a shell8 element.
%
%   Scontour = scontour_shell8(Node,Stress,Svalues) returns a matrix used
%   for plotting contours in a shell4 element.
%
%   Node       Node definitions           [x y z] (4 * 3)
%   Stress     Stress in nodes         (4 * 1) or (8 * 1)
%   Svalues    Values of contours          (nContour * 1)
%   Scontour   Coordinates of contours in GCS
%
%   See also PLOTSTRESSCONTOUR, PLOTSHELLFCONTOUR.

% Miche Jansen
% 2010
    


nPoints = 5;
[xi,eta] = meshgrid(linspace(-1,1,nPoints),linspace(-1,1,nPoints));

sigma = Stress(:);
if length(sigma)<=4
%coef bevat de coefficienten van de bilineaire interpolatie
coef = 0.25*[1 1 1 1; -1 1 1 -1;-1 -1 1 1;1 -1 1 -1]*sigma;

C = contourc(xi(1,:),eta(:,1),coef(1)+coef(2)*xi+coef(3)*eta+coef(4)*xi.*eta,Svalues);
teller = 1;
Scontour = [];
while teller < size(C,2)
   num = C(2,teller);
   xi = C(1,teller+1:teller+num);
   eta= C(2,teller+1:teller+num);
   xyz = zeros(length(xi),3);
   for indi=1:length(xi);
       N = sh_qs8(xi(indi),eta(indi));
       xyz(indi,:) = N.'*Node;
   end
   Scontour = [Scontour;C(1,teller) num 0; xyz ];
   teller = teller + num+1;
end

else % indien nodal solution wordt doorgegeven (in de 8 knopen)
%coef bevat de coefficienten van de interpolatie
coef=0.25*[-1 -1 -1 -1  2  2  2  2;
            0  0  0  0  0  2  0 -2;
            0  0  0  0 -2  0  2  0;
            1 -1  1 -1  0  0  0  0;
            1  1  1  1 -2  0 -2  0;
            1  1  1  1  0 -2  0 -2;
           -1 -1  1  1  2  0 -2  0;
           -1  1  1 -1  0 -2  0  2]*sigma;

C = contourc(xi(1,:),eta(:,1),coef(1)+coef(2)*xi+coef(3)*eta...
    +coef(4)*xi.*eta+coef(5)*xi.^2+coef(6)*eta.^2+coef(7)*xi.^2.*eta+coef(8)*xi.*eta.^2,Svalues);
teller = 1;
Scontour = [];
while teller < size(C,2)
   num = C(2,teller);
   xi = C(1,teller+1:teller+num);
   eta= C(2,teller+1:teller+num);
   xyz = zeros(length(xi),3);
   for indi=1:length(xi);
       N = sh_qs4(xi(indi),eta(indi));
       xyz(indi,:) = N.'*Node;
   end
   Scontour = [Scontour;C(1,teller) num 0; xyz ];
   teller = teller + num+1;
end
end

end

   