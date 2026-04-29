function Scontour = scontour_shell4(Node,Stress,Svalues)

%SCONTOUR_SHELL4   Matrix to plot contours in a shell4 element.
%
%   Scontour = scontour_shell4(Node,Stress,Svalues) returns a matrix used
%   for plotting contours in a shell4 element.
%
%   Node       Node definitions           [x y z] (4 * 3)
%   Stress     Stress in nodes                    (4 * 1)
%   Svalues    Values of contours          (nContour * 1)
%   Scontour   Coordinates of contours in GCS
%
%   See also PLOTSTRESSCONTOUR, PLOTSHELLFCONTOUR.

% Miche Jansen
% 2010
    

nPoints = 5;
[xi,eta] = meshgrid(linspace(-1,1,nPoints),linspace(-1,1,nPoints));

sigma = Stress(:);

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
       N = sh_qs4(xi(indi),eta(indi));
       xyz(indi,:) = N.'*Node;
   end
   Scontour = [Scontour;C(1,teller) num 0; xyz ];
   teller = teller + num+1;
end

end
