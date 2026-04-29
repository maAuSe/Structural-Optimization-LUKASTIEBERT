function [SdiagrLCS,loc,Extrema] = sdiagrlcs_beam(stype,Forces,DLoadLCS,L,Section,Points)

%SDIAGRLCS_BEAM   Stress diagram for a beam element in LCS.
%
%   [SdiagrLCS,loc,Extrema] = sdiagrlcs_beam(ftype,Forces,DLoadLCS,L,Points) 
%   computes the stresses at the specified points. The extreme values are
%   analytically determined.
%
%   stype      'snorm'      Normal stress due to normal force
%              'smomyt'     Normal stress due to bending moment  
%                           around the local y-direction at the top
%              'smomyb'     Normal stress due to bending moment  
%                           around the local y-direction at the bottom
%              'smomzt'     Normal stress due to bending moment  
%                           around the local z-direction at the top 
%              'smomzb'     Normal stress due to bending moment  
%                           around the local z-direction at the bottom
%              'smax'       Maximal normal stress (normal force and bending moment)
%              'smin'       Minimal normal stress (normal force and bending moment)
%   Forces     Element forces in LCS (beam convention) [N; Vy; Vz; T; My; Mz](12 * 1)
%   DLoadLCS   Distributed loads in LCS [n1localX; n1localY; n1localZ; ...](6 * 1)
%   Points     Points in the local coordinate system (1 * nPoints)
%   Section    Section definition         [A ky kz Ixx Iyy Izz yt yb zt zb]
%   SdiagrLCS  Stresses at the points (1 * nPoints)
%   loc        Locations of the extreme values (nValues * 1)
%   Extrema    Extreme values (nValues * 1)
%
%   See also SDIAGRGCS_BEAM.

% David Dooms
% November 2008

% PREPROCESSING
Points=Points(:).';

switch lower(stype)
   case 'snorm'
      A=[(DLoadLCS(1)-DLoadLCS(4))*L/2;  -DLoadLCS(1)*L;   Forces(1)]/Section(1);
   case 'smomyt'
      A=[(DLoadLCS(3)-DLoadLCS(6))*L^2/6;  -DLoadLCS(3)*L^2/2;   Forces(3)*L;  Forces(5)]/(Section(5)+1e-250)*Section(9);
   case 'smomyb'
      A=-[(DLoadLCS(3)-DLoadLCS(6))*L^2/6;  -DLoadLCS(3)*L^2/2;   Forces(3)*L;  Forces(5)]/(Section(5)+1e-250)*Section(10);
   case 'smomzt'
      A=[(DLoadLCS(2)-DLoadLCS(5))*L^2/6;  -DLoadLCS(2)*L^2/2;   Forces(2)*L;  Forces(6)]/(Section(6)+1e-250)*Section(7);
   case 'smomzb'
      A=-[(DLoadLCS(2)-DLoadLCS(5))*L^2/6;  -DLoadLCS(2)*L^2/2;   Forces(2)*L;  Forces(6)]/(Section(6)+1e-250)*Section(8);
   case {'smax','smin'}
      A(:,1)=[0; (DLoadLCS(1)-DLoadLCS(4))*L/2;  -DLoadLCS(1)*L;   Forces(1)]/Section(1) ...
        +[(DLoadLCS(3)-DLoadLCS(6))*L^2/6;  -DLoadLCS(3)*L^2/2;   Forces(3)*L;  Forces(5)]/(Section(5)+1e-250)*Section(9) ...
        +[(DLoadLCS(2)-DLoadLCS(5))*L^2/6;  -DLoadLCS(2)*L^2/2;   Forces(2)*L;  Forces(6)]/(Section(6)+1e-250)*Section(7);
      A(:,2)=[0; (DLoadLCS(1)-DLoadLCS(4))*L/2;  -DLoadLCS(1)*L;   Forces(1)]/Section(1) ...
        +[(DLoadLCS(3)-DLoadLCS(6))*L^2/6;  -DLoadLCS(3)*L^2/2;   Forces(3)*L;  Forces(5)]/(Section(5)+1e-250)*Section(9) ...
        -[(DLoadLCS(2)-DLoadLCS(5))*L^2/6;  -DLoadLCS(2)*L^2/2;   Forces(2)*L;  Forces(6)]/(Section(6)+1e-250)*Section(8);
      A(:,3)=[0; (DLoadLCS(1)-DLoadLCS(4))*L/2;  -DLoadLCS(1)*L;   Forces(1)]/Section(1) ...
        -[(DLoadLCS(3)-DLoadLCS(6))*L^2/6;  -DLoadLCS(3)*L^2/2;   Forces(3)*L;  Forces(5)]/(Section(5)+1e-250)*Section(10) ...
        +[(DLoadLCS(2)-DLoadLCS(5))*L^2/6;  -DLoadLCS(2)*L^2/2;   Forces(2)*L;  Forces(6)]/(Section(6)+1e-250)*Section(7);
      A(:,4)=[0; (DLoadLCS(1)-DLoadLCS(4))*L/2;  -DLoadLCS(1)*L;   Forces(1)]/Section(1) ...
        -[(DLoadLCS(3)-DLoadLCS(6))*L^2/6;  -DLoadLCS(3)*L^2/2;   Forces(3)*L;  Forces(5)]/(Section(5)+1e-250)*Section(10) ...
        -[(DLoadLCS(2)-DLoadLCS(5))*L^2/6;  -DLoadLCS(2)*L^2/2;   Forces(2)*L;  Forces(6)]/(Section(6)+1e-250)*Section(8);
      for iCombin=1:4
          SdiagrLCS(iCombin,:)=polyval(A(:,iCombin),Points);
      end
   otherwise
      error('Unknown stress.')
end


switch lower(stype)
   case {'snorm','smomyt','smomyb','smomzt','smomzb'}
      SdiagrLCS=polyval(A,Points);
      
      if (length(A)-find(A~=0,1))>1
         loc=roots(polyder(A));
         loc=[0; 1; loc(loc>0 & loc<1)];   
      else
         loc=[0; 1];
      end
      
      Extrema=polyval(A,loc);     
   case 'smax'
      [SdiagrLCS,I]=max(SdiagrLCS);
      loc=[];
      Extrema=[];
      for iCombin=1:4
          if any(I==iCombin)
              if (length(A(:,iCombin))-find(A(:,iCombin)~=0,1))>1
                  temploc=roots(polyder(A(:,iCombin)));
                  temploc=temploc(temploc>0 & temploc<1);
                  if length(temploc)>0
                  iright=find(Points>temploc(1),1);
                  if I(iright)==iCombin && I(iright-1)==iCombin
                      loc=[loc; temploc(1)];
                      Extrema=[Extrema; polyval(A(:,iCombin),temploc(1))]; 
                  end
                  end
                  if length(temploc)==2
                  iright=find(Points>temploc(2),1);
                  if I(iright)==iCombin && I(iright-1)==iCombin
                      loc=[loc; temploc(2)];
                      Extrema=[Extrema; polyval(A(:,iCombin),temploc(2))]; 
                  end
                  end
              end
          end
      end
      loc=[0; 1; loc];
      Extrema=[SdiagrLCS(1,1); SdiagrLCS(1,end); Extrema];
   case 'smin'
      [SdiagrLCS,I]=min(SdiagrLCS);
      loc=[];
      Extrema=[];
      for iCombin=1:4
          if any(I==iCombin)
              if (length(A(:,iCombin))-find(A(:,iCombin)~=0,1))>1
                  temploc=roots(polyder(A(:,iCombin)));
                  temploc=temploc(temploc>0 & temploc<1);
                  if length(temploc)>0
                  iright=find(Points>temploc(1),1);
                  if I(iright)==iCombin && I(iright-1)==iCombin
                      loc=[loc; temploc(1)];
                      Extrema=[Extrema; polyval(A(:,iCombin),temploc(1))]; 
                  end
                  end
                  if length(temploc)==2
                  iright=find(Points>temploc(2),1);
                  if I(iright)==iCombin && I(iright-1)==iCombin
                      loc=[loc; temploc(2)];
                      Extrema=[Extrema; polyval(A(:,iCombin),temploc(2))]; 
                  end
                  end
              end
          end
      end
      loc=[0; 1; loc];
      Extrema=[SdiagrLCS(1,1); SdiagrLCS(1,end); Extrema];
end

