function DispScal=animdisp(Nodes,Elements,Types,DOF,U,varargin)

%ANIMDISP   Animate the displacements.
%
%   DispScal=animdisp(Nodes,Elements,Types,DOF,U)
%   animates the displacements. 
%
%   Nodes      Node definitions        [NodID x y z]
%   Elements   Element definitions     [EltID TypID SecID MatID n1 n2 ...]
%   Types      Element type definitions  {TypID EltName Option1 ... }
%   DOF        Degrees of freedom  (nDOF * 1)
%   U          Displacements (nDOF * nSteps)
%   DispScal   Displacement scaling
%
%   ANIMDISP(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%   'DispScal'     Displacement scaling.  Default: 'auto'.
%   'Handle'       Plots in the axis with this handle.  Default: current axis.
%   'Fps'          Frames per second.  Default: 12.
%   'CreateMovie'  Saves the movie in the userdata of the axis of the figure.  
%                  Use getmovie to get the movie from the axis. Default: 'off'.
%   'Counter'      Displays the number of the frame for transient displacements.
%                  Default: 'on'.
%   Additional parameters are redirected to the PLOTDISP function which plots 
%   the individual frames of the movie.
%
%   See also GETMOVIE, PLOTDISP.

% David Dooms
% March 2008

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
  warning('ANIMDISP not supported by Octave.');
  return
end									   
% PARAMETERS
[DispScal,paramlist]=cutparam('DispScal','auto',varargin);
[DispMax,paramlist]=cutparam('DispMax','off',paramlist);
[haxis,paramlist]=cutparam('Handle',gca,paramlist);
[Fps,paramlist]=cutparam('Fps',12,paramlist);
[Counter,paramlist]=cutparam('Counter','on',paramlist);
[CreateMovie,paramlist]=cutparam('CreateMovie','off',paramlist);

% PREPARE FIGURE
drawnow;

if strcmpi(DispScal,'auto')
    [DispScal,h,Ax,Ay,Az,B]=plotdisp(Nodes,Elements,Types,DOF,max(abs(U),[],2),paramlist{:},'Handle',haxis,'DispMax',DispMax);
else
    [DispScal,h,Ax,Ay,Az,B]=plotdisp(Nodes,Elements,Types,DOF,max(abs(U),[],2),paramlist{:},'DispScal',DispScal,'Handle',haxis,'DispMax',DispMax);
end

limaxesX(1,:)=get(haxis,'XLim');
limaxesY(1,:)=get(haxis,'YLim');
limaxesZ(1,:)=get(haxis,'ZLim');

X=-DispScal*Ax*max(abs(U),[],2)+B(:,1);
Y=-DispScal*Ay*max(abs(U),[],2)+B(:,2);
Z=-DispScal*Az*max(abs(U),[],2)+B(:,3);
set(h.displacements,'XData',X,'YData',Y,'ZData',Z);

limaxesX(2,:)=get(haxis,'XLim');
limaxesY(2,:)=get(haxis,'YLim');
limaxesZ(2,:)=get(haxis,'ZLim');

limaxesX=[min(limaxesX(:,1)) max(limaxesX(:,2))];
limaxesY=[min(limaxesY(:,1)) max(limaxesY(:,2))];
limaxesZ=[min(limaxesZ(:,1)) max(limaxesZ(:,2))];

set(haxis,'XLim',limaxesX,'YLim',limaxesY,'ZLim',limaxesZ);
set(haxis,'XLimMode','manual','YLimMode','manual','ZLimMode','manual');

hlabel=annotation('TextBox','Position',[0 0 1 1],'LineStyle','none','HorizontalAlignment','center','VerticalAlignment','bottom','String','0','FontSize',get(0,'defaulttextfontsize'),'Visible','off');

if strcmpi(get(haxis,'NextPlot'),'replace')
    set(get(haxis,'Parent'),'NextPlot','replace');
end

% INITIALIZE TIMER
if size(U,2)==1
    takstoexecute=Inf;
else
    takstoexecute=size(U,2)-1;
    if strcmpi(Counter,'on')
        set(hlabel,'Visible','on');
    end
end

t=timer('ExecutionMode','fixedRate','Period',round(1000/Fps)/1000,'TasksToExecute',takstoexecute,'BusyMode','queue');
set(haxis,'DeleteFcn',{@StopTimer,t});

% PLOT FIRST FRAME
tic;
event=struct([]);
TimerAnim(t,event,haxis,h.displacements,Ax,Ay,Az,B,U,DispScal,hlabel,0)

% INITIALIZE MOVIE
if strcmpi(CreateMovie,'on')
    hmov=get(haxis,'Parent');
    UserData.mov=getframe(hmov);
    UserData.busy=true;
    set(haxis,'UserData',UserData);
else
    hmov=0;
end

% DEFINE TIMER
t.TimerFcn={@TimerAnim,haxis,h.displacements,Ax,Ay,Az,B,U,DispScal,hlabel,hmov};
t.StopFcn={@StopAnim,haxis};
elap=toc;
if elap < 1/Fps
    remain=round(1000*(1/Fps-elap))/1000;
else
    remain=0;
end
t.StartDelay=max([0.01; remain]);

start(t);

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout < 1
    clear('DispScal');
end

end

%-------------------------------------------------------------------------------
% IF AXIS ARE DELETED
function StopTimer(obj,event,t);

stop(t);

end

%-------------------------------------------------------------------------------
% TIMER FUNCTION
function TimerAnim(obj,event,haxis,hdispl,Ax,Ay,Az,B,U,DispScal,hlabel,hmov);

iFrame=get(obj,'TasksExecuted');
try
if size(U,2)==1                           % mode shapes or complex displacements
    if (hmov~=0) && (iFrame<41)
        UserData=get(haxis,'UserData');
        UserData.mov(iFrame)=getframe(hmov);
        if iFrame==40
            UserData.busy=false;
        end
        set(haxis,'UserData',UserData);
    end
    iFrame=rem(iFrame,40);
    Uplot=real(U*exp(i*(iFrame*pi/20)));
    X=DispScal*Ax*Uplot+B(:,1);
    Y=DispScal*Ay*Uplot+B(:,2);
    Z=DispScal*Az*Uplot+B(:,3);
    set(hdispl,'XData',X,'YData',Y,'ZData',Z);

else                                      % transient displacements
    if hmov~=0
        UserData=get(haxis,'UserData');
        UserData.mov(iFrame)=getframe(hmov);
        set(haxis,'UserData',UserData);
    end
    X=DispScal*Ax*U(:,iFrame+1)+B(:,1);
    Y=DispScal*Ay*U(:,iFrame+1)+B(:,2);
    Z=DispScal*Az*U(:,iFrame+1)+B(:,3);
    set(hdispl,'XData',X,'YData',Y,'ZData',Z);
    set(hlabel,'String',sprintf('%d/%d',iFrame+1,size(U,2)));
    if (hmov~=0) && (iFrame==(size(U,2)-1))
        UserData=get(haxis,'UserData');
        UserData.mov(iFrame+1)=getframe(hmov);
        UserData.busy=false;
        set(haxis,'UserData',UserData);
    end
end
end

end

%-------------------------------------------------------------------------------
% IF TIMER IS STOPPED
function StopAnim(obj,event,haxis);

if ishandle(haxis)
    set(haxis,'DeleteFcn',{});
end

delete(obj)

end