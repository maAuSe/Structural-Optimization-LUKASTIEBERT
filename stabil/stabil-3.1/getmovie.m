function mov=getmovie(h)

%GETMOVIE   Get the movie from a figure where an animation has been played.     
%
%   mov=getmovie(h)
%   gets the movie from the userdata of the axis of a figure where an animation
%   has been played. In order to save the movie in the userdata the animation 
%   should have been played using animdisp(...,'CreateMovie','on'). This 
%   function blocks the command prompt until the movie has become available. 
%
%   h          Axis handle.
%   mov        Structured array with movie frames.
%
%   The movie can be played with movie(gcf,mov).
%
%   See also ANIMDISP, MOVIE.

% David Dooms
% March 2008

UserData=get(h,'UserData');

if ~isfield(UserData,'mov') 
    error('There is no movie attached to this figure. Use animdisp(...,''CreateMovie'',''on'') to attach the movie.') ;
end

while UserData.busy
    pause(1);
    UserData=get(h,'UserData');
end

mov=UserData.mov;
