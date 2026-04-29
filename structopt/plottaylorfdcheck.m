function hplt = plottaylorfdcheck(h,r0,r)

%PLOTTAYLORFDCHECK   Plot Taylor finite difference check results.
%
%   hplt = PLOTTAYLORFDCHECK(h,r0,r) plots the Taylor remainder
%   r_i = |f_i(x+hp) - f_i(x) - h(df_i/dx)p| as a function of the step size h, 
%   together with a reference O(h^2) curve.
%
%   INPUT ARGUMENTS
%
%   h         Vector of step sizes (1 * nFDstep).
%   r0        Taylor remainder for the objective (1 * nFDstep).
%   r         Taylor remainder for the constraints (m * nFDstep).
%
%   OUTPUT ARGUMENTS
%
%   hplt      Structure containing handles to the generated plot objects:
%               hplt.axes   Handle to the current axes (GCA).
%               hplt.oh2    Handle to the reference O(h^2) curve.
%               hplt.f0     Handle to the objective remainder curve.
%               hplt.f      Handles to the constraint remainder curves.
%               hplt.label  Handle to the 'O(h^2)' text label.

% Mattias Schevenels
% January 2026

h = h(:).';         

C  = max([r0(1); r(:,1)]);                       
Oh2 = (h/h(1)).^2*C;

hplt.axes = gca;

hplt.oh2 = loglog(h,Oh2,'--','Color',[0.7 0.7 0.7],'LineWidth',2);
hold('on');
hplt.f0 = loglog(h,r0,'k');
hplt.f = loglog(h,r);

set(hplt.axes,'XDir','reverse','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off');

xlabel('$h$','Interpreter','latex');
ylabel('$\left|f_i(\mathbf{x}+h\mathbf{p})-f_i(\mathbf{x})-h\left(\frac{\partial f_i}{\partial\mathbf{x}}\right)\mathbf{p}\right|$','Interpreter','latex');

hplt.label = text(h(end),Oh2(end),'\ $O(h^2)$','Interpreter','latex','FontSize',hplt.axes.FontSize*hplt.axes.LabelFontSizeMultiplier);

if nargout==0
  clear('hplt');
end