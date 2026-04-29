function h=plotgcs(lref,h)

h.xaxis=plot3([0 lref],[0 0],[0 0],'r','Clipping','off');
h.xaxislabel=text(lref,0,0,'x','HorizontalAlignment','center','VerticalAlignment','bottom','Color',[1 0 0]);
h.yaxis=plot3([0 0],[0 lref],[0 0],'g','Clipping','off');
h.yaxislabel=text(0,lref,0,'y','HorizontalAlignment','center','VerticalAlignment','bottom','Color',[0 1 0]);
h.zaxis=plot3([0 0],[0 0],[0 lref],'b','Clipping','off');
h.zaxislabel=text(0,0,lref,'z','HorizontalAlignment','center','VerticalAlignment','bottom','Color',[0 0 1]);
