function lref=reflength(Nodes)

deltanodes=Nodes(:,2:4)-repmat(mean(Nodes(:,2:4)),size(Nodes,1),1);
for iNode=1:size(Nodes,1)
    distance(iNode,:)=norm(deltanodes(iNode,:));
end
indnodes=find(distance <= (mean(distance)+1.84*std(distance)));
lref=0.10*mean(max(Nodes(indnodes,2:4))-min(Nodes(indnodes,2:4)));
