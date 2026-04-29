function [SeGCS,SeLCS,vLCS] = se_truss(Node,Section,Material,UeGCS,Options,gcs)

SeGCS = nan(72,size(UeGCS,2));
SeLCS = nan(72,size(UeGCS,2));
if nargout > 2
    t = trans_truss(Node);
    vLCS = t.';
    vLCS = vLCS(:).';
end

end