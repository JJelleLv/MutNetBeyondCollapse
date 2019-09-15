% 
%function [data, y] = coco_test_HH(prob, data, u) %#ok<INUSL>
%test function for Hopf-Hopf bifurcation (co-dim 2)
%
function [data,y] = coco_test_FH(~, data, ~)
y=NaN;
if data.xdim>=3
    la=eig(data.ep_Jx);
    [~,ndx]=sort(abs(real(la)));
    leadla=la(ndx(1:3));
    if sum(imag(leadla)==0)==1
        y=real(leadla(3)-leadla(1));%not really continuous test
    end

end
