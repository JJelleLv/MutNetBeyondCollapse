% 
%function [data, y] = coco_test_HH(prob, data, u) %#ok<INUSL>
%test function for Hopf-Hopf bifurcation (co-dim 2)
%
function [data,y] = coco_test_HH(~, data, ~)
y=NaN;
if data.xdim>=4
    la=eig(data.ep_Jx);
    [~,ndx]=sort(abs(real(la)));
    leadla=la(ndx(1:4));
    if sum(imag(leadla)==0)==0 %all 4 imaginary
        y=real(leadla(4)-leadla(1));%not really continuous test
    end

end
