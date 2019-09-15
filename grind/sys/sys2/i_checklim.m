function lim=i_checklim(oldlim,s)
lim=i_checkstr(s);
if isempty(lim)||any(isnan(lim))
    if ~isempty(oldlim)
        lim=oldlim;
    else
        lim=[0 10];
    end

end

if (length(lim)==1)&&(lim>0)
    lim=[0 lim];
end

if (length(lim)==1)&&(lim==0)
    lim=oldlim;
end

if lim(1)>lim(2)
    lim=[lim(2) lim(1)];
end

