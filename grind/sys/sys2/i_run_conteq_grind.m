function  i_run_conteq_grind(obj)
global g_grind;
settings=obj.eval_settings;
N0=settings.N0;
i_keep(N0);
apar=settings.allpars{settings.freepars(1)};
currval=evalin('base',apar);
if isempty(obj.settings.parranges{1})
    settings.parranges{1}=[-Inf Inf];
end

parrange=settings.parranges{1};
maxstep=settings.coco.cont.h_max;
minstep=settings.coco.cont.h_min;
maxpnts=settings.coco.cont.PtMX;

if isnan(parrange(1))||parrange(1)<currval
    [cont]=continueeq(apar, currval, parrange(1), maxstep, minstep,maxpnts);
else
    cont.Y=[];
    cont.eig=[];
    cont.parvalues=[];
end
if isnan(parrange(2))||parrange(2)>currval
    [cont1]=continueeq(apar, currval, parrange(2),maxstep,minstep,maxpnts);
    if ~isoctave&&verLessThan('matlab','8.4.0') %R2014b
        %irritating that flipdim is renamed..
        cont.Y = cat(3,flipdim(cont.Y,3), cont1.Y);   %#ok<DFLIPDIM>
        cont.eig = cat(3,flipdim(cont.eig,3), cont1.eig);  %#ok<DFLIPDIM>
    else
        cont.Y = cat(3,flip(cont.Y,3), cont1.Y);
        cont.eig = cat(3,flip(cont.eig,3), cont1.eig);
    end

    cont.parvalues = [flipud(cont.parvalues); cont1.parvalues];
end

obj.run.pars=settings.allpars(settings.freepars(1));
if ~isempty(obj.run.Y)
    obj.run.Y(:, :, end + 1) = NaN;
    obj.run.eig(:, :, end + 1) = NaN;
    obj.run.parvalues(end + 1, :) = NaN;
    lenY = size(obj.run.Y, 3);
    obj.run.Y = cat(3,obj.run.Y,cont.Y);
    obj.run.eig = cat(3, obj.run.eig, cont.eig);
    obj.run.parvalues =  [obj.run.parvalues ; cont.parvalues];
else
    lenY = 0;
    obj.run.Y = cont.Y;
    obj.run.eig  = cont.eig;
    obj.run.parvalues =  cont.parvalues;
end

obj.runs(end).endndx=size(obj.run.parvalues,1);
eig=permute(cont.eig,[3,2,1]);
if g_grind.solver.isdiffer
    [domeig,imax]=max(1-abs(eig),[],2);
else
    [domeig,imax]=max(real(eig),[],2);
end

ibif=find(diff(domeig>0)~=0);
for i=1:length(ibif)
    if abs(domeig(ibif(i)+1))< abs(domeig(ibif(i)))
        j=ibif(i)+1;
    else
        j=ibif(i);
    end

    if abs(imag(eig(j,imax(j))))>0
        %addpoint(obj,labname,labnr,runid,runndx,ic)
        obj.addpoint('HB',[],[],lenY+j);
    else
        obj.addpoint('SN',[],[],lenY+j);
    end

end


i_keep(N0);
end
%% continueq
function [cont, foldexpected] = continueeq(apar, start, nend, maxstep, minstep, maxpnts)
global g_grind;
oldN0 = i_initvar;
if start<nend
    h1=1;
else
    h1=-1;
end
onstart;
[N0, eqfound] = findeq('Display','off');
if ~eqfound
    error('GRIND:conteq:NoEquilibrium','Cannot find equilibrium');
end

%nmax = 200;
dim = length(N0);
cont.par=apar;
cont.Y = NaN * ones(1, dim,maxpnts);
if g_grind.solver.haslags
    cont.eig = NaN * ones(1, 10*dim, maxpnts);
else
    cont.eig = NaN * ones(1, dim, maxpnts);
end
cont.parvalues = zeros(maxpnts, 1);
p = evalin('base', apar);
oldpar = p;
%[N0, eqfound] = findeq(0);

foldexpected = 0;
[~, eigenval] = i_eigen(1,g_grind.solver.iters);
[isstable1, issaddle1, isspiral1] =  i_stability(eigenval, g_grind.solver.isdiffer, g_grind.solver.haslags);
s = 'Continuing ';
if issaddle1 && ~isspiral1
    s = [s 'a saddle equilibrium from %s=%0.5g till %0.5g\n'];
else
    if isstable1
        s = [s 'a stable '];
    else
        s = [s 'an unstable '];
    end

    if isspiral1
        s = [s 'spiral from %s=%0.5g till %0.5g\n'];
    else
        s = [s 'node from %s=%0.5g till %0.5g\n'];
    end

end

fprintf(s, apar, start, nend);
cont.Y(1, :, 1) = transpose(N0);
N1 = N0;
cont.eig(1, :, 1) = transpose(eigenval);
cont.parvalues(1) = p;
i_waitbar(1);
n = 2;
while eqfound && (n < maxpnts) && p~=nend
    onstart;
    p=p+h1*maxstep;
    if p*h1>nend*h1
        p=nend;
    end

    [N0,eqfound,isstable,eigenval] = dostep(apar,p,n,N1);
    if xor(isstable1,isstable)||~eqfound
        i_keep(N1);
        p1=p-h1*maxstep;
        biffound=false;
        N02=N0;
        while ~biffound&&(n < maxpnts)
            p1=p1+h1*minstep;
            [N2,eqfound2,isstable2,eigenval2] = dostep(apar,p1,n,N02);
            if eqfound2
                cont.Y(1, :, n) = transpose(N2);
                cont.eig(1, :, n) = transpose(eigenval2);
                cont.parvalues(n) = p1;
                n = n + 1;
                i_keep(2 * N02 - N2);
                N02 = N2;
                biffound=xor(isstable2,isstable1);
            else
                if n>1
                    cont.Y(1, :, n) = cont.Y(1, :, n-1);
                    eigval=cont.eig(1, :, n-1);
                    if ~g_grind.solver.isdiffer&&isstable1
                        eigval(eigval==max(eigval))=0.01;
                    end

                    cont.eig(1, :, n) = eigval;
                    cont.parvalues(n) = p1;
                    n=n+1;
                end

                biffound=true;
            end

        end

    end

    if eqfound
        cont.Y(1, :, n) = transpose(N0);
        cont.eig(1, :, n) = transpose(eigenval);
        cont.parvalues(n) = p;
        n = n + 1;
        i_keep(2 * N0 - N1);
        N1 = N0;
    end

    %     else
    %         foldexpected = true;
    %         fprintf('Cannot find equilibrium anymore at %s = %0.5g, Fold bifurcation?\n',apar,p);
    % end

end
multassignin('base', apar, oldpar);
i_keep(oldN0);
cont.Y = cont.Y(1, :, 1:n - 1);
cont.eig = cont.eig(1, :, 1:n - 1);
cont.parvalues = cont.parvalues(1:n - 1, 1);
end
%%
function [N1,eqfound,isstable,eigenval,n] = dostep(apar,pnew,n,N1)
global g_grind;
multassignin('base', apar, pnew);
[N0, eqfound] = findeq('Display','off');
if eqfound
    [~, eigenval] = i_eigen(1,g_grind.solver.iters);
    isstable =  i_stability(eigenval, g_grind.solver.isdiffer, g_grind.solver.haslags);
%     cont.Y(1, :, n) = N0.';
%     cont.eig(1, :, n) = eigenval';
%     cont.parvalues(n) = pnew;
    n = n + 1;
    i_keep(2 * N0 - N1);
    N1 = N0;
else
    isstable=false;
    eigenval=nan(size(N1));
end

end

