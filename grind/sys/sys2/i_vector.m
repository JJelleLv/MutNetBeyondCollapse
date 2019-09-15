%i_vector internal function to create a matrix with vectors
%method=0 - normal
%method=1 - relative
%method=2 - accel (not working)
%method=4 - quasi steady state
function [X, Y, Vect, Vals] = i_vector(npoints, iX, Xaxis, iY, Yaxis, N0, method)
global g_grind;
if nargin<7
   method=0;
end

%Eps=1e-4;
t=0;
if abs(Xaxis(1)) < 0.00001*abs(Xaxis(2))
   Xaxis(1) = 0.00001*abs(Xaxis(2));
end
if abs(Yaxis(1)) < 0.00001*abs(Yaxis(2))
   Yaxis(1) = 0.00001*abs(Yaxis(2));
end
if numel(npoints)==1
    npoints=[npoints;npoints];
end

if (nargin < 6)||isempty(N0)
   N = i_initvar;
else
   N=N0;
end

if ismember(method,[0 1])&&(iX.isvar||iX.ispar)&&(iY.isvar||iY.ispar)
    [X,Y]=meshgrid(linspace(Xaxis(1),Xaxis(2),npoints(1)),linspace(Yaxis(1),Yaxis(2),npoints(2)));
    X=transpose(X);
    Y=transpose(Y);
    if ~isempty(iX.transform)
        if nargin(iX.transform.invfun)==2
            X1=iX.transform.invfun(X,X);
        else
            X1=iX.transform.invfun(X);
        end
    else
        X1=X;
    end
    if ~isempty(iY.transform)
        if nargin(iY.transform.invfun)==2
            Y1=iY.transform.invfun(Y,X1);
        else
            Y1=iY.transform.invfun(Y);
        end
    else
        Y1=Y;
    end
    N0 = reshape(repmat(transpose(N),prod(npoints),1),npoints(1),npoints(2),length(N));
    if ~iX.ispar
       N0(:,:,iX.no)=X1;
    end

    if ~iY.ispar
       N0(:,:,iY.no)=Y1;
    end

    if ~(iX.ispar||iY.ispar)
       Vect = i_runsinglestep(t + 0.00001, N0, true);
    else
       if iX.ispar&&iY.ispar
           P0(:,:,1)=X1;
           P0(:,:,2)=Y1;
           parno=[iX.no iY.no];
       elseif iX.ispar
           P0=X;
           parno=iX.no;
       elseif iY.ispar
           P0=Y1;
           parno=iY.no;
       end
       Vect = i_runsinglestep_withpars(t + 0.00001, N0, true ,parno,P0);   
    end

    if g_grind.solver.isdiffer
        Vect=Vect-N0;
    end

    if method==1
        if iX.isvar
          Vect(:,:,iX.no)=Vect(:,:,iX.no)./X1;
        end

        if iY.isvar
           Vect(:,:,iY.no)=Vect(:,:,iY.no)./Y1;
        end

    end

    return;
end

i_runsinglestep([],[],true);

if length(N)>2&&method==4
    %find partial equilibrium if there are more than 2 variables
    mask=true(size(N));
    mask(iX.no)=false;
    mask(iY.no)=false;
    nxtra=sum(mask);
    method=0;
else
    nxtra=0;
end

isstable=true;
oldX = i_getparvar(iX);
oldY = i_getparvar(iY);
if iX.isext
   oldactX=g_grind.externvars{iX.no}.options.active;
   g_grind.externvars{iX.no}.options.active=0;
end
if iY.isext
   oldactY=g_grind.externvars{iY.no}.options.active;
   g_grind.externvars{iY.no}.options.active=0;
end

try
   nvar = g_grind.statevars.dim;
   Vect = zeros(npoints(1), npoints(2), nvar);
   X = zeros(npoints(1), npoints(2));
   Y = zeros(npoints(1), npoints(2));
   N2=[];
   if nxtra>0
       Vals = zeros(npoints(1), npoints(2), nxtra);
   else
       Vals=[];
   end

   incrY = (Yaxis(2) - Yaxis(1)) / (npoints(2)-1);
   incrX = (Xaxis(2) - Xaxis(1)) / (npoints(1)-1);
   isdiffer = g_grind.solver.isdiffer;
   for y1 = 1:npoints(2)
      py = (y1 - 1) * incrY + Yaxis(1);
      if ~isempty(N2)
          N=N2;
      end

      N=i_setparvar(iY,N,py);
      if ~iY.isvar
          i_runsinglestep([],[],true);
      end
      for x1 = 1:npoints(1)
         px = (x1 - 1) * incrX + Xaxis(1);
         N=i_setparvar(iX,N,px);
         if ~iX.isvar
            i_runsinglestep([],[],true);
         end
         if nxtra>0
             [N1,isstable1]=findequilibrium(N,mask);
             j=1;
             while ~isstable1&&j<20
                 N(mask)=rand(sum(mask),1)*10;
                [N1,isstable1]=findequilibrium(N,mask);
                j=j+1;
             end

             if ~isstable1
                 isstable=isstable1;
             end

             if x1==1
                 N2=N1;
             end

             Vals(x1,y1,:)=N1(mask);
             N=N1;
             hwait=i_waitbar(1);
             if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
                 Vect=[];
                 return;
             end
         end

         Nres = transpose(i_runsinglestep(t + 0.00001, transpose(N), false));
         %recalculation of rednoise
         X(x1, y1) = px;
         Y(x1, y1) = py;
         if isdiffer
            Nres = Nres - N;
         end
         %for v = 1:nvar
         switch method
             case 0
                 Vect(x1, y1, :) = Nres;
             case 1
                 Vect(x1, y1, :) = Nres./N;
             case 2
               % N1=N;
               % N2=N1+Nres*Eps;
               eeps=1e-5;
%               accel=zeros(size(Nres));
               Nres2=i_runsinglestep(t + 0.00001, transpose(N+eeps*Nres),false);
               accel=(Nres2-Nres)./Nres;
%                 N3=N2+Nres2*Eps;
%                 RelS1=(log(N2)-log(N1))/Eps;
%                 RelS2=(log(N3)-log(N2))/Eps;
%                 accel=(RelS2-RelS1)./(N2-N1);
%                  accel(accel>1)=1;
%                  accel(accel<-1)=-1;
                Vect(x1, y1, :) = accel;
         end

      end

   end

   N=i_setparvar(iY,N,oldY);
   i_setparvar(iX,N,oldX);
   if ~isstable
       warning('grind:null','Some of the quasi steady states are unstable, the results can upredictable, try other initial values');
   end

if iX.isext
   g_grind.externvars{iX.no}.options.active=oldactX;
end
if iY.isext
   g_grind.externvars{iY.no}.options.active=oldactY;
end
   
catch err
 %  err=lasterror;
   N=i_setparvar(iY,N,oldY);
   i_setparvar(iX,N,oldX);
 if iX.isext
   g_grind.externvars{iX.no}.options.active=oldactX;
end
if iY.isext
   g_grind.externvars{iY.no}.options.active=oldactY;
end
  rethrow(err);
end


function [N0,isstable]=findequilibrium(N,mask)
global g_grind;
N0=findeq('Mask',mask,'N0',N,'Keep','no','Display','off');
jac= i_eigen(isempty(g_grind.syms.Jacobian),g_grind.solver.iters,N0);
jac1=jac(mask,mask);
eigenval=eig(jac1);
isstable=i_stability(eigenval, g_grind.solver.isdiffer);

