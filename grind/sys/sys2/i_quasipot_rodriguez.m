function res=i_quasipot_rodriguez(args)
% calculation of potential and curl for a two dimensional model. 
% Method Pablo Rodriguez Sanchez
global g_grind
%symbolic Jacobians for speed
try
    enterjac('-sym')
catch err
    if strcmp(err.identifier,'grind:enterjac:vector')
        disp('No symbolic Jacobian');
    else
        rethrow(err);
    end
end
%get the default axes for the phase plane
Xaxis=g_grind.xaxis.lim;
Yaxis=g_grind.yaxis.lim;
iX=i_getno(g_grind.xaxis.var);
iY=i_getno(g_grind.yaxis.var);

if length(args.npoints)==1
    args.npoints=[args.npoints args.npoints];
end
%generate two dimensions for the plot
[X,Y]=meshgrid(linspace(Xaxis(1),Xaxis(2),args.npoints(2)),linspace(Yaxis(1),Yaxis(2),args.npoints(1)));
%if the model has more dimensions some tricks are needed
N0=i_initvar;
N0 = repmat(N0,1,numel(X));
N0(iX.no,:)=X(:).';
N0(iY.no,:)=Y(:).';
%function handles to evaluate the ode or the Jacobians
han=i_getodehandle('Jacobian');
han2=i_getodehandle('normal');
%determine all Jacobians
if strcmpi(g_grind.solver.opt.Vectorized,'on');
    %very fast evaluation of the model is vectorized:
    Jacs=han(0,N0);
    N1=han2(0,N0);
else
    %somewhat slower if the Jacobians are not vectorized:
    Jacs=zeros(size(N0,1),size(N0,1),args.npoints(1)*args.npoints(2));
    N1=zeros(size(N0));
    for i=1:size(N0,2)
        Jacs(:,:,i)=han(0,N0(:,i));
        N1(:,i)=han2(0,N0(:,i));
    end
end

%calculate the potential and error for subsequent points
V=zeros(args.npoints);
Jacs=reshape(Jacs,size(Jacs,1),size(Jacs,2),args.npoints(1),args.npoints(2));
N1=reshape(N1.',args.npoints(1),args.npoints(2),size(N1,1));
N0=reshape(N0.',args.npoints(1),args.npoints(2),size(N0,1));
err=V;
for i=1:args.npoints(1)
    for j=1:args.npoints(2)
        %choose the closest point, V(1,1) = 0
        if ~(j==1&&i==1)
            if j==1
                %take x direction for first y
                i1=i-1;
                j1=j;
            else
                %y direction
                i1=i;
                j1=j-1;
            end
            Jsymm=0.5*(Jacs(:,:,i,j)+Jacs(:,:,i,j).');
            Jskew=0.5*(Jacs(:,:,i,j)-Jacs(:,:,i,j).');
            diffN0=N0(i,j,:)-N0(i1,j1,:);
            diffN0=diffN0(:); %make vertical vector
            NN1=N1(i,j,:);
            deltaV=-NN1(:).'*diffN0-0.5*diffN0.'*Jsymm*diffN0;
            V(i,j)=V(i1,j1)+deltaV;
            %err(i,j)=norm(Jskew*diffN0); %correct?
            err(i,j)=norm(Jskew*diffN0)/(norm(Jskew*diffN0)+norm(Jsymm*diffN0)); %correct?
        end
    end
end
pot=V;

if nargout>0
    res.X=X;
    res.Y=Y;
    res.pot=pot;
    res.curl=err;
end


