% CREATE_NETWORK - create different kinds of networks
% create a matrix with ones for connections and zeros for no connection, using different rules:
%    -exponential (e) (created by building a network, adding nodes one by
%    one, creating a number of edges (determined by PAR) between the added node
%    and existing nodes (connection with same node twice is allowed!)
%    -scale-free (s) (same as exponential, but with a preference (higher chance)
%    to connect with a node that has many connections already
%    -lattice 4 neighbors (l)
%    -random (r) pure random, no check for unconnected nodes
%    -random-regular (rr) random but with a fixed degree (all nodes have
%    degree=par
%
%  Usage:
%    [A,ndx]=CREATE_NETWORK(N, TYPE, PAR) N= number of nodes, type is kind of network ('e','s','l' etc)
%       PAR = a specific parameter (see below). The function returns A (sparse matrix with connections (use FULL to convert
%       to full matrix) and ndx, which  are the indices>0 of the matrix A
%
%    PAR: for exponential and scale-free PAR = the number of legs (= number of new connections)
%         for lattice PAR = 1 for periodic boundaries (and 0 for non periodic boundaries)
%         for random PAR = probability of a connection
%
%  See also:
%     SPARSE2NDX, NDX2SPARSE
%
function [A, ndxs] = create_network(n, type, par1)

if nargin < 1
    prompt={'Number of nodes','Kind of network (exponential (e), scale-free (s), lattice (l) or random (r)','Extra parameter (number of legs/probability)'};
    def={'100','e','1'};
    dlgTitle = 'Create Network';
    lineNo = 1;
    answer = inputdlg(prompt, dlgTitle, lineNo, def);
    n = str2double(answer{1});
    type = answer{2};
    par1 = str2double(answer{3});
elseif nargin < 2
    type = 'e';
    disp('exponential network'); %default network
end
random_order = 0; % switch to add the nodes in random order or not
if ~ischar(type)
    error('GRIND:create_network:UnknownType','Network type "%g" not supported', type);
end
if strncmpi(type, 'rr', 2)
    %%random regular network
    A=randomregular(n,par1);
    if nargout>1
       ndxs=sparse2ndx(A);
    end
elseif strncmpi(type, 'e', 1)
    %%%% EXPONENTIAL NETWORK %%%%
    
    %%%%% Gilarranz & Bascompte 2011 %%%%
    % Exponential networks such as the one displayed in Fig. 1c also
    % allow long-distance interactions, but as can be seen by the varia-
    % bility in node sizes, they show a more heterogeneous range of node?s
    % degree (Fortuna et al., 2006).
    % These networks were constructed following the algorithm by
    % Barabasi and Albert (1999) with random attachment:
    
    %%%%% Barabasi and Albert (1999)  %%%%%
    % Model A keeps the growing character of the network, but preferential
    % attachment is eliminated by assuming that a new vertex is connected
    % with equal probability to any vertex in the system [that is,
    % P(k)= const= 1/(m0+t-1)]. Such a model leads to P(k)= exp(-bk)
    
    if nargin < 3
        legs = 1; %default number of legs is one
    else
        legs = round(par1); %par1 are the number of legs
        if legs < 1
            error('GRIND:create_network:TooFewLegs','The number of new connections ("legs") should be at least 1');
        end
    end
    if random_order
        %shuffle the indexes of nodes (1:n in random order)
        rands = rand(n, 1); 
        [x, ndx] = sort(rands); %ndx are the nodes in random order
    else
        ndx = transpose(1:n);
    end
    ndx1 = ndx;
    %the nodes in random order are added one by one
    %now set for each node the number of nodes that should be already in the network
    maxndx = transpose(0:n - 1);
    ndxs = [];
    %draw the nodes to connect
    %     for j = 1:legs
    %         %choose a random node that is already in the network (don't mind double connections)
    %         ndx2 = ndx1(floor(rand(n, 1) .* maxndx) + 1);
    %         %save the subscripts of the connections
    %         ndxs = [ndxs; [ndx1(2:n), ndx2(2:n)]]; %no check for doubles
    %     end
    
    for i=1:n
        % number of legs for node 'n'
        legs=(rand(1,1)<par1-floor(par1))+floor(par1);
        for j = 1:legs
            %choose a random node that is already in the network (don't mind double connections)
            
            %determine with which node a connection will be made
            conn=floor(rand(1, 1) .* maxndx(i)) + 1;
            %fill ndx2
            ndx2(i,1) = ndx1(conn);
            %save the subscripts of the connections
            ndxs = [ndxs; [ndx1(2:i), ndx2(2:i)]]; %no check for doubles
        end
        
    end
    
    %fill the square matrix with connections
    A=ndx2sparse(ndxs);
    [~,ndx]=sortrows(sum(A,2));
    A=A(ndx,ndx);
    
elseif strncmpi(type, 's', 1)
    
    %%%% SCALE-FREE NETWORK %%%%
    % is similar as exponential, but with weighted selection
    
    %%%%% Gilarranz & Bascompte 2011 %%%%
    % Finally, scale-free networks (Fig. 1d) are extremely heteroge-
    % neous (Kininmonth et al., 2010). The majority of nodes have one
    % or a few links, but a few nodes are extremely well-connected.
    % The degree distribution follows a power law (Barabasi and Albert, 1999).
    
    %%%%% Barabasi and Albert (1999)  %%%%%
    % We next show that a model based on these
    % two ingredients naturally leads to the observed
    % scale-invariant distribution. To incorporate
    % the growing character of the network,
    % starting with a small number (m0) of vertices,
    % at every time step we add a new vertex with
    % m(<=m0) edges that link the new vertex to m
    % different vertices already present in the system.
    % To incorporate preferential attachment,
    % we assume that the probability P that a new
    % vertex will be connected to vertex i depends
    % on the connectivity ki of that vertex, so that
    % P(ki)=ki/sum(kj). After t time steps, the
    % model leads to a random network with t+ m0
    % vertices and mt edges.
    
    if nargin < 3
        par1 = 1; %default number of legs is one
    else
        %par1 are the number of legs
        if  floor(par1) < 1
            error('GRIND:create_network:TooFewLegs','The number of new connections ("legs") should be at least 1');
        end
    end
    if random_order
        %shuffle the indexes of nodes (1:n in random order)
        rands = rand(n, 1); 
        [x, ndx] = sort(rands); %ndx are the nodes in random order
    else
        ndx = transpose(1:n);
    end
    %the nodes in random order are added one by one
    %now set for each node the number of nodes that should be already in the network
    counts = zeros(n, 1);  % count the links
    %   A = sparse(0);
    ndxs = zeros(n * floor(par1)-(2*floor(par1)-1), 2); %was n*legs
    k = 1;
    %draw the nodes to connect
    for i = 2:n
        if (i == 2)  %first node is added
            ndxs(k, :) = transpose(ndx([i, i - 1])); %ndx could be shuffled
            k = k + 1;
            counts([i, i - 1]) = counts([i, i - 1]) + 1;
        else
            legs=(rand(1,1)<par1-floor(par1))+floor(par1);
            P = counts ./ sum(counts);
            P = P(1:i - 1);
            for j = 1:legs % could be made more efficient if length(added) > legs
                added = [];
                while isempty(added)
                    r = rand(i - 1, 1);
                    added = find(r < P);
                end
                anndx = [i, added(floor(rand(1) .* length(added)) + 1)]; %if more than one added,
                %choose at random one of these
                %(THIS IS DIFFERENT FROM LUISJO, he chooses allways the first);
                counts(anndx) = counts(anndx) + 1;
                ndxs(k, :) = transpose(ndx(anndx));
                k = k + 1;
            end
        end
    end
    %fill the square matrix with connections
    A=ndx2sparse(ndxs);
    [~,ndx]=sortrows(sum(A,2));
    A=A(ndx,ndx);
    
elseif strncmpi(type, 'r', 1)
    
    %%%%% RANDOM NETWORK %%%%%
    %pure random no check if nodes are all connected
    if nargin < 3
        prob = 0.5; %default probability of a connection = 0.5
    else
        prob = par1; %par1 is the probability
    end
    if prob<=1/n
        error('grind:createnetwork','Too small probablity to have at least one connection for each node');
    end        
    found_valid=false;
    i=0;
    maxtime=60; %try and check for one minute (could be short for big sparsely connected network)
    tic;
    while ~found_valid
        A = triu(rand(n)); %we only need one triangle as the matrix should be symmetric
        A(logical(speye(size(A)))) = 0;
        A = sparse(A > 1 - prob);  %we have a connection with a probability of prob (1 - p > rand)
        A =A+transpose(triu(A)); %make symmetric matrix
        i=i+1;
        found_valid=toc>maxtime||check_conn(A);
    end
    if toc>maxtime
        error('grind:createnetwork','Could not find a fully connected network after many trials');
    end
    if nargout > 1
        ndxs = sparse2ndx(A);
    end
    
elseif strncmpi(type, 'l', 1)
    %%%%% REGULAR LATTICE (4 neighbor) %%%%%
    % n can now be a vector [sizeX,sizeY]
    % if not n is assumed to be the total number of cells of a square matrix of [sqrt(n),sqrt(n)]
    if length(n) == 1
        siz = [round(sqrt(n)), round(sqrt(n))];
        if sqrt(n) - siz(1) > 0.01
            warning('GRIND:create_network:square','Matrix size assumed to be square: changed to [%d,%d] ',siz);
        end
    else
        siz = n;
    end
    if nargin < 3
        isround = 1; %default periodic boundaries
    else
        isround = par1; %par1 is the probability
    end
    ntot=prod(siz);
    ndxs = zeros(4*ntot,2);
    k=1;
    for i = 1:siz(1)
        for j = 1:siz(2)
            ii = sub2ind(siz, i, j);
            ndxs(k:k+3,1)=ii;
            ndxs(k,2)= sub2ind(siz, i, getndx(j - 1, siz(2), isround));
            ndxs(k+1,2)=sub2ind(siz, i, getndx(j + 1, siz(2), isround));
            ndxs(k+2,2)=sub2ind(siz, getndx(i + 1, siz(1), isround), j);
            ndxs(k+3,2)=sub2ind(siz, getndx(i - 1, siz(1), isround), j);
            %neighbors getndx takes the periodic boundaries into account
            k=k+4;
        end
    end
    A=ndx2sparse(ndxs);
    % A=ndx2sparse(ndxs(ndxs(:,1)>0,:));
else
    error('GRIND:create_network:UnknownType','Network type "%s" not supported', type);
end

function i = getndx(i, s1, isround)
if isround
    if i <= 0
        i = i + s1;
    elseif i > s1
        i = i - s1;
    end
else
    if i <= 0
        i = 1;
    elseif i > s1
        i = s1;
    end
end

function M=randomregular(nnodes,legs)
%Random regular network, a fast algorithm
degree=zeros(nnodes,1);
%if the required legs is larger than the half it is faster to remove the
%legs
if legs>nnodes/2
    legs=nnodes-legs;
    neg=true;
else
    neg=false;
end
%maximum 100 iterations before give up to get exactly the same degree
%this happens when we want high connectivity
numiters=0;

while any(degree~=legs)&&numiters<100
    M=sparse(nnodes,nnodes);
    degree=zeros(nnodes,1);
    found=false;
    iter=0;
    while ~found
        %select the fee nodes (i.e. nodes with less than legs nodes)
        freenodes=find(degree<legs);
        n=length(freenodes);
        if n>1
            %get two different random nodes
            i1=floor(rand(1)*n)+1;
            i2=floor(rand(1)*n)+1;
            if i1==i2
                ndx=[];
                i=1;
                while isempty(ndx)&&i<50
                    i2=floor(rand(50,1)*n)+1;
                    ndx=find(i2~=i1,1);
                    i=i+1;
                end
            else
                ndx=1;
            end
            i2=i2(ndx);
            n1=freenodes(i1);
            n2=freenodes(i2);
            if M(n1,n2) % problem if M is already true, possibly it is impossible to create a random network
                iter=iter+1;
            end
            M(n1,n2)=true;
            M(n2,n1)=true;
            degree=sum(M,2);
        else
            found=true;
        end
        if iter>100
            found=true;
        end
    end
    numiters=numiters+1;
end
if neg
    M=~M;
end


