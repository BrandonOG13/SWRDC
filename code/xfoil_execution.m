function [actions_str]=xfoil_execution(REref,AoAref,itersXFOIL,Ncrit,saveNoiseDataName)
%This function generates the boundary layer displacement database for the
%noise model

%Read REref and AoAref data
no_RE=length(REref);
no_aoa=length(AoAref);

%Construct command list for XFOIL
[actions_str,modairfoil]=xfoil_commandlist(REref,AoAref,no_RE,no_aoa,itersXFOIL,Ncrit);
disp('Please wait until XFOIL computations are finished')
 
%Run XFOIL
system([pwd '\xfoil.exe < actions.txt'])
disp('Please press ENTER key when XFOIL computations are finished')
pause
%Delete XFOIL commands file (actions.txt) & input airfoil profile
delete(actions_str)
delete(modairfoil)

%Read XFOIL output and collect relevant data
P = zeros(no_RE,no_aoa);
S = P;
for i = 1:no_RE
    for j = 1:no_aoa
        filename = ['bl_dstar\',sprintf('Re_%i_AoA_%i.txt',REref(i),AoAref(j))];
        A = importdata(filename);
        point = find(A.data(:,2)==1);
        P(i,j) = A.data(point(2),5);
        S(i,j) = A.data(point(1),5);
        clear A
    end
end

%Attempt to replace unconverged cases by interpolation
Pn=(P-min(min(P)))/(max(max(P))-min(min(P)));
Sn=(S-min(min(S)))/(max(max(S))-min(min(S)));
A=((abs(cat(1,zeros(1,size(P,2)),diff(Pn,1,1)))>0.09)&(abs(cat(1,diff(Pn,1,1),zeros(1,size(P,2))))>0.09));
B=((abs(cat(1,zeros(1,size(S,2)),diff(Sn,1,1)))>0.09)&(abs(cat(1,diff(Sn,1,1),zeros(1,size(S,2))))>0.09));
P(A)=NaN;
S(B)=NaN;

%Check for NaN
if any(any(isnan(P)))
    P=inpaint_nans(P,0);
end
if any(any(isnan(S)))
    S=inpaint_nans(S,0);
end

fid = fopen(saveNoiseDataName,'w');

for n = 1:no_RE
fprintf(fid,'%g ',P(n,:));

fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

for m = 1:no_RE
fprintf(fid,'%g ',S(m,:));

fprintf(fid,'\n');
end
fclose(fid);
disp('Airfoil noise database creation complete')
end

%% This function constructs the command list for all reynolds number &
%  angle of attack instances in XFOIL
function [actions_str,modairfoil]=xfoil_commandlist(REref,AoAref,no_RE,no_aoa,itersXFOIL,Ncrit)

modairfoil='airfoilXFOIL.txt';
actions_str='actions.txt';
fid = fopen(actions_str,'wt+');
fprintf(fid,'LOAD %s\n',modairfoil);
fprintf(fid,'%s\n','tempName');
actions1={'PANE';'MDES';'FILT 10.0';'EXEC';'';'PANE';'GDES';'CADD';...
    '';'';'';'UNIT';'EXEC';'';'PANEL';'OPER';'VPAR';sprintf('N %i',Ncrit);'';sprintf('ITER %i',itersXFOIL);...
    'VISC';num2str(REref(1))};

aoa_0=find(abs(AoAref)==min(abs(AoAref)));

lenTotal=no_RE*((2*aoa_0)+4+(2*(no_aoa-aoa_0)));

c=0;
actions2=cell(lenTotal,1);
for i=1:1:no_RE
    c=c+1;
    actions2(c)={sprintf('Re %i',REref(i))};
    for k=1:1:aoa_0
        c=c+1;
        alfa=AoAref(aoa_0+1-k);
        actions2(c)={sprintf('ALFA %i',alfa)};
        c=c+1;
        actions2(c)={cat(2,'DUMP ','bl_dstar\',sprintf('Re_%i_AoA_%i.txt',REref(i),alfa))};
    end
    c=c+1;
    actions2(c)={'INIT'};
    c=c+1;
    actions2(c)={sprintf('ALFA %i',AoAref(aoa_0))};
    for j=(aoa_0+1):1:no_aoa
        c=c+1;
        actions2(c)={sprintf('ALFA %i',AoAref(j))};
        c=c+1;
        actions2(c)={cat(2,'DUMP ','bl_dstar\',sprintf('Re_%i_AoA_%i.txt',REref(i),AoAref(j)))};
    end
    c=c+1;
    actions2(c)={'INIT'};
end
actions=[actions1;actions2;{'';'QUIT'}];
fprintf(fid,'%s\n',actions{:});
fclose(fid);

end

%% FROM MATLAB FILE EXCHANGE Interpolates (& extrapolates) NaN elements in a 2d array.
function B=inpaint_nans(A,method)
% INPAINT_NANS: in-paints over nans in an array
% usage: B=INPAINT_NANS(A)          % default method
% usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use
%       for the interpolation.) All methods are capable
%       of extrapolation, some are better than others.
%       There are also speed differences, as well as
%       accuracy differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor.
%       method  3 uses a better plate equation,
%                 but may be much slower and uses
%                 more memory.
%       method  4 uses a spring metaphor.
%       method  5 is an 8 neighbor average, with no
%                 rationale behind it compared to the
%                 other methods. I do not recommend
%                 its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a
%         linear system in the case of only a few
%         NaNs in a large array.
%         Extrapolation behavior is linear.
%
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts
%         of the array which do not have any contact with
%         NaNs. Uses a least squares approach, but it
%         does not modify known values.
%         In the case of small arrays, this method is
%         quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%
%         Note: method 2 has problems in 1-d, so this
%         method is disabled for vector inputs.
%
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result
%         in more accurate interpolations, at some cost
%         in speed.
%
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero)
%         connect each node with every neighbor
%         (horizontally, vertically and diagonally)
%         Since each node tries to be like its neighbors,
%         extrapolation is as a constant function where
%         this is consistent with the neighboring nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element.
%         This method is NOT recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1);
%  z0 = exp(x+y);
%  znan = z0;
%  znan(20:50,40:70) = NaN;
%  znan(30:90,5:10) = NaN;
%  znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06


% I always need to know which elements are NaN,
% and what size the array is for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
    method = 0;
elseif ~ismember(method,0:5)
    error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
    case 0
        % The same as method == 1, except only work on those
        % elements which are NaN, or at least touch a NaN.
        
        % is it 1-d or 2-d?
        if (m == 1) || (n == 1)
            % really a 1-d case
            work_list = nan_list(:,1);
            work_list = unique([work_list;work_list - 1;work_list + 1]);
            work_list(work_list <= 1) = [];
            work_list(work_list >= nm) = [];
            nw = numel(work_list);
            
            u = (1:nw)';
            fda = sparse(repmat(u,1,3),bsxfun(@plus,work_list,-1:1), ...
                repmat([1 -2 1],nw,1),nw,nm);
        else
            % a 2-d case
            
            % horizontal and vertical neighbors only
            talks_to = [-1 0;0 -1;1 0;0 1];
            neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
            
            % list of all nodes we have identified
            all_list=[nan_list;neighbors_list];
            
            % generate sparse array with second partials on row
            % variable for each element in either list, but only
            % for those nodes which have a row index > 1 or < n
            L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
            nl=length(L);
            if nl>0
                fda=sparse(repmat(all_list(L,1),1,3), ...
                    repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
                    repmat([1 -2 1],nl,1),nm,nm);
            else
                fda=spalloc(n*m,n*m,size(all_list,1)*5);
            end
            
            % 2nd partials on column index
            L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
            nl=length(L);
            if nl>0
                fda=fda+sparse(repmat(all_list(L,1),1,3), ...
                    repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
                    repmat([1 -2 1],nl,1),nm,nm);
            end
        end
        
        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);
        k=find(any(fda(:,nan_list(:,1)),2));
        
        % and solve...
        B=A;
        B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
        
    case 1
        % least squares approach with del^2. Build system
        % for every array element as an unknown, and then
        % eliminate those which are knowns.
        
        % Build sparse matrix approximating del^2 for
        % every element in A.
        
        % is it 1-d or 2-d?
        if (m == 1) || (n == 1)
            % a 1-d case
            u = (1:(nm-2))';
            fda = sparse(repmat(u,1,3),bsxfun(@plus,u,0:2), ...
                repmat([1 -2 1],nm-2,1),nm-2,nm);
        else
            % a 2-d case
            
            % Compute finite difference for second partials
            % on row variable first
            [i,j]=ndgrid(2:(n-1),1:m);
            ind=i(:)+(j(:)-1)*n;
            np=(n-2)*m;
            fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
                repmat([1 -2 1],np,1),n*m,n*m);
            
            % now second partials on column variable
            [i,j]=ndgrid(1:n,2:(m-1));
            ind=i(:)+(j(:)-1)*n;
            np=n*(m-2);
            fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
                repmat([1 -2 1],np,1),nm,nm);
        end
        
        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);
        k=find(any(fda(:,nan_list),2));
        
        % and solve...
        B=A;
        B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
        
    case 2
        % Direct solve for del^2 BVP across holes
        
        % generate sparse array with second partials on row
        % variable for each nan element, only for those nodes
        % which have a row index > 1 or < n
        
        % is it 1-d or 2-d?
        if (m == 1) || (n == 1)
            % really just a 1-d case
            error('Method 2 has problems for vector input. Please use another method.')
            
        else
            % a 2-d case
            L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
            nl=length(L);
            if nl>0
                fda=sparse(repmat(nan_list(L,1),1,3), ...
                    repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
                    repmat([1 -2 1],nl,1),n*m,n*m);
            else
                fda=spalloc(n*m,n*m,size(nan_list,1)*5);
            end
            
            % 2nd partials on column index
            L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
            nl=length(L);
            if nl>0
                fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
                    repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
                    repmat([1 -2 1],nl,1),n*m,n*m);
            end
            
            % fix boundary conditions at extreme corners
            % of the array in case there were nans there
            if ismember(1,nan_list(:,1))
                fda(1,[1 2 n+1])=[-2 1 1];
            end
            if ismember(n,nan_list(:,1))
                fda(n,[n, n-1,n+n])=[-2 1 1];
            end
            if ismember(nm-n+1,nan_list(:,1))
                fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
            end
            if ismember(nm,nan_list(:,1))
                fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
            end
            
            % eliminate knowns
            rhs=-fda(:,known_list)*A(known_list);
            
            % and solve...
            B=A;
            k=nan_list(:,1);
            B(k)=fda(k,k)\rhs(k);
            
        end
        
    case 3
        % The same as method == 0, except uses del^4 as the
        % interpolating operator.
        
        % del^4 template of neighbors
        talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
            0 1;0 2;1 -1;1 0;1 1;2 0];
        neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
        
        % list of all nodes we have identified
        all_list=[nan_list;neighbors_list];
        
        % generate sparse array with del^4, but only
        % for those nodes which have a row & column index
        % >= 3 or <= n-2
        L = find( (all_list(:,2) >= 3) & ...
            (all_list(:,2) <= (n-2)) & ...
            (all_list(:,3) >= 3) & ...
            (all_list(:,3) <= (m-2)));
        nl=length(L);
        if nl>0
            % do the entire template at once
            fda=sparse(repmat(all_list(L,1),1,13), ...
                repmat(all_list(L,1),1,13) + ...
                repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
                repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
        else
            fda=spalloc(n*m,n*m,size(all_list,1)*5);
        end
        
        % on the boundaries, reduce the order around the edges
        L = find((((all_list(:,2) == 2) | ...
            (all_list(:,2) == (n-1))) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1))) | ...
            (((all_list(:,3) == 2) | ...
            (all_list(:,3) == (m-1))) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1))));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(all_list(L,1),1,5), ...
                repmat(all_list(L,1),1,5) + ...
                repmat([-n,-1,0,+1,n],nl,1), ...
                repmat([1 1 -4 1 1],nl,1),nm,nm);
        end
        
        L = find( ((all_list(:,2) == 1) | ...
            (all_list(:,2) == n)) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1)));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(all_list(L,1),1,3), ...
                repmat(all_list(L,1),1,3) + ...
                repmat([-n,0,n],nl,1), ...
                repmat([1 -2 1],nl,1),nm,nm);
        end
        
        L = find( ((all_list(:,3) == 1) | ...
            (all_list(:,3) == m)) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1)));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(all_list(L,1),1,3), ...
                repmat(all_list(L,1),1,3) + ...
                repmat([-1,0,1],nl,1), ...
                repmat([1 -2 1],nl,1),nm,nm);
        end
        
        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);
        k=find(any(fda(:,nan_list(:,1)),2));
        
        % and solve...
        B=A;
        B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
        
    case 4
        % Spring analogy
        % interpolating operator.
        
        % list of all springs between a node and a horizontal
        % or vertical neighbor
        hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
        hv_springs=[];
        for i=1:4
            hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
            k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
            hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
        end
        
        % delete replicate springs
        hv_springs=unique(sort(hv_springs,2),'rows');
        
        % build sparse matrix of connections, springs
        % connecting diagonal neighbors are weaker than
        % the horizontal and vertical springs
        nhv=size(hv_springs,1);
        springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
            repmat([1 -1],nhv,1),nhv,nm);
        
        % eliminate knowns
        rhs=-springs(:,known_list)*A(known_list);
        
        % and solve...
        B=A;
        B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;
        
    case 5
        % Average of 8 nearest neighbors
        
        % generate sparse array to average 8 nearest neighbors
        % for each nan element, be careful around edges
        fda=spalloc(n*m,n*m,size(nan_list,1)*9);
        
        % -1,-1
        L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % 0,-1
        L = find(nan_list(:,3) > 1);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % +1,-1
        L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % -1,0
        L = find(nan_list(:,2) > 1);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % +1,0
        L = find(nan_list(:,2) < n);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % -1,+1
        L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % 0,+1
        L = find(nan_list(:,3) < m);
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % +1,+1
        L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
        nl=length(L);
        if nl>0
            fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
                repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
                repmat([1 -1],nl,1),n*m,n*m);
        end
        
        % eliminate knowns
        rhs=-fda(:,known_list)*A(known_list);
        
        % and solve...
        B=A;
        k=nan_list(:,1);
        B(k)=fda(k,k)\rhs(k);
        
end

% all done, make sure that B is the same shape as
% A was when we came in.
B=reshape(B,n,m);
end

% ====================================================
%      end of main function
% ====================================================
% ====================================================
%      begin subfunctions
% ====================================================
function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans
%   themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element
%      nan_list(i,2) == row index of i'th nan element
%      nan_list(i,3) == column index of i'th nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%
%      For example, talks_to = [-1 0;0 -1;1 0;0 1]
%      means that each node talks only to its immediate
%      neighbors horizontally and vertically.
%
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
    % use the definition of a neighbor in talks_to
    nan_count=size(nan_list,1);
    talk_count=size(talks_to,1);
    
    nn=zeros(nan_count*talk_count,2);
    j=[1,nan_count];
    for i=1:talk_count
        nn(j(1):j(2),:)=nan_list(:,2:3) + ...
            repmat(talks_to(i,:),nan_count,1);
        j=j+nan_count;
    end
    
    % drop those nodes which fall outside the bounds of the
    % original array
    L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m);
    nn(L,:)=[];
    
    % form the same format 3 column array as nan_list
    neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
    
    % delete replicates in the neighbors list
    neighbors_list=unique(neighbors_list,'rows');
    
    % and delete those which are also in the list of NaNs.
    neighbors_list=setdiff(neighbors_list,nan_list,'rows');
    
else
    neighbors_list=[];
end
end