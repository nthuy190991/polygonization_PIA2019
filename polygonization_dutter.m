function [vertices,points,S,Psub]= polygonization_dutter(P,GP,Vspec,Mspec,nb_pts,verbose)
% This script was written based on the work of Dutter [1] and was used in
% the our work [2].
% Note: The two algorithms only work for 3 levels of shapes (rectangle;
% L-, T- or Z-shape; and U-shape).
% The best way to understand correctly how these algorithms work is to read
% the respective references. This script doesn't contain the full details.
% There are places in this script where a page number is provided. The
% readers are suggested to refer to such pages in [1] for fully understanding.
%
% Reference:
% [1] Dutter, M. (2007). "Generalization of building footprints derived
%     from high resolution remote sensing data", Institut für Photogrammetrie
%     und Fernerkundung, Technische Universität Wien
% [2] T. H. Nguyen et al. (2019). "Unsupervised Automatic Building Extraction
%     Using Active Contour Model on Unregistered Optical Imagery and Airborne
%     LiDAR Data," Int. Arch. Photogramm. Remote Sens. Spatial Inf. Sci.,
%     XLII-2/W16, 181-188. DOI: 10.5194/isprs-archives-XLII-2-W16-181-2019
%
% The algorithm will result in a number of split-points. The points between
% each pair of split-points is a segment.
% The index of split-points in P will be saved in S.
% The subsets of points in the segments will be saved in Psub.
%  (Note: Psub is a struct, each struct element contains points in 1 segment)
%
% Inputs:
% - P: set of points to be polygonized (Mx2)
% - GP: guided points (Nx2) only used to compute MBR - Improved Polygonization
%       Set GP=[] or GP=P to perform the original algorithm by Dutter.
% User-defined parameters:
%  - Vspec: variation parameter (cf. page 46)
%  - Mspec: minimal length of edge (cf. page 46)
% nb_pts: number of points for each segment of the resulting polygon,
%         Otherwise set nb_pts=0
% verbose: verbose mode flag

% Outputs:
% - vertices: the vertices of the resulting polygon
% - points: the resulting polygon with each segment is resampled with nb_pts
% - S: set of split-point indices (Nx2)
% - Psub: set of subsets of points

%% Credits
% If you use any code in this script, cite our paper:
%
%   T. H. Nguyen et al., 2019, "Unsupervised Automatic Building Extraction Using
%   Active Contour Model on Unregistered Optical Imagery and Airborne LiDAR Data,"
%   Int. Arch. Photogramm. Remote Sens. Spatial Inf. Sci., XLII-2/W16, 181-188.
%   DOI: 10.5194/isprs-archives-XLII-2-W16-181-2019
%
% Written by T. H. Nguyen, Université Laval, 2019-2020.
%   For any other questions/issues/bug reports, please open an issue on the
%   Issues tracker (on the Github page).
%
% Last modified on April 30, 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN OF MAIN PROGRAM %
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine of the Minimum Bounding Rectangle (MBR)
% Original algo.: use the MBR of the input points
% Improved algo.: use the MBR of the projected LiDAR points (more reliable)
if isempty(GP)
    mbr = minBoundingBox([P(:,1)'; P(:,2)']); % mbr: 2x4
else
    mbr = minBoundingBox([GP(:,1)'; GP(:,2)']); % mbr: 2x4
end

% level-of-shape flags:
Level1=1; % rectangle
Level2=0; % L-, T- or Z-shape
Level3=0; % U-shape

if verbose
    close all
    figure
    plot(P(:,1),P(:,2),'b.','MarkerSize',10)
    hold on, grid on, axis equal
    plot([mbr(1,1:4) mbr(1,1)],[mbr(2,1:4) mbr(2,1)],'k--')
end

%% Level 1: rectangular (cf. page 49)
% At this level, 4 split-points are determined.
% Split-points s(i), i=1:4 are the points nearest to the MBR corner t(i).

S1=[]; % split-point index (at level 1), i.e. S1=IndexOf(s(i)) in P, i=1:4

% First, calculate distance from the MBR corners to every points in P
for i=1:4
    for j=1:size(P,1)
        d(j)=dist(mbr(:,i)',P(j,:));
    end
    [~,idx]=min(d);
    S1=[S1; idx];
end

[S1,I]=sort(S1); % sort S1 according to point indice (ascending order)
M0=mbr(:,I);     % sort MBR vertices according to S1

M1=[M0 M0(:,1)]; % M1 is the matrix of edges. Here, there are 4 edges:
% [M0(1,:),M0(2,:)], [M0(2,:),M0(3,:)], [M0(3,:),M0(4,:)],
% and [M0(4,:),M0(1,:)]
% ** Note: above the edges are written as [A,B] where
% A=(xA,yA)' and B=(xB,yB)' are the edge's end-points. **

% Stock subsets of points into a struct
Psub1={}; % struct of subsets of point (at level 1)
% Psub1(1,:) contains the points between split-point 1 and 2
% Psub1(2,:) contains the points between split-point 2 and 3
% Psub1(3,:) contains the points between split-point 3 and 4
% Psub1(4,:) contains the points between split-point 4 and 1

for i=1:4
    if i<4
        Psub1(i,:)={P(S1(i):S1(i+1),:)};
    else
        Psub1(i,:)={[P(S1(i):end,:); P(1:S1(1),:)]};
    end
end

if verbose
    figure
    plot(P(:,1), P(:,2), 'b.')
    hold on, grid on, axis equal
    plot([mbr(1,1:4) mbr(1,1)],[mbr(2,1:4) mbr(2,1)],'k--')
    cl=rand(4,3);
    for i=1:4
        % text(P(S1(i),1), P(S1(i),2), ['s_' num2str(i)]) % split-points
        p=Psub1{i};
        plot(p(:,1),p(:,2),'.','Color',cl(i,:),'MarkerSize',10)
        % plot([M1(1,i) M1(1,i+1)],[M1(2,i) M1(2,i+1)],'-','Color',cl(i))
        % text(mean([M1(1,i),M1(1,i+1)]),mean([M1(2,i), M1(2,i+1)]),...
        %    ['m_' num2str(i)]) % MBR edges
        
        % text(M1(1,i),M1(2,i),['t_' num2str(i)]) % MBR vertices
    end
    title('Level 1 result')
end

% Check variance of distances D from subset points to the corresponding MBR edge
% V(i), i=1:4 are the stds of the distances from each subset points to each edge
V=zeros(1,4);
for i=1:4
    p=Psub1{i};
    D=zeros(1,size(p,1));
    for j=1:size(p,1)
        D(j)=dist_to_line(p(j,:),[M1(:,i) M1(:,i+1)]);
    end
    V(i)=std(D);
end

% If all of V(i)<=Vspec, then Level 1 is accepted.
% if not (i.e. any of V(i)>Vspec), then proceed to Level 2.
if all(V<=Vspec)
    Level2=0;
    Level3=0;
else
    Level2=1;
    Level3=1;
end

if verbose
    disp('Level1: DONE')
    disp(['Level2=', num2str(Level2)])
end

%% Level 2: L-, T- or Z-shape (cf. page 52)
% At level 2, the level-1 polygon (rectangular) is "refined" by dividing one
% or more segments into smaller segments. Depending on how many and which
% segments are divided, the resulting polygon has L-, T-, Z-shape.

% If a split-point s(i) (from level 1) has a large distance from its
% corresponding MBR-corner t(i), the combination of subsets Psub(i-1) and Psub(i)
% shoud be divided into four new subsets.
% Psub(i-1): the subset containing the points between split-point i-1 and i
% Psub(i)  : the subset containing the points between split-point i and i+1

stop_flag=0;
S=S1; M=M1; Psub=Psub1;
iter=1;
while ~stop_flag
    disp(['iter=' num2str(iter)])
    if Level2
        
        % First, calculate the distance d(s(i),t(i)), i=1:4
        % where s(i): split-points, t(i): MBR corners.
        dist_si_ti=zeros(1,length(S));
        for i=1:length(S)
            dist_si_ti(i)=dist(P(S(i),:),M(:,i)');
        end
        idx=find(dist_si_ti>Mspec);
        
        if isempty(idx)
            Level2=0; % if no segment has dist_si_ti>Mspec, reject Level 2
            % M2=M;
        else
            if verbose
                disp(['Level 2: number of level1-segments to be divided=' num2str(length(idx))])
            end
            
            % Split each segment with distance d(s(i),t(i))>Mspec into 3 new split-points:
            % sb (b for backward), sf (f for forward) and sm (m for middle)
            sb=zeros(length(idx),2);
            sf=zeros(length(idx),2);
            h=zeros(length(idx),2);
            sm=zeros(length(idx),2);
            
            % Determine (sb, sm, sf) for each segment having dist_si_ti>Mspec
            for j=1:length(idx)
                
                % Backward split-point (sb) - See Equation (3.4) page 55
                if idx(j)==1
                    B=[4 1];
                else
                    B=[idx(j)-1 idx(j)]; % B is used to get the BACKWARD edge index
                end
                
                if idx(j)==4
                    F=[4 1];
                else
                    F=[idx(j) idx(j)+1]; % F is used to get the FORWARD edge index
                end
                
                % For example, if the split-point s(2) is to be divided,
                % Backward: refer to edge (1,2) (i.e. edge between split-points 1 and 2)
                % Forward:  refer to edge (2,3) (i.e. edge between split-points 2 and 3)
                
                % Calculate distance from points in the level1-segment to the backward edge
                d_backward=zeros(1,length(Psub{B(1)}));
                for i=1:length(Psub{B(1)})
                    d_backward(i)=dist_to_line(Psub{B(1)}(i,:),[M(:,B(1)) M(:,B(2))]);
                end
                idx_sb=max(find(d_backward<2)); % the constant 2 meters chosen by Dutter
                sb(j,:)=Psub{B(1)}(idx_sb,:);
                
                % Forward split-point (sf) - See Equation (3.4) page 55
                d_forward=zeros(1,length(Psub{F(1)}));
                for i=1:length(Psub{F(1)})
                    d_forward(i)=dist_to_line(Psub{F(1)}(i,:),[M(:,F(1)) M(:,F(2))]);
                end
                idx_sf=min(find(d_forward<2)); % the constant 2 meters chosen by Dutter
                sf(j,:)=Psub{F(1)}(idx_sf,:);
                
                % Rotate sb and sf with origin of t(i) with the positive x-axis
                % going through t(i-1) and positive y-axis through t(i+1)
                if idx(j)==1
                    h(j,:)=parallel_intersection(sf(j,:)',[M(:,4) M(:,idx(j))],...
                        sb(j,:)',[M(:,idx(j)) M(:,idx(j)+1)],0);
                else
                    h(j,:)=parallel_intersection(sf(j,:)',[M(:,idx(j)-1) M(:,idx(j))],...
                        sb(j,:)',[M(:,idx(j)) M(:,idx(j)+1)],0);
                end
                
                % Middle split-point (sm) - Equation (3.6) page 54
                Ptmp=[Psub{B(1)}; Psub{idx(j)}];
                d_middle=zeros(1,size(Ptmp,1));
                for i=1:size(Ptmp,1)
                    d_middle(i)=dist(h(j,:),Ptmp(i,:));
                end
                [~,idx_mid]=min(d_middle);
                sm(j,:)=Ptmp(idx_mid,:);
            end
            
            % if verbose
            %     % display(sb)
            %     % display(sf)
            %     for j=1:length(idx)
            %         plot(sb(j,1),sb(j,2),'rx')
            %         text(sb(j,1)+rand,sb(j,2)+rand,['s_{' num2str(idx(j)) 'b}'])
            %         plot(sf(j,1),sf(j,2),'rx')
            %         text(sf(j,1)+rand,sf(j,2)+rand,['s_{' num2str(idx(j)) 'f}'])
            %         plot(h(j,1), h(j,2), 'r+')
            %         text(h(j,1), h(j,2), ['h_' num2str(idx(j))])
            % 
            %         plot(sm(j,1),sm(j,2),'rd')
            %         text(sm(j,1)+rand,sm(j,2)+rand,['s_{' num2str(idx(j)) 'm}'])
            %     end
            % end
            
            % Add new split-points into S2
            % S2: split-point index (at level 2)
            S2=[];
            for j=1:length(idx)
                idx_sb=find(P(:,1)==sb(j,1) & P(:,2)==sb(j,2));
                idx_sf=find(P(:,1)==sf(j,1) & P(:,2)==sf(j,2));
                idx_sm=find(P(:,1)==sm(j,1) & P(:,2)==sm(j,2));
                
                % set of split-points
                S2b=S(1:idx(j)-1);
                S2f=S(idx(j)+1:end);
                S2=[S2; S2b; idx_sb; idx_sm; idx_sf; S2f];
            end
            S2=unique(S2);
            
            % set of MBR edges to be compared with the new split edges
            M2=repmat(M0,[1,length(S2)]);
            
            % Stock subsets of points into a struct
            Psub2={}; % struct of subsets of point (at level 2)
            for i=1:length(S2)
                if i<length(S2)
                    if S2(i)>S2(i+1)
                        Psub2(i,:)={[P(S2(i):end,:); P(1:S2(i+1),:)]};
                    else
                        Psub2(i,:)={P(S2(i):S2(i+1),:)};
                    end
                else %i=length(S2)
                    if S2(i)>S2(1) %S2(end)>S2(1)
                        Psub2(i,:)={[P(S2(i):end,:); P(1:S2(1),:)]};
                    else
                        Psub2(i,:)={P(S2(i):S2(1),:)};
                    end
                end
            end
            
            if verbose
                % display(S2)
                % display(M)
                % display(M2)
                % display(Psub2)
                figure
                plot(P(:,1), P(:,2), 'b.')
                hold on, grid on, axis equal
                plot([mbr(1,1:4) mbr(1,1)],[mbr(2,1:4) mbr(2,1)],'k--')
                %             cl=rand(length(S2)+1,3);
                cl=rand(length(S2)+1,3);
                for i=1:length(S2)
                    % text(P(S2(i),1), P(S2(i),2), ['s_{' num2str(i) '}']) % split-points
                    p=Psub2{i};
                    plot(p(:,1),p(:,2),'.','Color',cl(i,:),'MarkerSize',10)
                end
                title('Level 2 result')
            end
            
            % Check variance of distances D from subset points to the corresponding MBR edge
            V2=zeros(1,length(S2));
            for i=1:length(S2)
                p=Psub2{i};
                D=zeros(1,size(p,1));
                for j=1:size(p,1)
                    D(j)=dist_to_line(p(j,:),[M2(:,i) M2(:,i+1)]);
                end
                V2(i)=std(D);
            end
            
            % If all of V(i)<=Vspec, then Level 2 is accepted, otherwise proceed to Level 3.
            if all(V2<=Vspec)
                Level3=0;
            else
                Level3=1;
                S=S2;
                % M=M2;
                V=V2;
                Psub=Psub2;
                % if (verbose)
                %     display(V2)
                %     display(M0)
                %     display(M1)
                %     display(M2)
                % end
            end
        end %end of if ~isempty(idx) /line88
        
        if verbose
            disp('Level2: DONE')
            disp(['Level3=', num2str(Level3)])
        end
    end
    
    %% Level 3: U-shape (cf. page 55)
    if Level3
        
        % If Vi of subset Pi (from level 2)>Vspec, the aim of level 3 is to
        % replace a split-point s(i) with 5 new split-points (si, si1, si2, si3, si4)
        % i.e. split one segment into 5 smaller segments
        idx2=find(V>Vspec);
        if verbose
            disp('U-building')
            disp(['Number of level-2 split-points to be divided=' num2str(length(idx2))])
            disp(['Level-2-segments to be divided: ' num2str(idx2)])
        end
        
        % split the subset with variance>Vspec into 5
        si=zeros(length(idx2),2);
        si1=zeros(length(idx2),2);
        si2=zeros(length(idx2),2);
        si3=zeros(length(idx2),2);
        si4=zeros(length(idx2),2);
        si5=zeros(length(idx2),2); % I use si5 for the s(i+1)
        
        h1=zeros(length(idx2),2); % h1 is used to determine si1
        h4=zeros(length(idx2),2); % h4 is used to determine si4
        for i=1:length(idx2)
            p=Psub{idx2(i)};
            si(i,:)=p(1,:);
            si5(i,:)=p(end,:);
        end
        
        for i=1:length(idx2)
            % recalculate Di (i.e. distance from Pi to either MBR edge m(1) or m(2))
            k=idx2(i);
            p=Psub{k};
            n=size(p,1); % count of points in Psub(i)
            D_left=zeros(1,n);
            D_right=zeros(1,n);
            
            for j=1:n
                w=2-mod(k,2);
                if w==1
                    D_left(j) =dist_to_line(p(j,:),[M0(:,1) M0(:,2)]);
                    D_right(j)=dist_to_line(p(j,:),[M0(:,3) M0(:,4)]);
                else
                    D_left(j) =dist_to_line(p(j,:),[M0(:,2) M0(:,3)]);
                    D_right(j)=dist_to_line(p(j,:),[M0(:,4) M0(:,1)]);
                end
            end
            if mean(D_left)<mean(D_right)
                D=D_left;
            else
                D=D_right;
            end
            
            % Histogram of D(i) in 6 bin, in the order of ascending d
            [N,hist_edges] = histcounts(D,6);
            % if there are too few D for 6 bins, compute only 4 bins
            if N(end)<2
                [N,hist_edges] = histcounts(D,4);
            end
            
            % if verbose
            %     figure
            %     subplot(121), plot(D), title('Distances D(i)')
            %     subplot(122), hist(D,4), title('Histogram of D(i) in 6 bins')
            % end
            
            % only consider the last bin C6 which contains the largest values
            U=N(end)/n;
            
            Uref=0; %0.11; % Reference value (Dutter prefers an Uref of 0.11)
            % An U(i) smaller value than Uref means that the data is of poor quality
            % or the geometry is too complex for this method.
            % I think it is okay without Uref.
            
            if U<Uref % Ineligible for an U-shape
                
                disp(['Segment ' num2str(k) ', Uref=' num2str(Uref) ' not satisfied.'])
                disp('Reason: poor quality or geometry too complex.')
                disp('This segment should be manual edited')
                
            else % eligible for an U-shape
                
                % The first and last occurence of a point in C6 are the new
                % split-points si2 and si3, page 59
                si2(i,:)=p(find(D>=hist_edges(end-1) & D<=hist_edges(end),1),:);
                si3(i,:)=p(find(D>=hist_edges(end-1) & D<=hist_edges(end),1,'last'),:);
                
                % Find si1, page 59
                
                % h1 is the intersection between the line going through si2
                % parallel to m(i-1) with the line going to si and parallel to m(i)
                if k>1
                    h1(i,:)=parallel_intersection(si2(i,:)',[M2(:,k-1) M2(:,k)],...
                        si(i,:)', [M2(:,k) M2(:,k+1)],0);
                else
                    h1(i,:)=parallel_intersection(si2(i,:)',[M2(:,k+1) M2(:,k+2)],...
                        si(i,:)', [M2(:,k) M2(:,k+1)],0);
                end
                
                d_i1=zeros(1,n);
                for j=1:size(p,1)
                    d_i1(j)=dist(h1(i,:),p(j,:));
                end
                [~,idx_i1]=min(d_i1);
                si1(i,:)=p(idx_i1,:); % si1 is the points in the segment closest to h1
                
                % Find si4, page 59
                if k>1
                    h4(i,:)=parallel_intersection(si3(i,:)',[M2(:,k-1) M2(:,k)],si5(i,:)',[M2(:,k) M2(:,k+1)],0);
                else
                    h4(i,:)=parallel_intersection(si3(i,:)',[M2(:,k+1) M2(:,k+2)],si5(i,:)',[M2(:,k) M2(:,k+1)],0);
                end
                d_i4=zeros(1,n);
                for j=1:size(p,1)
                    d_i4(j)=dist(h4(i,:),p(j,:));
                end
                [~,idx_i4]=min(d_i4);
                si4(i,:)=p(idx_i4,:); % si4 is the points in the segment closest to h4
            end
            % plot(si1(i,1),si1(i,2),'m*')
            % plot(si2(i,1),si2(i,2),'b*')
            % plot(si3(i,1),si3(i,2),'k*')
            % plot(si4(i,1),si4(i,2),'g*')
        end
        
        % Add new split-points into S3
        % S3: split-point index (at level 3)
        S3=[];
        for j=1:length(idx2)
            idx_si =find(P(:,1)==si(j,1) & P(:,2)==si(j,2));
            
            if all(si2(j,:)~=[0 0])
                idx_si1=find(P(:,1)==si1(j,1) & P(:,2)==si1(j,2));
                idx_si2=find(P(:,1)==si2(j,1) & P(:,2)==si2(j,2));
                idx_si3=find(P(:,1)==si3(j,1) & P(:,2)==si3(j,2));
                idx_si4=find(P(:,1)==si4(j,1) & P(:,2)==si4(j,2));
            else % case U<Uref
                idx_si1=[];
                idx_si2=[];
                idx_si3=[];
                idx_si4=[];
            end
            % set of split-points
            S3b=S(1:idx2(j)-1);
            S3f=S(idx2(j)+1:end);
            S3=[S3; S3b; idx_si; sort([idx_si1, idx_si2, idx_si3, idx_si4])'; S3f];
        end
        S3=unique(S3);
        
        % Stock subsets of points into a struct
        Psub3={}; % struct of subsets of point (at level 3)
        for i=1:length(S3)
            if i<length(S3)
                if S3(i)>S3(i+1)
                    Psub3(i,:)={[P(S3(i):end,:); P(1:S3(i+1),:)]};
                else
                    Psub3(i,:)={P(S3(i):S3(i+1),:)};
                end
            else %i=length(S3)
                if S3(i)>S3(1) %S3(end)>S2(1)
                    Psub3(i,:)={[P(S3(i):end,:); P(1:S3(1),:)]};
                else
                    Psub3(i,:)={P(S3(i):S3(1),:)};
                end
            end
        end
        
        % set of MBR edges to be compared with the new split edges
        M3=repmat(M0,[1,length(S3)]);
        
        if verbose
            % display(S2)
            % display(M)
            % display(M2)
            % display(Psub2)
            figure
            plot(P(:,1), P(:,2), 'b.')
            hold on, grid on, axis equal
            plot([mbr(1,1:4) mbr(1,1)],[mbr(2,1:4) mbr(2,1)],'k--')
            cl=rand(length(S3)+1,3);
            for i=1:length(S3)
                % text(P(S3(i),1), P(S3(i),2), ['s_{' num2str(i) '}']) % split-points
                p=Psub3{i};
                plot(p(:,1),p(:,2),'.','Color',cl(i,:),'MarkerSize',10)
            end
            title('Level 3 result')
        end
    end
    
    
    %% Saving results
    if Level3
        S=S3;
        M=M3;
        Psub=Psub3;
    elseif Level2
        S=S2;
        M=M2;
        Psub=Psub2;
    else
        S=S1;
        M=M1;
        Psub=Psub1;
    end
    % display(S)
    % display(M)
    
    %% Compute the position of edges (each being parallel to one MBR edge)
    % The resulting edge is distanced from the corresponding MBR edge.
    % Such a distance (denoted by e) is computed by the median of points to the MBR edge.
    e=zeros(1,length(S));
    for i=1:length(S)
        p=Psub{i};
        D=zeros(1,size(p,1));
        for j=1:size(p,1)
            D(j)=dist_to_line(p(j,:),[M(:,i) M(:,i+1)]);
        end
        e(i)=median(D);
    end
    % display(e)
    
    % Decide between MBR_edge(i)+e(i) or MBR_edge(i)-e(i)
    edges={};
    for i=1:length(S)
        [p1,q1]=makeParallelLine([M(:,i) M(:,i+1)], e(i));
        [p2,q2]=makeParallelLine([M(:,i) M(:,i+1)], -e(i));
        
        % Get the line that has lower average distance to points of the edge
        p=Psub{i};
        D_left=zeros(1,size(p,1));
        D_right=zeros(1,size(p,1));
        for j=1:size(p,1)
            D_left(j)=dist_to_line(p(j,:),[p1,q1]);
            D_right(j)=dist_to_line(p(j,:),[p2,q2]);
        end
        if (mean(D_left)<mean(D_right))
            edges(i,:)={real(p1), real(q1)};
        else
            edges(i,:)={real(p2), real(q2)};
        end
    end
    
    vertices=zeros(2,length(S));
    for i=1:length(S)
        % display([edges{i,1} edges{i,2}])
        if i<length(S)
            vertices(:,i)=line_intersect([edges{i,1} edges{i,2}],[edges{i+1,1} edges{i+1,2}]);
        else
            vertices(:,i)=line_intersect([edges{i,1} edges{i,2}],[edges{1,1} edges{1,2}]);
        end
    end
    vertices=[vertices(:,end) vertices];
    % display(vertices)
    
    if verbose
        figure
        hold on, grid on, axis equal
        cl=rand(length(S)+1,3);
        for i=1:length(S)
            p=Psub{i};
            plot(p(:,1),p(:,2),'.','Color',cl(i,:))
            text(mean(p(:,1)),mean(p(:,2)),num2str(i),'Color',cl(i,:))
            plot([vertices(1,i) vertices(1,i+1)],[vertices(2,i) vertices(2,i+1)],'Color',cl(i,:));
        end
    end
    
    % Check variance of distances D from subset points to the corresponding MBR edge
    M3=vertices;
    V3=zeros(1,length(S));
    for i=1:length(S)
        p=Psub{i};
        D=zeros(1,size(p,1));
        for j=1:size(p,1)
            D(j)=dist_to_line(p(j,:),[M3(:,i) M3(:,i+1)]);
        end
        V3(i)=std(D);
    end
    
    iter=iter+1;
    
    if all(V3<=Vspec) || (iter>2) % only run twice
        stop_flag=1;
    else
        stop_flag=0;
        M=M3; V=V3; 
    end
    
end


if verbose
    figure
    p1=plot(P(:,1), P(:,2), 'b.');
    hold on, grid on, axis equal
    p2=plot([mbr(1,1:4) mbr(1,1)],[mbr(2,1:4) mbr(2,1)],'k--');
    for i=1:length(S)
        p3=plot([edges{i,1}(1) edges{i,2}(1)],[edges{i,1}(2) edges{i,2}(2)], 'g--');
        % text(mean([edges{i,1}(1),edges{i,2}(1)]),mean([edges{i,1}(2),edges{i,2}(2)]),num2str(i))
    end
    p4=plot(vertices(1,:),vertices(2,:),'r*');
    p5=plot([vertices(1,:) vertices(1,1)],[vertices(2,:) vertices(2,1)],'r');
    title('Final result (verbose mode)')
    legend([p1,p2,p3,p4,p5],{'Input points','MBR','Parallel edges','Final vertices','Final edges'})
end

figure
plot(P(:,1), P(:,2), 'b.', 'MarkerSize', 10)
hold on, grid on, axis equal
plot([mbr(1,1:4) mbr(1,1)],[mbr(2,1:4) mbr(2,1)],'k--','LineWidth',2)
plot(vertices(1,:),vertices(2,:),'r*','MarkerSize', 10)
plot([vertices(1,:) vertices(1,1)],[vertices(2,:) vertices(2,1)],'r','LineWidth',2)
title('Final result')
l=legend('Input points','MBR','Final vertices','Final edges');
l.FontSize=14;

%% Resampling each resulting segment (given nb_pts>0)
points=[];
if nb_pts>0
    for i=1:length(S)
        pts=zeros(nb_pts,2);
        if i<length(S)
            for j=1:2
                pts(:,j)=linspace(vertices(j,i),vertices(j,i+1),nb_pts)';
            end
        else
            for j=1:2
                pts(:,j)=linspace(vertices(j,i),vertices(j,1),nb_pts)';
            end
        end
        points=[points; pts];
    end
    if nb_pts<2
        disp('Warning: number of points for resampling should be greater than 2, preferably > 20')
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%
% END OF MAIN PROGRAM %
%%%%%%%%%%%%%%%%%%%%%%%

% Sub-functions

%% Euclid distance between two points P and Q
function d=dist(p,q)
if all(size(p)==size(q))
    d=norm(p-q);
else
    d=norm(p'-q);
end
end

%% Distance d(M,l) from point M to line l (given 2 end-points P and Q)
% Denote alpha: angle between PM and PQ
% cos(alpha)=vector(PM)*vector(PQ)/(norm(vector(PM))*norm(vector(PQ)))
% d(M,l)=norm(vector(PM))*sin(alpha)
function d=dist_to_line(m,l)
p=l(:,1)';
q=l(:,2)';
d = norm(m-p)*sin(acos(((m-p)*(q-p)')/(norm(m-p)*norm(q-p))));
end

%% Rotation
function mat_rotated=myRotate(origin,angle,mat)
R= [cos(angle) sin(angle);
    -sin(angle) cos(angle)];
mat_centered = mat-origin;
mat_rotated_centered = R*mat_centered;
mat_rotated = mat_rotated_centered+origin;
end

%% Intersection of two lines l1 and l2
% l1=[xA,xB;yA,yB] and l2=[xC,xD;yC,yD] are given by their end-points
function m=line_intersect(l1,l2)
% linear polynomial fitting
p1 = polyfit([l1(1,1) l1(1,2)],[l1(2,1) l1(2,2)],1);
p2 = polyfit([l2(1,1) l2(1,2)],[l2(2,1) l2(2,2)],1);

% calculate intersection
x_intersect = fzero(@(x) polyval(p1-p2,x),1);
y_intersect = polyval(p1,x_intersect);

m=[x_intersect y_intersect]';
end

%% Intersection between two lines (m1,//l1) and (m2,//l2)
% Given two points m1 and m2, h is the intersection between the line through
% m1 and parallel to l1 (denoted //l1) and the line through m2 and //l2.
% In this algorithm, l1 and l2 should be perpendicular.
function h=parallel_intersection(m1,l1,m2,l2,verbose)
m=line_intersect(l1,l2);

angle=atan((l1(2,1)-l1(2,2))/(l1(1,1)-l1(1,2)));  % atan((yA-yB)/(xA-xB))

m1_rotated=myRotate(m,angle,m1);
m2_rotated=myRotate(m,angle,m2);

tmp=myRotate(m,-angle,[m2_rotated(1),m1_rotated(2)]');
h=tmp';

if verbose
    figure
    hold on
    plot(m1(1),m1(2),'r.')
    plot(m2(1),m2(2),'b.')
    text(m1(1),m1(2),'m1')
    text(m2(1),m2(2),'m2')
    
    plot(m(1),m(2),'g.')
    plot([l1(1,1),l1(1,2)],[l1(2,1),l1(2,2)],'r-')
    plot([l2(1,1),l2(1,2)],[l2(2,1),l2(2,2)],'b-')
    plot(h(1),h(2),'g*')
    axis equal
    grid on
end
end

%% Make a line that are parallel to a baseline with a distance d
function [p,q]=makeParallelLine(baseline, d)
p0=baseline(:,1);
q0=baseline(:,2);

angle=atan((q0(2)-p0(2))/(q0(1)-p0(1)));
q0_prime=myRotate(p0,angle,q0);
q_prime=zeros(2,1);
q_prime(1)=q0_prime(1);
q_prime(2)=q0_prime(2)+d;
q=myRotate(p0,-angle,q_prime);

p0_prime=myRotate(q0,angle,p0);
p_prime=zeros(2,1);
p_prime(1)=p0_prime(1);
p_prime(2)=p0_prime(2)+d;
p=myRotate(q0,-angle,p_prime);
end

%%%%%%%
% END %
%%%%%%%
