% The function provides the paths, lengths, and ends of all the path forks referred to the vertex di.
% Inputs:
    %di- the vertex di;
    %bw- the final pore network;
    %pois- the collections (vertices) of inlets, outlets, and intersectional points & blocks of the pore network;
    %blocks- the intersectional blocks of the pore network.
% Outputs:
    %outpath- the fork paths of the vertex di;
    %outlpath- the lengths of the fork paths of the vertex di;
    %branchends- the vertex number of the fork ends referred to all the vertices in the pore network.
function [outpath,outlpath,branchends]=neiberdis1(di,bw,pois,blocks)
[mm,nn] = size(bw);
nnn= length(pois(:,1));
if di<nnn-length(blocks)+0.1
    x=pois(di,1);% x coordinate of di vertex
    y=pois(di,2);% y coordinate of di vertex
    k=1;
    for adi= 1:3
        for adj= 1:3
            % skip the pixels exceeded the image boundaries
            if (x+adi-2-mm>10^(-6) || x+adi-2-1<-10^(-6) || y+adj-2-nn>10^(-6) || y+adj-2-1<-10^(-6))
                continue;
            end
            % ignore the pixel itself
            if (adi==2 && adj==2)
                continue;
            end
            if abs(bw(x+adi-2,y+adj-2)-1)<10^(-10)% search for the starts of the path forks.
                % paths- a vector saves the start pixels of the path forks referred to the vertex di
                paths(k,:)=[x+adi-2 y+adj-2];
                k=k+1;
            end
        end
    end
    npath= k-1; % number of path forks
end
if di>nnn-length(blocks)+0.1
    % figure out the contours of the blocks
    bl= blocks{di-nnn+length(blocks)};
    bl1=[bl(:,1)-1 bl(:,2)-1];
    bl2=[bl(:,1)+1 bl(:,2)-1];
    bl3=[bl(:,1)-1 bl(:,2)+1];
    bl4=[bl(:,1)+1 bl(:,2)+1];
    bl=union(bl,bl1,'rows');
    bl=union(bl,bl2,'rows');
    bl=union(bl,bl3,'rows');
    bl=union(bl,bl4,'rows');
    abl=setdiff(bl,blocks{di-nnn+length(blocks)},'rows');% abl-a contour of a block
    k=1;
    % search the start pixels of the path forks referred to the intersectional blocks
    for i=1:length(abl)
        % skip the pixels exceeded the image boundaries
        if (abl(i,1)-mm>10^(-6) || abl(i,1)-1<-10^(-6) ||...
                abl(i,2)-nn>10^(-6) || abl(i,2)-1<-10^(-6))
            continue;
        end
        if abs(bw(abl(i,1),abl(i,2))-1)<10^(-10)
            paths(k,:)=[abl(i,1),abl(i,2)];
            k=k+1;
        end
    end
    npath= k-1;% number of path forks
end
% initiation of the path matrix
maxl= max(nn,mm);
outpath= zeros(maxl,2*npath);
outlpath=zeros(npath,1);
branchends=zeros(npath,1);
% determine the paths referred to di one by one
for k=1:npath
    flag= 0;% a flag parameter used to break the circulates
    path= zeros(maxl,2);
    % the first pixels of all the paths linked to di are the vertex di itself
    path(1,1)= pois(di,1);
    path(1,2)= pois(di,2);
    % the second pixels of all the paths are the start points of all the path forks
    path(2,1)= paths(k,1);
    path(2,2)= paths(k,2);
    % if a second point is linked to any other vertex in the vector pois, break the circulates by letting flag=1
    for ii=1:nnn
        if( (path(2,1)-pois(ii,1))^2+(path(2,2)-pois(ii,2))^2<0.1 )
            branchends(k)= ii;
            outpath(:,2*k-1)=path(:,1);
            outpath(:,2*k)=path(:,2);
            path1=path;
            path1(all(path1==0,2),:)= []; % delete the zero rows
            nnn1=length(path1);
            sumlength=0;
            for kk=1:nnn1-1
                diffd= ((path1(kk,1)-path1(kk+1,1))^2.0+(path1(kk,2)-path1(kk+1,2))^2.0)^0.5;
                sumlength= sumlength+diffd; % length of a path fork
            end
            outlpath(k)=sumlength;
            flag= 1;
            break;
        end
    end
    if(abs(flag-1)<10^(-10))
        continue;
    end
    % the third pixels of all the paths are the start points of all the path forks
    for adi= 1:3
        for adj= 1:3
            if (adi== 2 && adj== 2)
                continue;
            end
            % use the term (adi+adj)/2 to divide the 8-neighbourhood pixels into two groups: 4-neighbourhood pixels and diagonal pixels.
            % calculate the 4-neighbourhood pixels first by skipping the diagonal pixels
            if (mod(adi+adj,2)==1)
                continue;
            end
            nsx=path(2,1)+adi-2;
            nsy=path(2,2)+adj-2;
            % avoid exceding the image edges
            if (nsx-mm>10^(-6) || nsx-1<-10^(-6) || nsy-nn>10^(-6) || nsy-1<-10^(-6))
                continue;
            end
            % avoid tracing back
            if ( (nsx-path(1,1))^2+(nsy-path(1,2))^2<2.1 )
                continue;
            end
            if (bw( nsx, nsy )== 1)
                path(3,1)= nsx;
                path(3,2)= nsy;
            end
        end
    end
    for adi= 1:3
        for adj= 1:3
            if (adi== 2 && adj== 2)
                continue;
            end
            % calculate the diagonal pixels by skipping the 4-neighbourhood pixels
            if (mod(adi+adj,2)==0)
                continue;
            end
            nsx=path(2,1)+adi-2;
            nsy=path(2,2)+adj-2;
            % avoid exceding the image edges
            if (nsx-mm>10^(-6) || nsx-1<-10^(-6) || nsy-nn>10^(-6) || nsy-1<-10^(-6))
                continue;
            end
            % avoid tracing back
            if ( (nsx-path(1,1))^2+(nsy-path(1,2))^2<2.1 )
                continue;
            end
            if (bw( nsx, nsy )== 1)
                path(3,1)= nsx;
                path(3,2)= nsy;
            end
        end
    end
    % if a third point is linked to any other vertex in the vector pois, break the circulates by letting flag=1
    for ii=1:nnn
        if( (path(3,1)-pois(ii,1))^2+(path(3,2)-pois(ii,2))^2<0.1 )
            branchends(k)= ii;
            outpath(:,2*k-1)=path(:,1);
            outpath(:,2*k)=path(:,2);
            path1=path;
            path1(all(path1==0,2),:)= [];
            nnn1=length(path1);
            sumlength=0;
            for kk=1:nnn1-1
                diffd= ((path1(kk,1)-path1(kk+1,1))^2.0+(path1(kk,2)-path1(kk+1,2))^2.0)^0.5;
                sumlength= sumlength+diffd;
            end
            outlpath(k)=sumlength;
            flag= 1;
            break;
        end
    end
    if(abs(flag-1)<10^(-10))
        continue;
    end
    for i=1:maxl-3 % the maximum length of a preset path segment is set to be the image edge length
        flag=0;
        for ii=1:nnn
            if( (path(i+2,1)-pois(ii,1))^2+(path(i+2,2)-pois(ii,2))^2<2.1 )
                flag=1;
                branchends(k)= ii;
                path(i+3,1)=pois(ii,1);
                path(i+3,2)=pois(ii,2);
                break;
            end
            
        end
        if(abs(flag-1)<10^(-10))
            break;
        end
        flagg=1;
        for adi= 1:3
            for adj= 1:3
                if (adi== 2 && adj== 2)
                    continue;
                end
                % calculate the 4-neighbourhood pixels first by skipping the diagonal pixels
                if (mod(adi+adj,2)==1)
                    continue;
                end
                nsx=path(i+2,1)+adi-2;
                nsy=path(i+2,2)+adj-2;
                % avoid exceding the image edges
                if (nsx-mm>10^(-6) || nsx-1<-10^(-6) ||...
                        nsy-nn>10^(-6) || nsy-1<-10^(-6))
                    continue;
                end
                % avoid tracing back
                if ( (nsx== path(i+2,1) && nsy== path(i+2,2)) ||...
                        (nsx== path(i+1,1) && nsy== path(i+1,2)) ||...
                        (nsx== path(i,1) && nsy== path(i,2)) )
                    continue;
                end
                if (bw( nsx, nsy )== 1)
                    path(i+3,1)= nsx;
                    path(i+3,2)= nsy;
                    flagg=0;
                    break;
                end
            end
        end
        for adi= 1:3
            for adj= 1:3
                if (adi== 2 && adj== 2)
                    continue;
                end
                % calculate the diagonal pixels by skipping the 4-neighbourhood pixels
                if (mod(adi+adj,2)==0)
                    continue;
                end
                nsx=path(i+2,1)+adi-2;
                nsy=path(i+2,2)+adj-2;
                % avoid exceding the image edges
                if (nsx-mm>10^(-6) || nsx-1<-10^(-6) || nsy-nn>10^(-6) || nsy-1<-10^(-6))
                    continue;
                end
                % avoid tracing back
                if ( (nsx== path(i+2,1) && nsy== path(i+2,2)) ||...
                        (nsx== path(i+1,1) && nsy== path(i+1,2)) ||...
                        (nsx== path(i,1) && nsy== path(i,2)) )
                    continue;
                end
                if (bw( nsx, nsy )== 1)
                    path(i+3,1)= nsx;
                    path(i+3,2)= nsy;
                    flagg=0;
                    break;
                end
            end
        end
        if(flagg==1)
            break;
        end
    end
    outpath(:,2*k-1)=path(:,1);
    outpath(:,2*k)=path(:,2);
    path1=path;
    path1(all(path1==0,2),:)= [];
    nnn1=length(path1);
    sumlength=0;
    for kk=1:nnn1-1
        diffd= ((path1(kk,1)-path1(kk+1,1))^2.0+(path1(kk,2)-path1(kk+1,2))^2.0)^0.5;
        sumlength= sumlength+diffd;
    end
    outlpath(k)=sumlength;
end
end