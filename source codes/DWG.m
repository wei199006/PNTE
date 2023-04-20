% The function segments the pore space from a cast thin section image.
% Inputs:
    % bw- the final pore network;
    % pois- the collections (vertices） of inlets, outlets, and intersectional points & blocks of the pore network.
% Outputs:
    % dis- the distance-weighted graph;
    % tnote- the matrix of the intermediate vertices used in the Floyd-Warshall algorithm;
    % tpath- the fork paths of each vertex in the pore network;
    % tbns- the vertex number of the fork ends referred to all the vertices in the pore network.
function [dis,tnote,tpath,tbns]=DWG(bw,pois)
[pm,pn]=size(pois);
dis=1./zeros(pm,pm); % the distance-weighted graph
tnote=ones(pm,pm); % the pore path matrix
for ti=1:pm
    tnote(:,ti)=ti; % initiate tnote.
end
for ti=1:pm
    tnote(ti,ti)= NaN;% NaN means that no intermediate vertex exists between the vertex ti and itself.
    dis(ti,ti)=0;% let the distance from a vertex to itself being zero.
end
for di=1:pm % di- the vertex di
    % if there are some the intersectional blocks, the function used is slightly different.
    % path1- the fork paths of the vertex di
    % lpath1- the lengths of the fork paths of the vertex di
    % bns1- the vertex number of the fork ends referred to all the vertices in the pore network
    if exist('blocksxy','var')
        [path1,lpath1,bns1]=neiberdis1(di,bw,pois,blocks);
    else
        [path1,lpath1,bns1]=neiberdis2(di,bw,pois);
    end
    [path,lpath,bns]=drpath(path1,lpath1,bns1);
    tpath{di}=path;
    tbns{di}=bns;
    for dj=1:length(bns)
        if (bns(dj)~=0)
            dis(di,bns(dj))=lpath(dj); % generate the distance-weighted graph
        end
    end
end
for ti=1:pm % let the distance from ti to itself being zero
    dis(ti,ti)=0;
end
end