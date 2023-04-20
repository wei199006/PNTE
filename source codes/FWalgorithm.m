% The function is Floyd-Warshall algorithm for solving all pairs shortest path problems.
% Inputs:
    % dis- the distance-weighted graph;
    % tnote- the matrix of the intermediate vertices used in the Floyd-Warshall algorithm;
    % pois- the collections (vertices） of inlets, outlets, and intersectional points & blocks of the pore network.
% Outputs:
    % dis- the distance-weighted graph after applying the Floyd-Warshall algorithm;
    % tnote- the matrix of the intermediate vertices after applying the Floyd-Warshall algorithm.
function [dis,tnote]=FWalgorithm(dis,tnote,pois)
[pm,pn]=size(pois);
disnew=dis;
for bm=1:pm
    for bi= 1:pm
        for bj=1:pm
            if bi==bj
                continue;
            end
            % Eqn.(9)
            if (dis(bi,bm)+dis(bm,bj)-dis(bi,bj)< -10^(-10))
                disnew(bi,bj)= dis(bi,bm)+dis(bm,bj);
                tnote(bi,bj)= tnote(bi,bm);
            end
        end
    end
    dis=disnew;
end
end