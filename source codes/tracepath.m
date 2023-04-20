% The function figures out the trace and path between two vertices.
% Inputs:
    % startp- start pixel;
    % endp- end pixel;
    % tnote- the matrix of the intermediate vertices used in the Floyd-Warshall algorithm;
    % tpath- the fork paths of each vertex in the pore network;
    % tbns- the vertex number of the fork ends referred to all the vertices in the pore network;
    % mm- the total number of rows of the final pore network image.
% Outputs:
    % trace- the distance-weighted graph;
    % spath- the matrix of the intermediate vertices used in the Floyd-Warshall algorithm.
function [trace,spath]=tracepath(startp,endp,tnote,tpath,tbns,mm)
% build up the trace vector
trace(1)= endp;
trace(2)= tnote(endp,startp);
for i=1:5*mm-2
    if ( abs(trace(i+1)-startp)<10^(-10) )
        break;
    end
    trace(i+2)= tnote(trace(i+1),startp);
    if ( isnan(trace(i+2)) )
        trace(i+2)=startp;
        break;
    end
end
% delete the zero elements.
trace=trace(trace>0);
spath=zeros(1,2);
% calculate the path from the trace..
for i=1:length(trace)-1
    st= tpath{trace(i)}; % st- the referred to 
    sb= tbns{trace(i)}; % sb- the referred to 
    for j=1:length(sb)
        if (abs(sb(j)-trace(i+1))<10^(-6))
            se= st(:,2*j-1:2*j);
            se(any(se,2)==0,:)=[]; % delete the zeros
            spath=[spath;se];
            spath(any(spath,2)==0,:)=[]; % delete the zeros.
        end
    end
end
spath = unique(spath,'rows','stable');  
end