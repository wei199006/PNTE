% The function traces the shortest paths.
% Inputs:
    % mm- the total number of rows of the final pore network image;
    % sp- the inlets of the pore network;
    % pi- the outlets referred to each pi;
    % bw- the final pore network image;
    % tnote- the matrix of the intermediate vertices used in the Floyd-Warshall algorithm;
    % tpath- the fork paths of each vertex in the pore network;
    % tbns- the vertex number of the fork ends referred to all the vertices in the pore network;
    % nt- the indices of the minimum tortuosity among all pi.
% Outputs:
    % nettrace- the indices of the paths of all pi;
    % netpath- the paths of all pi.
function [nettrace,netpath]=SPtrace(mm,sp,pi,bw,tnote,tpath,tbns,nt)
for i=1:length(sp)
    startp=i;
    if pi(i)>10^(-6)
        endp=length(sp)+pi(i);
        [trace,spath]= tracepath(startp,endp,tnote,tpath,tbns,mm);
        nettrace{i}=trace;
        netpath{i}=spath;
    end
    if abs(pi(i))<10^(-6)
        nettrace{i}= NaN;
        netpath{i}= NaN;
    end
end
end