% The function deals with the cases that more than one path exist between two vertices.
% The function is designed to figure out the shortest path among all the paths linking the two vertices.
% Inputs:
    %path- the fork paths of the vertex di;
    %lpath- the lengths of the fork paths of the vertex di;
    %bns- the vertex number of the fork ends referred to all the vertices in the pore network
% Outputs:
    %opath- the output fork paths of the vertex di;
    %olpath- the output lengths of the fork paths of the vertex di;
    %obns- the output vertex number of the fork ends referred to all the vertices in the pore networ.
function [opath,olpath,obns]=drpath(path,lpath,bns)
% check if more than one path exist between two vertices.
X = unique(bns);
if abs(length(bns)-length(X))>0.5
    for i=1:length(X)
        m{i}= find(bns==X(i));
    end
    for i= 1:length(m)
        if (length(m{i})-1.5>0)
            pl= lpath(m{i});
            sbns=m{i};
            mint= min(pl);
            mnum= find(pl==mint);
            if length(mnum)>1.5
                mnum=mnum(1);
            end
            obns(i)= bns(sbns(mnum));
            olpath(i)= lpath(sbns(mnum));
            opath(:,2*i-1:2*i)= path(:,2*sbns(mnum)-1:2*sbns(mnum));
        else
            obns(i)= bns(m{i});
            olpath(i)= lpath(m{i});
            opath(:,2*i-1:2*i)= path(:,2*m{i}-1:2*m{i});
        end
    end
else
    opath=path;
    olpath=lpath;
    obns=bns;
end
end