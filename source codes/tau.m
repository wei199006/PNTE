% The function segments the pore space from a cast thin section image.
% Inputs:
    % sp- the inlets of the pore network;
    % ep- the outlets of the pore network;
    % dis- the distance-weighted graph after applying the Floyd-Warshall algorithm;
    % dtype- =1 width of the image, =2 straight-line distance between the inlet and the outlet;
    % ttype- =1 average tortuosity among all pi, =2 minimum tortuosity among all pi;
    % qq- the column numbers of the cast thin section image.
% Outputs:
    % taupi- tortuosity of the shortest path pi;
    % taun- tortuosity of the pore network;
    % pi- the outlets referred to each pi;
    % nt- the indices of the minimum tortuosity among all pi.
function [taupi,taun,pi,nt]=tau(sp,ep,dis,dtype,ttype,qq)
nsp= length(sp); % number of inlets
nep= length(ep); % number of outlets
minv= 1./zeros(nsp,1); % A intermediate parameter used to acquired pi
pi= zeros(nsp,1);
taupi= 1./zeros(nsp,1); % tortuosity of pi
for i=1:nsp
    minv(i)= unique(min(dis(nsp+1:nsp+nep,i)));
    if minv(i)==Inf
        continue;
    end
    findmin= find(dis(nsp+1:nsp+nep,i)==minv(i));
    pi(i)=findmin(1);
    % calculate the tortuosity for each pi with the straight-line distance defined by Eqn. (3)
    if abs(dtype-1)<10^(-10)
        taupi(i)= dis(nsp+pi(i),i)/(qq-1);
    end
    if abs(dtype-2)<10^(-10)
        taupi(i)= dis(nsp+pi(i),i)/((sp(i,1)-ep(pi(i),1)).^2+(qq-1)^2).^0.5; % = Inf means not connected
    end
end
% Calculate tortuosity of pore network according to Eqn. (4)
if abs(ttype-1)<10^(-10)
    taun= mean(taupi(~isinf(taupi)));
    nt= find(taupi==min(taupi));
end
if abs(ttype-2)<10^(-10)
    taun= min(taupi);
    nt= find(taupi==min(taupi)); 
end
end