function [ dist ] = shortestPath( dist )
%SHORTESTPATH Summary of this function goes here
%   Detailed explanation goes here
N=size(dist,1);


for k=1:N
    for i=1:N
        for j=1:N
            if dist(i,j)>dist(i,k)+dist(k,j)
                dist(i,j)=dist(i,k)+dist(k,j);
            end
        end
    end
end

end

