close all
clear all
clc

fID=fopen('T0711.txt');
[pdbData] = textscan(fID,'%s %d %s %s %d %f %f %f %f %f %s');
fclose(fID);
T(:,1)=pdbData{1,6};
T(:,2)=pdbData{1,7};
T(:,3)=pdbData{1,8};
figure
plot3(T(:,1),T(:,2),T(:,3))
title('native structure')


distM1=pdist2(T,T);

contactM2=distM1;
for i=1:size(contactM2,1)
   for j=1:size(contactM2,1)
      if(contactM2(i,j) <=8)
          contactM2(i,j)=1;
      else
          contactM2(i,j)=0;
      end
   end 
end
   
distM3=distM1;
for i=1:size(distM3,1)
   for j=1:size(distM3,1)
      if(distM3(i,j) <=8)
          distM3(i,j)=6;
      else
          distM3(i,j)=Inf;
      end
   end 
end
%set diagnol value to 0;
for i=1:size(distM3,1)
    distM3(i,i)=0;
end  
%set +1/-1 diagnol value to 3.8
for i=1:size(distM3,1)-1
    distM3(i,i+1)=3.8;
    distM3(i+1,i)=3.8;
end
distM3=shortestPath(distM3);

mask=contactM2;
%set diagnol value to 0;
for i=1:size(mask,1)
    mask(i,i)=0;
end  
%set +1/-1 diagnol value to 3.8
for i=1:size(mask,1)-1
    mask(i,i+1)=3.8;
    mask(i+1,i)=3.8;
end
P=cmdscale(distM3);
figure
imagesc(distM3)
colorbar
title('heatmap of distMat reconstructing 1st model')
[~,Z]=procrustes(T,P(:,1:3));
figure
plot3(Z(:,1),Z(:,2),Z(:,3));
title('1st reconstruct structure')

iteration=3;
for iter=2:iteration
   
    temp=pdist2(P(:,1:3),P(:,1:3));
    for row=1:size(mask,1)
       for col=1:size(mask,1) 
            if((mask(row,col)==1&&temp(row,col)<8)||(mask(row,col)==0&&temp(row,col)>8))
                
            elseif(mask(row,col)==1&&temp(row,col)>8)
                temp(row,col)=6;
            elseif(mask(row,col)==0&&temp(row,col)<8)
                temp(row,col)=Inf;
            end
       end
    end
    %set diagnol value to 0;
    for i=1:size(temp,1)
        temp(i,i)=0;
    end  
    %set +1/-1 diagnol value to 3.8
    for i=1:size(temp,1)-1
        temp(i,i+1)=3.8;
        temp(i+1,i)=3.8;
    end
    temp=shortestPath(temp);
    P=cmdscale(temp);
    figure
    imagesc(temp);
    colorbar
    if iter==2
        title('heatmap of distMat reconstructing 2nd model')
    elseif iter==3
        title('heatmap of distMat reconstructing 3rd model')
    end
    [~,Z]=procrustes(T,P(:,1:3));
    figure
    plot3(Z(:,1),Z(:,2),Z(:,3));
    if iter==2
        title('2nd reconstruct structure')
    elseif iter==3
        title('3rd reconstruct structure')
    end

 
    
end
