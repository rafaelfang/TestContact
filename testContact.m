function [totalNumOfChange]=testContact(filename)


fID=fopen(filename);
[pdbData] = textscan(fID,'%s %d %s %s %d %f %f %f %f %f %s');
fclose(fID);
T(:,1)=pdbData{1,6};
T(:,2)=pdbData{1,7};
T(:,3)=pdbData{1,8};
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
          distM3(i,j)=12;
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

accuracyPercent=1;
if accuracyPercent==1
    
else
    noisyPercent=1-accuracyPercent;
    totalNoisyPoints=floor(noisyPercent*size(distM3,1)*size(distM3,1));
    numOfNoisyPoints=totalNoisyPoints;
    flag=zeros(size(distM3,1),size(distM3,1));
    while(numOfNoisyPoints>0)
        xCoord=randi(size(distM3,1));
        yCoord=randi(size(distM3,1));
        if(xCoord==yCoord||xCoord==yCoord+1||xCoord+1==yCoord||flag(xCoord,yCoord)==1)
            continue;
        else
            if(distM3(xCoord,yCoord)==12)
                distM3(xCoord,yCoord)=6;
                flag(xCoord,yCoord)=1;
                distM3(yCoord,xCoord)=6;
                flag(yCoord,xCoord)=1;
            else
                distM3(xCoord,yCoord)=12;
                flag(xCoord,yCoord)=1;
                distM3(yCoord,xCoord)=12;
                flag(yCoord,xCoord)=1;
            end
            
        end
        numOfNoisyPoints= numOfNoisyPoints-1;
        
    end
end

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
%P=mdscale(distM3,3);
first=char(pdbData{1,1});
second=pdbData{1,2};
third=char(pdbData{1,3});
forth=char(pdbData{1,4});
fifth=pdbData{1,5};
ninth=pdbData{1,9};
tenth=pdbData{1,10};
eleventh=char(pdbData{1,11});

mytable=table(first,second,third,forth,fifth,P(:,1),P(:,2),P(:,3),ninth,tenth,eleventh);
writetable(mytable,strcat('g1_',filename),'FileType','text','Delimiter','\t');

fileID=fopen(strcat('g1_',filename),'w');
for i=1:size(first,1)
    
   str=sprintf('ATOM  %5d %4s %3s   %3d    %8.3f%8.3f%8.3f%6.2f%6.2f           %s',...
                 second(i,:),third(i,:),forth(i,:),fifth(i,:),P(i,1),P(i,2),P(i,3),...
                    ninth(i,:),tenth(i,:),eleventh(i,:)); 
   fprintf(fileID,'%s\n',str);
end
fclose(fileID);




iteration=100;
candidateDistM=cell(iteration,1);
totalNumOfChange=zeros(iteration,1);
candidateDistM{1,1}=distM3;
for iter=2:iteration
    counter=0;
    temp=pdist2(P(:,1:3),P(:,1:3));
    for row=1:size(mask,1)
       for col=1:size(mask,1) 
            if((mask(row,col)==1&&temp(row,col)<8)||(mask(row,col)==0&&temp(row,col)>8))
                
            elseif(mask(row,col)==1&&temp(row,col)>8)
                counter=counter+1;
                temp(row,col)=6;
            elseif(mask(row,col)==0&&temp(row,col)<8)
                temp(row,col)=12;
                counter=counter+1;
            end
       end
    end
    totalNumOfChange(iter,1)=counter;
    %set diagnol value to 0;
    for i=1:size(temp,1)
        temp(i,i)=0;
    end  
    %set +1/-1 diagnol value to 3.8
    for i=1:size(temp,1)-1
        temp(i,i+1)=3.8;
        temp(i+1,i)=3.8;
    end
    candidateDistM{iter,1}=temp;
    P=cmdscale(temp);
    %P=mdscale(temp,3);
    first=char(pdbData{1,1});
    second=pdbData{1,2};
    third=char(pdbData{1,3});
    forth=char(pdbData{1,4});
    fifth=pdbData{1,5};
    ninth=pdbData{1,9};
    tenth=pdbData{1,10};
    eleventh=char(pdbData{1,11});

    mytable=table(first,second,third,forth,fifth,P(:,1),P(:,2),P(:,3),ninth,tenth,eleventh);
    writetable(mytable,strcat('g',num2str(iter),'_',filename),'FileType','text','Delimiter','\t');

    fileID=fopen(strcat('g',num2str(iter),'_',filename),'w');
    for i=1:size(first,1)

       str=sprintf('ATOM  %5d %4s %3s   %3d    %8.3f%8.3f%8.3f%6.2f%6.2f           %s',...
                     second(i,:),third(i,:),forth(i,:),fifth(i,:),P(i,1),P(i,2),P(i,3),...
                        ninth(i,:),tenth(i,:),eleventh(i,:)); 
       fprintf(fileID,'%s\n',str);
    end
    fclose(fileID);
    
end


end
