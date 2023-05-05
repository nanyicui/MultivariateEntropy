function H=binless_entropy_multivariate(spikeTrain,D)

%Computes differential entropy using binless estimation method
%D is the dimension, e.g. D=2 is joint differential entropy


C=size(spikeTrain,2);%C is number of cells, for matrices the number of columns is number of cells              

indices=ones(1,size(spikeTrain,2));   %for later use
    
S=D*pi^(D/2)/(gamma(D/2+1));

if iscell(spikeTrain)==1  %this is only for the MUA data with 16 cells that contain structs
    
  for c=1:C
       
    x{c}=((spikeTrain{c}.spikeTimes)); %rearranging data into a cell structure
                                                                           
    x{c}=unique_no_sort(x{c});% a code from Caitlin Bever 2007 MIT that does what is said on the tin
    
    xsize(c)=size(x{1,c},2);

  end
  
elseif ismatrix(spikeTrain)==1    %most of my simulated data that normally are matrices
    
  for c=1:C

    x{c}=(spikeTrain(:,c));           
    
    x{c}=unique_no_sort(x{c});
    
    xsize(c)=size(x{1,c},2);

  end
  
  
end

ISImatrix=zeros((sum(xsize)-D*C+1),D*C);

                                      %constructing an empty matrix 

  for n=1:length(ISImatrix)
    
      
      for m=1:C
          
        ISImatrix(n,((m-1)*D+1):m*D)=[x{1,m}(indices(m):(indices(m)+D-1))];
        
        %fill the nth row of empty matrix with D consecutive ISIs from each cell
        %in the order of ascending cell number.

      end
      
      for i=1:C
        
        xsum(i)=sum(x{1,i}(1:(indices(i)+D-1)));
        
        %xsum is the sum of from the 1st ISI upto the (indices(i)+D)th ISIs in 
        %each cell.

      end
      
        [Y,I]=min(xsum,[],2);
        
        %find the minimum xsum of all cells
        %this will then be registered in 'indices'.

        indices(I)=indices(I)+1;
    
        %increment the indecies that corresponds to cell number with smallest 
        %xsum that had been filled into ISImatrix.
        %so this means the ISImatrix is updated everytime a cell fired.
        %ISImatrix row number n would represent the timeline of all spiking 
        %events

        clear xsum I

        if sum((xsize-indices)==(D-1))==1

            break 
        
        %stop the code when it has reached the end of any cell
            
        end
  end

ISImatrix=ISImatrix(sum(ISImatrix(:,:),2)~=0,:); %eliminate zeros for safety

  for k=1:length(ISImatrix)
  

      %disp(k/length(ISImatrix)*100)     
      %disp('percent')
   
    [index(k,:),euclidist(k,:)]=knnsearch(ISImatrix,ISImatrix(k,:),'k',2);
        
        %find nearest neighbour, and calculate nearest euclidean distance
  end

euclidist=euclidist(:,2);
euclidist=euclidist(euclidist~=0); %eliminate zeros

H=D/sum(xsize)*sum(log2(euclidist))+log2(S*sum(xsize)/D)+0.5772156649/log(2);

return
