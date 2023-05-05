function [H, distan]=binless_entropy_array(SpikeTrain, D)

%Input spikeTrain is an (N x 1) or (1 x N) array
%Computes differential entropy using binless estimation method
%D is the dimension, e.g. D=2 is joint differential entropy

S=D*pi^(D/2)/(gamma(D/2+1));%S changes as Dimension changes

SpikeTrain=unique_no_sort(SpikeTrain);%deletes repetitions which results in log(0),
                                                     %and re-odering to original order

ISIarray=diff(SpikeTrain);

N=length(ISIarray); %building D-dimensional matrix

ISIlist=zeros((N-D+1),D);

for n=1:(N-D+1)
    
    ISIlist(n,:)=ISIarray(n:(n+D-1)); 
    %list of D-dimensional points
end

for i=1:(N-D+1)
    
    [index(i,:),dista(i,:)]=knnsearch(ISIlist,ISIlist(i,:),'k',2);
    %find nearest neighbour, and calculate nearest euclidean distance
end

dista_=dista(:,2);
distan=dista_(dista_~=0); %eliminate zeros

H=D/length(SpikeTrain)*sum(log2(distan))+log2(S*(length(SpikeTrain))/D)+0.5772156649/log(2);



return