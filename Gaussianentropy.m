%code for Gaussian distribution and discrete entropy

function Hgauss=Gaussianentropy(isilist,k)
%k is the denominator for determining no. of bins
%isilist=randn(nspikes,1);
nbins=length(isilist)/k;
[nisi,isibinned] = hist(isilist,nbins);

%figure(2); hist(isilist,8);
pisi = nisi./length(isilist);

%figure(1);
%plot(isibinned,pisi,'-o');

pisi=pisi(pisi~=0);
plogp=zeros(1,length(pisi));

   for o=1:size(pisi,2)
   plogp(o)=pisi(o).* log2(pisi(o)/((max(isilist)-min(isilist))/nbins));
   Hgauss=-sum(plogp);
   end 
end
   
   