%% search neighbourhood of specifilied center node according to the input search
% criteria.
% Note: only work in 2d cases.
% By Guofeng Cao
function IndexClose = searchnodes2(inGrid,gridSpecs,index,srchSpecs,nCloseMax)
% index is the result of geoeas2ixigiz 

nDim=2;
vDim = gridSpecs(:,1);
nRow = vDim(1);
nCol = vDim(2);
dPixelSize = gridSpecs(:,3)';
vRadMax = (srchSpecs(2:3));
vRad = max(vRadMax)./dPixelSize;
nRadius = max(vRad);
nCount = 1;

xyz=indgeoeas2ixiyiz(index,vDim,nDim,1);
nLen=size(inGrid,1);

IndexClose = zeros(nCloseMax,1);
%ValClose = zeros(nCloseMax,1);
nRadius = min(max(nRow-xyz(1),nCol-xyz(2)),nRadius);

for i = 1:nRadius
   for j = -i:i
          vTemp = [xyz(1)+j,xyz(2)-i ];
          if vTemp(1)>nRow||vTemp(1)<=0||vTemp(2)>nCol||vTemp(2)<=0
              continue;
          end
          index = (vTemp(2)-1)*nCol+vTemp(1);
         if ~isnan(inGrid(index))
            IndexClose(nCount) = index;
%             ValClose(nCount) = inGrid(index);
            nCount = nCount+1;
            if nCount > nCloseMax
                IndexClose = IndexClose(1:nCount-1);
%                 ValClose = ValClose(1:nCount-1);
                return;
            end
         end
   end
   for j = -i+1:i
          vTemp = [xyz(1)+i,xyz(2)+j];
          if vTemp(1)>nRow||vTemp(1)<=0||vTemp(2)>nCol||vTemp(2)<=0
              continue;
          end
          index = (vTemp(2)-1)*nCol+vTemp(1);

         if ~isnan(inGrid(index))
            IndexClose(nCount) = index;
%             ValClose(nCount) = inGrid(index);
            nCount = nCount+1;
            if nCount > nCloseMax
               IndexClose = IndexClose(1:nCount-1);
%                ValClose = ValClose(1:nCount-1);
               return;
            end
         end
   end
   for j = i-1:-1:-i
          vTemp = [xyz(1)+j,xyz(2)+i];
          if vTemp(1)>nRow||vTemp(1)<=0||vTemp(2)>nCol||vTemp(2)<=0
              continue;
          end
          index = (vTemp(2)-1)*nCol+vTemp(1);
         if ~isnan(inGrid(index))
            IndexClose(nCount) = index;
%             ValClose(nCount) = inGrid(index);
            nCount = nCount+1;
            if nCount > nCloseMax
               IndexClose = IndexClose(1:nCount-1);
%                ValClose = ValClose(1:nCount-1);
               return;
            end
         end
   end
   for j = i-1:-1:-i+1
         vTemp = [xyz(1)-i,xyz(2)+j];
         if vTemp(1)>nRow||vTemp(1)<=0||vTemp(2)>nCol||vTemp(2)<=0
              continue;
          end
         index = (vTemp(2)-1)*nCol+vTemp(1);
         if ~isnan(inGrid(index))
            IndexClose(nCount) = index;
%             ValClose(nCount) = inGrid(index);
            nCount = nCount+1;
            if nCount > nCloseMax
               IndexClose = IndexClose(1:nCount-1);
%                ValClose = ValClose(1:nCount-1);
               return;
            end
         end
   end
end
IndexClose = IndexClose(1:nCount-1);
% ValClose = ValClose(1:nCount-1);
return;
