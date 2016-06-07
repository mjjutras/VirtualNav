function col = storeColor(x, y, xStore, yStore)
  minIndex = 1;
  minDistance = mDistance(x, y, xStore(1), yStore(1));
  for i=2:length(xStore)
    d = mDistance(x,y,xStore(i),yStore(i));
    if d < minDistance
      minIndex = i;
      minDistance = d; 
    end
  end

  if minIndex == 1
    col = 'r';
  elseif minIndex == 2
    col = 'g';
  elseif minIndex == 3
    col = 'b';
  elseif minIndex == 4
    col = 'y';
  end
end

function d = mDistance(x, y, xStore, yStore)
    d = abs(x-xStore) + abs(y-yStore);
end
