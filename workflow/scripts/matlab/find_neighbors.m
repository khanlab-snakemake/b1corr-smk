function [nbrs] = find_neighbors(faces,vertices)
% find neighboring relations between vertices in a surface

surf.nverts = size(vertices,1);
surf.nfaces = size(faces,1);
surf.faces = faces;
surf.coords = vertices;

num_nbrs=zeros(surf.nverts,1);
for i=1:surf.nfaces
    num_nbrs(surf.faces(i,:))=num_nbrs(surf.faces(i,:))+1;
end
max_num_nbrs=max(num_nbrs);

nbrs=zeros(surf.nverts,max_num_nbrs);
for i=1:surf.nfaces
    for j=1:3
        vcur = surf.faces(i,j);
        for k=1:3
            if (j ~= k)
                vnbr = surf.faces(i,k);
                if find(nbrs(vcur,:)==vnbr)
                else
                    n_nbr = min(find(nbrs(vcur,:) == 0));
                    nbrs(vcur,n_nbr) = vnbr;
                end
            end
        end
    end
end

return