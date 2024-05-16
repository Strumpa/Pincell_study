function indices = index_mat_to_vec(nvol,nsurf)

indpos=@(i,j) max(i,j).*(max(i,j)-1)./2+min(i,j) ;

indices = zeros(45,3) ;

for i=1:(nvol+nsurf) 
    for j=1:(nvol+nsurf)
        indices(indpos(i,j),1) = j ;
        indices(indpos(i,j),2) = i ;
        indices(indpos(i,j),3) = indpos(i,j) ;
    end
end