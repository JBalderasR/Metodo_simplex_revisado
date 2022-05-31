function pivotcolum(S::Any,i::Int64)
    col=size(S)[2]
    N=S[:,1:col-1]
    M=[N S[:,col]]
    M[i,:]=M[i,:]./M[i,col]
    for k in 1:(size(M)[1])
        #for j in 1:(size(M)[2])
            if k!=i
            #M[k,j]=-round(M[k,j]-M[k,col].*M[i,j],digits=14)
            M[k,:]=M[k,:]-M[k,col].*M[i,:]
            end
        #end
    end
    M[i,col]=0.0
    return M
end
