function elegir_vsale(M::Any)
    nu_vbasicas=size(M)[1]-1
    cocientes=[]
    indexx=[]
    for i=1:nu_vbasicas
        if M[i,size(M)[2]]>0
            cocientes=vcat(cocientes,( M[i,size(M)[2]-1] )/ M[i,size(M)[2]])
            indexx=push!(indexx,i)
        end
    end
    if length(cocientes)==0
        return 0
    else
        if length(cocientes)==1
            return indexx[1]
        else
           return argmin(cocientes)
        end
    end
end
