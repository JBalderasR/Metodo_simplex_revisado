
function actualizacion_vbasicas(vbasicas_anterior::Any,posicion::Int64,v_entra::Int64)
    long=length(vbasicas_anterior)
    vbasicas_siguiente=zeros(Int64,1,long)
    for i=1:long
        if i==posicion
        vbasicas_siguiente[i]=v_entra
        else
        vbasicas_siguiente[i]=vbasicas_anterior[i]
        end
    end
    return vbasicas_siguiente
end
