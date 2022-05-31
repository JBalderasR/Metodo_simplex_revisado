#Función que construye un vector que guarda la posicion de las va. basicas
# A partir de la matriz "A" 

function createv_basicas(A::Any)
        n=size(A)[1] #encuentra el número de filas de una matriz
        ncol=size(A)[2]
        v_basicas=zeros(Int64,1,n)
        j=1
        while j<=ncol
            if (sum(A[:,j])==1.0) & (minimum(A[:,j])==0.0) & (maximum(A[:,j])==1.0)
                i=1
                posicion=argmax(A[:,j])
                acum=0.0
                    while  (A[i,j]==0.0) & (i<=n)
                        acum=acum+0.0
                        i=i+1
                    end
                    if  acum==0.0
                    v_basicas[posicion]=j
                    end
            end
            j=j+1
        end
        return v_basicas
end
