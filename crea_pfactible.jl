#Función que construye y evalua el punto factible de un probl. de PL
#A partir de la posicion de las variables basicas(v_basicas)
#y a partir de los terminos independientes(b)
#ENTRADA: Numero de variables del problema (reales,holgura,artificiales)
#,vector de posiciones de variables basicas
#vector de actuales terminos independientes.
#SALIDA: Vector que representa el punto factible.
function create_pfactible(n_variables::Int64,v_b::Any,b::Any)
    x=zeros(Float64,1,n_variables)
    for i=1:length(v_b)
         for k=1:n_variables
                if k==v_b[i]
                x[k]=b[i]
                end
         end
    end
    return x
end
