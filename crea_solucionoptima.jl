function solucion_optima(CR_L::Any)
        if minimum(CR_L)>=0
            #println("La solucion factible es optima")
            return 1
        else
            #println("La solución factible no es optima")
            return 0
        end
end
