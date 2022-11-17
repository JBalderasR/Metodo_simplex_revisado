include("./crea_pfactible.jl")
include("./crea_vbasicas.jl")
include("./crea_vlibres.jl")
include("./crea_solucionoptima.jl")
include("./crea_eleccionvsale.jl")
include("./crea_actualizacionvbasicas.jl")
include("./crea_pivotematrix.jl")
include("./delete_vartificial.jl")

using LinearAlgebra,DataFrames,Printf
function MSR(b,c_a,A,v_artificiales,costos)
   @printf("\t Método Simplex Revisado\n")
   no_col= size(A)[2] #(m+h+a)
   lista_vbasicas=[]
   v_basicas = createv_basicas(A)
   A0 = [A b; c_a' 0]
   n= size(v_basicas)[2]
   A0_0=A[1:n,1:length(costos)]
   c_b=c_a[v_basicas]'
   menos_z = (-1)c_b'*b
   R0 = [Matrix(1.0I,n,n) zeros(n,1) b zeros(n,1) ; zeros(1,n) 1.0 menos_z 0.0]
   Inversa_B=R0[1:n,1:n]
   R0[n+1,1:n]=-c_b'*Inversa_B
   solucion=0
   solucion_acotada=1
   solucion2=0
   #inicio iteraciones
   R = []
   l_vbasicas=[]
   l_vlibres=[]
   x=[]
   costosreducidos=[]
   lista_ventra=[]
   R=[]
   z=[]
   iteraciones=0
   parar=0
   fase1=1.0
   v_libres0=createv_libres(v_basicas,no_col)
   lista_vbasicas=push!(lista_vbasicas,v_basicas)
   l_vlibres=push!(l_vlibres,v_libres0)
   k=0
   LK=[]

   while (solucion==0) & (solucion_acotada==1)&(parar==0)&(k<17)
      LK=push!(LK,k)
      x0=create_pfactible(no_col,v_basicas,R0[1:n,size(R0)[2]-1])
      z=push!(z,round.(-R0[n+1,n+2],digits=14))
      x=push!(x,round.(x0,digits=14))
      L0=A[:,v_libres0]
      C_reducidosL=R0[n+1,1:n]'*L0
      costosreducidos=push!(costosreducidos,round.(C_reducidosL,digits=14))
      solucion=solucion_optima(C_reducidosL)#cambiar nombre por óptima
      posicion_e=argmin(C_reducidosL')
      v_libre_e=v_libres0[posicion_e]#Revisión de los signos de los costos red
      lista_ventra=push!(lista_ventra,v_libre_e)
      c_reducidoe=C_reducidosL[posicion_e]
      Ae=R0[1:n,1:n]*A[1:n,v_libre_e]
      Ae_1=Ae;push!(Ae_1,c_reducidoe)
      R0[1:n+1,n+3]=Ae_1
      R0
      R=push!(R,round.(R0,digits=14))

       if (length(v_artificiales)==0)&((round.(-R0[n+1,n+2],digits=12))==0)
          A0 = [A0_0 b; costos' 0]
          no_col= size(A0)[2]-1
          Inversa_B=R0[1:n,1:n]
          c_b=costos[v_basicas]
          z1= (-c_b*Inversa_B)*b
          R0=[R0[1:n,1:n+3]; -c_b*Inversa_B R0[n+1,n+1] z1 R0[n+1,n+2]]
          x0=create_pfactible(no_col,v_basicas,R0[1:n,size(R0)[2]-1])
          z=push!(z, round.(-R0[n+1,n+2],digits=15))
          lista_vbasicas = push!(lista_vbasicas,v_basicas)
          x=push!(x,round.(x0,digits=14))
          v_libres0=createv_libres(v_basicas,no_col)
          l_vlibres=push!(l_vlibres,v_libres0)
          L0=A0[:,v_libres0]
          C_reducidosL=R0[n+1,1:n+1]'*L0
          costosreducidos=push!(costosreducidos,round.(C_reducidosL,digits=14))
          solucion=solucion_optima(C_reducidosL)#cambiar nombre por óptima
          posicion_e=argmin(C_reducidosL')
          v_libre_e=v_libres0[posicion_e]#Revisión de los signos de los costos red
          lista_ventra=push!(lista_ventra,v_libre_e)
          c_reducidoe=C_reducidosL[posicion_e]
          Ae=R0[1:n,1:n]*A[1:n,v_libre_e]
          Ae_1=Ae;push!(Ae_1,c_reducidoe)
          R0[1:n+1,n+3]=Ae_1
          fase1=1.0
          R=push!(R,round.(R0,digits=14))
          k1=k
          k=k1
          LK=push!(LK,k)
      end

      posicion_vsale0=elegir_vsale(R0,lista_ventra)
      if posicion_vsale0==0
         solucion_acotada=0 #si la variable solucion_acotada=0 y cuando es igual a 1 es no acotada
      else
          solucion_acotada=1
      end
     if   minimum(costosreducidos[k+1])>0
         @printf("\t Solución óptimo no acotado\n")
         parar=1
     end
      vbasica_sale0=v_basicas[posicion_vsale0]
      v_artificiales=delete_vartificiales(vbasica_sale0,v_artificiales)
      v_basicas = actualizacion_vbasicas(v_basicas,posicion_vsale0,v_libre_e)
      v_libres0=createv_libres(v_basicas,no_col)
      R1=Array{Float64,2}(undef,n+1,n+3)
      lista_vbasicas=push!(lista_vbasicas,v_basicas)
      l_vlibres=push!(l_vlibres,v_libres0)
      R1 = pivotcolum(R0,posicion_vsale0)
      R0=R1
      L0=A0[:,v_libres0]
      C_reducidosL=R0[n+1,1:n+1]'*L0
      solucion=solucion_optima(C_reducidosL)
      if (z[k+1]==0.0) & (length(v_artificiales)==0.0)
         solucion2=0
      else
         solucion2=1
         solucion=0
      end
      if (solucion==1)&(fase1==1.0) & (length(v_artificiales)==0.0)
         R=push!(R,round.(R0,digits=14))
         k+=1
         LK=push!(LK,k)
         z=push!(z,round.(-R0[n+1,n+2],digits=14))
         x0=create_pfactible(no_col,v_basicas,R0[1:n,size(R0)[2]-1])
         x=push!(x,x0)
         costosreducidos=push!(costosreducidos,round.(C_reducidosL,digits=14))
         parar=1
      end
      k+=1
   end

      ET1=[]
      for i in LK
      etiqueta_columna=[ "$i" " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " b" "a.e"]
      ET1=push!(ET1,etiqueta_columna)
      end
      ET=[]
      Basicas=[]
      for k in 1:length(lista_vbasicas)
         etiqueta_fila=["x_"*"$(lista_vbasicas[k][1])" ;"x_"*"$(lista_vbasicas[k][2])";"x_"*"$(lista_vbasicas[k][3])" ;"x_"*"$(lista_vbasicas[k][4])";"x_"*"$(lista_vbasicas[k][5])" ;"x_"*"$(lista_vbasicas[k][6])";
         "x_"*"$(lista_vbasicas[k][7])" ;"x_"*"$(lista_vbasicas[k][8])";"x_"*"$(lista_vbasicas[k][9])" ;"x_"*"$(lista_vbasicas[k][10])";"x_"*"$(lista_vbasicas[k][11])" ;"x_"*"$(lista_vbasicas[k][12])"; "x_"*"$(lista_vbasicas[k][13])" ;"x_"*"$(lista_vbasicas[k][14])";"-z"]
         ET=push!(ET,etiqueta_fila)
         VB=["x_"*"$(lista_vbasicas[k][1])"  "x_"*"$(lista_vbasicas[k][2])"]
         Basicas=push!(Basicas,VB)
      end


      VL=Array{String}[]
      E=String[]
      for j in 1:length(l_vlibres)
         for i in l_vlibres[j]
            E=push!(E,"x_"*"$(i)")
            if length(E)==length(l_vlibres[j])
               VL=push!(VL,E)
               E=[]
            end
         end
      end

      Basicas=Array{String}[]
      VB=String[]
      for j in 1:length(lista_vbasicas)
         for i in lista_vbasicas[j]
            VB=push!(VB,"x_"*"$(i)")
            if length(VB)==length(lista_vbasicas[j])
               Basicas=push!(Basicas,VB)
               VB=[]
            end
         end
      end

      Tablas=[]
      for i in 1:length(R)
      df_R=[ET1[i]; ET[i] R[i]]
       Tablas=push!(Tablas,df_R)
      end

      for i in 1:length(R)-1
      print( "\n" ,"\n",DataFrame(Tablas[i],[" ","c_1","c_2","c_3","c_4","c_5","c_6","c_7","c_8","c_9","c_10","c_11","c_12","c_13","c_14","c_15","c_16","c_17"] ) )
      print( "\n" ,"\n","Solución factible z = ",z[i] )
      print( "\n" ,"Punto factible x = ",x[i] )
      print( "\n" ,"Variables libres = ",VL[i] )
      print( "\n" ,"Variables básicas = ", Basicas[i] )
      end

      print( "\n" ,"\n",DataFrame(Tablas[k+1],[" ","c_1","c_2","c_3","c_4","c_5","c_6","c_7","c_8","c_9","c_10","c_11","c_12","c_13","c_14","c_15","c_16","c_17"]  ))
      print( "\n" ,"\n","Solución óptima z* = ",z[k+1] )
      print( "\n" ,"Punto óptimo x* = ",x[k+1] )
      print( "\n" ,"Variables libres = ",VL[k+1] )
      print( "\n" ,"Variables básicas = ", Basicas[k+1] )
end


include("./crea_pfactible.jl")
include("./crea_vbasicas.jl")
include("./crea_vlibres.jl")
include("./crea_solucionoptima.jl")
include("./crea_eleccionvsale.jl")
include("./crea_actualizacionvbasicas.jl")
include("./crea_pivotematrix.jl")
include("./delete_vartificial.jl")

using LinearAlgebra,DataFrames,Printf

function MSR(b,c_a,A,v_artificiales,costos)
   @printf("\t Método Simplex Revisado\n")
   
   no_col= size(A)[2]
   lista_vbasicas=[]
   v_basicas = createv_basicas(A)
   A0 = [A b; c_a' 0]
   n= size(v_basicas)[2]
   A0_0=A[1:n,1:n+2]
   c_b=c_a[v_basicas]'
   menos_z = (-1)c_b'*b
   R0 = [Matrix(1.0I,n,n) zeros(n,1) b zeros(n,1) ; zeros(1,n) 1.0 menos_z 0.0]
   Inversa_B=R0[1:n,1:n]
   R0[n+1,1:n]=-c_b'*Inversa_B
   solucion=0
   solucion_acotada=1
   solucion2=0
 
   R = []
   l_vbasicas=[]
   l_vlibres=[]
   x=[]
   costosreducidos=[]
   R=[]
   z=[]
   iteraciones=0
   parar=0
   v_libres0=createv_libres(v_basicas,no_col)
   lista_vbasicas=push!(lista_vbasicas,v_basicas)
   l_vlibres=push!(l_vlibres,v_libres0)
   k=0
   LK=[]
   while (solucion==0) & (solucion_acotada==1)&(parar==0)
      LK=push!(LK,k)
      x0=create_pfactible(no_col,v_basicas,R0[1:n,size(R0)[2]-1])
      z=push!(z,round.(-R0[n+1,n+2],digits=14))
      x=push!(x,round.(x0,digits=14))
      L0=A[:,v_libres0]
      C_reducidosL=R0[n+1,1:n]'*L0
      costosreducidos=push!(costosreducidos,round.(C_reducidosL,digits=14))
      solucion=solucion_optima(C_reducidosL)
      posicion_e=argmin(C_reducidosL')
      v_libre_e=v_libres0[posicion_e]
      c_reducidoe=C_reducidosL[posicion_e]
      Ae=R0[1:n,1:n]*A[1:n,v_libre_e]
      Ae_1=Ae;push!(Ae_1,c_reducidoe)
      R0[1:n+1,n+3]=Ae_1
      R0
      R=push!(R,round.(R0,digits=14))

       if length(v_artificiales)==0
          A0 = [A0_0 b; costos' 0]
          no_col= size(A0)[2]-1
          Inversa_B=R0[1:n,1:n]
          c_b=costos[v_basicas]
          z1= (-c_b*Inversa_B)*b
          R0=[R0[1:n,1:n+3]; -c_b*Inversa_B R0[n+1,n+1] z1 R0[n+1,n+2]]
          x0=create_pfactible(no_col,v_basicas,R0[1:n,size(R0)[2]-1])
          z=push!(z, round.(-R0[n+1,n+2],digits=14))
          lista_vbasicas = push!(lista_vbasicas,v_basicas)
          x=push!(x,round.(x0,digits=14))
          v_libres0=createv_libres(v_basicas,no_col)
          l_vlibres=push!(l_vlibres,v_libres0)
          L0=A0[:,v_libres0]
          C_reducidosL=R0[n+1,1:n+1]'*L0
          costosreducidos=push!(costosreducidos,round.(C_reducidosL,digits=14))
          solucion=solucion_optima(C_reducidosL)
          posicion_e=argmin(C_reducidosL')
          v_libre_e=v_libres0[posicion_e]
          c_reducidoe=C_reducidosL[posicion_e]
          Ae=R0[1:n,1:n]*A[1:n,v_libre_e]
          Ae_1=Ae;push!(Ae_1,c_reducidoe)
          R0[1:n+1,n+3]=Ae_1
          R=push!(R,round.(R0,digits=14))
          k1=k
          k=k1
          LK=push!(LK,k)
      end

      posicion_vsale0=elegir_vsale(R0)
      vbasica_sale0=v_basicas[posicion_vsale0]
      if posicion_vsale0==0
         solucion_acotada=0 
      else
          solucion_acotada=1
      end

      v_artificiales=delete_vartificiales(vbasica_sale0,v_artificiales)
      v_basicas = actualizacion_vbasicas(v_basicas,posicion_vsale0,v_libre_e)
      v_libres0=createv_libres(v_basicas,no_col)
      R1=Array{Float64,2}(undef,n+1,n+3)
      lista_vbasicas=push!(lista_vbasicas,v_basicas)
      l_vlibres=push!(l_vlibres,v_libres0)
      R1 = pivotcolum(R0,posicion_vsale0)
      R0=R1
      L0=A0[:,v_libres0]
      C_reducidosL=R0[n+1,1:n+1]'*L0
      solucion=solucion_optima(C_reducidosL)
      if (z[k+1]==0.0) & (length(v_artificiales)==0.0)
         solucion2=0
      else
         solucion2=1
         solucion=0
      end
      if (solucion==1)&(solucion2==0)
         R=push!(R,round.(R0,digits=14))
         k+=1
         LK=push!(LK,k)
         z=push!(z,round.(-R0[n+1,n+2],digits=14))
         x0=create_pfactible(no_col,v_basicas,R0[1:n,size(R0)[2]-1])
         x=push!(x,x0)
         costosreducidos=push!(costosreducidos,round.(C_reducidosL,digits=14))
         parar=1
      end
      k+=1
   end


   ET1=[]
   for i in LK
   etiqueta_columna=[ " $i" " " " " " " "  b" "a.e"]
   ET1=push!(ET1,etiqueta_columna)
   end
   ET=[]
   Basicas=[]
   for k in 1:length(lista_vbasicas)
      etiqueta_fila=["x_"*"$(lista_vbasicas[k][1])" ;"x_"*"$(lista_vbasicas[k][2])"; "-z"]
      ET=push!(ET,etiqueta_fila)
      VB=["x_"*"$(lista_vbasicas[k][1])"  "x_"*"$(lista_vbasicas[k][2])"]
      Basicas=push!(Basicas,VB)
   end


   VL=Array{String}[]
   E=String[]
   for j in 1:length(l_vlibres)
      for i in l_vlibres[j]
         E=push!(E,"x_"*"$(i)")
         if length(E)==length(l_vlibres[j])
            VL=push!(VL,E)
            E=[]
         end
      end
   end

   Basicas=Array{String}[]
   VB=String[]
   for j in 1:length(lista_vbasicas)
      for i in lista_vbasicas[j]
         VB=push!(VB,"x_"*"$(i)")
         if length(VB)==length(lista_vbasicas[j])
            Basicas=push!(Basicas,VB)
            VB=[]
         end
      end
   end

   Tablas=[]
   for i in 1:length(R)
   df_R=[ET1[i]; ET[i] R[i]]
    Tablas=push!(Tablas,df_R)
   end

   for i in 1:length(R)-1
   print( "\n" ,"\n",DataFrame(Tablas[i],[" ","c_1","c_2","c_3","c_4","c_5"] ) )
   print( "\n" ,"\n","Solución factible z = ",z[i] )
   print( "\n" ,"Punto factible x = ",x[i] )
   print( "\n" ,"Variables libres = ",VL[i] )
   print( "\n" ,"Variables básicas = ", Basicas[i] )
   end

   print( "\n" ,"\n",DataFrame(Tablas[k+1],[" ","c_1","c_2","c_3","c_4","c_5"]  ))
   print( "\n" ,"\n","Solución óptima z* = ",z[k+1] )
   print( "\n" ,"Punto óptimo x* = ",x[k+1] )
   print( "\n" ,"Variables libres = ",VL[k+1] )
   print( "\n" ,"Variables básicas = ", Basicas[k+1] )
end
