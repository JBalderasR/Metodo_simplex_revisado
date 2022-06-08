# Metodo_simplex_revisado
Se implementa el método simplex revisado (método de la gran M) de Programación Lineal en el lenguaje Julia. Usa las librerías DataFrames y LinearAlgebra.

Entrada:

  A: Matriz de coeficientes.

  b: vector de términos independientes.

  v_artificiales: vector de etiquetas de las variables artificiales.

  c_a: vector de coeficientes de las variables artificiales.

  costos: vector de costos de la función objetivo.

Salida:

Por cada iteración se tienen la tabla simplex, solución factible, punto factible, vector de variables libres, vector de variables básicas.

__________________________________________________________________________________________________________________________________________________________________

Ejemplo: Solución óptima acotada.

Entrada:

julia> b=[12;14]
2-element Vector{Int64}:
 12
  4

julia> c_a=[0; 0; 0; 0; 1]

5-element Vector{Int64}:
 0
 0
 0
 0
 1
 julia> A=[1 1 -1 0 1;5 -2 0 1 0]
 
2×5 Matrix{Int64}:
 1   1  -1  0  1
 5  -2   0  1  0

julia> v_artificiales=[5]
1-element Vector{Int64}:
 5

julia> costos=[4;3;0;0]
4-element Vector{Int64}:
 4
 3
 0
 0

 Se carga el siguiente comando para ejecutar el método MSR: MSR(b, c_a, A,v_artificiales, costos)
 
 Salida:
 
 julia> MSR(b,c_a,A,v_artificiales,costos)
 
         Método Simplex Revisado


4×6 DataFrame

 Row │      c_1   c_2  c_3  c_4    c_5  
 
     │ Any  Any   Any  Any  Any    Any  
     
─────┼──────────────────────────────────

   1 │  0                     b    a.e
   
   2 │ x_5  1.0   0.0  0.0  12.0   1.0
   
   3 │ x_4  0.0   1.0  0.0  4.0    5.0
   
   4 │ -z   -1.0  0.0  1.0  -12.0  -1.0

Solución factible z = 12.0

Punto factible x = [0.0 0.0 0.0 4.0 12.0]

Variables libres = ["x_1", "x_2", "x_3"]

Variables básicas = ["x_5", "x_4"]

4×6 DataFrame

 Row │      c_1   c_2   c_3  c_4    c_5  
 
     │ Any  Any   Any   Any  Any    Any  
     
─────┼───────────────────────────────────

   1 │  1                      b    a.e
   
   2 │ x_5  1.0   -0.2  0.0  11.2   1.4
   
   3 │ x_1  0.0   0.2   0.0  0.8    -0.4
   
   4 │ -z   -1.0  0.2   1.0  -11.2  -1.4

Solución factible z = 11.2

Punto factible x = [0.8 0.0 0.0 0.0 11.2]

Variables libres = ["x_2", "x_3", "x_4"]

Variables básicas = ["x_5", "x_1"]

4×6 DataFrame

 Row │      c_1       c_2        c_3  c_4  c_5  
 
     │ Any  Any       Any        Any  Any  Any   
     
─────┼───────────────────────────────────────────────

   1 │  2                               b  a.e
   
   2 │ x_2  0.714286  -0.142857  0.0  8.0  -0.714286
   
   3 │ x_1  0.285714  0.142857   0.0  4.0  -0.285714
   
   4 │ -z   0.0       0.0        1.0  0.0  0.0

Solución factible z = -0.0

Punto factible x = [4.0 8.0 0.0 0.0 0.0]

Variables libres = ["x_3", "x_4", "x_5"]

Variables básicas = ["x_2", "x_1"]

4×6 DataFrame

 Row │      c_1       c_2        c_3  c_4    c_5   
 
     │ Any  Any       Any        Any  Any    Any  
     
─────┼─────────────────────────────────────────────────

   1 │  2                               b    a.e
   
   2 │ x_2  0.714286  -0.142857  0.0  8.0    -0.142857
   
   3 │ x_1  0.285714  0.142857   0.0  4.0    0.142857
   
   4 │ -z   -3.28571  -0.142857  1.0  -40.0  -0.142857

Solución factible z = 40.0

Punto factible x = [4.0 8.0 0.0 0.0]

Variables libres = ["x_3", "x_4"]

Variables básicas = ["x_2", "x_1"]

4×6 DataFrame

 Row │      c_1   c_2  c_3  c_4    c_5 
 
     │ Any  Any   Any  Any  Any    Any 
     
─────┼─────────────────────────────────

   1 │  3                     b    a.e
   
   2 │ x_2  1.0   0.0  0.0  12.0   0.0
   
   3 │ x_4  2.0   1.0  0.0  28.0   0.0
   
   4 │ -z   -3.0  0.0  1.0  -36.0  0.0

Solución óptima z* = 36.0

Punto óptimo x* = [0.0 12.0 0.0 28.0]

Variables libres = ["x_1", "x_3"]

Variables básicas = ["x_2", "x_4"]
