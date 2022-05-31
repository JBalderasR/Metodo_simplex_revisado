
function createv_libres(v::Any,n::Int64)
    C_libre = convert(Vector,1:n)
    for  i in 1:length(v)
        deleteat!(C_libre, (j for j in eachindex(C_libre) if C_libre[j] == v[i]))
    end
    return C_libre
end

#deleteat!(C_libre, (j for j in eachindex(C_libre) if C_libre[j] == v[i]))
