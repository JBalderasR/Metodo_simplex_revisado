

function delete_vartificiales(vbs::Int64,va::Any)
  deleteat!(va, (j for j in eachindex(va) if va[j] == vbs))
return va
end
