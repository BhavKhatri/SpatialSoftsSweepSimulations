#B = Allele Variance Matrix(x)
function allele_var_mtx(x)  #x has to be a vector
  l=length(x)
  B=zeros(l,l)
  for j in 1:l, k in 1:l
     if j==k
        B[j,k]=x[j]*(1-x[j])
      elseif k>j #only fill the upper triangle part
        B[j,k]=-x[j]*x[k]
      end
   end
  #using LinearAlgebra
  B=triu(B,1)'+B #diagonal-symmetric!
  return B
 end
