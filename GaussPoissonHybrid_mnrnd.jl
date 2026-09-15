using Distributions #MvNormal
using StatsBase #isposdef
using LinearAlgebra #eigen, Hermitian
using PoissonRandom
include("Allele_Var_Mtx.jl")

function GaussPoissonHybrid_mnrnd(N,x)

 #Hybrid Gaussian-Poisson approx to multinomial random number generator
 #N should be an integer
 #x should be a vector of frequencies/probabilities with sum = 1
 if length(findall(x.<0)) > 0
   println("Error: frequency vector must not have negative entries")
 end

 if sum(x) !=1
    x = x/sum(x)
 end

      if N<1e3
    #The approximation requires that the threshold used nth<<N — nth will
    #typically be about 10 so if N<1e3 seems reasonable to just use exact
    #multinomial random numbers — for N<1e3 multinomial runs fast
         n = rand(Multinomial(Int(N),x))
      else
         #Threshold to use Poisson appprox
          nthr = Float64(10)
          indnot0 = x.!=0
          n = zeros(size(x))

          #find all entries which have expected number in next generation <=gamma
          indsmall = N*x .<= nthr

          if sum(indsmall.==1) == 0 #i.e. there are no "small" entries
                  nsum = 0
          else
             #Instead treating these "multinomially" as they are small compared to
              #N: draw from Poisson distribution instead
             n[indsmall .& indnot0] = pois_rand.(N.*x[indsmall .& indnot0])
             #Sum the number of "small" alleles
              nsum = sum(n[indsmall .& indnot0])
              #index ~indsmall are those alleles with expected number >nthr & by def !=0
           end

           NN = N-nsum

           #Convert frequencies of "large" alleles to frequencies within subset of
            #large alleles and also take only on simplex

           indlarge  = findall(N*x .> nthr)

           if size(indlarge,1) > 1  #more than one large allele & there has to be at least 1!
              #xx should be on the simplex of the subset of alleles which are "large"
              xx = x[indlarge[1:end-1]]/sum(x[indlarge])
              #So xx does not sum to one

              #Draw random number from multivariate normal distribution with
              #variance matrix of allele frequencies given by B

              B = allele_var_mtx(xx)/NN

              nn = NN*rand(MvNormal(xx,B))


              while sum(nn)>NN || isempty(findall(nn.<0)) == 0 #So rarely negative entries can be produced
                      nn = NN*rand(MvNormal(xx,B))
              end

         #The Gaussian approx will produce non-integer nn for the large
         #subset=> round them

                 nn = round.(nn)
                 n[indlarge[1:end-1]] = nn

                 #The total sum of the large subset should equal NN
                   nlast = NN - sum(nn)
                   n[indlarge[end]] = nlast

         elseif size(indlarge,1) == 1
            n[indlarge[1]] = NN # Changed by KC, 07/07/2022

         end
      end #if N <1e3
      return n
end #function
