using MAT
using PoissonRandom
using Distributions
using LinearAlgebra
using StatsBase


include("SummaryStats_Spatial.jl")



N=4e8
s=0.2
# s=0.015
hdom=0.5
D=200

t=250
#number of months per generation: choice of 1 reflects Gambiae, a choice of gamma=348 would be human generation time of 29 years.
gamma=1
samplingyear = 2013 #choose calendar year that sequences are sampled

#I want torigin to be such that the final number of generations at the final sampling time 
#in samplingyear the number of generations passed since torigin is t (where torigin is a calendar year)
# =>  t = (samplingyear - torigin)*12/gamma
yearorigin = samplingyear-gamma*t/12


println("N = ",round(N;sigdigits=4))
println("s = ",round(s;sigdigits = 3))
println("D = ",round(D;sigdigits=4))
println("year of origin = ",round(yearorigin;sigdigits=5))


A=1e7
Ar=2

# p=0.05
Gaussian=0

μ=2e-9

θ=2*N*μ #1.0


R=6378.137#km radius at the equator

#rectangular shape so these shape variables not needed
lon_shape=[]
lat_shape=[]
#sample locations will be determined in SummaryStats_Spatial.jl as a subset of sample vertices not continous latitude and longitude coordinstes like real data
latlon=[]
latlontime=[]

# randomspatialsampling=0
nspatialsites = 50
nsample = 2218

#Create an Nsamp=Nsampt which is nsample chromosomes distributed across nspatialsites
q = 1/nspatialsites*ones(nspatialsites); #vector of the probability of each spatial site being selected is equal
Nsamp = rand(Multinomial(nsample,q));
Nsampt = Nsamp;
uniquespatialindx = collect(Int,1:nspatialsites);



println("size(Nsamp) = ", size(Nsamp))

pltfigs=1

# @time SummaryStats_Spatial(N::Real, A::Real, θ::Real, s::Real, D::Real, p::Real, Ar::Real, Gaussian::Real, xwt_target, Nsamp, mint::Int, maxt::Int,lonsites,latsites,lon_shape,lat_shape,pltfigs)
@time η, Rg, H, NeiD, xmutsamp, Xt, tgen, ηmean, ηstd, n, Nsampt, Nsamp, nx, Ademe, xsamp, ysamp  = SummaryStats_Spatial(N::Real, A::Real, θ::Real, s::Real, D::Real, hdom::Real,
                                                                                                                            Ar::Real, Gaussian::Real, Nsamp, Nsampt, 
                                                                                                                            yearorigin::Real, gamma,latlon, latlontime,
                                                                                                                            uniquespatialindx,lon_shape,lat_shape,
                                                                                                                            pltfigs);
                                                                                    #    SummaryStats_Spatial(N::Real, Area::Real, twoNmu::Real, s::Real, D::Real, p::Real, Ar::Real, Gaussian::Real, Nsamp, Nsampt, torigin::Real, latlon, latlontime,uniquespatialindx,lon_shape,lat_shape,pltfigs)
println("")
println("[η, Rg, H, NeiD, xmutsamp, ηmean, ηstd] = ", round.([η, Rg, H, NeiD, xmutsamp, ηmean, ηstd]; sigdigits = 3))

results = Dict("tgen"=>tgen,
                 "eta"=>η,"Rg"=>Rg,"H"=>H,"NeiD"=>NeiD,"xmut"=>xmutsamp,
                 "eta_mean"=>ηmean, "eta_std"=>ηstd,
                 "N"=>N,"Area"=>A,"twoNmu"=>θ,"s"=>s,"D"=>D,"yearorigin"=>yearorigin,
                 "n"=>n,"Nsampt" => Nsampt, "Nsamp" => Nsamp, "nx"=>nx, "Ademe"=>Ademe,
                 "xsamp"=>xsamp,"ysamp"=>ysamp)

           

matwrite("Gambiae_Results.mat",results;compress=true)