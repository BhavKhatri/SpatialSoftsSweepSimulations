using MAT

include("SummaryStats_Spatial.jl")

Species = "Human"



N= 2*5.8767e+05 #Number of chromosomes 

#From the data the range of origins is about 1000km, so in about Δt = 500 generations 
#and since origins are not overlapping significantly this means distance covered is
# Δr = vΔt => v = Δr/Δt = 1000/500 = 2km/gen => v = 2sqrt(Ds) = 2 => Ds≈1 m^2/gen^2
#So selection coefficient estimates of 14010 are between 0.01 and 0.1 so chose 0.05 and so D = 20 to start exploration
#this is too small a diffusion constant, since if D is small we get classic Fisher waves 
#but D=2000 and s=0.03 gives thinly spread fisher waves as we seem to see in data. But then v = 2*sqrt(D*s) ≈15km/gen and so Δt = Δr/v = 1000/15 ≈ 65 generations ! 
#But this is only true of the first origin and so is the size of this origin. Also this doesn't take account of the time of establishment and so in 65 generations 
#and an establishment time of 1/2s = 16 generations, so 

s=0.082794
hdom = 1
D=3820


μ=6e-8  #1.2e-7  #human bp mutation rate is μ0=1.2e-8 so for 5 SNPs we have μ = 5μ0 = 6e-8

# t=250
vars = matread("LCT/n_array_latlon_Liebert2017_Afr_3SNPs.mat")
gamma = vars["gamma"]; #number of months per generation
println("γ = ", gamma," months per generation")
#I want yearorigin to be such that the final number of generations at the final sampling time 
#in 2013 the number of generations passed since yearorigin is t (where yearorigin is a calendar year)
# =>  t = (2013 - yearorigin)*12/gamma
maxsamplingyear =2010
# yearorigin = maxsamplingyear-gamma*t/12
yearorigin = -7000#-8050.4 
# samplingtime = samplingtime .+ (2000 - yearorigin)*12/gamma

println("N = ",round(N;sigdigits=4))
println("s = ",round(s;sigdigits = 3))
println("D = ",round(D;sigdigits=4))
println("yearorigin = ",round(yearorigin;sigdigits=5))



#these next two are irrelevant, since they are determined by the shape file
A=NaN
Ar=NaN

Gaussian=0


θ=2*N*μ #1.0

R=6378.137#km radius at the equator


#read in the shape outline of Gambiae range
vars = matread("JoshsAfricaMapDemes/africa_shape_outline.mat")
lon_shape = vars["lon"]
lat_shape = vars["lat"]

randomspatialsampling=0
nspatialsites = 50
nsample = 2218


if randomspatialsampling==0

    #read in the sample locations
    vars = matread("LCT/n_array_latlon_Liebert2017_Afr_3SNPs.mat")
    latlon = vars["unique_latlon"]
    latlontime = vars["unique_latlontime"];
    Nsamp = vars["Nsamp"];
    Nsampt = vars["Nsampt"];
    # Nsamp = vars["Nsamp"];
    gamma = vars["gamma"];
    uniquespatialindx = vars["uniquespatialindx"];
    uniquespatialindx = Int.(uniquespatialindx);


else #draw random space time sites but all at the same time! And distribute nsample chromosomes across them evenly

    vars = matread("LCT/n_array_latlon_Liebert2017_Afr_3SNPs.mat")
    gamma = vars["gamma"]; #number of months per generation


    latlon = zeros(nspatialsites,2)

    minlon = minimum(lon_shape);
    minlat = minimum(lat_shape);


    #So now minlon =minlat=0
    maxlon = maximum(lon_shape);
    maxlat = maximum(lat_shape);
    
    polygonGambiae = Point.(vec(lon_shape),vec(lat_shape));

    #now draw pairs of locations with 0 <= lon <= maxlon and 0 <= lat <= maxlat and then test if they are in the Gambaie range.
        
    let np=1

        while np<=nspatialsites

            Lon = rand(Uniform(minlon,maxlon))
            Lat = rand(Uniform(minlat,maxlat))

            P = Point.(Lon,Lat)

            if isinside(P,polygonGambiae; allowonedge=true)
                latlon[np,1] = Lat
                latlon[np,2] = Lon
                np += 1
            end

        end

    end

    #last sampling time is appromximately dec 2012, so let's say 13 years after 2000 => 
    maxt = (maxsamplingyear-2000)*12 #in months
    # maxt = maxt/gamma #in generations

    println("maxt = ",maxt, " months")
    times = maxt*ones(nspatialsites);
    latlontime  = hcat(latlon,times); #we'll be sampling all at the same time the last time point of the simulation

    q = 1/nspatialsites*ones(nspatialsites); #vector of the probability of each spatial site being selected is equal
    Nsamp = rand(Multinomial(nsample,q));

    Nsampt = Nsamp;

    uniquespatialindx = collect(Int,1:nspatialsites);







end


println("size(Nsamp) = ", size(Nsamp))

pltfigs=1
@time η, Rg, H, NeiD, xmutsamp, Xt, tgen, ηmean, ηstd, n, Nsampt, Nsamp, nx, Ademe, xsamp, ysamp = SummaryStats_Spatial(N::Real, A::Real, θ::Real, s::Real, D::Real, hdom::Real, Ar::Real, Gaussian::Real, Nsamp, Nsampt, yearorigin::Real, gamma, latlon, latlontime,uniquespatialindx,lon_shape,lat_shape,pltfigs,Species);

println("")
println("[η, Rg, H, NeiD, xmutsamp, ηmean, ηstd] = ", round.([η, Rg, H, NeiD, xmutsamp, ηmean, ηstd]; sigdigits = 3))


results = Dict("tgen"=>tgen,
                 "eta"=>η,"Rg"=>Rg,"H"=>H,"NeiD"=>NeiD,"xmut"=>xmutsamp,
                 "eta_mean"=>ηmean, "eta_std"=>ηstd,
                 "N"=>N,"Area"=>A,"twoNmu"=>θ,"s"=>s,"D"=>D,"yearorigin"=>yearorigin,
                 "n"=>n,"Nsampt" => Nsampt, "Nsamp" => Nsamp, "nx"=>nx, "Ademe"=>Ademe,
                 "xsamp"=>xsamp,"ysamp"=>ysamp)
           

matwrite("Human_Results.mat",results;compress=true)