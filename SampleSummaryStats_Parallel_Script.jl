using Distributions
using Random
using MAT


include("SampleSummaryStats_Parallel.jl")


#returns args[i] if it was supplied and non-empty, otherwise [] (so each positional
#arg can be left blank, e.g. julia script.jl "" "1.2e-7", to skip it while still
#supplying a later one)
function arg_or_empty(args,i)
    return (length(args)>=i && !isempty(args[i])) ? args[i] : []
end

ContinueSim_arg = arg_or_empty(ARGS,1)
ContinueSim = isempty(ContinueSim_arg) ? [] : parse(Int128,ContinueSim_arg)

#optional 2nd arg: mutation rate override, passed in from the PBS submission script.
#If not supplied (or empty), the default μ set below in the Species block is used.
mu_arg = arg_or_empty(ARGS,2)

#optional 3rd arg: number of replicates M, passed in from the PBS submission script.
#If not supplied (or empty), the default M set below is used (handy for quick local
#testing with a small M without that value reaching the cluster's submitted jobs).
M_arg = arg_or_empty(ARGS,3)

println("ContinueSim = ",ContinueSim)
println("mu_arg = ",mu_arg)
println("M_arg = ",M_arg)

PPC=0 #set to 1 to do posterior predictive check (PPC), where posterior samples are read in from file 
#and the maxmin values below ignored, otherwise set to 0 for normal behaviour to draw params between limits shown
#in SampleSummaryStats_Parallel 
#the file specified when PPC=1 is read in which contains posterior samples from Manolo's NN-ABC analysis
#a new directory is also created for storing the simulation results for each posteror sample set.
PPC_path = "Manolo_Analysis/Human_LCT_FieldSampling/NNABC_Rerun"
PPC_filename = "PosteriorSamples_mu=1.2e-7_h=1.mat" #ONly used if PPC=1

#Species switch: determines which priors/constants/data files are used below. 
#Set to either "AGambiae" (Anopheles gambiae, Vgsc locus) or "Human" (lactase persistence, LCT locus).
Species = "Human"
# Species = "AGambiae"

if Species=="AGambiae"

    #Priors used for Anopheles Gambiae
    Nminmax = [1e7,5e9] #if μ=2e-9 then Nμ = [0.02 to 200]
    sminmax = [2e-2,2.0]
    Dminmax = [2.5,500.0]
    Tminmax = [1940.0, 2000.0] #calendar years

    #Anopheles mutation rate: Rashid et al, 202x Sci Reps.
    μ=2e-9 # mu0 = 1e-9

    #Anopheles dominance at Vgsc locus
    hdom=1/2

    #read in the shape outline of the Gambiae range
    rangepath = "JoshsAfricaMapDemes/GambiaeRange.mat"

    #Anopheles sample-location data file(s), used further below
    samplelocationpath = "Phase2_data/Phase2_haplotype_SpaceTime_location.mat"

elseif Species=="Human"

    #Priors used for humans in Africa for LCT locus
    Nminmax = [5e4,5e7] #if μ=2e-9 then Nμ = [0.02 to 200]
    sminmax = [1e-3,1.0]
    Dminmax = [50.0,20000.0]
    Tminmax = [-13000.0, -3000.0] #calendar years

    #Human mutation rate per base-pair per generation
    μ=1.2e-7 #6e-8 #μ0 = 1.2e-8 (Scally & Durbin, 2012, Nat. Rev. Gen) & Long et al 2012, Nature.

    #Human dominance coefficient for lactase persistence
    hdom=1

    #Human generation time 29 years taken from Langergraber, K. E. et al. 
    #Generation times in wild chimpanzees and gorillas suggest earlier 
    #divergence times in great ape and human evolution. Proc. Natl. Acad. Sci. 109, 15716–15721 (2012).
    #This value of gamma = 29*12=348 is stored in the data file read in below

    #read in the shape outline of Africa (used as the human sampling range)
    rangepath = "JoshsAfricaMapDemes/africa_shape_outline.mat"

    #Human sample-location data file, used further below
    samplelocationpath = "LCT/n_array_latlon_Liebert2017_Afr_3SNPs.mat"

else

    error("Species must be \"AGambiae\" or \"Human\", got: ", Species)

end

#override the default mutation rate above if one was passed in via ARGS[2]
if !isempty(mu_arg)
    μ = parse(Float64,mu_arg)
end
println("μ = ",μ)

vars = matread(rangepath)
lon_shape = vars["lon"];
lat_shape = vars["lat"];

R=6378.137#km radius at the equator

randomspatialsampling=0

#Note that none of the code this script calls needs to read these files again, since pltfigs=0 
#in SampleSummaryStats_Parallel.jl, which calls SummaryStat_Spatial.jl with this value
if randomspatialsampling==0

    #read in the sample locations (samplelocationpath set above by the Species switch)
    vars = matread(samplelocationpath)
    # lonsites = vars["unique_lon"];
    # latsites = vars["unique_lat"];
    latlon = vars["unique_latlon"]
    latlontime = vars["unique_latlontime"];
    Nsamp = vars["Nsamp"];
    Nsampt = vars["Nsampt"];
    # Nsamp = vars["Nsamp"];
    gamma = vars["gamma"];
    uniquespatialindx = vars["uniquespatialindx"];
    uniquespatialindx = Int.(uniquespatialindx);


else #draw random space time sites but all at the same time! And distribute nsample chromosomes across them evenly

    #only used for random sampling so not for real data for Anopheles or Human
    #it is used to create a latlontime array which is the same type as for real data.
    if Species=="AGambiae"
        maxsamplingyear = 2013
    else
        maxsamplingyear = 2010 
    end

    
    nspatialsites = 500
    nsample = 22180

    vars = matread(samplelocationpath)
    gamma = vars["gamma"]; #number of months per generation



    #ok how do I draw random spatial sites within sub Saharan Africa

    latlon = zeros(nspatialsites,2)
    # latlontime = zeros(nspatialsites,3)
    minlon = minimum(lon_shape)
    minlat = minimum(lat_shape)
    maxlon = maximum(lon_shape)
    maxlat = maximum(lat_shape)


  
    polygonGambiae = Point.(vec(lon_shape),vec(lat_shape))
    # polygonGambiae = Point.(vec(x_shape),vec(y_shape))

    #now draw pairs of locations with 0 <= lon <= maxlon and 0 <= lat <= maxlat and then test if they are in the Gambaie range.
        
    let np=1

        while np<=nspatialsites

            Φ = rand(Uniform(minlon,maxlon))
            Θ = rand(Uniform(minlat,maxlat))


            P = Point.(Φ,Θ)
            
            if isinside(P,polygonGambiae; allowonedge=true)

                latlon[np,1] = Θ
                latlon[np,2] = Φ
                np += 1
            end

        end

    end

    #the 3rd column of latlontime is that sampling time from jan 2000 
    #last sampling time is appromximately dec 2012, so let's say 13 years after 2000 => 
    maxt = (maxsamplingyear-2000)*12 #in months
    maxt = maxt/gamma #in generations #so gamma is the number of months per generation 
    #and 1/gamma the number of generations per month
    #Note that this is the time from 2000, so not the actual maxt in generations 
    #and we add in the time of origin later using the following commented code is executed in SummaryStat_Spatial.jl:
        # samplingtime = samplingtime .+ (2000 - torigin)*12/gamma

        # maxt = maximum(samplingtime)


    println("maxt = ",maxt)
    times = maxt*ones(nspatialsites)
    latlontime  = hcat(latlon,times) #we'll be sampling all at the same time the last time point of the simulation

    q = 1/nspatialsites*ones(nspatialsites) #vector of the probability of each spatial site being selected is equal
    Nsamp = rand(Multinomial(nsample,q))

    Nsampt = Nsamp

    uniquespatialindx = collect(Int,1:nspatialsites)







end






Aminmax = [1e4,1e4] #this is ignored if being used with a shape file (as area set by shape file)
θminmax = [μ,μ] 

Gaussian=0
Ratio=2
# Nsamp = Int(5000)
# xmut_target = 0.78
# xwt_target=1-xmut_target
xwt_target = []


if Species=="AGambiae"
    prescript = string("AnophelesGambiaeRangeSubSaharanAfrica_h=",hdom)
else #Species=="Human"
    prescript = string("HumanLCT_AfricaRange_h=",hdom) # 0.8"
end

M=3#256

#override the default M above if one was passed in via ARGS[3]
if !isempty(M_arg)
    M = parse(Int,M_arg)
end
println("M = ",M)


if PPC ==0
    
    println("pwd() = ",pwd())

    if !isdir("SampleSummaryStats")
        mkdir("SampleSummaryStats")
    end
    cd("SampleSummaryStats")


    if !isempty(latlontime)

        if randomspatialsampling!=0
            dirname = string(prescript,"_RandomSampling_N=",Nminmax[1],"_",Nminmax[2],"_D=",Dminmax[1],"_",Dminmax[2],
                                "_s=",sminmax[1],"_",sminmax[2],"_Tmin=",Tminmax[1],"_Tmax=",Tminmax[2],"_nspatialsites=",nspatialsites,"_nsample=",nsample,"_mu=",μ)

        else

            dirname = string(prescript,"_N=",Nminmax[1],"_",Nminmax[2],"_D=",Dminmax[1],"_",Dminmax[2],
                                "_s=",sminmax[1],"_",sminmax[2],"_Tmin=",Tminmax[1],"_Tmax=",Tminmax[2],"_mu=",μ)
        end

    else


        if Aminmax[1]==Aminmax[2] 
            dirname = string("N=",Nminmax[1],"_",Nminmax[2],"_D=",Dminmax[1],"_",Dminmax[2],"_s=",sminmax[1],"_",sminmax[2],"Area=",Aminmax[1],"_Nsamp=",Nsamp)

        else
            dirname = string("N=",Nminmax[1],"_",Nminmax[2],"_D=",Dminmax[1],"_",Dminmax[2],"_s=",sminmax[1],"_",sminmax[2],"Area=",Aminmax[1],"_",Aminmax[2],"_Nsamp=",Nsamp)
        end

    end





    if !isdir(dirname)
        mkdir(dirname)
    end

    cd(dirname)



else
    #cd to posterior directory and create necessary directories
    cd(PPC_path)
    # cd("Manolo_Analysis/Human_LCT_FieldSampling/NNABC_Rerun")


end



N,s,D, yearorigin, A, θ, η, Rg, H, NeiD, xmut, ηmean, ηstd = SampleSummaryStats_Parallel(Nminmax::Vector{Float64}, Aminmax::Vector{Float64}, θminmax::Vector{Float64}, 
                                                                                        sminmax::Vector{Float64}, Dminmax::Vector{Float64}, Tminmax::Vector{Float64}, hdom,
                                                                                        Ratio::Real, Gaussian::Real, Nsamp, Nsampt, M::Integer, latlon, latlontime, 
                                                                                        uniquespatialindx, lon_shape, lat_shape, ContinueSim, PPC, gamma, Species,PPC_filename)


 



cd("../..")


