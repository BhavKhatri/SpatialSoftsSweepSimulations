using Distributions
using Random
using MAT
using CSV
using DataFrames
using Dates



include("SummaryStats_Spatial.jl")




function SampleSummaryStats_Parallel(Nminmax::Vector{Float64}, Aminmax::Vector{Float64}, θminmax::Vector{Float64}, sminmax::Vector{Float64}, Dminmax::Vector{Float64}, Tminmax::Vector{Float64}, hdom,
                                                                                Ratio::Real, Gaussian::Real, Nsamp, Nsampt, M::Integer, latlon, latlontime, uniquespatialindx, lon_shape, lat_shape, ContinueSim, PPC, gamma, Species::String="AGambiae", PPC_filename)
 
    #M is the number of replicates
    #Each replicate has parameters drawn from a log-uniform distribution

    #if Aminmax have equal elements then fix area to A = Aminmax[1]=Aminmix[2]
    #if θminmax are equal then assume value to be a fixed mutation rate μ = θminmax[1]=θminmax[2] and then θ[k] = 2*N[k]*μ

    if PPC==0

        N = zeros(M)
        s = zeros(M)
        D = zeros(M)        
        yearorigin = zeros(M)
    else
        #readin values from file
        # var = matread("PosteriorSamples_mu=1.2e-7_NewPrior_h=1.mat")
        var = matread(PPC_filename)
        N = var["N"]
        s = var["s"]
        D = var["D"]
        yearorigin = var["t0"] #this is the number of years before 2013 the last sample date in Phase2_data
        
        M = length(N)

    end

    A = zeros(M)
    θ = zeros(M)
    NDeme = zeros(M)
    ND = zeros(M)

    k=1

    Area=0

    if !isempty(lon_shape)
        R=6378.137#km radius at the equator

        x,y = LongLat2km(lon_shape,lat_shape,R)

        polygonAfrica = Point.(vec(x),vec(y)) #create polygon object from x and y positions defining outline/shape of Africa

        #Calculate area from polyarea in Luxor and shape file of gambiae range
        Area = polyarea(polygonAfrica)
        println("Area from range shape file A = ",round(Area;sigdigits=2))


        
    end

    η = zeros(M)
    Rg = zeros(M)
    H = zeros(M)
    NeiD = zeros(M)
    xmut = zeros(M)
    ηmean = zeros(M)
    ηstd = zeros(M)


    exclude = [false for i=1:M]

    if isempty(ContinueSim)

        if PPC==0
            while k<=M 
                println("k = ", k)

                #drawnumbers
                N[k], s[k], D[k], yearorigin[k], θ[k], ND[k], NDeme[k], μ, A[k] = DrawRandomParameters(Nminmax,sminmax,Dminmax,Tminmax,Aminmax,Area,θminmax)

                #check

                if NDeme[k]>=10 && ND[k] <5e5 && ND[k]>1
                    k=k+1
                end



            
            end

            


            jj=1
            dirname = "Individual_samples"
            # dirname = string("Individual_samples_M=",M)
            if !isdir(dirname)
                mkdir(dirname)
            else
                dirnamenew = string(dirname,"_",jj)
                while isdir(dirnamenew)
                    jj+=1
                    dirnamenew = string(dirname,"_",jj)
                end

                dirname=dirnamenew

                mkdir(dirname)
                

            end
            
            cd(dirname)


        
            priorparams = Dict("N"=>N,"s"=>s, "D"=>D,"yearorigin"=>yearorigin,"theta"=>θ,
                        "Area"=>A,"ND"=>ND,"NDeme"=>NDeme,"mu"=>μ)



            matwrite("PriorParams.mat",priorparams;compress=true)
            CSV.write("PriorParams.csv", DataFrame(priorparams)) #DataFrames are essentially like Matlab tables and this converts non-arrays (scalars) 
            #to arrays the same length as other arrays, but with same entry, so all the arrays have to be the same length. 
            #The csv files columns are these arrays with labels given by those defined in the Dict


        else

            for k=1:M
                

                if any(isnan.([N[k],s[k],D[k],yearorigin[k]])) #somtimes Manolo's posterior samples using NN regression give NaNs (particularly for yearorigin)
                    exclude[k] = true
                    θ[k] = NaN
                    ND[k] = NaN
                    NDeme[k] = NaN
                    A[k] = NaN
                else
                    _, _, _, _, θ[k], ND[k], NDeme[k], μ, A[k] = DrawRandomParameters(N[k],s[k],D[k],yearorigin[k],Aminmax,Area,θminmax)
                end

            end
            #calculate remaining paramters from the posterior in file for posterior predictive check
            dirname = "PosteriorPredictiveCheck_NN_mu=1.2e-7"
                
            if !isdir(dirname)
                mkdir(dirname)
            end

            cd(dirname)
        end






    else #read in params file 

        println("ContinueSim = ",ContinueSim)

        if ContinueSim==0
            strname = "Individual_samples"
        else
            strname = string("Individual_samples_",ContinueSim)
        end

        cd(strname)

        var = matread("PriorParams.mat")
        N = var["N"]
        s = var["s"]
        D = var["D"]
        yearorigin = var["yearorigin"]
        θ = var["theta"]
        A = var["Area"]
        ND = var["ND"]
        NDeme = var["NDeme"]
        μ = var["mu"]

        files = filter(x -> occursin("Sum", x), readdir()) #this should return all saved mat files, starting with "SummaryStat.." in array files

        # m = zeros(length(files))

        overwrite=0 #read in each replicate and overwrite — normally set this equal to zero, to read-in only

        for k in eachindex(files)

            ind1 = findfirst('=', files[k])
            ind2 = findfirst('.', files[k])


            mk = parse(Int,files[k][ind1+1:ind2-1]) #equivalent of str2num in matlab, but you need to specify type of number
            exclude[mk] = true
            # m[k] = mk
            vark = matread(files[k])

            η[mk] = vark["eta"]
            Rg[mk] = vark["Rg"]
            H[mk] = vark["H"]
            NeiD[mk] = vark["NeiD"]
            xmut[mk] = vark["xmut"]
            ηmean[mk] = vark["eta_mean"]
            ηstd[mk] = vark["eta_std"]


        end

    end

    println("M = ",M)

    starttime = [0.0 for i=1:M]
    rerun_starttime = [0.0 for i=1:M]
    endtime = [0.0 for i=1:M]

    initialtime = now()

    
    Threads.@threads for i=1:M
    # for i=1:M

        Δt = now() - initialtime
        starttime[i] = Δt.value/1000 #time in seconds

        filename =  string("SummaryStats_sample_i=",i,".mat")
    
        if !exclude[i] & !isfile(filename)

            NN = N[i]
            DD = D[i]
            ss = s[i]
            yyearorigin = yearorigin[i]
            Area = A[i]
            Θ = θ[i]
            
            NND = ND[i]
            NNDeme  = NDeme[i]

            pltfigs=0

            println("")
            println("New simulation")
            ηη, RRg, HH, NNeiD, xxmut, Xt, tgen, ηηmean, ηηstd, n, Nsampti, Nsampi,nx, Ademe = SummaryStats_Spatial(NN, Area, Θ, ss, DD, hdom, Ratio, Gaussian, Nsamp, Nsampt, yyearorigin, gamma, latlon,latlontime,uniquespatialindx,lon_shape,lat_shape,pltfigs,Species)

            println("η = ",ηη)
            summarystats = [ηη,ηηmean,ηηstd,RRg,HH,NNeiD,xxmut]

            # while isnan(ηη) || (ηη==0)
            while any(isnan.(summarystats)) || (ηη==0)

                println("Summary stats = ",summarystats)
                
                Δt = now() - initialtime
                rerun_starttime[i] = Δt.value/1000 #time in seconds

                if PPC==0 #redraw parameters, otherwise if PPC=1 then keep the same parameters — need to be careful if summary stats produced keep being NaNs then need to change this

                    uu = 0

                    while uu==0
                        NN, ss, DD, yyearorigin, Θ, NND, NNDeme, μ = DrawRandomParameters(Nminmax,sminmax,Dminmax,Tminmax,Aminmax,Area,θminmax)

                        if NNDeme>=10 && NND <5e5 && NND>1
                        uu=1
                        end
                    end
            
                end

                # NN = rand(LogUniform(Nminmax[1],Nminmax[2]))
                # ss = rand(LogUniform(sminmax[1],sminmax[2]))
                # DD = rand(LogUniform(Dminmax[1],Dminmax[2]))

                # if θminmax[1]==θminmax[2]
                #     μ = θminmax[1]
                #     Θ = 2*NN*μ
                # else
                #     Θ = rand(LogUniform(θminmax[1],θminmax[2]))
                # end
                

                
                println("")
                println("New simulation: rerun")
                
                ηη, RRg, HH, NNeiD, xxmut, Xt, tgen, ηηmean, ηηstd, n, Nsampti, Nsampi, nx, Ademe = SummaryStats_Spatial(NN, Area, Θ, ss, DD, hdom, Ratio, Gaussian, Nsamp,Nsampt, yyearorigin, gamma, latlon,latlontime,uniquespatialindx,lon_shape,lat_shape,pltfigs,Species)
  
                summarystats = [ηη,ηηmean,ηηstd,RRg,HH,NNeiD,xxmut]

                if PPC==0
                    # summarystats = [ηη,ηηmean,ηηstd,RRg,HH,NNeiD,xxmut]
                    if !(any(isnan.(summarystats)) || (ηη==0))
                    # if !isnan(ηη) || (ηη==0)

                        N[i] = NN
                        s[i] = ss
                        D[i] = DD
                        yearorigin[i] = yyearorigin
                        θ[i] = Θ
                        ND[i] = NND
                        NDeme[i] = NNDeme
                        println("Success after rerun")

                        priorparams = Dict("N"=>N,"s"=>s, "D"=>D,"yearorigin"=>yearorigin,"theta"=>θ,"Ademe"=>Ademe,
                        "Area"=>A,"ND"=>ND,"NDeme"=>NDeme,"mu"=>μ)
        
                        wait_until_file_closed("PriorParams.mat")
                        # while isopen("PriorParams.mat")
                        #     sleep(10)
                        # end

                        matwrite("PriorParams.mat",priorparams;compress=true)

                        wait_until_file_closed("PriorParams.csv")
                        # while isopen("PriorParams.csv")
                        #     sleep(10)
                        # end

                        CSV.write("PriorParams.csv", DataFrame(priorparams))
                    end
            
                end

            end


            println("i = ",i)
            println("time in generations from introdction of insecticides to last sample time = ",tgen)
            println("η = ",ηη)
            η[i]  = ηη
            Rg[i] = RRg
            H[i] = HH
            NeiD[i] = NNeiD
            xmut[i] = xxmut
            ηmean[i] = ηηmean
            ηstd[i] = ηηstd


            Δt = now() - initialtime
            endtime[i] = Δt.value/1000 #time in seconds
            
 

            results = Dict("tgen"=>tgen, "NDeme"=>NDeme[i],"ND"=>ND[i], 
                 "eta"=>η[i],"Rg"=>Rg[i],"H"=>H[i],"NeiD"=>NeiD[i],"xmut"=>xmut[i],
                 "eta_mean"=>ηmean[i], "eta_std"=>ηstd[i],
                 "N"=>NN,"Area"=>Area,"twoNmu"=>Θ,"s"=>ss,"D"=>DD,"yearorigin"=>yyearorigin,
                 "n"=>n,"Nsampt" => Nsampti, "Nsamp" => Nsampti, "nx"=>nx, "Ademe"=>Ademe,
                 "starttime"=>starttime[i],"rerun_starttime"=>rerun_starttime[i],"endtime"=>endtime[i])

           

        
    
            wait_until_file_closed(filename)


            matwrite(filename,results;compress=true)
 
            



        else

            println("Excluding replicate ",i," as this has already been completed")

        end
            


          
            
    end

    timings = Dict("starttime"=>starttime,"rerun_starttime"=>rerun_starttime,"endtime"=>endtime)
    matwrite("timings.mat",timings;compress=true)

    cd("../")

    return N,s,D, yearorigin, A, θ, η, Rg, H, NeiD, xmut, ηmean, ηstd

end


function DrawRandomParameters(Nminmax,sminmax,Dminmax,Tminmax,Aminmax,Area,θminmax)

    if length(Nminmax)==2

        N = rand(LogUniform(Nminmax[1],Nminmax[2]))
        s = rand(LogUniform(sminmax[1],sminmax[2]))
        D = rand(LogUniform(Dminmax[1],Dminmax[2]))
        # yearorigin = Int(round(rand(Uniform(Tminmax[1],Tminmax[2]))))
        yearorigin = rand(Uniform(Tminmax[1],Tminmax[2])) #yearorigin — in years — is a continous variable: 
                                                    #in the simulation code the condition to end is that t>max(yearorigin)+1, 
                                                    #where t is an integer index for generation time tgen = t-1

    else#only one element to each of below so just calculate deme size and number

        N = Nminmax
        s = sminmax
        D = Dminmax
        yearorigin = Tminmax


    end
    p=0.05
    Ademe = 4*D/p #square demes, where we assume p*w^2 = 4D is the MSD in one generation and w is the deme width


    if θminmax[1]==θminmax[2]
        μ = θminmax[1]
        θ = 2*N*μ
    else
        θ = rand(LogUniform(θminmax[1],θminmax[2]))
    end


    if Area==0

        if Aminmax[1]==Aminmax[2]
            A = Aminmax[1]
        else
            A = rand(LogUniform(Aminmax[1],Aminmax[2]))
        end

    
        NND = A/Ademe; 
        println("1st estimate of number of demes: Area/Ademe =",NND)
        ny = round(Int,sqrt(NND/Ratio))
        nx = Int(Ratio*ny)
    
        println("nx=",nx)
        println("ny=",ny)
        ND = nx*ny
        println("Total number of demes after correcting for specified Aspect Ratio = ",ND)
        AA = Ademe*ND
        println("Modified total area (km^2) =",AA)
        println("Requested area (km^2) = ", A[k])

        if ND > 0
    
            NDeme = round(Int, N / ND) #Population of mosquitos in a Deme 
            println("Population size per deme NDeme = ", NDeme)
            println("Actual total population size ND*NDeme = ", ND*NDeme)

        else
            NDeme = 0
        
        end
    
    else    

       
        #Approx estimate of number of demes and population size per deme
        ND = Area/Ademe
        NDeme = N/ND

        println("Initial estimate of number of demes ND = ",round(ND))
        println("Initial estimate of population size per deme  = ",round(NDeme))

        A=Area

    end



    return N, s, D, yearorigin, θ, ND, NDeme, μ, A



    



end


function wait_until_file_closed(path; check_interval=10)
    while true
        try
            open(path, "a") do _ end  # try opening for append
            break  # success means file is not locked
        catch e
            sleep(check_interval)  # wait before retrying
        end
    end
end