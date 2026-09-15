
using PoissonRandom
using Distributions
using LinearAlgebra
# using SparseArrays
using StatsBase
# using DelimitedFiles
using Luxor
using Plots
using Contour
using Interpolations
#using LaTeXStrings
using LazySets
# using GeometryOps

# Force Julia to use StatsBase's sample function by default
import StatsBase: sample

import Interpolations: Line

include("NumberofOrigins_Spatial_ComplexHaplotypes_SqDemes_mnmigration.jl") #this includes functions to calculate x, y position from latitude
include("GaussPoissonHybrid_mnrnd.jl")

@doc """
This runs the spatial soft sweep code outputting the final state in Xt, which is an array of an array where Xt[i] is a complex vector which lists all the haplotypes in the ith deme. 
Each entry is a complex number h+im*f whose real part is the haplotype identifier h and the imaginary part it's frequency f. The rest of the entries are NaNs.

This code first converts this to an array of an array X, where X[i] holds a fixed length vector of length nhaps+1, 
where nhaps is the total number of origins and +1 for wt, and each entry is the frequency of that origin (i.e. X[i][k] is the frequency of the (k+1)th origin/wt (k=1 is wt))

This is then sampled either 1) by randomly selecting demes, or 2) if specified by the list of sample locations in lon and lat 

Input:
- N: population size of chromosomes (if simulating diploid organism with additive fitness effects divide by 2),
- Area: Total area (this is ignored if shape outline specified in lon and lat aren't empty),
- theta: population-scaled mutation rate (2Nμ),
- s: selection coefficient,
- D: diffusion constant [km^2 per generation]
- Ar: aspect ratio (width / height) (this is ignored if shape outline specified in lon and lat aren't empty)
- Gaussian: 1 (using Gaussian approximation) or 0 (multinomial distribution)
- Nsamp: number of chromosomes/haploid individuals sampled at each spatial location 
- Nsampt: number of chromosomes/haploid individuals sampled at each space-time location
- yearorigin: the calendar year of origin of selective pressure
- gamma number of months per generation of organism
- latlon a nsites x 2 array containing the latitude and longitude (in this order) of the sample locations 
- latlontime a nsites x 3 array containing the latitude and longitude and sampling time —in months— (in this order) of the space-time sample locations
- uniquespatialindx: index of length Nsampt into Nsamp array indicating the unique spatial locations
- lon_shape and lat_shape contain arrays defining the species range in longitude and latitude
- pltfigs: plot figure of haplotype frequency maps (option below to uncomment section to plot pie charts of frequencies of haplotypes, but currently locations are not automated and need to be configured in file read-in)

Output:
- η: number of independent origins across all demes in sample size Nsamp
- Rg: average radius of gyration over all origins
- H: average heterozygosity/Simpsons diversity per deme
- NeiD: average (log) distance between all pairs of demes
- xmut: frequency of mutant
- Xt: frequency array of an array with complex haplotypes/frequency information in each deme
- tgen: final generation time
- ηmean: mean number of origins over all sample sites
- ηstd: standard deviation of origins over all sample sites 
- n: nsites x nhaps array of numbers of of haplotype — n[i,j] = number of jth haplotype sample in ith spatial location
- Nsamp: number of chromosomes/haploid individuals sampled at each spatial location 
- Nsampt: number of chromosomes/haploid individuals sampled at each space-time location 
- nx: Number of demes in x direction
- Ademe: Area of each deme
- NDeme: Number of haploid individuals/chromosomes per deme
- xx: x-position of sample sites in km
- yy: y-position of sample sites in km

"""
function SummaryStats_Spatial(N::Real, Area::Real, theta::Real, s::Real, D::Real, hdom::Real, Ar::Real, Gaussian::Real, Nsamp, Nsampt, yearorigin::Real, gamma, latlon, latlontime,uniquespatialindx,lon_shape,lat_shape,pltfigs,Species::String="AGambiae")

    p=0.05 # percentage that migrate out of deme per generation — keep fixed at p=0.05 unless you have good reason!
    Ademe = 4*D/p #square demes, where we assume p*w^2 = 4D is the MSD in one generation and w is the deme width
    
    R=6378.137#km radius at the equator — defined here (not inside the if-block below) so it's always available, including when lon_shape is empty
    println("Radius of Earth at the equator R = ",R," km")

    if !isempty(lon_shape)

        
        #calculate x,y positions of shape outline in lon and lat

        ND, nx, ny, x_shape, y_shape, minx, miny = CalculateDemes(Ademe, Area, Ar, lon_shape,lat_shape,R)


        wdeme = sqrt(Ademe)
        #Create a mask array for Africa/sub-Saharan Africa from shape outline input in x and y in kms
        polygonGambiae = Point.(vec(x_shape),vec(y_shape)) #create polygon object from x and y positions defining outline/shape of Africa
        xp = 0:wdeme:(nx-1)*wdeme #vector of deme positions/points in x-direction
        yp = 0:wdeme:(ny-1)*wdeme #vector of deme positions/points in y-direction
        demepoints = vec(Point.(xp',yp)) #calling Point in this way, where one argument is row vector and the other column 
        #— like an outer product — gives a *vector* of deme positions: size(demepoints) = nx*ny
        DemesInAfrica = [isinside(p,polygonGambiae; allowonedge=true) for p in demepoints] #size(DemesInAfrica) = nx*ny — are these column order and which way?
        #Both demepoints and DemesInAfrica are column-ordered with increasing y-values, then increasing x values:
        #e.g. demepoints = (0,0) , (0,wdeme), (0,2*wdeme), ... (0, (nx-1)*wdeme), (wdeme,0), (wdeme, wdeme), (wdeme, 2*wdeme),...
        # => my arbirary integer grid points for demes will ba mapped onto space in this way


        #Update DemesInAfrica for any demes that have zero nearest neighbours (e.g. on peninsula)
        _, _, DemesInAfrica = CalcNearestNeighbours(nx,ny,DemesInAfrica)
        

        if pltfigs==1
            #Add total outline of Africa — not just Gambiae range
            vars = matread("JoshsAfricaMapDemes/africa_shape_outline.mat")
            xAfr = vars["x"]
            yAfr = vars["y"]

            xAfr = xAfr .- minx
            yAfr = yAfr .- miny
            polygonAfrica = Point.(vec(xAfr),vec(yAfr))

        end

        lon = latlontime[:,2]
        lat = latlontime[:,1]
        samplingtime = latlontime[:,3] #this is the sampling time in months from Jan 2000
        
      
        #if yearorigin is a year date e.g. yearorigin = 1940 then we need to add this to all times, and convert to generations
        #where gamma is the number of months per generation 
        samplingtime = (samplingtime .+ (2000 - yearorigin)*12)/gamma

        maxt = maximum(samplingtime)
        println("maxt = ",maxt)

        samplingtime = Int.(round.(samplingtime)) #round to give an integer generation

        println("")
        println("max simulation time (corresponding to last date of sample) = ",maxt, " generations")
        println("measured in generations from the year inputted that insecticides were introduced (yearorigin = ",yearorigin,")")
        println("")
       



        #Calculate positions of space-time sample sites in kms
        x,y = LongLat2km(lon,lat,R)

        #Simulations performed with offset in x and y below from the 
        x = x .- minx
        y = y .- miny

        samplepoints = Point.(vec(x),vec(y))
        SamplesInAfrica = [isinside(p,polygonGambiae; allowonedge=true) for p in samplepoints]

        println("Number of space-time samples sites before Gambiae range filter = ", length(x))
        
        x = x[SamplesInAfrica]
        y = y[SamplesInAfrica]
        
        samplingtime = samplingtime[SamplesInAfrica]
        
        
        Nsampt = Nsampt[SamplesInAfrica] #so I filter Nsampt based on whether the locations are within species range but not Nsamp!
                
        uniquespatialindx = uniquespatialindx[SamplesInAfrica]
        #So this is a list of indices into the array which has unique spatial locations for which Nsamp is the number of chromosomes sampled at each.
        


        #Calculate the positions of each spatial sample site
        lon = latlon[:,2]
        lat = latlon[:,1]

        #Calculate positions of spatial sample sites in kms - which will be used to calculate the radius of gyration of origins 
        xx,yy = LongLat2km(lon,lat,R)
        xx = xx .- minx
        yy = yy .- miny

        samplepoints = Point.(vec(xx),vec(yy))
        SamplesInAfrica = [isinside(p,polygonGambiae; allowonedge=true) for p in samplepoints]

        println("length(SamplesInAfrica)) = ",length(SamplesInAfrica))

        println("length(Nsamp)) = ",length(Nsamp))
    

        #Return the closest deme number for each space-time sample site located at x and y
        spacetimesample_deme_number = Calculate_Closest_Demes_For_SampleLocations(x,y,demepoints,ny,DemesInAfrica,wdeme)

        uniquesamplingtime = sort(unique(samplingtime))
        numbuniquetimepoints = length(uniquesamplingtime)

        println("Number of (space-time) samples sites after Gambiae range filter = ", length(x))

   


    else #simulation on rectangular space (no shape file / real sample locations available): 
         #Nsamp/Nsampt (as passed in) already give the number of chromosomes to sample at each of nspatialsites sites — 
         #site *locations* are randomly chosen here (as distinct demes), and every site is sampled at the same time 
         #(the final generation of the simulation, from yearorigin and gamma)

        ND, nx, ny = CalculateDemes(Ademe, Area, Ar,lon_shape,lat_shape,R)

        wdeme = sqrt(Ademe)
        DemesInAfrica = fill(true,ND) #no shape file to exclude any demes, so every deme is available/in-range

        #deme centre positions, using the same column-major (y varies fastest, then x) convention as the shape-file branch above
        xp = 0:wdeme:(nx-1)*wdeme
        yp = 0:wdeme:(ny-1)*wdeme
        demepoints = vec(Point.(xp',yp))

        #Rectangular boundary polygon (padded by half a deme-width so it encloses deme *areas*, not just centres) — 
        #used by PlotHaplotypes below in place of the shape-file outline, since there isn't one in rectangular mode.
        #PlotHaplotypes snaps out-of-bounds contour points to the *nearest polygon vertex*, so the boundary needs many 
        #points along each edge (not just the 4 corners) for that snapping to spread evenly along the edges.
        xlo, xhi = -wdeme/2, (nx-1)*wdeme+wdeme/2
        ylo, yhi = -wdeme/2, (ny-1)*wdeme+wdeme/2
        nbpts = max(4*(nx+ny), 500) #number of boundary points per edge, scaled to grid size
        rect_x = vcat(range(xlo,xhi,length=nbpts), fill(xhi,nbpts), range(xhi,xlo,length=nbpts), fill(xlo,nbpts))
        rect_y = vcat(fill(ylo,nbpts), range(ylo,yhi,length=nbpts), fill(yhi,nbpts), range(yhi,ylo,length=nbpts))
        polygonGambiae = Point.(rect_x, rect_y)

        #Nsamp (as passed in) gives the number of chromosomes to sample at each of nspatialsites sites.
        #Site *locations* are chosen here: nspatialsites distinct demes, picked uniformly at random without replacement.
        nspatialsites = length(Nsamp)
        if nspatialsites > ND
            error("nspatialsites ($nspatialsites, from length(Nsamp)) exceeds the number of available demes ($ND) — reduce nspatialsites or increase Area / decrease Ademe.")
        end
        siteidx = sample(1:ND, nspatialsites, replace=false) #distinct randomly-chosen deme locations

        xx = [demepoints[k].x for k in siteidx]
        yy = [demepoints[k].y for k in siteidx]

        #All sites sampled at the same time: the final generation of the simulation, computed from yearorigin and gamma
        #(using year 2000 as the reference point, consistent with the shape-file branch above, where sampling times are measured in months from Jan 2000)
        maxt = ((2013 - yearorigin)*12)/gamma
        println("maxt = ",maxt)
        samplingtime = fill(Int(round(maxt)), nspatialsites)

        Nsampt = Nsamp #each site is sampled only once (at the single shared time above), so space-time sample counts equal the per-site counts
        uniquespatialindx = collect(1:nspatialsites) #each space-time sample maps 1:1 to its own unique spatial site (no repeated time-points per site)

        #Sample sites coincide exactly with deme centres, so the "closest deme" for each site is just that deme itself
        spacetimesample_deme_number = siteidx

        uniquesamplingtime = sort(unique(samplingtime))
        numbuniquetimepoints = length(uniquesamplingtime)

        maxnumsamples=1

        println("Number of randomly-placed spatial sample sites (no shape file) = ", nspatialsites)

        if pltfigs==1
            #No wider geographic context available in rectangular mode, so reuse the same rectangle for the background layer
            polygonAfrica = polygonGambiae
        end

    end


    # xmut_target = 0.072669 #[]# if empty this ensures that simulations end at t=maxt, otherwise they end when this mutant frequency is reached.
    xmut_target = []

    println("")
    println("Starting main spatial soft sweeps simulation")
    @time Xt, ncomplexhap, nx, Ademe, NDeme, t, truefreq = NumberOfOrigins_Spatial(N, Area, theta, s, D, hdom, Ar, Gaussian, maxt, xmut_target,lon_shape,lat_shape,Nsamp,samplingtime,Nsampt,uniquespatialindx,spacetimesample_deme_number,uniquesamplingtime);
    println("Time to run main spatial soft sweeps simulation")

    #Remove from ncomplexhap empty entries which should correspond to sites outside Gambiae range, as they haven't been accessed by entries in uniquespatialindx
    ind_empty = findall(.!isempty.(ncomplexhap))
    ncomplexhap = ncomplexhap[ind_empty]
    Nsamp = Nsamp[ind_empty]
    nspatialsites = length(Nsamp)
    println("Nsamp = ",Nsamp)
    xx = xx[ind_empty]
    yy = yy[ind_empty]
    println("length(xx) = ",length(xx))
    println("length(yy) = ",length(yy))

    samplepoints = Point.(vec(xx),vec(yy)) #samplepoints will only be used for plotting  

    ND = length(Xt)
    ny = Int(ND/nx);

    println("nx = ",nx)
    println("ny = ", ny)

    #Create CartesianIndices over deme array of size nx x ny 
    CartInd = CartesianIndices((1:ny,1:nx))


    #Need to transform format of X to X[i][j,k], where size(X[i])[1] is the number of different haplotypes 
    #and j and k refer to haplotype number and time point respectively
    #for case of only tracking last the last time point this is X[i][j], where for each deme X[i] is an array of length
    #maxh, containing the frequency of jth allele in ith deme, so that in general X[i][j] can be zero as certain 
    #alleles j don't currently exist in that deme
    #maxh below is defined to include the WT

    #find total number of haplotypes that have existed
    maxh = 0.0;

    for i in eachindex(Xt)

        if DemesInAfrica[i]
            indNotNan  = findall(x->!isnan(x),Xt[i])
            
            if isempty(indNotNan)
                println("indNotNan = ",indNotNaN)
                println("Xt[i] = ",Xt[i])
                println("i = ",i)
                prinln("ND = ",ND)
                println("NDeme = ",NDeme)
                println("D = ",D)
                println("N = ",N)
            end


            maxh = max(maxh,maximum(real(Xt[i][indNotNan])))

    
        end
    end



    maxh = Int(maxh+1); #+1 to include the wild type: hmarker=0.0 is the wild type
    # maxt = size(Xt[1])[2]

    println("maxh = ",maxh)
    # println("maxt = ",maxt)


    X = [fill(0.0,maxh) for iy=1:ny, ix=1:nx] #current generation
    
    Xhap = [fill(0.0,ny,nx) for ii=1:maxh]
    
    
    for i in CartInd
        
        # Get linearindex
        k = LinearIndices(CartInd)[i]


        if DemesInAfrica[k]

            hmarker = real(Xt[i][:])
            freq = imag(Xt[i][:])


            indNotNan  = findall(x->(!isnan(x)),hmarker) #this includes WT

            h = hmarker[indNotNan]
            f = freq[indNotNan]


            for j in eachindex(indNotNan)

                X[i][Int(h[j]+1)] = f[j] #hmarker=0.0 is the WT and so for indexing strating with 1, we need hmarker+1 for correct index
               
                Xhap[Int(h[j]+1)][i] = f[j]

            end


        end
    end

    maxhn = 0.0;

    for i in eachindex(ncomplexhap)

        
        indNotNan  = findall(x->!isnan(x),ncomplexhap[i])
        
        maxhn = max(maxhn,maximum(real(ncomplexhap[i][indNotNan])))

    
        
    end

    maxhn = Int(maxhn+1); #+1 to include the wild type: hmarker=0.0 is the wild type
    println("maxhn = ",maxhn)
    println("maxh = ",maxh)
    #Construct n: a numspacetimesamplessites x nhaps arrayfrom ncomplexhap instead of nt, where length of ncomplexhap is the number of sample sites in data
    n = zeros(nspatialsites,maxhn)
    xhapsamp = zeros(nspatialsites,maxhn)
    ntotal = 0

    for kk in eachindex(ncomplexhap) #ncomplexhap is a 1D array of length the number of unique spatial sample sites

        h = real(ncomplexhap[kk])
        nsamp = imag(ncomplexhap[kk])

        if !isempty(nsamp)

            for j in eachindex(nsamp)
                n[kk, Int(h[j]+1)] = n[kk, Int(h[j]+1)] + nsamp[j]  # for a given sample point in a deme there may be 
                #contributions of the same haplotype from different time-points => sum over them (n array is initialised to be zero)

                ntotal = ntotal + nsamp[j]
            end

            if Nsamp[kk]!=0
                xhapsamp[kk,:] = n[kk, :]./Nsamp[kk] 
            end

        end          

    end



    #Check all chromosomes accounted for:
    println("n total chromosomes (from nt) = ", ntotal)
    println("sum(n) = ",sum(n))
    println("this should both be = ",sum(Nsampt)," the total number of chromosomes samples across all sample sites in Phase 2 including FX/SX haplotypes")
    println("and also equal to sum(Nsamp) = ",sum(Nsamp))

    println("size(n) = ",size(n))
    println("this should be (nspatialsites=",nspatialsites,")x(maxh=",maxh,")")

    #Remove any haplotypes that don't have more than 2 chromosomes across all samples
    #n is number of demes x number of haplotypes x max number samples, so want to sums over 1st dimension
    #This filter only applies for Anopheles gambiae: haplotypes there are determined by clustering of gene 
    #sequences, so a minimum-copy-count threshold is needed to distinguish a real haplotype from clustering noise. 
    #For Human data, haplotypes are called directly from known SNPs, so no such threshold is needed/applied.
    if Species=="AGambiae"

        nthreshold = 2
        nHarray = sum(n,dims=1) #total number of each haplotype across all spatial sample sites

        # println("nHArray = ",nHarray)

        indthr = vec(nHarray .> nthreshold)

        # println("indthr = ", indthr)
        n = n[:,indthr] #only keep those haplotypes that across all sample sites have more than nthreshold chromosomes (this also removes any zero columns)

        println("Number of chromosomes in simulation (after removing haplotypes that are sampled < ", nthreshold+1," times each): ",sum(n))
        println("Number of chromosomes removed:",sum(Nsamp)-sum(n))
        #N.B. Nsamp does not change which is deliberate! Since in the real data these are removed guard against rare haplotypes clusters due to sequencing errors and so for simplicity 
        #we effectively we reassign these as wild type. The sample mutant frequency from real data is calculated in this way.

    elseif Species=="Human"

        println("Species==\"Human\": skipping the nthreshold haplotype-count filter (not applicable — human haplotypes are called from known SNPs, not sequence clustering)")

    else

        error("Species must be \"AGambiae\" or \"Human\", got: ", Species)

    end

    println("size(n) = ",size(n))


#Check xhapsamp[kk,j] the frequency of the jth haplotype in the kkth unique spatial sample location is correctly constucted
    println("sum(xhapsamp,dims=2) = ",sum(xhapsamp,dims=2))
    



    nhaps,  = size(X[1])
    println("maxh = ",nhaps)
    println(length(Xhap))
    # println("maxt = ",T)




    #Calculate each summary statistic in turn 

    println("")
    println("Calculating number of independent origins η")
    ##Total number of origins η
    @time η, xmutsamp = NumberOfSampledOrigins(n,Nsamp)


    
    println("")
    println("Calculating mean and standard deviation of independent origins η over spatial sample sites")
    ##Mean number of origins η
    @time ηmean, ηstd = MeanVarianceNumberOfSampledOrigins(n)


    println("")
    println("Calculating Radius of gyration origins Rg")
    ##RMSD of origins — i.e. calculate rmsd of each origin and then average over all origins
    @time Rg = RadiusOfGyration(n,xx,yy)


    println("")
    println("Calculating average heterozygosity per spatial sample sites H")
    ##Average heterozygosity per deme (~Simpson's diversity)
    @time H = AveHeterozygosity(n,Nsamp) #N.B. using Ksite which is a ND x maxnumbsamplesperdeme array, since in this version n is number samples cummulated over all time points at a given samples site


    
    println("")
    println("Calculating average Nei's genetic distance between all pairs of spatial sample sites NeiD")
    # ##Average pair-wise Nei's "genetic" or haplotype distance between all demes 
    @time NeiD = AveNeiD(n,Nsamp)



    


    if pltfigs==1    
        println("")   
        println("Plotting haplotypes")

        xthr=0.01

        # pltdemes,skip, col = PlotHaplotypes(Xhap, xthr,xp,yp,wdeme,polygonAfrica,polygonGambiae,samplepoints)
        pltdemes,skip, col = PlotHaplotypes(Xhap, xthr,xp,yp,wdeme,polygonAfrica,polygonGambiae,samplepoints,nx,ny)

        println("skip = ", skip)

        if skip==0
 
            title!(string("η = ",η,"; Rg = ",round(Rg;digits=1),"kms; H = ",round(H;sigdigits=2),"; NeiD = ",round(NeiD;sigdigits=2),
                            ";  x_samp = ",round(xmutsamp;sigdigits=2),";  x_true = ",round(truefreq[end];sigdigits=2),";  xwt_true = ",round(1-truefreq[end];sigdigits=3)),titlefontsize=6)

            display(pltdemes)

            #Uncomment for pie chart plots on map
            # plot_geographic_pies!(xx, yy, vec(xoffset), vec(yoffset), Int.(n), 
            #                     rgb_matrix; #note that I could probably use col directly, by modifying function called, since it just converts it to RGB array!
            #                     marginfactor=2,
            #                     scalefactor=20,
            #                     transparency_alphaalpha=0.7,
            #                     haplotypelabels=nothing,
            #                     plot_original=false,
            #                     show_lines=true)

            savefig("HaplotypeMap.svg")

        end
    end


    tgen = t-1

    return η, Rg, H, NeiD, xmutsamp, Xt, tgen, ηmean, ηstd, n, Nsampt, Nsamp, nx, Ademe, xx, yy

end



function PlotHaplotypes(xhap,xthr,xp,yp,wdeme,polygonAfrica,polygonGambiae,samplepoints,nx,ny)
   

    #Determine the km-to-points scale from Plots.jl's own DEFAULT canvas size (not an explicit size= override — 
    #that was clobbering Plots' automatic title/margin layout) combined with the plotted domain's actual km extent. 
    #With aspect_ratio=:equal, whichever axis (width or height) is more constrained by the canvas determines the 
    #actual data-to-pixel scale that gets rendered, so we take the smaller of the two implied px/km ratios.
    #This is what keeps msize (further down) consistent with deme spacing (touching circles) for any domain 
    #shape/extent, since it always reflects the canvas Plots.jl is actually going to draw, not an assumed one.
    dpi_plt = 300
    defaultsize = Plots.default(:size) #(width_px, height_px) that Plots.jl will use for this plot
    xdom = [p.x for p in polygonAfrica]
    ydom = [p.y for p in polygonAfrica]
    domwidth = maximum(xdom)-minimum(xdom)
    domheight = maximum(ydom)-minimum(ydom)
    pxperkm = min(defaultsize[1]/domwidth, defaultsize[2]/domheight)

    pltdemes = plot(Tuple.(polygonAfrica), color = :gray, legend=false, grid=false, aspect_ratio=:equal, framestyle = :box, dpi=dpi_plt, tickfontsize=16)
    
    skip=0

    plot!(Tuple.(polygonGambiae), color = :black, legend=false, grid=false, aspect_ratio=:equal)



    
    nhaps = length(xhap)
    

    nht = 0
    for h=2:nhaps
        if !isempty(findall(x-> x>=xthr,xhap[h]))
            nht+=1
        end
    end

    println("number haplotypes with frequency somewhere >xthr = ",nht)

    if nht==1
        col = get(ColorSchemes.Spectral, range(0.0, 1.0, length=2))
        # col = cols[1]
    else
        col = get(ColorSchemes.Spectral, range(0.0, 1.0, length=nht)) #so number of colours is number of haplptypes exc WT
    end
    

    println("Number of Haplotypes above xthr (= ",xthr,") : ",nht)

    if nht>0
        skip=0
    else
        skip=1
    end

    nnh = 1
    for h = 2:nhaps #I think I need to turn into 2D array for each haplotype..


        if !isempty(findall(x-> x>=xthr,xhap[h]))
            




            xq = [q[1] for q in polygonGambiae]
            yq = [q[2] for q in polygonGambiae]


            #uncomment this for loop to plot contours using Contour.jl package
            # c = contours(xp,yp,xhap[h]',[xthr, 0.9])
            c = contours(xp,yp,xhap[h]',[0.1, 0.9])
            LW = [1,2]

            ii=1
            for cl in Contour.levels(c)
                lvl = Contour.level(cl)
 

                lw = LW[ii]
                ii+=1


                for line in lines(cl)
                    xc,yc = coordinates(line)

                    # plot!(xc,yc, color=col[nnh], linewidth=lw)

                    # xc = moving_average(xc, 3)
                    # yc = moving_average(yc, 3)
                   
                    #Uncomment below for interpolation and smoothing of contours (upto lines with "xc = rc.*cos.(θc) .+ xcm" and"yc = rc.*sin.(θc) .+ ycm")
                    #(But note it doesn't work well with rectangular domain)
                    # xcm = mean(xc)
                    # ycm = mean(yc)

                    # rc = sqrt.((xc .- xcm).^2 .+ (yc .- ycm).^2)
                    # θc = angle.((xc .- xcm) .+ im*(yc .- ycm))

                    # indθ = sortperm(θc)
                    # θc = θc[indθ]
                    # rc = rc[indθ]

                    # #add points at exactly ±π to ensure smoothed contour closes



                    # #Distance in angle from first point to -π
                    # #-π -θc[1] = -(π + θc[1]) so below is absolute angle not signed
                    # δθ1 = π + θc[1]
                    # #Distance in angle of last point to π
                    # δθ2 = π - θc[end]

                    # #need to account for if either δθ1 or δθ2 are zero

                    
                    # #Linear interpolation 
                    # rπ = rc[1] + (rc[end]-rc[1])/(δθ1+δθ2) * δθ1


                    # #To avoid duplicate knot warnings avoid identical values of angle (presumably -π and π)
                    # θc = [-π;θc;π]
                    # rc = [rπ;rc;rπ]

       

                    # #Apparently, this produces duplicate values of θc, so remove them
                    # θc, rc = deduplicate_data(θc, rc)
      
                    # # plot!(xc,yc, color=col[nnh], linewidth=2lw, ls=:dash)

                    # #Interpolate
                    # θi = collect(-π:π/500:π)

                    # ri = zeros(size(θi))
                    # idx = trues(size(θi))      

                    # intf = linear_interpolation(vec(θc), vec(rc), extrapolation_bc=Line())

                    # ri[idx] = intf(θi[idx])

                    # #Smooth again
                    # rc = moving_average(ri, 10)
                    # rπ = 1/2*(rc[1] + rc[end])
                    # rc[1] = rπ
                    # rc[end] = rπ

                    # θc = θi
                    # # rc = ri #effectively remove smoothing, so comment out for smoothing

                    # xc = rc.*cos.(θc) .+ xcm
                    # yc = rc.*sin.(θc) .+ ycm


                    #Determine if any of these points are outside Gambiae range and then snap to that if so
                    pointsc = Point.(vec(xc),vec(yc))
                    

                    for k in eachindex(pointsc)
                        if !isinside(pointsc[k],polygonGambiae; allowonedge=true)
                            #find the closest point on Gambie range
                            
                            ρ = (xq .- xc[k]).^2 + (yq .- yc[k]).^2
                            indmin = findall(ρ.==minimum(ρ))
                            indmin = indmin[1]
                            xc[k] = only(xq[indmin])
                            yc[k] = only(yq[indmin])
                        end
   
                    end


                    plot!(xc,yc, color=col[nnh], linewidth=lw)
                   


                end
            end


            #markersize (in points) computed from wdeme (deme diameter, km) using the SAME pxperkm/dpi_plt scale 
            #factor used to set the canvas size above — this is what keeps adjacent deme circles just-touching 
            #regardless of the plotted domain's extent or shape (previously this was a fixed empirical constant 
            #tied to D alone, which only worked because every plot happened to share the same Africa-map canvas).
            #NOTE: Plots.jl's markersize/canvas-size relationship can vary slightly by backend — if circles overlap 
            #or leave small gaps, adjust the msizescale multiplier below (it only needs tuning once, not per-domain).
            msizescale = 1.9
            points_per_km = pxperkm*72/dpi_plt
            msize = msizescale*wdeme*points_per_km


            indthr = findall(x-> x>=xthr,xhap[h])

            if !isempty(indthr)

                xxhap = zeros(length(indthr))
                yyhap = zeros(length(indthr))
                dd=1
                for d in indthr
                    yyhap[dd] = (d[1]-1)*wdeme
                    xxhap[dd] = (d[2]-1)*wdeme
                    dd+=1

                end

                happoints = Point.(vec(xxhap),vec(yyhap))
                # println("happoints = ",happoints)
                scatter!(Tuple.(happoints), ms=msize, ma = 0.3, markerstrokewidth=0.0, color=col[nnh]) #, markershape=:square
                # scatter(p,Tuple.(happoints), ms=msize, ma = 0.3, markerstrokewidth=0.0, color=col[nnh]) #, markershape=:square
            end





            nnh+=1
        end
    end


    #Plot where sample points are — doing here so they come after haplotype plotting and can be seen
    scatter!(Tuple.(samplepoints), ms=5, ma = 0.6, markerstrokewidth=0.0, markershape=:star5, color = :red)

  
    return pltdemes, skip, col

end

# Manual deduplication
function deduplicate_data(θc, rc)
    # Create a dictionary to collect all rc values for each θc
    θ_to_r = Dict{Float64, Vector{Float64}}()
    
    for (θ, r) in zip(θc, rc)
        if haskey(θ_to_r, θ)
            push!(θ_to_r[θ], r)
        else
            θ_to_r[θ] = [r]
        end
    end
    
    # Sort by angle and average duplicate r values
    θc_unique = sort(collect(keys(θ_to_r)))
    rc_unique = [mean(θ_to_r[θ]) for θ in θc_unique]
    
    return θc_unique, rc_unique
end



function DiscreteContour(nx,ny,wdeme,f, thr)
    

    #xp and yp are a list of x and y positions, respectively, of all demes
    #f is an array that contains the frequency at each deme position


    #find points for which f>thr
    ind = findall(x->x>=thr,f) #returns a list of Cartesian indices that obey this condition

    xc = []
    yc = []

    for i in ind
        #each i is one of these Cartesian indices
        iy = i[1]
        ix = i[2]

        #nearest neighbour index
        nnsubind=[[iy,ix-1],[iy-1,ix],[iy+1,ix],[iy,ix+1]]

        #for each calculate how many of these nearest neighbours are not part of ind
        nnn=0

        for kk=1:4

            iix = nnsubind[kk][2]
            iiy = nnsubind[kk][1]

            if (0 < iiy <=ny) & (0 < iix <= nx)
                if f[iiy,iix]<thr
                    nnn += 1
                end
            end

        end


        if nnn!=0 #i.e. there is at least one nn that is outside f>=thr
            #then this deme must be at the perimeter
            
            xc = [xc;(ix-1)*wdeme]
            yc = [yc;(iy-1)*wdeme]

        end

    end

    @show(xc)
    @show(yc)
    println("typeof(xc) = ",typeof(xc))
    println("")

 





    

    if !isempty(xc) | !isempty(yc)



        #calculate convex hull of points (i.e. the most exterior set of points that are convex)
        #package being used is LazySets.jl and it needs a vector of 2 element vectors!

        xyc = zeros(length(xc))
        #convert xc and yc to xyc a vector of 2-element vectors
        xyc = [[only(xc[i]), only(yc[i])] for i in eachindex(xc)]


        println("typeof(xyc) = ",typeof(xyc))
        println("")

        hull = convex_hull(xyc)

        xc = zeros(length(hull))
        yc = zeros(length(hull))

        for i in eachindex(hull)
            xc[i] = hull[i][1]
            yc[i] = hull[i][2]

        end

        ##Smooth points

        #Approximation of centre of mass: doesn't need to be exact I just want the origin 
        #for the polar co-ords to be encompassed within the contour so that -π<=θ<=π
        xcm = mean(xc)
        ycm = mean(yc)

        rc = sqrt.((xc .- xcm).^2 .+ (yc .- ycm).^2)
        θc = angle.((xc .- xcm) .+ im*(yc .- ycm))

        indθ = sortperm(θc)
        θc = θc[indθ]
        rc = rc[indθ]

        xc = rc.*cos.(θc) .+ xcm
        yc = rc.*sin.(θc) .+ ycm

        # @show[θc]
        # println("size(θc) = ", size(θc))
        # println("")


        #smooth
        rcs = moving_average(rc, 5) 

        # rr = 1/2*(rc[1]+rc[end])
        # rc[1] = rr
        # rc[end] = rr

        #add points at exactly ±π to ensure smoothed contour closes
        #Distance in angle from first point to -π
        δθ1 = π + θc[1]
        #Distance in angle of last point to π
        δθ2 = π - θc[end]

        #Linear interpolation 
        rπ = rcs[1] + (rcs[end]-rcs[1])/(δθ1+δθ2) * δθ1

        θcs = [-π;θc;π]
        rcs = [rπ;rc;rπ]

        xcs = rcs.*cos.(θcs) .+ xcm
        ycs = rcs.*sin.(θcs) .+ ycm

    else
        xcs=[]
        ycs=[]

    end






    return xc, yc, xcs, ycs



end






function DemePosition2D(k,nx,Ademe)

    wdeme = sqrt(Ademe)

    ky = Int(ceil((k)/nx))
    kx = Int(k - (ky-1)*nx)

    #Uncomment for zero starting indices
    #ky=ky-1;
    #kx=kx-1;

    x = wdeme*kx
    y = wdeme*ky

    return x,y

end




function NumberOfSampledOrigins(n,Nsamp)

    #n is a (number of sample sites)x(number of haplotypes)
    
    nHap = sum(n,dims=1) #Sum of sampled individuals for each haplotype(column) across all demes 

    Numbsamp = sum(Nsamp) #total number of sampled individuals before removing low copy number haplotypes
    # println("Number of total chromosomes sampled inc FX and SX = ",Numbsamp)
    nHap = nHap[2:end] #take just the mutants excluding WT
    # println("Number of individuals of each haploytpe exc WT: nHap = ",nHap)

    η = count(a->a>0, nHap) #n:number of origins
    println("η = ",η)
    # println("Number of mutants = ",sum(nHap))
    x  = sum(nHap)/Numbsamp #sample frequency of mutants
    println("xmut = ",x)


   
    return η, x
end

function MeanVarianceNumberOfSampledOrigins(n)

    #n is a (number of sample sites)x(number of haplotypes)   

    sz = size(n)
    nspatialsites = sz[1]
    ηdeme = zeros(nspatialsites)


    for k=1:nspatialsites
        nHap = n[k,:] #take the list of haplotypes for the kth spatial sampling site

        nHap = nHap[2:end] #take just the mutants excluding WT
        # println("Number of individuals of each haploytpe exc WT: nHap = ",nHap)

        ηdeme[k] = count(a->a>0, nHap) #n:number of origins
        # println("η = ",η)
        # println("Number of mutants = ",sum(nHap))
        # x  = sum(nHap)/Numbsamp #sample frequency of mutants



    end

    #Calculate mean and variance of number of origins 
    ind = ηdeme.!=0

    if !isempty(ind)
        ηsample = ηdeme[ind]

        ηmean = mean(ηsample)
        ηstd =  std(ηsample)
    else #there are no demes that have non-zero sampling of mutant haplotypes
        ηmean = 0
        ηstd = 0
    end

    #Check size
    # println("length(ηsample) = ",length(ηsample))
    println("<η> = ",ηmean)
    println("sqrt(<<η^2>>) = ",ηstd)


    return ηmean, ηstd
end



function RadiusOfGyration(n,x,y)

    # ND,nhaps = size(n)
    # ND,nhaps,nsamps = size(n)
    sz = size(n)
    ND = sz[1]
    nhaps = sz[2]

    # if length(sz)==3
    #     nsamps = sz[3]

    # else
    #     nsamps = 1
    # end

    if length(sz)==3
        nsamps = sz[3]
        numtimepoints = 1

    elseif length(sz)==4
        nsamps = sz[3]
        numtimepoints = sz[4]
    else
        nsamps = 1
        numtimepoints = 1
    end


    Rg=0.0
    # print("Rg=")
    # println(Rg)
    # println("Hello")

    nzerohaps = 0

    for h=2:nhaps

        # print("k=")
        # println(k)

        xcm = 0.0
        ycm = 0.0
        M=0.0

        #Calculate the centre of mass of hth origin
        for k=1:ND


            # x,y = DemePosition2D(k,nx,Ademe)
            # nxy[ky,kx] = n[j,k]
            # println(nxy[ky,kx])
            xcm = xcm + n[k,h]*x[k]
            ycm = ycm + n[k,h]*y[k]
            M += n[k,h] #M counts the total number of individuals for a given haplotype



        end 
        
        Rgk = 0.0
        if M!=0         
            xcm = xcm/M
            ycm = ycm/M
            nzerohaps += 1 #why do I have this? Because many of the haplotypes aren't sampled, so this counts the number of non-zero haplotypes
            
            #Calculate the radius of gyration 
        
            
            # print("Rgk=")
            # println(Rgk)
            for k=1:ND
                # println(j)
                # x,y = DemePosition2D(k,nx,Ademe)
                # nxy[ky,kx] = n[j,k]
                rsqd = (x[k]-xcm)^2 + (y[k]-ycm)^2
                Rgk += n[k,h]*rsqd
                # println(Rg)

            end        
        
            Rgk = sqrt(Rgk/M)
        end

        # print("Rgk=")
        # println(Rgk)
        
        Rg += Rgk
        # print("Rg=")
        # println(Rg)


    end

    if nzerohaps==0
        Rg = 0
    else
        Rg = Rg/nzerohaps #Normalise by number of non-zero sampled haplotypes
    end

    println("Rg = ",Rg)



    return Rg
end








function AveHeterozygosity(n,K)

    #K should be Nsamp (and not Nsampt)
    # ND = length(K) #ND is effectively nspatialsites

    #Now size of K is ND x1 
    # n is  (ND x nhaps) array..
    #in this function it doesn't care about where these samples are so I can just treat the extra samples 
    #as extra demes and reshape the arrays, so that calculation below works

    #Now K can also be ND x maxnumsamples x numuniquetimepoints — while n is still ND x maxnumhaplotypes x maxnumsamples
    #want to reshape K as above:
    #Can first reshape KK to be ND x newmaxnumsamples
    
    println("size(K) = ", size(K))
    println("size(n) = ", size(n))

    sz = size(K)

    # println("size(K) = ",sz)

    # if length(sz)>1
    #     nsamps = size(K)[2]
    # else
    #     nsamps = 1
    # end

    if length(sz)==1
        nsamps = 1        
    elseif length(sz) ==2
        nsamps = size(K)[2]        
    end

    # println("nsamps = ", nsamps)

    if nsamps>1
        Ktemp = K[:,1]
        ntemp = n[:,:,1]

        for k=2:nsamps
            Ktemp = [Ktemp;K[:,k]]
            ntemp = [ntemp;n[:,:,k]]
        end

        K = Ktemp
        n = ntemp


    end


    #find index to all demes that have more than 1 sample/individuals
    ind_nonzero = findall(x->(x>1),K)#don't need this anymore, since K is Nsamp not an NDx1 array

    #n is a (ND x nhaps) array and K is a (ND x 1) array of samples sizes in each deme, 
    #so I want to calculate n[j](n[j]-1)/K[j]/(K[j]-1) — but I can do this direcly by

    nn = n[ind_nonzero,:]
    KK = K[ind_nonzero]

    println("size(K) = ", size(KK))
    println("size(n) = ", size(nn))

    NND = length(KK)

    H = nn.*(nn.-1)./(KK)./(KK.-1) 
    #where this calculates element by element n*(n-1) for all origins and all demes 
    #N.B. that as KK is NND x 1 and nn is NND x nhaps, then we can use a dot multiplication or divide
    #such nn./KK means each nn[i,:]/KK[i]
    
    H = 1-sum(H)/NND #want to calculate sum over origins in each deme and then average over all demes*samples/deme, 
    #which is a sum scaled by ND: sum(H) sums all elements irrespective of dimension of array
    
    
    if isnan(H)
        println("H = ", H)
        println("ND = ",NND)
        println("K = ",K)

    end

    println("H = ", H)

    return H

end





function AveNeiD(x,K)

    #This is very slow

    #Average NeiD between all pairs of space-time sample sites:
    #if in space-time sample site k there are nh haplotypes and their number is n[h] for hth haplotype then we want to calculate 
    #NeiD[k,j] = sum_h n[h,k]n[h,j]/K[k]/K[j]/ {sqrt(sum_h(n[h,k].*n[h,k]/K[k]/K[k])*sum(n[h,j].*n[h,j]/K[j]/K[j]))}
    #where k and j are different sampling locations
    #Ave NeiD is then sum_k,j NeiD[k,j]/Npairs
    #where Npairs is the number of pairs of locations

    #n is  (ND x nhaps x maxnumsamples x numbuniquetimepoints) array.

    #K is ND x maxnumsamples x numuniquetimepoints 
    

    #Challenge is to have a sum over all pairs of demes, samples within demes at a given time point and then over time points


    #Now size of K is ND x maxnumsamples, where nsamps represents multiple samples from the same deme 
    #Similarly n is now  (ND x nhaps x maxnumsamples) array..
    #in this function it doesn't care about where these samples are so I can just treat the extra samples 

    # println("")
    #as extra demes and reshape the arrays, so that calculation below works
    
    # indK = findall(x->(x.>0),K)
    # println("length(indK) = ", length(indK)) #these aren't used!
    println("size(K) = ", size(K))
    println("size(n) = ", size(x))

    sz = size(x)

    ND = sz[1]
    nhaps = sz[2]


    npairs = 0

    # # HH = zeros(ND,nsamp,numbuniquetimepoints)
    ρ=0.0

    for k=1:ND
        for kk=1:ND
            if k!=kk

                Z1 = 0.0
                Z2 = 0.0
                ρρ = 0.0
                
                # println("K = ",K[k,j,t])
                # println("KK = ",K[kk,jj,tt])
                # K1 = only(K[k,j,t])
                # K2 = only(K[kk,jj,tt])
                K1 = K[k]
                K2 = K[kk]
                

                if K1>0 && K2>0
                    # println("K1 = ",K1)
                    # println("K2 = ",K2)
                    npairs += 1

                    # println("npairs = ",npairs)

                    for h=1:nhaps

                    

                        Z1 += x[k,h]^2
                        Z2 += x[kk,h]^2

                        ρρ += x[k,h]*x[kk,h]

                        

                    end     
                    
                    ρ += ρρ/sqrt(Z1*Z2)

                end


            end

        end
    end

    # npairs = npairs/2  #scale by 1/2 to correct for overcounting pairs
    
    # println("npairs = ",npairs)

    ρ = ρ/npairs




    # ρ = ρ/Npairs #Calculate average over all pairs of demes 
    println("ρ = ", ρ)
    # # NeiD = NeiD/(ND^2-ND) #Calculate average over all pairs of demes 

    NeiD = -log(ρ)
    println("Nei's Genetic distance = ",NeiD)

    # if isinf(NeiD) || NeiD<0 || isnan(NeiD)
    #     # println(K)
    #     # println(n)
    #     println("Nei's Genetic distance = ",NeiD)
    #     println("Correlation ρ=",ρ)
    #     # println("Number of demes = ", NDsamp)
    #     println("Number of space-time samples location pairs = ", npairs)
    # end

    return NeiD

end




function AveNeiD_optimised(x, K)

    #Should input K as Nsamp
    # Collect all valid samples where K > 0
    samples = findall(x -> x > 0, K) #this is a vector of CartesianIndices, each the same dimension as K
    M = length(samples)
    M < 2 && return NaN  # Not enough samples to form pairs


    sz =  size(x)

    ND = sz[1]
    nhaps = sz[2]
    # nsamps = sz[3]
    # numbuniquetimepoints = sz[4]

    if length(sz)==3
        nsamps = sz[3]
        numtimepoints = 1

    elseif length(sz)==4
        nsamps = sz[3]
        numtimepoints = sz[4]
    else
        nsamps = 1
        numtimepoints = 1
    end



    # ND, nhaps, nsamps, numbuniquetimepoints = size(x)

    

    # Precompute vectors v and their squared norms z for each valid sample
    X = Matrix{Float64}(undef, nhaps, M)
    z = Vector{Float64}(undef, M)
    for (α, s) in enumerate(samples) #syntax  for (index,value) in enumerate(array)
                                     #so as samples is an array of CartesianIndices, α is the ath index of this array and s is the ath CartesianIndex 
        @views xx = x[s, :] #essentially a list of sample frequencies of each haplotype
        X[:, α] .= xx #this vector is the list of sample frequencies of each haplotype for the αth sample location x_hα
        z[α] = xx'*xx #sum(abs2, x) #z[α] is effectively the (self-correlation)^2 of haplotype type frequencies for the αth location


    end

    # Check the X matrix has columns that sum to 1
    # println("sum(X,dims=1) = ", sum(X,dims=1))

    # Compute covariance matrix and pairwise terms
    C = X' * X  # Covariance matrix of dot products #V is a matrix nhapsxM, where X[h,α] is the sample frequency of the hth haplotype in αth location
                # so C is M x M so C[α,α']  = sum_h X[h,α]X[h,α'] is the unnormalised correlation/covariance of haplotype structure between α and α' location   
    z_sqrt = sqrt.(z)  #this is the sqrt(self-correlation) in each sample location
    Z = z_sqrt .* z_sqrt' #where Z = z_sqrt * z_sqrt' is MxM matrix the same size as C and Z[α,α'] = sqrt(z[α])*sqrt(z[α'])
    S = C ./ Z # Normalize by product of norms,  

    # # Calculate total sum of off-diagonal elements and average
    # println("")
    # println("sum(S) = ",sum(S))
    # println("M = ",M)
    total = sum(S) - M  # Subtract diagonal elements (each is 1.0)
    npairs = M * (M - 1) #it is M*(M-1)/2 for the upper diagonal, so *2 for all which is what above calculates..
    # println("number of pairs = ", npairs)
    ρ = total / npairs  # Average over unordered pairs
    println("Correlation ρ = ", ρ)
    NeiD = -log(ρ)
    # println("NeiD = -ln(ρ) = ", NeiD)


    # # Handle numerical issues
    # if isinf(NeiD) || NeiD < 0 || isnan(NeiD)
    #     println("Nei's Genetic distance = ", NeiD)
    #     println("Correlation ρ = ", ρ)
    #     println("Number of space-time sample pairs = ", npairs ÷ 2)
    # end

    return NeiD
end









function moving_average(A::AbstractArray, m::Int)
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m÷2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    return out
end





function Calculate_Closest_Demes_For_SampleLocations(x,y,demepoints,ny,DemesInAfrica,wdeme)

    #x,y are x and y coords of sample sites in Africa 
    #demepoints are a list of demes positions as a Points array
    #DemesInAfrica is a ND length vector of true positive values of whether that deme is in Gambiae range
    #Nsamp is a vector of numbers of chromosome samples in each sample site in x,y 
    #sampletime is the times at which each of the locations in x,y were sampled

    #this is very slow for large ND
    #surely I can just loop through each sample location (x[i],y[i]) calculate x[i]/wdeme
    
    println("length(x) = ",length(x))
    # println("length(samplingtime) = ",length(samplingtime))
    # println("samplingtime = ",samplingtime)
    
    #List all demes in Gambiae range, by x, y pairs
    xy = zeros(length(findall(DemesInAfrica)),2)
    println("size(xy) = ", size(xy))

    kxy = zeros(Int,length(x))

    ii=1    

    for k in eachindex(DemesInAfrica)

        # Get linearindex
        # k = LinearIndices(CartInd)[i]

        if DemesInAfrica[k]

            xy[ii,1] = demepoints[k].x
            xy[ii,2] = demepoints[k].y

            ii+=1
        end

    end


    println("")



    for i in eachindex(x)    #loop over all sample points #I don't know which of these correspond to the same sample location

            

        ρ = sqrt.( (xy[:,1] .- x[i]).^2 + (xy[:,2] .- y[i]).^2 )


        indmin = findall(x->x.==minimum(ρ),ρ)

        

        #x and y position of closest deme
        xpi = only(xy[indmin,1])
        ypi = only(xy[indmin,2])

   
        indx = round(Int,xpi/wdeme) + 1
        indy = round(Int,ypi/wdeme) + 1


        # #Find linear index this closest deme corresponds to (for the ith sample site and jth time point)
        kxy[i] = indy+(indx-1)*ny



    end

    return kxy


end


function plot_geographic_pies!(x::Vector{T}, y::Vector{T}, 
                               xoffset::Vector{T}, yoffset::Vector{T},
                               data::Matrix{Int},  # nsites × nhaplotypes
                               colors::Matrix{Float64};  # nhaplotypes × 3 RGB
                               marginfactor::Real=2,
                               scalefactor::Real=20,
                               transparency_alpha::Real=0.7,
                               haplotypelabels::Union{Nothing, Vector{String}}=nothing,
                               plot_original::Bool=true,
                               show_lines::Bool=true) where T <: Real
    
    nsites = length(x)
    nhaplotypes = size(data, 2) #this includes WT
    
    # Validate inputs
    @assert length(y) == nsites "x and y must have same length"
    @assert length(xoffset) == nsites "xoffset must match x"
    @assert length(yoffset) == nsites "yoffset must match y"
    @assert size(data, 1) == nsites "data rows must match number of sites"
    @assert size(colors, 1) >= nhaplotypes "colors must have at least nhaplotypes rows"
    
    # Calculate pie radii
    piesums = sum(data, dims=2)[:, 1]  # Total at each site
    pieradii = sqrt.(piesums) .* scalefactor
    
    # Plot lines from original positions to pie centers
    if show_lines
        for i in 1:nsites
            if piesums[i] > 0
                xx = x[i] + marginfactor * pieradii[i] * xoffset[i]
                yy = y[i] + marginfactor * pieradii[i] * yoffset[i]

                #Draw line emanating from centre but starting from edge of pie chart
                #circle:
                #Calculate x and y co-ordinates of this point
                Z = sqrt((xx-x[i])^2 + (yy-y[i])^2) # => unit vector along this direction is u = 1/Z [x-xx, y-yy] and position vector from [xx,yy] is r*u
                r = pieradii[i]
                xedge = xx + r/Z*(x[i]-xx)
                yedge = yy + r/Z*(y[i]-yy)
                
                # Draw line
                plot!([x[i], xedge], [y[i], yedge], 
                      color=:black, linewidth=1, linestyle=:solid, 
                      label=false)
            end
        end
    end
    
    # Plot original positions
    if plot_original
        scatter!(x, y, color=:red, marker=:pentagon, markersize=6, 
                 label=false, markeralpha=0.8)
    end
    
    # Plot pie charts
    for i in 1:nsites
        if piesums[i] == 0
            continue
        end
        
        # Calculate frequencies
        frequencies = data[i, :] ./ piesums[i]
        
        # Calculate pie center
        xx = x[i] + marginfactor * pieradii[i] * xoffset[i]
        yy = y[i] + marginfactor * pieradii[i] * yoffset[i]
        
        # Plot pie chart
        plot_pie_at_location!(xx, yy, frequencies, pieradii[i], 
                              colors[1:nhaplotypes, :], i, transparency_alpha)
    end
    
    # Add legend
    # if !isnothing(haplotypelabels)
    #     create_legend!(colors[1:nhaplotypes, :], haplotypelabels, transparency_alpha)
    # else
    #     # Create default labels
    #     default_labels = ["WT"; ["H$j" for j in 1:nhaplotypes-1]]
    #     create_legend!(colors[1:nhaplotypes, :], default_labels, transparency_alpha)
    # end

    create_legend!(colors, haplotypelabels, transparency_alpha)
    
    return pieradii
end

"""
    plot_pie_at_location!(center_x, center_y, frequencies, radius, colors, site_number, transparency_alpha)
    
Plot a single pie chart at the specified location.
"""
function plot_pie_at_location!(center_x::Real, center_y::Real, 
                               frequencies::Vector{Float64}, 
                               radius::Real, 
                               colors::Matrix{Float64},
                               site_number::Int,
                               transparency_alpha::Real)
    
    # Skip if no data
    if sum(frequencies) ≈ 0
        return
    end
    
    # Calculate wedge angles
    percentages = frequencies ./ sum(frequencies)
    cumpercentages = cumsum(percentages)
    start_angles = [0.0; cumpercentages[1:end-1]] .* 2π
    end_angles = cumpercentages .* 2π
    
    # Plot each wedge
    for k in 1:length(frequencies)
        if frequencies[k] ≈ 0
            continue
        end
        
        # Create wedge with rotation (matching MATLAB: -theta + π/2)
        wedge = create_wedge(center_x, center_y, radius, 
                             start_angles[k], end_angles[k], 
                             colors[k, :], transparency_alpha)
        
        # Plot the wedge
        plot!(wedge[1, :], wedge[2, :], 
              fill=true, fillcolor=RGBA(colors[k, 1], colors[k, 2], colors[k, 3], transparency_alpha),
              linewidth=0.25, linecolor=RGB(0.5, 0.5, 0.5),
              label=false)
    end
    
    # Optional: Add site number text (uncomment if needed)
    # annotate!(center_x, center_y, text(string(site_number), 10, :black))
end

"""
    create_wedge(center_x, center_y, radius, start_angle, end_angle, color, transparency_alpha)
    
Create wedge vertices for a pie slice.
"""
function create_wedge(center_x::Real, center_y::Real, 
                      radius::Real, start_angle::Real, end_angle::Real,
                      color::Vector{Float64}, transparency_alpha::Real)
    
    # Number of points for smooth arc
    npoints = max(50, ceil(Int, (end_angle - start_angle) * 50 / (2π)))
    
    # Generate angles with rotation (matching MATLAB: -theta + π/2)
    angles = range(start_angle, stop=end_angle, length=npoints)
    transformed_angles = -angles .+ π/2  # Apply rotation
    
    # Create wedge vertices
    wedge_x = vcat([center_x], 
                   center_x .+ radius .* cos.(transformed_angles), 
                   [center_x])
    wedge_y = vcat([center_y], 
                   center_y .+ radius .* sin.(transformed_angles), 
                   [center_y])
    
    return [wedge_x wedge_y]'
end

"""
    create_legend!(colors, haplotypelabels, transparency_alpha)
    
Add a legend to the current plot with colored patches, matching MATLAB behavior.
"""
function create_legend!(colors::Matrix{Float64}, 
                        haplotypelabels::Union{Nothing, Vector{String}}=nothing,
                        transparency_alpha::Real=0.7)
    
    n = size(colors, 1)
    
    # Create default labels if none provided
    if haplotypelabels === nothing
        haplotypelabels = Vector{String}(undef, n)
        for i in 1:n
            if i == 1
                haplotypelabels[i] = "WT"
            else
                haplotypelabels[i] = "H$(i-1)"
            end
        end
    end
    
    # Create dummy scatter plots for legend entries
    # Using NaN coordinates so they don't appear in the plot
    for i in 1:n
        scatter!([NaN], [NaN], 
                 markershape=:rect,  # Square marker for legend
                 markercolor=RGBA(colors[i, 1], colors[i, 2], colors[i, 3], transparency_alpha),
                 markersize=8,
                 markerstrokewidth=0.5,
                 markerstrokecolor=RGB(0.5, 0.5, 0.5),
                 label=haplotypelabels[i])
    end
    
    # Return the labels in case they were generated
    return haplotypelabels
end

