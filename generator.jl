using PyPlot
using Polyhedra
using CDDLib
using Convex
using SCS
using Distributions

set_default_solver(SCSSolver(verbose=0))
const ϵ = 1e-24
const δ = 1e-14

function generate_polyeder( dim, no_planes ) # dany Ax≥b
    while true
        A=randn(no_planes, dim)
        b=randn(no_planes)

        x0=rand(dim)*2-ones(dim)
        for i in range(1,no_planes)
            if A[i,:]⋅x0 < b[i]
                A[i,:] *= -1
                b[i] *= -1
            end
        end

        bounded = true
        x = Variable(dim)
        for i=1:2*dim
            if i>dim
                p = maximize(x[i-dim])
            else
                p = minimize(x[i])
            end
            p.constraints += A*x >= b
            solve!(p)
            if !(p.status == :Optimal)
                bounded = false
                break
            end
            # if p.optval*p.optval < 1  # TODO mozno odkomentovat
            #     bounded = false
            #     break
            # end
        end
        if !bounded
            continue
        end

        VRep=collect(Polyhedra.points(polyhedron(hrep(-A,-b), CDDLibrary())))
        vertices=zeros(length(VRep), dim)
        i=1
        for point in VRep
            for j in range(1,length(point))
                vertices[i,j]=point[j]
            end
            i+=1
        end

        # odstran zbytocne rovnice z H reprezentacie
        bounding_set = []
        for i=1:no_planes
            bounds=false
            for j=1:size(vertices, 1)
                if abs(A[i,:]⋅vertices[j,:] - b[i]) <ϵ
                    bounds=true
                    break
                end
            end
            if bounds
                append!(bounding_set,i)
            end
        end
        A_reduced = A[bounding_set,:]
        b_reduced = b[bounding_set]

        return (A_reduced,b_reduced, vertices, x0)
    end
end

function find_MVEE(Fx, supp_ini, γ=4, eff=1-1e-9, it_max=10^12, t_max=30) # pouzitim REX algoritmu
    n = size(Fx,1)
    m = size(Fx,2)
    L = min(n, γ*m)

    supp = supp_ini
    w = zeros(n)
    w[supp] = 1 / length(supp)
    M = (sqrt.(w[supp]) .* Fx[supp,:])'*((sqrt.(w[supp]) .* Fx[supp,:]))
    @show det(M)
    for n_iter = 1:it_max
        K = length(supp)
        w_supp = w[supp]
        M_inv=inv(M)
        d_fun = ((Fx * (chol( (M_inv+M_inv')/2 ))' ).^2) * ones(m) / m
        ord = reverse(sortperm(d_fun))[1:L]
        if d_fun[ord[1]] < 1/eff
            break
        end
        lx_vec = shuffle(ord)
        kx_vec = shuffle(supp)

        n_iter += 1
        kb = supp[ findmin(d_fun[supp])[2] ]
        lb = ord[1]
        v = [kb; lb]
        cv = Fx[v, :] * (M \ Fx[v, :]')
        α = 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + δ)
        α = min(w[kb], α)
        w[kb] -= α
        w[lb] += α
        M += α * ((Fx[lb,:])*(Fx[lb,:]') - (Fx[kb,:])*(Fx[kb, :]'))
        M = (M+M')/2    # pre istotu
        for l = 1:L
            lx = lx_vec[l]
            Alx = Fx[lx, :]*Fx[lx, :]'
            for k = 1:K
                kx = kx_vec[k]
                v = [kx; lx]
                cv = Fx[v, :] * (M \ Fx[v, :]')
                α = 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + ϵ)
                α = min(w[kx], max(-w[lx], α))
                wkx_temp = w[kx] - α
                wlx_temp = w[lx] + α
                if ((w[kb] >= δ) || (wkx_temp < δ) || (wlx_temp < δ))
                    w[kx] = wkx_temp
                    w[lx] = wlx_temp
                    M += α * (Alx - (Fx[kx,:])*(Fx[kx,:]'))
                    M = (M+M')/2    # pre istotu
                end
            end
        end
        supp = (1:n)[find(λ -> (λ>δ), w)]
    end

    # vypocitaj MVEE
    reg = Fx[:, 2:m]
    Z = zeros(m-1)
    H0 = zeros(m-1,m-1)
    for i=1:n
        Z += w[i]*reg[i,:]
    end
    for i=1:n
        H0 += w[i]*(reg[i,:]-Z)*(reg[i,:]-Z)'
    end
    H=inv(H0)/(m-1)
    return (H,Z)
end

function generate_on_sphere( dim )
    x_sym = randn( dim )
    return (x_sym/norm(x_sym))
end

function generate_in_MVEE( H )
    return H*( generate_on_sphere(size(H,1))*rand(Uniform(0,1))^(1/size(H,1) ) )
end

function gibbs(x, A, b)
    ϵ=1e-14
    for j = 1:size(A,2)
        lb = -10^14
        ub = 10^14
        for l=1:size(A,1)
            if A[l,j] > ϵ
                lb = max(lb, (b[l] -sum(A[l,:]⋅x) +A[l,j]*x[j] ) /A[l,j])
            elseif A[l,j] < -ϵ
                ub = min(ub, (b[l] -sum(A[l,:]⋅x) +A[l,j]*x[j] ) /A[l,j])
            end
        end
        x[j] = rand(Uniform(lb, ub))
    end
    return x
end

no_setups = 10
setup_size = 10^5

times=zeros(no_setups,3)
generations=zeros(no_setups)
print("\nProgram started\n")
for setup=2:no_setups  # inicializacia testu
    dimension=setup
    # dimension=2^(setup)
    X=zeros(setup_size, dimension) # zoznam vygenerovanych bodov - nie je nutny
    (A,b,vertices,x0)=generate_polyeder(dimension, dimension*10)
    print("Polyhedra generated:  ")
    @show dimension, size(A,1), size(vertices,1)

    # # REX generator
    # starttime=time()
    # ff = [ones(size(vertices,1)) vertices]
    # qq = randperm(size(vertices,1))[1:dimension+4]
    # (H,Z)=find_MVEE(ff, qq)
    # # (H,Z)=find_MVEE([ones(size(vertices,1)) vertices], randperm(size(vertices,1))[1:dimension+4])
    # print("MVEE found\n")
    # for i=1:setup_size
    #     X[i,:]=generate_in_MVEE(H)
    #     count=1
    #     while any(x ->(x<δ), A*X[i,:]-b)
    #         count+=1
    #         if count==10^7
    #             break
    #         end
    #         X[i,:]=generate_in_MVEE(H)
    #     end
    #     if count==10^7
    #         @show setup
    #         generations[setup] = 10^8
    #         break
    #     end
    #     generations[setup] += count
    #     # @show count
    # end
    # endtime=time()
    # times[setup,1]=endtime-starttime
    # generations[setup] /= setup_size

    # Hit-and-Run generator
    burn=100
    starttime=time()
    w=zeros(size(A,1))
    x=deepcopy(x0)
    error()
    for i=(1-burn):setup_size
        D = generate_on_sphere(size(A,2))

        for j=1:size(A,1)
            w[j]=(A[j,:]⋅x-b[j])/(A[j,:]⋅D)
        end
        @show w

        w_neg= filter(y->(y<0), w)
        w_pos= filter(y->(y>0), w)
        lb=0
        ub=0
        if length(w_neg)>0
            lb = maximum(w_neg)
        end
        if length(w_pos)>0
            ub = minimum(w_pos)
        end
        if lb==ub
            @show w
            error()
        end
        dist=rand(Uniform(lb,ub))
        x += dist*D

        if i>0
            X[i,:]=x
        end
    end
    # TODO shuffle(X[i,:])
    endtime=time()
    times[setup,2]=endtime-starttime

    # Gibbs generator
    starttime=time()
    burn=100
    x_next=deepcopy(x0)
    for i=(1-burn):setup_size
        x=x_next
        x_next = gibbs(x, A, b)
        if i>0
            X[i,:]=x
        end
    end
    # TODO shuffle(X[i,:])
    endtime=time()
    times[setup,3]=endtime-starttime
    print("end of setup\n")
end

print("average time: ", mean(times), "\n")

# scatter(1:no_setups, times[:,1], label="REX")
scatter(1:no_setups, 1000*times[:,2], label="Hit-and-Run")
scatter(1:no_setups, 1000*times[:,3], label="Gibbs")
# xlabel("log dimension")
ylabel("time runned [ms]")

# scatter(1:no_setups, generations, label="REX")
# ylabel("pocet pokusov")
legend()
# close()
