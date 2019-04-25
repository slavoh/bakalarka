using PyPlot
using Polyhedra
using CDDLib
# using Distributions
using Convex
using SCS


# zatial na mensej skale
no_setups = 3
setup_size = 10^3
set_default_solver(SCSSolver(verbose=0))

function generate_polyeder( dim, no_planes ) # dany Ax≥b
    # count=0
    while true
        # count += 1
        # if mod(count,1000)==0
        #     @show count
        # end
        A=randn(no_planes, dim) # TODO
        b=randn(no_planes)      # TODO prenormovat?

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
            # if p.optval*p.optval < 1  # TODO
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
                if (A[i,:]⋅vertices[j,:] - b[i])^2 <ϵ
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

const ϵ = 1e-24
const δ = 1e-14
function find_MVEE(Fx, supp_ini, γ=4, eff=1-1e-9, it_max=Inf, t_max=30) # pouzitim REX algoritmu
    n = size(Fx,1)
    m = size(Fx,2)
    eff_inv = 1/eff
    n_iter = 0
    L = min(n, γ*m)
    lx_vec = zeros(L)
    index = 1:n
    one = ones(m)

    supp = supp_ini
    K = length(supp)
    Fx_supp = Fx[supp, :]
    w = zeros(n)
    w[supp] = 1 / K
    w_supp = w[supp]

    # M=zeros(m, m)
    # for i in supp
    #     M += Fx[i,:]*Fx[i,:]'
    # end
    # M /= size(supp,1)
    # if findmin(M'-M)[1]>10e-10
    #     print("ERROR1")
    # end
    M = (sqrt.(w_supp) .* Fx_supp)'*((sqrt.(w_supp) .* Fx_supp))
    M_inv=inv(M)
    @show rank(M)
    d_fun = ((Fx * (chol((M_inv+M_inv')/2))' ).^2) * one / m
    ord = reverse(sortperm(d_fun))[1:L]
    lx_vec = shuffle(ord)
    kx_vec = shuffle(supp)

    print("zaciatok cyklu\n")
    while true
        n_iter += 1
        ord1 = findmin(d_fun[supp])[2]
        kb = supp[ord1]
        lb = ord[1]
        v = [kb; lb]
        cv = Fx[v, :] * (inv(M) * Fx[v, :]')     # pre istotu
        # cv = Fx[v, :] * (M \ Fx[v, :]')
        α = 0_5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + δ)
        α = min(w[kb], α)
        w[kb] -= α
        w[lb] += α
        M += α * ((Fx[lb,:])*(Fx[lb,:]') - (Fx[kb,:])*(Fx[kb, :]'))
        M = (M+M')/2    # pre istotu


        if (w[kb] < δ) # LBE je nulujuci
            print("nulujuci\n")
            for l = 1:L
                lx = lx_vec[l]
                Alx = Fx[lx, :]*Fx[lx, :]'
                for k = 1:K
                    kx = kx_vec[k]
                    v = [kx, lx]
                    print(rank(M))

                    cv = Fx[v, :] * (inv(M) * Fx[v, :]')     # pre istotu
                    # cv = Fx[v, :] * (M \ Fx[v, :]')

                    α = 0_5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + ϵ)
                    α = min(w[kx], max(-w[lx], α))

                    wkx_temp = w[kx] - α
                    wlx_temp = w[lx] + α
                    if ((wkx_temp < δ) || (wlx_temp < δ))
                        w[kx] = wkx_temp
                        w[lx] = wlx_temp
                        M += α * (Alx - (Fx[kx,:])*(Fx[kx,:]'))
                        M = (M+M')/2    # pre istotu
                    end
                end
            end
        else # LBE je nenulujuci
            print("nenulujuci\n")
            for l = 1:L
                lx = lx_vec[l]
                Alx = Fx[lx, :]*Fx[lx, :]'
                for k = 1:K
                    kx = kx_vec[k]
                    v = [kx; lx]
                    cv = Fx[v, :] * (inv(M) * Fx[v, :]')    # pre istotu
                    # cv = Fx[v, :] * (M \ Fx[v, :]')
                    α = 0_5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + δ)
                    α = min(w[kx], max(-w[lx], α))
                    w[kx] -= α
                    w[lx] += α
                    M += α * (Alx - (Fx[kx,:])*(Fx[kx,:]'))
                    M = (M+M')/2    # pre istotu
                end
            end
        end
        supp = index[find(λ -> (λ>δ), w)]
        K = length(supp)
        w_supp = w[supp]
        # print("\n")
        # @show size(M)
        @show rank(M)
        M_inv=inv(M)
        # @show rank(M_inv)
        if findmin(M'-M)[1]>10e-10    # pre istotu
            print("ERROR1")
        end
        d_fun = ((Fx * (chol( (M_inv+M_inv')/2 ))' ).^2) * one / m
        ord = reverse(sortperm(d_fun))[1:L]

        lx_vec = shuffle(ord)
        kx_vec = shuffle(supp)

        eff_act =  1 / d_fun[ord[1]]
        @show eff_act

        if ((d_fun[ord[1]] < eff_inv) || (n_iter >= it_max))
            break
        end
    end

    Phi_best = det(M)^(1/m)
    eff_best = 1/d_fun[ord[1]]

    Z=0
    H=0
    for i=1:n
        Z=Z+w[i]*z[i]
    end
    for i=1:n
        H=H+w[i]*(z[i]-Z)'*(z[i]-Z)
    end
    H=inv(H)/(m-1)
    return (H,Z)
end

function generate_on_sphere( dim )
    x_sym = randn( dim )
    return (x_sym/norm(x_sym))
end

function generate_in_MVEE( P )
    return P*( ( generate_on_sphere(size(P,2))*Uniform(0,1) )^(1/d) )
end

function generate_in_MVEE( P )
    x_sym = randn( size(P,2) )
    x_ball = (x_sym/norm(x_sym)*Uniform(0,1))^(1/d)
    return P*x_ball
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

print("\nProgram started\n")
times=zeros(no_setups,3)
for setup=2:no_setups
    # inicializacia testu
    dimension=2^(setup)
    X=zeros(setup_size, dimension) # zoznam vygenerovanych bodov - nie je nutny
    (A,b,vertices,x0)=generate_polyeder(dimension, dimension*10)
    print("Polyhedra generated\n")
    @show dimension

    # REX generator
    starttime=time()
    P=find_MVEE([ones(size(vertices,1)) vertices], randperm(size(vertices,1))[1:dimension+4])
    for i=1:setup_size
        X[i,:]=generate_in_MVEE(P)
        while any(x ->(x<0), A*X[i,:]-b)
            X[i,:]=generate_in_MVEE(P)
            # print("_") # na debuggovanie
        end
    end
    endtime=time()
    times[setup,1]=endtime-starttime

    # Hit-and-Run generator
    # burn=100
    # starttime=time()
    # w=zeros(size(A,1))
    # x=deepcopy(x0)
    # for i=(1-burn):setup_size
    #     D = generate_on_sphere(size(A,2))
    #
    #     for j=1:size(A,1)
    #         w[j]=(A[j,:]⋅x-b[j])/(A[j,:]⋅D)
    #     end
    #
    #     lb = maximum(filter(y->(y<0), w))
    #     ub = minimum(filter(y->(y>0), w))
    #     dist=rand(Uniform(lb,ub))
    #
    #     x += dist*D
    #
    #     if i>0
    #         X[i,:]=x
    #     end
    # end
    # TODO shuffle(X[i,:])
    # endtime=time()
    # times[setup,2]=endtime-starttime
    #
    # # Gibbs generator
    # starttime=time()
    # burn=100
    # x_next=deepcopy(x0)
    # for i=(1-burn):setup_size
    #     x=x_next
    #     x_next = gibbs(x, A, b)
    #     if i>0
    #         X[i,:]=x
    #     end
    # end
    # TODO shuffle(X[i,:])
    # endtime=time()
    # times[setup,3]=endtime-starttime
end

print("average time: ", mean(times), "\n")

scatter([1:no_setups], times[:,1], label="REX")
scatter([1:no_setups], times[:,2], label="Hit-and-Run")
scatter([1:no_setups], times[:,3], label="Gibbs")
xlabel("log dimension")
ylabel("time runned")
legend()
# close()
