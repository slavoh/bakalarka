using PyPlot
using Polyhedra
using CDDLib
using Convex
using Distributions


# zatial na mensej skale
no_setups = 9
setup_size = 10^3

function generate_polyeder( dim, no_planes ) # given by Ax≥b
    while true
        A=randn(no_planes, dim)
        b=randn(no_planes)

        c=randn(dim)
        for i in range(1,dim)
            if A[i,:]'*c < b[i]
                A[i,:] *= -1
                b[i] *= -1
            end
        end
        VRep=collect(Polyhedra.points(polyhedron(hrep(A,b), CDDLibrary())))
        if length(VRep)>0
            # odstran zbytocne rovnice z H reprezentacie
            # count=0
            # bounding = zeros(no_planes)
            # for i in range(0, no_planes)
            #     if minimum(A*vertices[i]-b)=0
            #         bounding[i] = 1
            #         break
            #     end
            # end
            # TODO remove bounding==0

            vertices=randn(length(VRep), dim)
            i=1
            for point in VRep
                for j in range(1,length(point))
                    vertices[i,j]=point[j]
                end
                i+=1
            end

            return (A,b, vertices)
        end
    end
end

function find_MVEE(Fx, supp_ini, ver=1, γ=4, eff=1-1e-9, it_max=Inf, t_max=30) # using REX algorithm
    δ = 1e-14
    ϵ = 1e-24
    n = size(Fx,1)
    m = size(Fx,2)
    eff_inv = 1/eff
    n_iter = 0
    L = min(n, γ*m)
    lx_vec = zeros(L)
    index = [1:n]
    one = ones(m)

    supp = supp_ini # TODO
    K = size(supp,1)
    Fx_supp = Fx[supp, :]
    w = zeros(n)
    w[supp] = 1 / size(supp,1)
    w_supp = w[supp]
    Q=sqrt.(w_supp)⋅Fx_supp
    M=Q'*Q
    # M = (sqrt.(w_supp) .* Fx_supp)'*((sqrt.(w_supp) .* Fx_supp))
    d_fun = ((Fx * (chol(pinv(M)))' ).^2) * one / m
    ord = reverse(sort(d_fun))
    lx_vec = shuffle(ord)[1:L]
    kx_vec = shuffle(supp)

    while true
        n_iter = n_iter + 1
        ord1 = findmin(d_fun[supp])[2]
        # print(ord1)
        kb = supp[ord1]
        lb = ord[1]
        v = [kb, lb]
        cv = Fx[v, :] * (M \ Fx[v, :]')
        α = 0_5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + δ)
        α = min(w[kb], α)
        w[kb] = w[kb] - α
        w[lb] = w[lb] + α
        M = M + α * ((Fx[lb,:])*(Fx[lb,:]') - (Fx[kb,:])*(Fx[kb, :]'))

        if ((w[kb] < δ) && (ver==1)) # LBE is nullifying and the version is 1
            for l = 1:L
                lx = lx_vec[l]
                Alx = (Fx[lx, :])*(Fx[lx, :]')
                for k = 1:K
                    kx = kx_vec[k]
                    v = [kx, lx]
                    cv = Fx[v, :] * (M \ Fx[v, ]')
                    α = 0_5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + ϵ)
                    α = min(w[kx], max(-w[lx], α))
                    wkx_temp = w[kx] - α
                    wlx_temp = w[lx] + α
                    if ((wkx_temp < δ) || (wlx_temp < δ))
                        w[kx] = wkx_temp
                        w[lx] = wlx_temp
                    end
                end
            end
        else # LBE is non-nullifying or the version is 0
            M = M + α * (Alx - (Fx[kx, :])*(Fx[kx, :]'))
            for l = 1:L
                lx = lx_vec[l]
                Alx = (Fx[lx, :])*(Fx[lx, :]')
                for k = 1:K
                    kx = kx_vec[k]
                    v = [kx, lx]
                    cv = Fx[v, :] * (M \ Fx[v, :]')
                    α = 0_5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + δ)
                    α = min(w[kx], max(-w[lx], α))
                    w[kx] = w[kx] - α
                    w[lx] = w[lx] + α
                    M = M + α * (Alx - (Fx[kx,:])*(Fx[kx,:]'))
                end
            end
        end

        supp = index[x -> (x>δ), w]
        K = size(supp,1)
        w_supp = w[supp]
        d_fun = ((Fx * (chol(pinv(M)))' ).^2) * one / m
        ord = reverse(sort(d_fun))[1:L]

        lx_vec = shuffle(ord)
        kx_vec = shuffle(supp)

        eff_act =  1 / d_fun[ord[1]]
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
    H=pinv(H)/(m-1)
    return (H,Z)
end

function generate_in_MVEE( P )
    x_sym=MvNormal( size( P, 2), 1 )
    x_ball=(x_sym/norm(x_sym)*Uniform(0,1))^(1/d)
    return P*x_ball
end

function gibbsrD(x,j, A, b)
    lb=-10^9
    ub=10^9
    ϵ=1e-14
    for l=1:size(A,1)
        if(A[1,j]>ϵ)
            ub=min(ub, (b[l]-sum(A[l,length(x)-j]*x[length(x)-j])) / A[l,j] )
        elseif A[l,j] < ϵ
            lb=min(ub, (b[l]-sum(A[l,length(x)-j]*x[length(x)-j])) / A[l,j] )
        end
    end
    return rand(Uniform(lb, ub), 1, 1 )
end

times=zeros(no_setups,3)
for setup=2:no_setups # initiate setup
    dimension=2^(setup)
    X=zeros(setup_size, dimension) # zoznam vygenerovanych bodov - nie je nutny
    (A,b,vertices)=generate_polyeder(dimension, dimension*2)
    # VRep=Polyhedra.points(polyhedron(SimpleHRepresentation(A,b), CDDLibrary()))
    # VRep=collect(Polyhedra.points(polyhedron(hrep(A,b), CDDLibrary())))

    # REX generate
    starttime=time()
    P=find_MVEE(vertices, [])
    for i=1:setup_size
        X[i,:]=generate_in_MVEE(P)
        while any(x ->(x<0), A*X[i,:]-b)
            X[i,:]=generate_in_MVEE(P)
            # print("_") # na debuggovanie
        end
    end
    endtime=time()
    times[setup,1]=endtime-starttime

    # Hit-and-Run generate
    burn=100
    starttime=time()
    for i=(1-burn):setup_size
        x=randn(dimension)
        D= # generate U_ball TODO
        for j=1:dimension
            w[j]=(A[j,:]'*x-b[j])/(A[j,:]'*D)
        end
        for j=1:dimension
            x_min = Variable(dimension)
            prob_min = maximize(x_min, [x_min<0])
            prob_max = maximize(x_max, [A*x >= b])
            solve!(prob_min)
            solve!(prob_max)
            w=Uniform(prob_min.optval, prob_max.optval)
            x += w*D
            if i>0
                X[i,:]=x
            end
        end
        if all(λ ->(λ<=0),A*x-b)
            break
        end
    end
    endtime=time()
    times[setup,2]=endtime-starttime

    # Gibbs generate
    starttime=time()
    burn=100
    for i=2:(setup_size+burn)
        X[i,:] = X[i-1,:]
        for j=1:size(A,2)
            X[i, j] = gibbsrD(X[i,:], j, A, b)
        end
    end
    endtime=time()
    times[setup,3]=endtime-starttime
end

print("average time: ", mean(times), "\n")

scatter([1:no_setups], times[:,1], label="REX")
scatter([1:no_setups], times[:,2], label="Hit-and-Run")
scatter([1:no_setups], times[:,3], label="Gibbs")
xlabel("log dimension")
ylabel("time runned")
legend()
