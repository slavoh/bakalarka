using PyPlot

# zatial na mensej skale
no_setups = 3
setup_size = 10^3

function generate_polyeder( dA, db ) # given by Ax≥b
    A=randn(dA, db)
    b=randn(db)
    # TODO dopisat podla clanku may1982

    return (A,b)
end

function find_MVEE( A, b, γ, eff_target ) # using REX algorithm
    P=randn(size(b,1), size(b,1))

    w=randn( size(b,1) )
    w=w/(w'*w)
    while eff_target < eff
        # LBE step using w
        k=argmin gw cez nosic
        l=argmax gw cez 𝔛
        α=argmax Φ(M( w + α2(e_l-e_k) ))

        w=w+α(e_l-e_k)

        # set support
        supp=zeros(0)
        for ...
            push!(supp, ...)
        end
        K=size(supp,1)

        # set greedy
        L=min(ceiling(Int, γ*m), n)
        Sgreedy= L najvacsich prvkov gw

        # subspace step
        K_perm=shuffle(supp)
        L_perm=shuffle(Sgreedy)
        for k=1:K
            for l=1:L
                α=argmax Φ(M( w + α2(e_l-e_k) ))
                if α==w[K_perm[k]] or α==-w[L_perm[l]]
                    w[K_perm[k]] = w[K_perm[k]] - α
                    w[L_perm[l]] = w[L_perm[l]] + α
                end
            end
        end
    end

    return P
end

function generate_in_MVEE( P )
    x_sym=randn( size( P, 2) )    # TODO dat d-rozmerne normalne rozdelenie
    x_ball=(x_sym/sqrt(x_sym'*x_sym))^(1/d)
    return P*x_ball
end

function is_in_polyhedra(A, x, b) # check whether Ax≥b
    c=A*x-b
    for i=1:size(c,1)
        if(c[i]<0)
            return false
        end
    end
    return true
end

function gibbsrD(x,j, A, b)
    lb=-∞
    ub=∞
    for i=1:size(A,1)
        if(A[1][j]>ϵ)
            ub=min(c(ub, (b[l]-sum(A[l,-j]*x[-j])) / A[l][j] ))
        elseif A[l][j] < ϵ
            lb=min(c(ub, (b[l]-sum(A[l,-j]*x[-j])) / A[l][j] ))
        end
    end
    return runif(1, min=lb, max=ub)
end

times=zeros(no_setups,3)
for setup=1:no_setups # initiate setup
    dimension=10^(setup)
    X=zeros(setup_size, dimension) # zoznam vygenerovanych bodov - nie je nutny
    (A,b)=generate_polyeder(dimension, dimension) # TODO mozno ine rozmery

    # REX generate
    starttime=time()
    P=find_MVEE(A, b)
    for i=1:setup_size
        X[i,:]=generate_in_MVEE(P)
        while (!is_in_polyhedra(A,X[i,:],b))
            X[i,:]=generate_in_MVEE(P)
            # print(".") # na debuggovanie
        end
    end
    endtime=time()
    times[setup,1]=endtime-starttime

    # Hit-and-Run generate
    starttime=time()
    for i=1:setup_size
        while true
            x=zeros(dimension)
            for j=1:dimension
                x[j]=randn(simplexmin(), simplexmax()) # TODO
            end
            if is_in_polyhedra(A,x,b)
                break
        end
    end
    endtime=time()
    times[setup,2]=endtime-starttime

    # Gibbs generate
    starttime=time()
    # TODO Gibbs
    ϵ=10^(-14)
    burn=100
    for i=2:(setup_size+burn)
        x[i,:] = x[i-1,:]
        for j=1:size(A,2)
            x[i, j] = gibbsrD(x[i,:], j)
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
