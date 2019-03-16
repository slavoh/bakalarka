using PyPlot

# zatial na mensej skale
no_setups = 9
setup_size = 10^3

function generate_polyeder( dA, db ) # given by Ax≥b
    is_bounded=false
    while !is_bounded
        C=randn(m+2*n, n)  #constraints matrix
        perm=shuffle([1:m+2n])
        B=zeros(n,n)
        s=zeros(n)
        for i=1:n
            B[i,:]=C[perm[i],:]
            s[i]=(B[i,:]')*(B[i,:])
        end
        Binv=pinv(B)
        V=Binv*(s')

        c=randn(n)
        y=c'*Binv

        for i=1:n
            if y[i]<0   # otoc znamienko rovnice
                C[i,:]=-C[i,:]
            end
        end

        is_bounded=true
        for int i=1:n
            if simplexmin[i] < -1 || simplexmax[i]>1    #TODO
                is_bounded=false
                break
            end
        end
    end


    A=randn(dA, db)
    b=randn(db)
    # TODO dopisat podla clanku may1982

    return (A,b)
end

function find_MVEE( A, b, γ, eff_target ) # using REX algorithm
    # function REX(Fx, supp.ini, ver=1, γ=4, eff=1-1e-9, it.max=Inf, t.max=30)
    δ = 1e-14
    ϵ = 1e-24
    n = size(Fx,1)
    m = size(Fx,2)
    eff.inv = 1/eff
    n.iter = 0
    L = min(n, γ*m)
    lx.vec = zeros(L)
    index = [1:n]
    one = ones(m)

    supp = supp.ini # TODO
    K = size(supp,1)
    Fx.supp = Fx[supp, :]
    w = zeros(n)
    w[supp] = 1 / size(supp,1)
    w.supp = w[supp]
    M = (sqrt(w.supp) * Fx.supp)'*((sqrt(w.supp) * Fx.supp))
    d.fun = ((Fx * (cholfakt(pinv(M))')^2)' * one / m
    ord = sort(d.fun) # decreasing
    lx.vec = shuffle(tail(ord,L))
    kx.vec = shuffle(supp)

    while true
        n.iter = n.iter + 1
        ord1 = which.min(d.fun[supp]) # TODO
        kb = supp[ord1]
        lb = ord[1]
        v = [kb, lb]
        cv = Fx[v, :] * (M \ Fx[v, :]')
        α = 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + δ)
        α = min(w[kb], α)
        w[kb] = w[kb] - α
        w[lb] = w[lb] + α
        M = M + α * ((Fx[lb,:])*(Fx[lb,:]') - (Fx[kb,:])*(Fx[kb, :]'))

        if ((w[kb] < δ) && (ver==1)) # LBE is nullifying and the version is 1
            for l = 1:L
                lx = lx.vec[l]
                Alx = (Fx[lx, :])*(Fx[lx, :]')
                for k = 1:K
                    kx = kx.vec[k]
                    v = [kx, lx]
                    cv = Fx[v, :] * (M \ Fx[v, ]')
                    α = 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + ϵ)
                    α = min(w[kx], max(-w[lx], α))
                    wkx.temp = w[kx] - α
                    wlx.temp = w[lx] + α
                    if ((wkx.temp < δ) || (wlx.temp < δ))
                        w[kx] = wkx.temp
                        w[lx] = wlx.temp
                        M = M + α * (Alx - (Fx[kx, :])*(Fx[kx, :]'))
                    end
                end
            end
        else # LBE is non-nullifying or the version is 0
            for l = 1:L
                lx = lx.vec[l]
                Alx = (Fx[lx, :])*(Fx[lx, :]')
                for k = 1:K
                    kx = kx.vec[k]
                    v = [kx, lx]
                    cv = Fx[v, :] * (M \ Fx[v, :]')
                    α = 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + δ)
                    α = min(w[kx], max(-w[lx], α))
                    w[kx] = w[kx] - α
                    w[lx] = w[lx] + α
                    M = M + α * (Alx - (Fx[kx,:])*(Fx[kx,:]'))
                end
            end
        end

        supp = index[w > δ] # TODO
        K = size(supp,1)
        w.supp = w[supp]
        d.fun = ((Fx * (cholfakt(pinv(M))')^2) * one / m
        ord.ind = (1:n)[d.fun >= -sort(-d.fun, partial=L)[L]]
        ord = ord.ind[order(d.fun[ord.ind], decreasing=TRUE)]
        # The two lines above can be replaced by simpler but usually
        # somewhat slower ord = order(d.fun, decreasing=TRUE)[1:L]
        lx.vec = sample(ord)
        kx.vec = sample(supp)

        eff.act =  1 / d.fun[ord[1]]
        if ((d.fun[ord[1]] < eff.inv) || (n.iter >= it.max))
            break
        end
    end

    Phi.best = det(M)^(1/m)
    eff.best = 1/d.fun[ord[1]]
    list(w.best=w, Phi.best=Phi.best, eff.best=eff.best, n.iter=n.iter, t.act=t.act)

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
            ub=min(ub, (b[l]-sum(A[l,-j]*x[-j])) / A[l][j] )
        elseif A[l][j] < ϵ
            lb=min(ub, (b[l]-sum(A[l,-j]*x[-j])) / A[l][j] )
        end
    end
    return rand(Uniform(lb, ub), 1, 1 )
end

times=zeros(no_setups,3)
for setup=1:no_setups # initiate setup
    dimension=2^(setup)
    X=zeros(setup_size, dimension) # zoznam vygenerovanych bodov - nie je nutny
    (A,b)=generate_polyeder(dimension, dimension) # TODO mozno ine rozmery
    # TODO sprav H reprezentaciu

    # REX generate
    starttime=time()
    P=find_MVEE(A, b )
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
