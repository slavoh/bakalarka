using PyPlot

# zatial na mensej skale
no_setups = 3
setup_size = 10^5

function generate_polyeder( dA, db ) # given by Ax≥b
    A=randn(dA, db)
    b=randn(db)
    # TODO dopisat generovanie nahodneho polyedru

    #print("generate_polyeder| A:", A,", b: ", b, "\n")
    return (A,b)
end

function find_MVEE( A, b ) # using REX algorithm
    P=randn(size(b,1), size(b,1))
    # TODO naprogramuj REX

    #print("find_MVEE| P:", P, "\n")
    return P
end

function generate_in_ball( d )
    x=randn( d )
    return x/(x'*x)
    # TODO preskalovat vzdialenost pomocou inverznej funkcie
end

function generate_in_MVEE( P )
    x=generate_in_ball( size( P, 2 ) )
    # TODO napisat vyraz pomocou P

    return P*x
end

function in_polyhedra(A, x, b) # check whether Ax≥b
    c=A*x-b
    for i=1:size(c,1)
        if(c[i]<0)
            return false
        end
    end
    return true
end

times=zeros(no_setups)
for setup=1:no_setups
    dimension=10^(setup)
    X=zeros(setup_size, dimension) # zoznam vygenerovanych bodov - nie je nutny

    (A,b)=generate_polyeder(dimension, dimension) # TODO mozno ine rozmery
    starttime=time()
    P=find_MVEE(A, b)
    for i=1:setup_size
        X[i,:]=generate_in_MVEE(P)
        while false # (!in_polyhedra(A,X[i,:],b))
            print(".") # na debuggovanie
            X[i,:]=generate_in_MVEE(P)
        end
    end
    endtime=time()
    times[setup]=endtime-starttime
end

print("average time: ", mean(times), "\n")
scatter(times, [1:no_setups])
xlabel("log dimension")
ylabel("time runned")


# scatter(randn(), grads, label="gradient computation")
# title(string("λ=", λ))
# xlabel("minibatch size")
# legend()
# savefig(string("/home/slavo/Desktop/project/graphs/timeplot_", lamb, "_", init) )
# close()
