##############################
#Exercise 1
##############################

function hello_world()
    println("hello world")
end

function hello(name::String)
    println("hello $name")
end

##############################
#Exercise 2
##############################

function f(n::Int)
    if n%2 == 0
        return div(n,2)
    else return 3*n+1
    end
end

function g(n::Int)
    k=0
    while n != 1
        n = f(n)
        k += 1
    end
    return k
end

function h()
    n=101
    d=Dict{Int,Int}() #Map the k to the n
        while true
            k=g(n)
            if haskey(d, k)
                return (d[k], n)
            end
            d[k]=n
            n +=1
        end
    end

############################################
#Exercise 3
############################################

function pascal_triangle(n::Int)
    for i=0:(n-1)
        print(" " ^(n-i))
        for k=0:i
            print(binomial(i,k), " ")
        end
        println(" ")
    end
end