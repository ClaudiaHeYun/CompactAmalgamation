using Oscar

# generates partitions whose substituted Schur polynomials have appropriate degrees
function partitions(v,n,m)
    if v[1] < v[2]
        v1 = v[2]
        v2 = v[1]
    else
        v1 = v[1]
        v2 = v[2]
    end

    list_p = [ Partition([1,1,1]) ]
    cnt = 1
    r = 1
    while r <= (2*(n-1)*m)/(2v1+v2)
        t = max(0,ceil(Int,((1-n)*m+(v1+v2)*r)/v2))
        while t <= ((n-1)*m-v1*r)/v2 && t <= r
            p = (t==0) ? Partition([r]) : Partition([r,t])
            push!(list_p, p)
            cnt += 1
            t += 1
        end
        r += 1
    end
    return (list_p,cnt)
end

# R is a univariate Laurent polynomial ring
function z_polynomial(R,n,m)
    p1 = R(0)
    p2 = R(0)
    for i in 0:n-1
        p1 += gen(R)^i
        p2 += gen(R)^(-i)
    end
    return p1^m*p2^m
end

function linear_form_to_vector(c)
    i = 1
    vec = coeff(c,1) * Vector{fmpz}(exponent_vector(c,1))
    while i<length(c)
        i += 1
        vec += coeff(c,i) * Vector{fmpz}(exponent_vector(c,i))
    end
    return AffineHyperplane(vec,-constant_coefficient(c))
end

function print_schur(parts, num)
    k = length(parts)
    nonzero = []
    for i = 1:k
        if num[i] != 0
            push!(nonzero,i)
        end
    end
    ans = ""
    for i = 1:length(nonzero)
        if num[nonzero[i]] == 1
            if i == length(nonzero)
                curly_bracket = replace(string("s_",parts[nonzero[i]]),"["=>"{","]"=>"}")
                ans = string(ans,curly_bracket)
            else
                curly_bracket = replace(string("s_",parts[nonzero[i]]," + "),"["=>"{","]"=>"}")
                ans = string(ans,curly_bracket)
            end
        else
            if i == length(nonzero)
                curly_bracket = replace(string(num[nonzero[i]],"s_",parts[nonzero[i]]),"["=>"{","]"=>"}")
                ans = string(ans,curly_bracket)
            else
                curly_bracket = replace(string(num[nonzero[i]],"s_",parts[nonzero[i]]," + "),"["=>"{","]"=>"}")
                ans = string(ans,curly_bracket)
            end
        end
    end

    return ans
end
    
function verify_solution(sol,v,n,m,N)
    R,z = LaurentPolynomialRing(QQ, "z")
    @show target = N*z_polynomial(R,n,m)
        
    parts, k = partitions(v,n,m)
    @show k;
    P = 0
    for i in 1:k
        P+= numerator(sol[i])*schur_polynomial(parts[i],3)
    end

    println("P = ",print_schur(parts, sol))

    Pv = P(z^(v[1]),z^(v[2]),z^(v[3]));
    @show Pv;

    print("dimension is ")
    @show P(1,1,1);
    
    Qx, x = PolynomialRing(QQ, "x");
    _, xi = NumberField(x^2 + x + 1, "xi");
    
    print("check 3rd root of unity ")
    @show P(xi,xi,xi);

    @show Pv == target;
end

function feasible_region(v,n,m,N)
    parts, k = partitions(v,n,m)

    R, a = PolynomialRing(ZZ, "a"=>1:k)
    LR, z = LaurentPolynomialRing(R, "z")
    
    target_z_poly = N*z_polynomial(LR,n,m)
    
    Pv = LR(0)
    for i in 1:k
        s = schur_polynomial(parts[i],3)
        sv = s(z^(v[1]),z^(v[2]),z^(v[3]))
        Pv += a[i]*sv
    end

    pz = Pv-target_z_poly

    min_d = AbstractAlgebra.Generic.trail_degree(pz)
    max_d = AbstractAlgebra.Generic.lead_degree(pz)

    coeffs_pz = map(x -> coeff(pz, x), min_d:max_d)
    coeffs_pz = filter(x -> x != 0, coeffs_pz)
    
    M = [ linear_form_to_vector(c) for c in coeffs_pz ]
    
    A = fill(0,(k,k))
            
    for i in 1:k
        A[i,i] = -1
    end

    b = fill(0,k)

    return Polyhedron((A,b),M)
end
    
function ilp_conj(v,n,m,N)
    P = feasible_region(v,n,m,N)
    _, k = partitions(v,n,m)
    obj = fill(1,k)
    return MixedIntegerLinearProgram(P,obj;k=0,convention = :min)
end
    
function lp_conj(v,n,m,N)
    P = feasible_region(v,n,m,N)
    _, k = partitions(v,n,m)
    obj = fill(1,k)
    return LinearProgram(P,obj;k=0,convention = :min)
end