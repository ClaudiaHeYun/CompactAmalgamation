using Oscar

function partitions_lex(k)
    list_p = [ Partition([1,1,1]) ]
    cnt = 1
    m = 1
    n = 0
    while cnt < k
        p = (n==0) ? Partition([m]) : Partition([m,n])
        push!(list_p, p)
        cnt += 1
        n += 1
        if n>m
            m += 1
            n = 0
        end
    end
    return list_p
end

# first k schur polynomials
function schur_lex(k)
    partitions = partitions_lex(k)
    return [schur_polynomial(partition,3) for partition in partitions]
end

function dot_prod(l1,l2)
    ans = 0
    for i = 1:length(l1)
        ans += l1[i]*l2[i]
    end
    return ans
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

function linear_form_to_vector(c)
    i = 1
    vec = coeff(c,1) * Vector{fmpz}(exponent_vector(c,1))
    while i<length(c)
        i += 1
        vec += coeff(c,i) * Vector{fmpz}(exponent_vector(c,i))
    end
    return LinearHyperplane(vec)
end

function display_solution(sol,k,v,w)
    sol_z = [numerator(a) for a in sol]
    sol_v = sol_z[1:k]
    sol_w = sol_z[k+1:2*k]
    
    partitions = partitions_lex(k)
    
    print("In the Schur basis, P = ")
    println(print_schur(partitions, sol_v))
    println("==================")
    print("In the Schur basis, Q = ")
    println(print_schur(partitions, sol_w))
    println("==================")
    
    list_s = schur_lex(k)
    
    @show P = dot_prod(sol_v,list_s)
    println("==================")
    @show Q = dot_prod(sol_w,list_s)
    println("==================")
    
    print("dimension is ")
    @show P(1,1,1);
    
    println("==================")
    
    LR,z = LaurentPolynomialRing(ZZ,"z")
    @show Pz = P(z^(v[1]),z^(v[2]),z^(v[3]));
end

function solution_hyperplanes(k,v,w)
    R, _ = PolynomialRing(ZZ, :a=>1:k, :b=>1:k)
    a = gens(R)[1:k]
    b = gens(R)[k+1:2*k]
    LR, z = LaurentPolynomialRing(R, "z")

    schur_list = schur_lex(k)
    lv = [s(z^(v[1]),z^(v[2]),z^(v[3])) for s in schur_list]
    lw = [s(z^(w[1]),z^(w[2]),z^(w[3])) for s in schur_list]
    
    pqz = LR(0)

    for i = 1:k
        pqz += a[i]*lv[i] - b[i]*lw[i]
    end
    
    min_d = AbstractAlgebra.Generic.trail_degree(pqz)
    max_d = AbstractAlgebra.Generic.lead_degree(pqz)
    
    coeffs_pqz = map(x -> coeff(pqz, x), min_d:max_d)
    coeffs_pqz = filter(x -> x != 0, coeffs_pqz)
    
    M = [ linear_form_to_vector(c) for c in coeffs_pqz ]

    return M
end

function solution_cone(k,v,w)
    M = solution_hyperplanes(k,v,w)
    C = cone_from_inequalities(-identity_matrix(ZZ,2*k),M)
    
    return C
end

function schur_deg_mod_3_lex(k)
    parts = partitions_lex(k)
    return [sum(part) % 3 == 0 ? 0 : 1 for part in parts]
end

function schur_dim_lex(k)
    schurs = schur_lex(k)
    return [s(1,1,1) for s in schurs]
end

function id_matrix(n)
    A = fill(0,(n,n))
            
    for i in 1:n
        A[i,i] = 1
    end
    return A
end

function inequalities(k;integer=false,upper_bound=10000)
    inj = -schur_deg_mod_3_lex(k)
    A = [-id_matrix(2*k);transpose(vcat(inj,fill(0,k)));transpose(vcat(fill(0,k),inj))]
            
    b = fill(0,2*k+2)
    b[2*k+1] = -1
    b[2*k+2] = -1

    if integer
        bound = fill(1,(1,2*k))
        A_bounded = [A; bound]
        b_bounded = [b; upper_bound]
        return A_bounded,b_bounded
    else
        return A,b
    end
end

function feasible_region(k,v,w; integer=false,upper_bound=10000)
    M = solution_hyperplanes(k,v,w)
    A,b = inequalities(k;integer,upper_bound)
    return Polyhedron((A,b),M)
end

function lp(k,v,w;integer=false, upper_bound=10000)
    obj = vcat(schur_dim_lex(k),fill(0,k))
    P = feasible_region(k,v,w; integer,upper_bound)
    if integer
        return MixedIntegerLinearProgram(P, obj; convention = :min)
    else
        return LinearProgram(P,obj;k=0,convention = :min)
    end
end

function verify_solution(sol,k,v,w)
    sol_z = [numerator(a) for a in sol]
    sol_v = sol_z[1:k]
    sol_w = sol_z[k+1:2*k]
    
    
    list_s = schur_lex(k)
    
    P = dot_prod(sol_v,list_s)
    Q = dot_prod(sol_w,list_s)

    LR, z = LaurentPolynomialRing(ZZ, "z")
    Pv = P(z^(v[1]),z^(v[2]),z^(v[3]));
    Qw = Q(z^(w[1]),z^(w[2]),z^(w[3]));

    print("dimension ");
    @show P(1,1,1);
    
    @show Pv == Qw;

    Qx, x = PolynomialRing(QQ, "x");
    _, xi = NumberField(x^2 + x + 1, "xi");
    
    print("check 3rd root of unity ");
    @show P(xi,xi,xi) == P(1,1,1);
    @show Q(xi,xi,xi) == Q(1,1,1);
end

function extract_vector(scip_sol_file)
    f = open(scip_sol_file,"r")
    t = read(f,String);
    close(f)

    re = r"x[0-9]+\s+([e0-9.-]+)"
    m = eachmatch(re,t)
    res = collect(m)
    vec = [round(Int64,parse(Float64,String(num[1]))) for num in res]
    return vec
end

# hydra: /opt/local/scip/8.0.3/bin/scip
# my laptop: /Users/yun/Documents/scipoptsuite-8.0.2/scip/bin/scip
function scip(filename,outputname)
    run(`/Users/yun/Documents/scipoptsuite-8.0.2/scip/bin/scip -f $filename -l $outputname -q`);
end

function feasible_scip(k,v,w)
    filename = string(v[1],v[2],"_",w[1],w[2],"_",k,".lp")
    outputname = string(v[1],v[2],"_",w[1],w[2],"_",k,"_output.txt")
    write_to_lp(k,v,w,filename)
    scip(filename,outputname)
    outputfile = open(outputname,"r")
    output = read(outputfile,String);
    close(outputfile)
    
    return occursin("optimal solution found",output)
end

function feasible(k,v,w; using_scip=false)
    if using_scip
        return feasible_scip(k,v,w)
    else
        return is_feasible(feasible_region(k,v,w))
    end
end

function find_min_k(max_k,v,w;min_k=1,using_scip=false)
    if feasible(max_k,v,w; using_scip)
        low = min_k - 1
        high = max_k + 1
        searching = true
        while searching
            k = floor(Int,(low+high)/2)
            @show k;
            if !feasible(k,v,w; using_scip)
                low = k
            else
                high = k
            end

            if high - low == 1
                searching = false
            end
        end

        print("minimal k is ",high)
        #return high
    else
        print("no solution found")
        #return -1
    end
end

# the following two functions exist because of a bug in OSCAR:
# save_lp does not work if the lp is infeasible
# however, determining feasibility is much faster in SCIP
function write_to_lp(k,v,w,filename)
    R, x = PolynomialRing(ZZ, "x"=>1:2*k)
    a = gens(R)[1:k]
    b = gens(R)[k+1:2*k]
    LR, z = LaurentPolynomialRing(R, "z")

    dim_list = schur_dim_lex(k)

    io = open(filename, "w")
    write(io, "MINIMIZE\n")
    ln = ""
    for i in 1:k
        ln = string(ln,"+",dim_list[i]," x",i," ")
    end
    ln = string("\tobj: ",ln,"\n")
    write(io, ln)
    write(io, "Subject To\n")
    for i in 1:2*k
        ln = string("\tie",i-1,": +1 x",i," >= 0\n")
        write(io, ln)
    end

    deg_list = schur_deg_mod_3_lex(k)
    
    inj_p = ""
    inj_q = ""
    for i in 1:k
        if deg_list[i] != 0
            inj_p = string(inj_p,"+1 x",i," ")
            inj_q = string(inj_q,"+1 x",i+k," ")
        end
    end
        
    inj_p = string("\tie",2*k,": ",inj_p,">= 1\n")
    write(io, inj_p)
    inj_q = string("\tie",2*k+1,": ",inj_q,">= 1\n")
    write(io, inj_q)

    schur_list = schur_lex(k)
    lv = [s(z^(v[1]),z^(v[2]),z^(v[3])) for s in schur_list]
    lw = [s(z^(w[1]),z^(w[2]),z^(w[3])) for s in schur_list]
    
    pqz = LR(0)

    for i = 1:k
        pqz += a[i]*lv[i] - b[i]*lw[i]
    end
    
    min_d = AbstractAlgebra.Generic.trail_degree(pqz)
    max_d = AbstractAlgebra.Generic.lead_degree(pqz)
    
    coeffs_pqz = map(x -> coeff(pqz, x), min_d:max_d)
    coeffs_pqz = filter(x -> x != 0, coeffs_pqz)
    
    for j in 1:length(coeffs_pqz)
        eq = coeffs_pqz[j]
        eq_s = ""
        for i in 1:length(eq)
            mon = monomial(eq,i)
            mon = replace(string(mon),"_"=>"","["=>"","]"=>"","{"=>"","}"=>"")
            c = coeff(eq,i)
            if c > 0
                eq_s = string(eq_s,"+",c," ",mon," ")
            else
                eq_s = string(eq_s,"-",-c," ",mon," ")
            end
        end
        ln = string("\teq",j-1,": ",eq_s,"= 0\n")
        write(io, ln)
    end
    
    write(io, "BOUNDS\n")
    for i in 1:2*k
        ln = string("\tx",i," free\n")
        write(io, ln)
    end
        
    write(io, "END")
    close(io)
end

function write_to_milp(k,v,w,filename)
    R, x = PolynomialRing(ZZ, "x"=>1:2*k)
    a = gens(R)[1:k]
    b = gens(R)[k+1:2*k]
    LR, z = LaurentPolynomialRing(R, "z")
    
    dim_list = schur_dim_lex(k)

    io = open(filename, "w")
    write(io, "MINIMIZE\n")
    ln = ""
    for i in 1:k
        ln = string(ln,"+",dim_list[i]," x",i," ")
    end
    ln = string("\tobj: ",ln,"\n")
    write(io, ln)
    write(io, "Subject To\n")
    for i in 1:2*k
        ln = string("\tie",i-1,": +1 x",i," >= 0\n")
        write(io, ln)
    end

    deg_list = schur_deg_mod_3_lex(k)
    
    inj_p = ""
    inj_q = ""
    for i in 1:k
        if deg_list[i] != 0
            inj_p = string(inj_p,"+1 x",i," ")
            inj_q = string(inj_q,"+1 x",i+k," ")
        end
    end
        
    inj_p = string("\tie",2*k,": ",inj_p,">= 1\n")
    write(io, inj_p)
    inj_q = string("\tie",2*k+1,": ",inj_q,">= 1\n")
    write(io, inj_q)

    schur_list = schur_lex(k)
    lv = [s(z^(v[1]),z^(v[2]),z^(v[3])) for s in schur_list]
    lw = [s(z^(w[1]),z^(w[2]),z^(w[3])) for s in schur_list]
    
    pqz = LR(0)

    for i = 1:k
        pqz += a[i]*lv[i] - b[i]*lw[i]
    end
    
    min_d = AbstractAlgebra.Generic.trail_degree(pqz)
    max_d = AbstractAlgebra.Generic.lead_degree(pqz)
    
    coeffs_pqz = map(x -> coeff(pqz, x), min_d:max_d)
    coeffs_pqz = filter(x -> x != 0, coeffs_pqz)
    
    for j in 1:length(coeffs_pqz)
        eq = coeffs_pqz[j]
        eq_s = ""
        for i in 1:length(eq)
            mon = monomial(eq,i)
            mon = replace(string(mon),"_"=>"","["=>"","]"=>"","{"=>"","}"=>"")
            c = coeff(eq,i)
            if c > 0
                eq_s = string(eq_s,"+",c," ",mon," ")
            else
                eq_s = string(eq_s,"-",-c," ",mon," ")
            end
        end
        ln = string("\teq",j-1,": ",eq_s,"= 0\n")
        write(io, ln)
    end
    
    write(io, "BOUNDS\n")
    for i in 1:2*k
        ln = string("\tx",i," free\n")
        write(io, ln)
    end
    
    write(io, "GENERAL\n")
    for i in 1:2*k
        ln = string("\tx",i,"\n")
        write(io, ln)
    end
        
    write(io, "END")
    close(io)
end