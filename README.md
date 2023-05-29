# Some thoughts and experiments on Bergman's compact amalgamation problem

This page contains the source code and explanations related to the computational results presented in the paper:  
Michael Joswig, Mario Kummer, Andreas Thom, and Claudia He Yun  
*Some thoughts and experiments on Bergman's compact amalgamation problem*  
ARXIV: [https://arxiv.org/abs/2304.08365](https://arxiv.org/abs/2304.08365)  
CODE: [https://mathrepo.mis.mpg.de/CompactAmalgamation/index.html](https://mathrepo.mis.mpg.de/CompactAmalgamation/index.html)  

### Abstract
We study the question whether copies of $S^1$ in $\mathrm{SU}(3)$ can be amalgamated in a compact group. This is the simplest instance of a fundamental open problem in the theory of compact groups raised by George Bergman in 1987. Considerable computational experiments suggest that the answer is positive in this case. We obtain a positive answer for a relaxed problem using theoretical considerations.

## Amalgamating representations of the circle group
Representations of the circle group $S^1 = \{z \in \mathbb{C} : \|z\|=1\}$ in $\mathrm{SU}(3)$ are given by $(a,b,c)\in \mathbb{Z}^3$ with $a+b+c=0$ and $\mathrm{gcd}(a,b,c)=1$, up to conjugation. Given two such representations $\psi_{a,b,c}$ and $\psi_{a',b',c'}$, we look for injective representations $f,g: \mathrm{SU}(3) \to \mathrm{U}(n)$ such that $f \circ \psi_{a,b,c} = g \circ \psi_{a',b',c'}$. Unitary representations of $\mathrm{SU}(3)$ are semi-simple. The irreducible representations have characters Schur polynomials $s_\lambda(x_1,x_2,x_3)$, where $\lambda$ is either the partition $(1,1,1)$ or a partition with at most two parts. Therefore, our goal is to find positive linear combinations $P,Q$ of suitable Schur polynomials such that $P(z^a,z^b,z^c) = Q(z^{a'},z^{b'},z^{c'})$. By Proposition 3.1, we can accomplish this by solving a particular integer linear program $(\mathrm{ILP}_k)$ on page 5. We set up our linear program using the computer algebra system `OSCAR`. We solve this linear program in `OSCAR` and `SCIP`.  

In Table 2 Column 4, we only record the dimension of the amalgamated representations. For explicit descriptions of these representations, see [this page](Table2.md).  

## Setting up computations

The Jupyter notebook `amalgamation.ipynb` contains a walkthrough of Section 3. We give a detailed explanation of Example 3.3 and show how we produced Tables 1 and 2. 

- Download the Jupyter notebook [amalgamation.ipynb](amalgamation.ipynb) and [source code](amalgamation.jl).
- Download and install [Julia](https://julialang.org/).
- Install the Julia package [OSCAR](https://docs.oscar-system.org/stable/).
- Download and install [SCIP](https://www.scipopt.org/) and configure it as follows.

### Configure `SCIP`
Install `SCIP` and record the path of the binary file. In the file amalgamation.jl, find the following function and replace `/Users/yun/Documents/scipoptsuite-8.0.2/scip/bin/scip` by your path.

```
function scip(filename,outputname)
    run(`/Users/yun/Documents/scipoptsuite-8.0.2/scip/bin/scip -f $filename -l $outputname -q`);
end
```

## Checking Conjecture 4.9

Let $$F_{a,b} = (1+z+\dots+z^{b-1})^a\cdot (1+z^{-1}+\dots+z^{-(b-1)})^a.$$ Let $v=(v_1,v_2,v_3) \in \mathbb{Z}^3$ be a vector such that $v_1+v_2+v_3=0$ and $\mathrm{gcd}(v_1,v_2,v_3)=1$.  

Conjecture 4.9. There are natural numbers $a_0,b_0 > 0$ such that for all $a \geq a_0$, all $b$ divisible by $b_0$ there is $N\in \mathbb{N}$ and a Schur positive symmetric polynomial $P$ in three variables such that $$N\cdot F_{a,b}=P_v.$$

To verify this conjecture, we pick values for $a,b$ and $N$ and search for appropriate $P$. We can set up a similar integer linear program by writing $$P=\lambda_1S_1+\dots+\lambda_kS_k,$$ where $\{S_1,\dots,S_k\}$ are Schur polynomials in some ordering. Then the linear program is given by matching coefficients of $P_v$ with $F_{a,b}$ while requiring $\lambda_i \geq 0$.  

Note that once $a$ and $b$ are chosen, they define a finite set of Schur polynomials that can appear in $P$ in the following sense. Let $(r,t)$ be a partition. Let $v=(v_1,v_2,v_3)$ as before, assuming $v_1\geq v_2 > 0$. Then $(s_{r,t})_v(z)$ has leading deg $rv_1+tv_2$ and trailing deg $-rv_1+(t-r)v_2$. In the meantime, the polynomial $F_{a,b}$ has leading degree $(b-1)a$ and trailing degree $(1-b)a$. Therefore, all Schur polynomials $s_{r,t}$ that appear in $P$ must satisfy $$rv_1+tv_2 \leq (b-1)a,$$ $$-rv_1+(t-r)v_2 \geq (1-b)a.$$ There are only finitely many such partitions. In our program, we use this finite set of partitions and order them lexicographically.  

All computations are done in `Julia` and `OSCAR`. In Table 3 Column 4, only dimensions of the representations given by $P$ are recorded. For explicit descriptions, see [this page](Table3.md).  

To reproduce these computations, download [conjecture_section4.ipynb](conjecture_section4.ipynb) and [source code](conjecture_section4.jl).  


Project page created: 19/05/2023

Project contributors: Michael Joswig, Mario Kummer, Andreas Thom and Claudia Yun

Software used: Julia (version 1.8.5), OSCAR (version 0.12.0), SCIP (version 8.0.2)

System setup used: 
- Hydra at the MPI MiS: 4x16-core Intel Xeon E7-8867 v3 CPU (3300 MHz) on Debian
GNU/Linux 5.10.149-2 (2022-10-21) x86_64
- MacBook Pro: 8-core Intel i7-6700HQ CPU (2600 MHz) on Darwin Kernel 21.6.0 (2022-09-29) x86_64 

Corresponding author of this page: Claudia Yun, [clyun@mis.mpg.de](clyun@mis.mpg.de)



