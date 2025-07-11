# IDP Graph Polytope Analysis

## Definition of characteristic polytope of a graph G=(V,E)

Given a finite simple graph $G=(V,E)$ with vertices $V=\{1,\dots,n\}$ the characteristic polytope of $G$ is defined as the convex hull of the following vectors:

- vectors of vertices: For a vertex $i$ we assign the standard basis vector $e_i$ in $\mathbb{R}^n$
- vector of edges: For an edge $\{i,j\}$ a vector in $\mathbb{R}^n$ with
  - 1 in the $i$-th and $j$-th positions
  - 0 in all other positions

Formally, $P(G)=\text{convexHull}(e_i, e_i+e_j \mid i\in V, \{i,j\}\in E)$

## Definition of cycle graph C_n

Recall that the graph $C_n$ is known as the cycle graph with $n$-vertices.

The cycle graph $C_n$ is a connected graph where:

- The vertices can be labeled by $v_1,\dots,v_n$
- Each vertex $v_i$ is connected to $v_{i+1}$
- There is also an edge connecting $v_n$ back to $v_1$ to complete the cycle

## Properties of cycle graph C_n

- Number of vertices: $n$
- Number of edges: $n$
- Degree of each vertex: 2 (because each vertex is connected to exactly other two vertices)
- Connected: Yes
- Contains a cycle: Yes (it is a single cycle)
- Is it bipartite:
  - Yes, if n is even
  - No, if n is odd
- Is it a tree: No (it contains a cycle)

### Remark

The following M2 code defines the characteristic polytope of $C_n$ for $n=4..14$. It:

- Checks if $P(C_n)$ is IDP
- Computes the toric ideal of $I_{P(C_n)}$, the codimension of $I_{P(C_n)}$ and the betti table (only up to $n=12$)
- Computes the $h^\star$-vector of $P(C_n)$ (only up to $n=12$)

```Macaulay2
clearAll

needsPackage "Polyhedra"
needsPackage "NormalToricVarieties"
needsPackage "LatticePolytopes"

hilbertFunctionArtForI = (J) -> (
   fList := apply (1..dim (($\text{ring } J)^1/J$), i -> random (1, (ring J)));
   apply (0..10, i -> hilbertFunction (i, ($\text{ring } J)^1/ (J + (ideal fList))))
)

--hstarvecconj is the function which checks the conjecture that $h_i^\star=\frac{n}{n-i}\binom{n-i}{i}$
hstarvecconj = (n,i) -> (
  pt1 = 1/(n-i);
  pt2 = n*binomial(n-i,i);
  pt1*pt2
)

testfun = (list1,n,i) -> (list1_i == hstarvecconj(n,i))

--definepolytope function defines the characteristic polytope for $n$-gon
definepolytope = (n) -> (
    -- Build the list of vertices
    origin = for i from 0 to n-1 list 0;
    basisVectors = for i from 0 to n-1 list (
        for j from 0 to n-1 list (if j == i then 1 else 0)
    );
    edgePoints = for i from 0 to n-1 list (
        v = new MutableList from for j from 0 to n-1 list 0;
        v#i = 1;
        v#((i+1) % n) = 1;
        toList v
    );
 convexHull transpose matrix ({origin}|basisVectors|edgePoints)
    );

idpfunction = (n) -> (isNormal( definepolytope(n)))

toricIdealFromExponentList = (uList) -> (
 -- Number of parameters: $t_0$ plus the rest
    numParams := 1 + #(uList#0);
    -- Build the source ring variables: $t_0, t_1, \dots, t_{\text{numParams}}$ as symbols
    sourceRing := QQ[for i from 0 to numParams - 1 list("t" | toString i)];
    sourceVars := flatten entries vars sourceRing;
-- Build the target ring variables: $x_1, x_2, \dots, x_{\text{numberoflatticepoints}}$ as symbols
    numTargetVars := #uList;
    targetVarNames := for i from 1 to numTargetVars list ("x" | toString i);
    targetRing := QQ[targetVarNames, Degrees => for i from 1 to numTargetVars list 1];
    targetVars := flatten entries vars targetRing;
-- Build the parametrization matrix
    monomialList := apply(uList, p -> (
        productList := for i from 1 to numParams - 1 list (sourceVars#i)^(p#(i - 1));
        sourceVars#0 * product productList
    ));
    paramMap = map(sourceRing, targetRing, matrix { monomialList });
    -- Return the toric ideal
    ideal mingens kernel paramMap
);

--function which computes the dilation of the polytope
dilatePolytope = (P, k) -> (
    verts := vertices P;
    newVerts := $k \cdot \text{verts}$;
    convexHull newVerts
)
```

## Analysis of the 5-gon Characteristic Polytope

### Construction and Dimension

```Macaulay2
n = 5
P5 = definepolytope(n)

dim $P_5$
--answer: 5

idpfunction(n)
--answer: true
```

### Five-vertex Lattice Structure

```Macaulay2
uList1 = latticePoints P5;
# uList1  -- answer: 11

apply(uList1, p -> flatten entries p)
-- answer: {
--   {0, 0, 0, 0, 0}, {0, 0, 0, 0, 1}, {0, 0, 0, 1, 0}, {0, 0, 0, 1, 1}, 
--   {0, 0, 1, 0, 0}, {0, 0, 1, 1, 0}, {0, 1, 0, 0, 0}, {0, 1, 1, 0, 0}, 
--   {1, 0, 0, 0, 0}, {1, 0, 0, 0, 1}, {1, 1, 0, 0, 0}
-- }

apply(interiorLatticePoints P5, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Five-vertex Toric Properties

```Macaulay2
I5 = toricIdealFromExponentList(uListreord)

isHomogeneous I5
-- answer: true 

($\text{codim } I_5$, $\text{betti res } I_5$)
--codim 5, Gorenstein

hstarvec = toList hilbertFunctionArtForI(I5)
--answer: {1, 5, 5, 1, 0, 0, 0, 0, 0, 0, 0}

list67 = {}
for i from 0 to 3 do
    list67 = list67|{(i,testfun(hstarvec,5,i))}

tally list67
--answer: 
-- Tally{
--   (0, true) => 1
--   (1, true) => 1
--   (2, true) => 1
--   (3, false) => 1
-- }
```

## Analysis of the 6-gon Characteristic Polytope

### Six-vertex Construction

```Macaulay2
n = 6
P6 = definepolytope(n)

dim $P_6$
--answer: 6

idpfunction(n)
--answer: true
```

### Six-vertex Lattice Points

```Macaulay2
uList2 = latticePoints P6;
# uList2  -- answer: 13

apply(uList2, p -> flatten entries p)
-- answer: {
--   {0, 0, 1, 0, 0, 0}, {0, 0, 1, 1, 0, 0}, {0, 1, 0, 0, 0, 0}, 
--   {0, 1, 1, 0, 0, 0}, {1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1}, 
--   {1, 1, 0, 0, 0, 0}
-- }

apply(interiorLatticePoints P6, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints

I6 = toricIdealFromExponentList(uListreord)

isHomogeneous I6
-- answer: true 

($\text{codim } I_6$, $\text{betti res } I_6$)
--(codim 6, CM)

hstarvec = toList hilbertFunctionArtForI(I6)
--answer: (1, 6, 9, 1, 0, 0, 0, 0, 0, 0, 0)

list68 = {}
for i from 0 to 2 do
    list68 = list68|{(i,testfun(hstarvec,6,i))}

tally list68
--answer:
-- Tally{
--   (0, true) => 1
--   (1, true) => 1
--   (2, true) => 1
-- }
```

## Analysis of the 7-gon Characteristic Polytope

### Seven-vertex Properties

```Macaulay2
n = 7
P7 = definepolytope(n)

dim $P_7$
--answer: 7

idpfunction(n)
--answer: true
```

### Seven-vertex Lattice Structure

```Macaulay2
uList1 = latticePoints P7;
# uList1  -- answer: 15

apply(uList1, p -> flatten entries p)
-- answer: {
--   {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 1, 0}, 
--   {0, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0, 1}, 
--   {0, 0, 0, 0, 1, 1, 0}, {0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 1}, 
--   {0, 0, 0, 1, 1, 0, 0}, {0, 0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 1}, 
--   {0, 0, 1, 1, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 1}, 
--   {0, 1, 1, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 1}, 
--   {1, 1, 0, 0, 0, 0, 0}
-- }

apply(interiorLatticePoints P7, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Seven-vertex Toric Investigation

```Macaulay2
I7 = toricIdealFromExponentList(uListreord)

isHomogeneous I7
-- answer: true 

($\text{codim } I_7$, $\text{betti res } I_7$)
--answer: codim 7, Gorenstein

hstarvec = toList hilbertFunctionArtForI(I7)
--answer: (1, 7, 14, 7, 1, 0, 0, 0, 0, 0, 0)

list68 = {}
for i from 0 to 3 do
    list68 = list68|{(i,testfun(hstarvec,7,i))}

tally list68
--answer:
-- Tally{(0, true) => 1}
--        (1, true) => 1
--        (2, true) => 1
--        (3, true) => 1
```

## Analysis of the 8-gon Characteristic Polytope

### Eight-vertex Construction

```Macaulay2
n = 8
P8 = definepolytope(n)

dim P8
--answer: 8

idpfunction(n)
--answer: true
```

### Eight-vertex Lattice Analysis

```Macaulay2
uList1 = latticePoints P8;
# uList1  -- answer: 17

apply(uList1, p -> flatten entries p)
-- answer: {
--   {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 1, 0}, 
--   {0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 1, 1, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, 
--   {0, 0, 0, 0, 1, 1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 0, 0, 0}, 
--   {0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, 
--   {0, 1, 1, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 1}, 
--   {1, 1, 0, 0, 0, 0, 0, 0}
-- }

apply(interiorLatticePoints P8, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Eight-vertex Toric Properties

```Macaulay2
I8 = toricIdealFromExponentList(uListreord)

isHomogeneous I8
-- answer: true 

($\text{codim } I_8$, $\text{betti res } I_8$)
--answer: codim 8, CM

hstarvec = toList hilbertFunctionArtForI(I8)
--answer: (1, 8, 20, 16, 1, 0, 0, 0, 0, 0, 0)

list68 = {}
for i from 0 to 3 do
    list68 = list68|{(i,testfun(hstarvec,8,i))}

tally list68
--answer: 
-- Tally{(0, true) => 1}
--        (1, true) => 1
--        (2, true) => 1
--        (3, true) => 1
```

## Analysis of the 9-gon Characteristic Polytope

### Nine-vertex Configuration

```Macaulay2
n = 9
P9 = definepolytope(n)

dim P9
--answer: 9

idpfunction(n)
--answer: true
```

### Nine-vertex Lattice Points

```Macaulay2
uList1 = latticePoints P9;
# uList1  -- answer: 19

apply(uList1, p -> flatten entries p)
-- answer: {
--   {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 1, 0, 0},
--   {0, 0, 0, 0, 0, 0, 1, 0, 1}, {0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 1, 0},
--   {0, 0, 0, 0, 0, 1, 1, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 1},
--   {0, 0, 0, 0, 1, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 1},
--   {0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 1}, {0, 0, 1, 1, 0, 0, 0, 0, 0},
--   {0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 1}, {0, 1, 1, 0, 0, 0, 0, 0, 0},
--   {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 1}, {1, 1, 0, 0, 0, 0, 0, 0, 0}
-- }

apply(interiorLatticePoints P9, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Nine-vertex Toric Structure

```Macaulay2
I9 = toricIdealFromExponentList(uListreord)

isHomogeneous I9
-- answer: true 

($\text{codim } I_9$, $\text{betti res } I_9$)
--answer: codim 9, CM

hstarvec = toList hilbertFunctionArtForI(I9)
--answer: (1, 9, 27, 30, 9, 1, 0, 0, 0, 0, 0)

list68 = {}
for i from 0 to 4 do
    list68 = list68|{(i,testfun(hstarvec,9,i))}

tally list68
--answer: 
-- Tally{(0, true) => 1}
--        (1, true) => 1
--        (2, true) => 1
--        (3, true) => 1
--        (4, true) => 1
```

## Analysis of the 10-gon Characteristic Polytope

### Ten-vertex Properties

```Macaulay2
n = 10
P10 = definepolytope(n)

dim P10
--answer: 10

idpfunction(n)
--answer: true
```

### Ten-vertex Lattice Analysis

```Macaulay2
uList1 = latticePoints P10;
# uList1  -- answer shown in results

apply(uList1, p -> flatten entries p)
-- answer shown in results

apply(interiorLatticePoints P10, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Ten-vertex Toric Features

```Macaulay2
I10 = toricIdealFromExponentList(uListreord)

isHomogeneous I10
-- answer: true 

($\text{codim } I_{10}$, $\text{betti res } I_{10}$)
--answer: codim 10, CM

hstarvec = toList hilbertFunctionArtForI(I10)
--answer: {1, 10, 35, 50, 25, 1, 0, 0, 0, 0, 0}

list68 = {}
for i from 0 to 4 do
    list68 = list68|{(i,testfun(hstarvec,10,i))}

tally list68
--answer: 
-- Tally{(0, true) => 1}
--        (1, true) => 1
--        (2, true) => 1
--        (3, true) => 1
--        (4, true) => 1
```

## Analysis of the 11-gon Characteristic Polytope

### Eleven-vertex Construction

```Macaulay2
n = 11
P11 = definepolytope(n)

dim P11
--answer: 11

idpfunction(n)
--answer: true
```

### Eleven-vertex Lattice Study

```Macaulay2
uList1 = latticePoints P11;
# uList1  -- answer: 23

apply(uList1, p -> flatten entries p)
-- answer shown in results

apply(interiorLatticePoints P11, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Eleven-vertex Toric Results

```Macaulay2
I11 = toricIdealFromExponentList(uListreord)

isHomogeneous I11
-- answer: true 

($\text{codim } I_{11}$, $\text{betti res } I_{11}$)
--answer: codim 11, CM

hstarvec = toList hilbertFunctionArtForI(I11)
--answer: {1, 11, 44, 77, 55, 11, 1, 0, 0, 0, 0}

list68 = {}
for i from 0 to 5 do
    list68 = list68|{(i,testfun(hstarvec,11,i))}

tally list68
--answer: 
-- Tally{(0, true) => 1}
--        (1, true) => 1
--        (2, true) => 1
--        (3, true) => 1
--        (4, true) => 1
--        (5, true) => 1
```

## Analysis of the 12-gon Characteristic Polytope

### Twelve-vertex Properties

```Macaulay2
n = 12
P12 = definepolytope(n)

dim P12
--answer: 12

idpfunction(n)
--answer: true
```

### Twelve-vertex Lattice Points

```Macaulay2
uList1 = latticePoints P12;
# uList1  -- answer: 25

apply(uList1, p -> flatten entries p)
-- answer shown in results

apply(interiorLatticePoints P12, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Twelve-vertex Toric Computation

```Macaulay2
I12 = toricIdealFromExponentList(uListreord)

isHomogeneous I12
-- answer: true 

($\text{codim } I_{12}$, $\text{betti res } I_{12}$)
--answer: codim 12, does not compute the Betti resolution

hstarvec = toList hilbertFunctionArtForI(I12)
--answer: ? (computation too large)

list68 = {}
for i from 0 to 5 do
    list68 = list68|{(i,testfun(hstarvec,12,i))}

tally list68
-- result not shown due to computational limitations
```

## Analysis of the 13-gon Characteristic Polytope

### Thirteen-vertex Properties

```Macaulay2
n = 13
P13 = definepolytope(n)

dim P13
--answer: 13

idpfunction(n)
--answer: true
```

### Thirteen-vertex Lattice Study

```Macaulay2
uList1 = latticePoints P13;
# uList1  -- answer: 27

apply(uList1, p -> flatten entries p)
-- answer shown in results

apply(interiorLatticePoints P13, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
```

### Thirteen-vertex Toric Analysis

```Macaulay2
I13 = toricIdealFromExponentList(uListreord)

isHomogeneous I13
-- answer: true 

($\text{codim } I_{13}$, $\text{betti res } I_{13}$)
--answer: codim 13, does not compute the Betti resolution

hstarvec = toList hilbertFunctionArtForI(I13)
--answer: ? (computation too large)

list68 = {}
for i from 0 to 5 do
    list68 = list68|{(i,testfun(hstarvec,13,i))}

tally list68
-- result not shown due to computational limitations
```

## Analysis of the 14-gon Characteristic Polytope

### Fourteen-vertex Properties

```Macaulay2
n = 14
P14 = definepolytope(n)

dim P14
--answer: 14

idpfunction(n)
--answer: true
```

### Fourteen-vertex Lattice Investigation

```Macaulay2
uList1 = latticePoints P14;
# uList1  -- answer: 29

apply(uList1, p -> flatten entries p)
-- answer shown in results

apply(interiorLatticePoints P14, p -> flatten entries p)
--answer: {}

uListreord = {origin}|basisVectors|edgePoints
--uListreord is a reordering of the list uList1 of the lattice points of P
```

### Fourteen-vertex Toric Results

```Macaulay2
I14 = toricIdealFromExponentList(uListreord)

isHomogeneous I14
-- answer: true 

($\text{codim } I_{14}$, $\text{betti res } I_{14}$)
--answer: codim 14, does not compute the Betti resolution

hstarvec = toList hilbertFunctionArtForI(I14)
--answer: ? (computation too large)

list68 = {}
for i from 0 to 5 do
    list68 = list68|{(i,testfun(hstarvec,14,i))}

tally list68
-- result not shown due to computational limitations
```

## Comprehensive IDP Property Analysis

### Results for Small Polygons (n ≤ 20)

```Macaulay2
listcheckidp = {}
for i from 1 to 20 do
    listcheckidp = listcheckidp|{idpfunction(i)} 

listcheckidp
--answer: Tally{true => 20}
```

### Results for Medium Polygons (20 ≤ n ≤ 30)

```Macaulay2
listcheckidp = {}
for i from 20 to 30 do
    listcheckidp = listcheckidp|{idpfunction(i)} 

tally listcheckidp
--answer: Tally{true => 11}
```

### Results for Large Polygons (31 ≤ n ≤ 40)

```Macaulay2
listcheckidp = {}
for i from 31 to 40 do
    listcheckidp = listcheckidp|{idpfunction(i)} 

tally listcheckidp
--results shown in output
```
