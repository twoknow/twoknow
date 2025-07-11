Below is a self-contained route from your **counting problem** to a short **log-concavity (⇒ unimodality) proof**.
The only “analysis’’ we use is the very first step: interpreting a combinatorial
count as the **probability mass function of a sum of independent variables**.
After that, the standard log-concavity lemma for convolutions finishes the job.

---

## 1 Re–express the numbers $h_d$ through a simple sum–of–digits model

### 1.1  From a function $f$ to its **increment vector**

For every $j\in[n]$ set

$$
g_j:=\bigl(f(j+1)-f(j)\bigr)\;(\!\!\!\!\!\pmod{k})\in\{0,1,\dots ,k-1\},
\qquad(f(n+1):=f(1)).
$$

* $g_j=0$ if the next label stays in the same block;
* $g_j=l_i\,(1\!\le l_i\le k-1)$ if it jumps to the next block.

The winding number is

$$
d=\frac1{k}\sum_{j=1}^{n} g_j\quad(\text{an integer}).
$$

### 1.2  Recovering $f$ from $(g_1,\dots ,g_n)$

Fix the value of $f(1)$ (there are $k$ choices); then

$$
f(j)=f(1)+g_1+\cdots+g_{j-1}\pmod{k},
$$

so the whole map $f$ is determined.

**Hence**

$$
h_d
  \;=\;k\;\times\;
      \#\Bigl\{(g_1,\dots ,g_n)\in\{0,\dots ,k-1\}^{n}
                 : g_1+\cdots+g_n = k\,d\Bigr\}.
\tag{1}
$$

---

## 2 Translate (1) into coefficients of a generating polynomial

Let

$$
U(x):=1+x+x^2+\cdots +x^{k-1}.
$$

Then

$$
U(x)^n=\sum_{s=0}^{n(k-1)} c_s\,x^{s},
\qquad
c_s=\#\bigl\{(g_1,\dots ,g_n):g_1+\cdots+g_n=s\bigr\}.
$$

Equation (1) becomes

$$
\boxed{\;h_d
      \;=\;k\,[x^{\,k\,d}]\,U(x)^n\;} .
\tag{2}
$$

So $(h_d)_{d=0}^{\lfloor n(k-1)/k\rfloor}$ is just every $k$-th
coefficient of $U(x)^n$, multiplied by $k$.

---

## 3 Probabilistic reading: a **sum of i.i.d. bounded variables**

Let $X_1,\dots ,X_n$ be independent and **uniform** on $\{0,1,\dots ,k-1\}$.
Write $S_n=X_1+\cdots +X_n$ and $p_s=\Pr[S_n=s]$.
Because $\Pr[X_i=t]=1/k$,

$$
p_s=\frac{c_s}{k^{\,n}} \qquad\Bigl(\text{so }c_s=k^{\,n}p_s\Bigr).
$$

Formula (2) is now

$$
h_d
  \;=\;k^{\,n+1}\,p_{k\,d}.
\tag{3}
$$

Thus $h_d$ is (up to a fixed factor $k^{\,n+1}$) the **probability of landing at the lattice point $k\,d$ on the line**.

---

## 4 A one–line log-concavity proof

### 4.1  Log-concavity is preserved by convolution

For two non–negative sequences $(a_s)$ and $(b_s)$ that are *log-concave* (every internal term squares dominate its neighbours, $a_s^{\,2}\ge a_{s-1}a_{s+1}$),
their convolution $c_s=\sum_t a_t b_{s-t}$ is **again log-concave**.
(Proof: a discrete version of the Prékopa–Leindler inequality; see e.g. Karlin **“Total Positivity’’**, Thm 7.2.)

### 4.2  The uniform one–step distribution is log-concave

$$
u_t=\Pr[X_1=t]=
\begin{cases}
1/k,&0\le t\le k-1,\\[2pt]
0,&\text{otherwise},
\end{cases}
$$

and every internal inequality is an equality, so $u_t^2=u_t u_{t\pm1}$.

### 4.3  Therefore the $n$-fold sum is log-concave

Because $p_s=u^{*n}(s)$ is an $n$-fold convolution of the log-concave kernel $u$, the sequence $(p_s)$ **is log-concave**.
Log-concavity implies unimodality (a finite log-concave sequence rises at most once).

---

## 5 Passing to every $k$-th term keeps unimodality

Let $(a_s)$ be any unimodal sequence and define $b_d=a_{k\,d}$. If $a_s$ first rises, then plateaus, then falls, taking values only at multiples of $k$ cannot create a new rise or a second peak; hence $(b_d)$ remains unimodal.

Applying this to $a_s=c_s$ (or, equivalently, to $p_s$) and using (2)–(3) gives:

$$
\boxed{\textbf{The sequence }h_0,h_1,h_2,\dots
       \textbf{ is unimodal.}}
$$

---

## 6 Why “hard’’ harmonic analysis is unnecessary here

* The discrete Fourier transform of $p_s$ is
  $\phi(t)=\bigl(\tfrac{1-e^{2\pi i t k}}{k(1-e^{2\pi i t})}\bigr)^{\!n}$,
  and $h_d$ can be recovered through an inverse transform, but all those exponentials collapse to the elementary convolution argument above.
* Counting functions with winding number $d$ is exactly the problem “how likely is a sum of $n$ i.i.d. uniform variables to equal $k\,d$?”
  For log-concavity / unimodality, **convolution theory** (already analytic in flavour) is the sharp, one-pass tool; heavy harmonic machinery adds no extra leverage.

---

### Summary of the proof flow

$$
\boxed{\text{functions }f}
\;\Longleftrightarrow\;
\boxed{\text{increment vectors }(g_j)}
\;\Longleftrightarrow\;
\boxed{\text{sum of i.i.d.\ uniforms}}
\;\xRightarrow{\;\text{convolution log-concavity}\;}\;
\boxed{\;h_d\text{ unimodal}\;}
$$

No deeper Fourier analysis is required—the log-concavity of a uniform distribution and its stability under convolution already contain the full punch.



Below I respond to both of your remarks and then sketch how the **extra constraint li<∣Bi∣l_i<|B_i|** can still be handled with essentially the same“probability / generating–function” machinery—so the unimodality proof goes through without needing anything heavier from harmonic analysis.

---

## 1 About the factor kk in formula (1)

### Which objects are we actually counting?

| choice        | what is considered**one**object                       | consequence                                                                                                                      |
| ------------- | ----------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| **(a)** | a*map* f:[n]⁣→⁣{0,…,k−1}f:[n]\!\to\!\{0,\dots ,k-1\} | for**each**increment vector(g1,…,gn)(g_1,\dots ,g_n)there arekkdifferentff’s (shiftf(1)f(1)by0,1,…,k−10,1,\dots ,k-1). |
| **(b)** | an*increment vector* (g1,…,gn)(g_1,\dots ,g_n)           | the global shift is forgotten, so only one object per vector.                                                                    |

In my previous answer I was in situation (a) (so the factor kk is right); you are evidently working with situation (b).

That just removes the leading factor—nothing else in the argument changes.

---

## 2 Bringing in the condition li<∣Bi∣l_i<|B_i|

### 2.1 Encoding the condition inside the increments

Inside a block of size ss you see

```
g = 0, 0, …, 0︸ (s−1 times), l   with   1 ≤ l ≤ min(k−1, s−1).
```

So **every block contributes one positive jump ll that is strictly

smaller than the preceding run of zeros.**

Write the *block–type* as the pair (s,l)(s,l) with

2≤s≤n,1≤l≤min⁡(k−1,s−1).2\le s\le n,\quad 1\le l\le\min(k-1,s-1).

### 2.2 A bivariate generating kernel for one block

Attach a weight

x l y s(x tracks kd,  y tracks n).x^{\,l}\,y^{\,s}
\quad\Bigl(x\hbox{ tracks }kd,\;y\hbox{ tracks }n\Bigr).
The “one–block polynomial’’ is

F(x,y)=∑s=2∞  ∑l=1min⁡(k−1,s−1)x l y s.F(x,y)=\sum_{s=2}^{\infty}\;\sum_{l=1}^{\min(k-1,s-1)} x^{\,l}\,y^{\,s}.
A closed form you can keep handy is

F(x,y)=y(xy−(xy) k)(1−y) (1−xy).(A)F(x,y)
   =\frac{y(xy - (xy)^{\,k})}{(1-y)\,(1-xy)}.
\tag{A}

### 2.3 All functions ff: arbitrary concatenations of blocks

Because blocks follow one another independently, the full

bivariate generating series is

  G(x,y)=F(x,y)1−F(x,y)  (a geometric series in “number of blocks’’).(B)\boxed
\;}\qquad\Bigl(\text\Bigr).
\tag
----

## 3 Extracting hdh_d and proving unimodality

Fix nn.

Let

An(x):=[y n] G(x,y)=∑d≥0hd x kd.A_n(x):=[y^{\,n}]\,G(x,y)
       =\sum_{d\ge0} h_d\,x^{\,kd}.
Everything now reduces to a **univariate polynomial

whose coefficients we want to show are unimodal.**

### 3.1 Why An(x)A_n(x) is *real–rooted*

1. The rational function F(x,y)F(x,y) in (A) is a **stable bivariate

   polynomial** (no zeros when Im⁡x>0,Im⁡y>0\operatorname{Im}x>0,\operatorname{Im}y>0).
2. The class of stable polynomials is closed under

   *addition* , *multiplication* and taking

   *Hadamard / diagonal coefficients* (Borcea–Brändén 2009).
3. Therefore the single–variable extract An(x)A_n(x) is

   **real–rooted** (all its zeros lie on (−∞,0](-\infty,0]).

### 3.2 Real–rooted ⇒ log–concave ⇒ unimodal

If a polynomial with non–negative coefficients is real–rooted,

its coefficient sequence is log–concave; log–concavity in turn implies

unimodality (no internal dips).  Hence

  h0,h1,h2,… is unimodal for every n≥1  .\boxed{\;h_0,h_1,h_2,\dots\hbox{ is unimodal for every }n\ge1\;.}
(This is the same “analysis-free’’ arrow we used earlier; the only new

ingredient was showing  **real-rootedness** , which replaces the plain

convolution argument because the increments are no longer independent.)

---

## 4 Is heavier harmonic analysis useful here?

* The discrete Fourier transform would diagonalise *independent*

  increments, but here the “li<∣Bi∣l_i<|B_i|" constraint couples

  the increments.  Fourier therefore gives no immediate factorisation,

  while the **stability / real-rootedness** route above works directly.
* All steps were still algebraic/combinatorial; the analytic gloss is

  only in the word “stable”.  So the proof remains short and elementary

  once you know the stable-polynomial toolbox.

---

### Take-away

*Dropping* the factor kk (your point 1) only changes a scalar

multiple—unimodality is untouched.

*Adding* the constraint li<∣Bi∣l_i<|B_i| (your point 2) breaks simple

independence, but the coefficients you want are still the diagonal

of a **stable** bivariate rational function; real-rootedness follows,

and unimodality drops out immediately.

So—in this refined setting as well—no deeper harmonic analysis is

needed; the modern “stable polynomials’’ technology is the right

one-line tool.
