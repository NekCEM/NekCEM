<TeXmacs|1.0.6>

<style|generic>

<\body>
  Here are Maxwell's Equations in 2D, for the TM case:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>H<rsub|x>>|<cell|=>|<cell|-\<partial\><rsub|y>E<rsub|z>,>>|<row|<cell|\<partial\><rsub|t>H<rsub|y>>|<cell|=>|<cell|\<partial\><rsub|x>E<rsub|z>,>>|<row|<cell|\<partial\><rsub|t>E<rsub|z>>|<cell|=>|<cell|\<partial\><rsub|x>H<rsub|y>-\<partial\><rsub|y>H<rsub|x>.>>>>
  </eqnarray*>

  We make the separation-of-variables ansatz

  <\equation>
    <label|eq:sov-ansatz>E<rsub|z>(x,y,t)=X(x)Y(x)T(t).
  </equation>

  Then, observe

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t><rsup|2>E<rsub|z>>|<cell|=>|<cell|\<partial\><rsub|t>\<partial\><rsub|x>H<rsub|y>-\<partial\><rsub|t>\<partial\><rsub|y>H<rsub|x>>>|<row|<cell|\<Leftrightarrow\>\<partial\><rsub|t><rsup|2>E<rsub|z>>|<cell|=>|<cell|\<partial\><rsub|x>\<partial\><rsub|t>H<rsub|y>-\<partial\><rsub|y>\<partial\><rsub|t>H<rsub|x>>>|<row|<cell|\<Rightarrow\>\<partial\><rsub|t><rsup|2>E<rsub|z>>|<cell|=>|<cell|\<partial\><rsub|x><rsup|2>E<rsub|z>+\<partial\><rsub|y><rsup|2>E<rsub|z>.>>>>
  </eqnarray*>

  We obtain from the ansatz (<reference|eq:sov-ansatz>)

  <\equation>
    <label|eq:sov-waveeq><frac|T<rprime|''>(t)|T(t)>=<frac|X<rprime|''>(x)|X(x)>+<frac|Y<rprime|''>(y)|Y(y)>.
  </equation>

  Thus any linear combination of <with|mode|math|E<rsub|z>> composed using
  the functions

  <\eqnarray*>
    <tformat|<table|<row|<cell|X(x)>|<cell|=>|<cell|A*exp(i\<alpha\>*x),>>|<row|<cell|Y(y)>|<cell|=>|<cell|C*exp(i\<beta\>y).>>|<row|<cell|(<reference|eq:sov-waveeq>)\<Rightarrow\><frac|T<rprime|''>(t)|T(t)>>|<cell|=>|<cell|(i\<alpha\>)<rsup|2>+(i\<beta\>)<rsup|2>=-\<alpha\><rsup|2>-\<beta\><rsup|2>>>|<row|<cell|\<rightsquigarrow\><space|1em>\<omega\>>|<cell|\<assign\>>|<cell|<sqrt|\<alpha\><rsup|2>+\<beta\><rsup|2>>>>|<row|<cell|\<Rightarrow\>T(t)>|<cell|=>|<cell|*exp(i\<omega\>t)>>>>
  </eqnarray*>

  is a solution of the PDE, where <with|mode|math|A,C\<in\>\<bbb-C\>> are
  constant.

  We assume the boundary conditions to be periodic in <with|mode|math|x> and
  PEC in <with|mode|math|y>:

  <\eqnarray*>
    <tformat|<table|<row|<cell|E<rsub|z>,\<b-H\><left|(>x-<frac|1|2>,y,t<right|)>>|<cell|=>|<cell|E<rsub|z>,\<b-H\><left|(>x+<frac|1|2>,y,t<right|)>,>>|<row|<cell|E<rsub|z>,H<rsub|y><left|(>x,y\<pm\><frac|1|2>,t<right|)>>|<cell|=>|<cell|0.>>>>
  </eqnarray*>

  Now consider the PEC boundary condition on <with|mode|math|E<rsub|z>> at
  <with|mode|math|y=1/2>. First, note that the same condition at
  <with|mode|math|y=-1/2> is redundant for symmetry reasons. Next, no single
  function of the form <with|mode|math|exp(i\<beta\>y)> has a zero. Thus let
  us consider linear combinations of the form

  <\equation*>
    Y(y)=C<rsup|+>exp(i\<beta\>y)+C<rsup|->exp(-i\<beta\>y).
  </equation*>

  We obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|0>|<cell|=>|<cell|C<rsup|+>*exp<left|(><frac|i\<beta\>|2><right|)>+C<rsup|->exp<left|(><frac|-i\<beta\>|2><right|)>>>|<row|<cell|\<Leftrightarrow\><frac|C<rsup|+>|C<rsup|->>exp(i\<beta\>)>|<cell|=>|<cell|-1.>>>>
  </eqnarray*>

  Two easy solutions of this equation are
  <with|mode|math|C<rsup|+>=C,C<rsup|->=-C, \<beta\>=2\<pi\>m> and
  <with|mode|math|C<rsup|+>=C,C<rsup|->=C,\<beta\>=\<pi\>(2m+1)>.

  The periodic boundary condition on <with|mode|math|E<rsub|z>> is satisfied
  if <with|mode|math|\<alpha\>=2\<pi\>n>. Combining <with|mode|math|A> and
  <with|mode|math|C> into a single factor <with|mode|math|A>, we obtain

  <\equation*>
    E<rsub|z>(x,y,t)=\<Re\><left|{>A*exp(i(\<alpha\>*x+\<omega\>t))<left|[>C<rsup|+>exp(i\<beta\>y)+C<rsup|\<um\>>exp(-i\<beta\>y)<right|]><right|}>
  </equation*>

  as a general solution.

  Next,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>H<rsub|x>>|<cell|=>|<cell|-\<partial\><rsub|y>E<rsub|z>=-i\<beta\><left|[>A*exp(i(\<alpha\>x+\<omega\>t))<left|[>C<rsup|+>exp(i\<beta\>y)-C<rsup|\<um\>>exp(-i\<beta\>y)<right|]><right|]>>>|<row|<cell|\<Rightarrow\>H<rsub|x>>|<cell|=>|<cell|-<frac|i\<beta\>|i\<omega\>><left|[>A*exp(i(\<alpha\>x+\<omega\>t))<left|[>C<rsup|+>exp(i\<beta\>y)-C<rsup|\<um\>>exp(-i\<beta\>y)<right|]><right|]>,>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>H<rsub|y>>|<cell|=>|<cell|\<partial\><rsub|x>E<rsub|z>=\<alpha\>i*<left|[>A*exp(i(\<alpha\>x+\<omega\>t))<left|[>C<rsup|+>exp(i\<beta\>y)+C<rsup|\<um\>>exp(-i\<beta\>y)<right|]>>>|<row|<cell|\<Rightarrow\>H<rsub|y>>|<cell|=>|<cell|<frac|i\<alpha\>|i\<omega\>><left|[>A*exp(i(\<alpha\>x+\<omega\>t))<left|[>C<rsup|+>exp(i\<beta\>y)+C<rsup|\<um\>>exp(-i\<beta\>y)<right|]>.>>>>
  </eqnarray*>

  It is not hard to see that the <with|mode|math|H>-fields also satisfy their
  respective boundary conditions.
</body>

<\initial>
  <\collection>
    <associate|page-type|letter>
    <associate|sfactor|4>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|eq:sov-ansatz|<tuple|1|?>>
    <associate|eq:sov-waveeq|<tuple|2|?>>
  </collection>
</references>