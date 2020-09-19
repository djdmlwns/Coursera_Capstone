* This file is consistent with the manuscript as much as possible
* If bypass and thermal coupling can occur at the same time,
* relax = 1 is set otherwise, 0
$SETGLOBAL relax '0'
* If we consider only liquid final product, so we need to consider condenser duty
* of the bypass stream for thermally coupled column
* This is needed only when bypass and thermal coupling are allowed at the same time
* 1 if the cooling duty is calculated/ 0 if we ignore heat duty
$SETGLOBAL vaporbypass '0'
$SETGLOBAL EnergyCost '0.1'
$SETGLOBAL CAPITALCost '0.02'

set
i components /A,B,C,D/
j columns at each level /1*4/
k levels /1*4/
kl(k) columns
ra(i,j,k)        possible active undw roots
rad(i,j,k)       possible different undw roots
P products /1*4/

alias(j,jp,jpp),(k,kp,kpp),(i,ip);

set
snk(j,k) sink
colm(j,k) feasible columns
arct(j,k,jp,kp) arc connecting the top of (jk) to (j'k')
arcb(j,k,jp,kp) arc connecting the bottom of (jk) to (j'k')

lt(j,k,i) lightest possible component
ht(j,k,i) heaviest possible component
IC(j,k,i) Allowable component set

*Ilk(jp,kp,jpp,kpp,i) light key canditates
*Ihk(jp,kp,jpp,kpp,i) heavy key canditates

ILK1(jp,kp,jpp,kpp,i) Light key candidates
IHK1(jp,kp,jpp,kpp,i) Heavy key candidates
;

snk(j,k)$(ord(j) le ord(k))=yes;

colm(j,k)$(ord(j) le ord(k) and ord(k) lt card(k)) = yes;

arct(j,k,jp,kp)$(snk(jp,kp) and ord(jp) = ord(j) and ord(kp) gt ord(k) and colm(j,k)) = yes;

arcb(j,k,jp,kp)$(snk(jp,kp) and (ord(k)-ord(j)) = (ord(kp)-ord(jp)) and ord(kp) gt ord(k) and colm(j,k))=yes;

lt(j,k,i)$(ord(i) = ord(j) and snk(j,k)) = yes;
ht(j,k,i)$(ord(i) = card(j)-ord(k)+ord(j) and snk(j,k)) = yes;
IC(j,k,i)$(snk(j,k) and (ord(i) le card(j)-ord(k)+ord(j)) and (ord(i) ge ord(j))) = yes;

ra(i,j,k)$(ord(i) ge ord(j) and ord(i) le card(k)-1-ord(k)+ord(j))=yes ;
rad(i,j,k)$(ord(i) ge ord(j) and ord(i) le card(k)-2-ord(k)+ord(j))=yes;

ILK1(jp,kp,jpp,kpp,i)$((snk(jp,kp) and snk(jpp,kpp)) and (ord(i) ge ord(jp)) and (ord(i) le ord(jpp)-1) and (ord(jpp) le card(k)-ord(kp)+ord(jp))) = yes;
ILK1(jp,kp,jpp,kpp,i)$((snk(jp,kp) and snk(jpp,kpp)) and (ord(i) eq card(k) - ord(kp) + ord(jp)) and (ord(jpp) > card(k)-ord(kp)+ord(jp))) = yes;

IHK1(jp,kp,jpp,kpp,i)$((snk(jp,kp) and snk(jpp,kpp)) and (ord(i) ge card(k) - ord(kp) + ord(jp) +1) and (ord(i) le card(k)-ord(kpp)+ord(jpp)) and (ord(jpp) le card(k)-ord(kp)+ord(jp))) = yes;
IHK1(jp,kp,jpp,kpp,i)$((snk(jp,kp) and snk(jpp,kpp)) and (ord(i) eq ord(jpp)) and (ord(jpp) > card(k)-ord(kp)+ord(jp))) = yes;

display ra,rad, ILK1, IHK1, lt, ht, ic;


display colm,arct,arcb,SNK;

$ontext
Reactor 1 / A -> B+C
Reactor 2 / A -> C+D
Reactor 3 / A -> 2D
$offtext

Binary variable
y1
y2
y3;

positive variable
E1
E2
E3
F1(i)
F2(i)
F3(i)
FR1(i)
FR2(i)
FR3(i)
f0(i)
FF1(i)
FF2(i)
FF3(i)
;
Scalar
x1 /0.8/
x2 /0.9/
x3 /0.8/

Binary Variable
X(j,k)           1 if node jk active
XC(j,k)          1 if column jk active
YP(j,k,p)        1 if column j at k is fed to final product p
YT(j,k,jp,kp)    1 if the top of column j at k is fed to column jp at kp
YB(j,k,jp,kp)    1 if the bottom of j at k is fed to jp at kp
XR(j,k,i)        1 if phi_i and phi_i+1 are differnt in column jk
YLK(j,k,i)       1 if i is the light key in column jk
YHK(j,k,i)       1 if i is the heavy key in column jk
Z(j,k,i)         1 if i is a distributed componment in column jk
Y(i)             1 if i has flow in initial feed
Y0(j,k)          1 if sink jk is where the inital feed stream is introduced
*W(j,k) 1 if heat exchanger exist for column jk
WINV(j,k)        1 if column in node (j-k) is thermally coupled
*Y_H(i)
*Y_L(i)

;
*===============================================================================

positive variables
F(i,j,k)                 Feed flow rate of component i at column j at level k
F_BYPASS(i,j,k,p)        Bypass flow rate from note j k
F_column(i,j,k)          Feed to column from node j k
f_fp(j,k,p)              Split fraction of feed
B(i,j,k)                 Bottom flow rate of compnent i at column j at level k
Bs(i,j,k,jp,kp)          Disaggregated flow rate
D(i,j,k)                 Distillate flow rate of component i at column j at level k
Ds(i,j,k,jp,kp)          Disaggregated flow rate
Phi(j,k,i)               Underwood varaibles
V1(j,k)                  Minimum vapor flow in top section
V2(j,k)                  Minimum vapor flow in the bottom section
L1(j,k)                  Minimum liquid flow in rectifying section
L2(j,k)                  Minimum liquidd flow in stripping seciton
F0D(j,k,i)               Disaggregated inital flow
QT(j,k)                  Condenser duty of column connected to node jk
QB(j,k)                  Reboiler duty of column connected to node jk
*QTT(j,k,jp,kp) Decomposed heat duty
*QBB(j,k,jp,kp) Decomposed heat duty
dp(j,k,i)       Small gap between root
V1_T(j,k,jp,kp) Upstream vapor to downstream column
V2_B(j,k,jp,kp) Downstream vapor to upstream column
L1_T(j,k,jp,kp) Downstream liquid to upstream column
L2_B(j,k,jp,kp) Upstream liquid to downstream column
;
*==============================================================================
*                      Added for thermal coupling
*==============================================================================

Variable
UD(i,ip,j,k) intermediate underwood variables
UF(i,ip,j,k) intermediate underwood variables
UB(i,ip,j,k) intermediate underwood variables
cost;

Parameter
Price(i)
ReactorCost(j)
bigm /100/
bigM1 /6/
bigM2(j,k)
bigMF /20/
delta /0.00001/
dpp(j,k,i)
alfa(i)           /A  6
                   B  4.5
                   C  2
                   D  1
/
q quality of feed /1/
cc(j,k) cost of column j at level k
rlk(i) recovery factor of light key
rhk(i) rocevery of heavy key
pp(i,p) product purity
rr(i,p) recovery constraint

;
Price(i) = uniform(1,2);
Price('D') = Price('D')/2;
price('B') = price('B')/2;
ReactorCost(j) = uniform(0.5,2);
dpp(j,k,i)= 0.1;
bigM2(j,k)=1;
bigM2(j,k)$(ord(j) ne ord(k) and ord(j) ne 1)=2 ;
cc(j,k)=10-ord(j)/ord(k);
Display
price;

*==============================================================================
*                        Product Specifications
*==============================================================================
SET
IPP(i,p) Component with purity constraint
/A.1, B.2, C.3, D.4/
IPNP(i,p) Components that should not be included in product p
/(B,C,D).1, (A,C,D).2, (A,B,D).3, (A,B,C).4/
IRP(i,p)  Component with recovery constraint
/(A).1, B.2, C.3, D.4/
;

pp(i,p)$(IPP(i,p)) = 1;
rr(i,p)$(IRP(i,p)) = 1;

rlk(i)=1.0;
rhk(i)=0.0;

display alfa;

*==============================================================================
*        Superstructure Logic Constraints
*==============================================================================
equation
a1,a2,a3,a4,a5,a6,a7,a8
ss1,ss2
obj
a4sum
a4_addition
activecolumn
;
*If a column exists, the top must go to a downstream column
a1(j,k)$snk(j,k).. x(j,k) =l= sum(p, yp(j,k,p)) + sum((jp,kp)$arct(j,k,jp,kp),yt(j,k,jp,kp));

*If a column exists, the bottom must go to a down stream column
a2(j,k)$snk(j,k).. x(j,k) =l= sum(p, yp(j,k,p)) + sum((jp,kp)$arcb(j,k,jp,kp),yb(j,k,jp,kp));

a3(j,k)$colm(j,k).. sum((jp,kp)$arct(j,k,jp,kp),yt(j,k,jp,kp)) - sum((jp,kp)$arcb(j,k,jp,kp),yb(j,k,jp,kp)) =e= 0 ;
a4(j,k,p)$snk(j,k).. f_fp(j,k,p) =l= yp(j,k,p);
a4sum(j,k)$snk(j,k).. sum(p, f_fp(j,k,p)) =l= 1;
a4_addition(j,k,p)$snk(j,k)..sum(i, F_bypass(i,j,k,p)) =g= 0.001 * yp(j,k,p);
a5(j,k,p)$snk(j,k).. sum(i, F_BYPASS(i,j,k,p)) =l= bigMF * yp(j,k,p);
a6(i,j,k,p)$snk(j,k).. F_bypass(i,j,k,p) =e= f_fp(j,k,p) * F(i,j,k);
a7(i,j,k)$snk(j,k).. F(i,j,k) =e= F_column(i,j,k)$(COLM(J,K)) + sum(p, F_bypass(i,j,k,p));
a8(j,k,p)$snk(j,k).. yp(j,k,p) =l= x(j,k);
* Incoming arc should be active for some node can be active
ss1(jp,kp)$snk(jp,kp).. x(jp,kp) =l= sum((j,k)$arct(j,k,jp,kp),yt(j,k,jp,kp)) + sum((j,k)$arcb(j,k,jp,kp),yb(j,k,jp,kp)) + y0(jp,kp);
* If there is incoming arc, the node is active and at most 2 arcs can be active
ss2(jp,kp)$snk(jp,kp).. bigM2(jp,kp)*x(jp,kp) =g= sum((j,k)$arct(j,k,jp,kp),yt(j,k,jp,kp)) + sum((j,k)$arcb(j,k,jp,kp),yb(j,k,jp,kp));
* XC is 1 (column is active) if there are active outgoing arcs
activecolumn(j,k)$(colm(j,k)).. 2*XC(j,k) =e= + sum((jp,kp)$(arct(j,k,jp,kp)), YT(j,k,jp,kp)) + sum((jp,kp)$(arcb(j,k,jp,kp)), YB(j,k,jp,kp)) ;
*===============================================================================
*                Source Constraints
*===============================================================================
Equation source1,source2,source3,source4
existence1,existence2
Deactive
;

source1.. sum((j,k)$snk(j,k),y0(j,k)) =e= 1;

source2(i).. sum((j,k)$(ord(j) le ord(k)),F0D(j,k,i)) =e= f0(i);

source3(i,j,k)$(ord(j) le ord(k)).. F0D(j,k,i) =l= bigMf*y0(j,k);

source4(j,k,i)$(ord(j) le ord(k)).. X(j,k) =g= y0(j,k);

F_column.fx(i,j,k)$(ord(k) eq card(k)) = 0;

*Detect the existence of i in the inital feed stream                                 s
existence1(i).. Y(i)*delta =l= f0(i);
existence2(i).. y(i)*bigmf =g= f0(i);
*Deactivation of nodes including component in IL or IH
Deactive(i,j,k)$(lt(j,k,i) or ht(j,k,i)).. x(j,k) =l= y(i);
*==============================================================================
*                Column Model
*==============================================================================
Equation mb1,mb2,mb3,
mb4,mb5,mb6,mb7,mb8,mb9,mb10
col1,col2,col3
recy2,recy3
purup,purlo,pur12,recovery_const
sharp1
sharp2
*======================================
* Underwood Equation
undw1,undw2,undw3,undw1a,undw2a
undw4,undw5,undw3a,undw7
undw1b,undw1c,undw1d,undw1e
undw2b,undw2c,undw2d,undw2e
undw3b,undw3c
phicut1, phicut2
key1,key2,key3,key4,key5,key6,key7,key8,key9
*======================================
key10,key11,key12,key13,key14
key15,key16
;

mb1(i,j,k)$(colm(j,k)).. F_column(i,j,k) =e= D(i,j,k) + B(i,j,k);
* feed flow to column
mb2(i,jp,kp)$snk(jp,kp)..  F(i,jp,kp) =e= sum((j,k)$arct(j,k,jp,kp),Ds(i,j,k,jp,kp))+sum((j,k)$arcb(j,k,jp,kp),Bs(i,j,k,jp,kp)) +F0D(jp,kp,i);

* distillate disaggregation
mb3(i,j,k,jp,kp)$arct(j,k,jp,kp)..Ds(i,j,k,jp,kp) =g= D(i,j,k) - bigMF*(1-yt(j,k,jp,kp));
mb4(i,j,k,jp,kp)$arct(j,k,jp,kp)..Ds(i,j,k,jp,kp) =l= D(i,j,k);
mb5(i,j,k,jp,kp)$arct(j,k,jp,kp)..Ds(i,j,k,jp,kp) =l= bigMF*yt(j,k,jp,kp);
* bottom disaggregation
mb6(i,j,k,jp,kp)$arcb(j,k,jp,kp)..Bs(i,j,k,jp,kp) =g= B(i,j,k) - bigMF*(1-yb(j,k,jp,kp));
mb7(i,j,k,jp,kp)$arcb(j,k,jp,kp)..Bs(i,j,k,jp,kp) =l= B(i,j,k);
mb8(i,j,k,jp,kp)$arcb(j,k,jp,kp)..Bs(i,j,k,jp,kp) =l= bigMF*yb(j,k,jp,kp);
*Minimum flow rate for activation of node
mb9(jp,kp)$snk(jp,kp)..sum(i,F(i,jp,kp)) =g= 0.0001*x(jp,kp);
*Minimum flow rate for postulated components
mb10(i,jp,kp)$((lt(jp,kp,i) or ht(jp,kp,i)) and snk(jp,kp)).. F(i,jp,kp) =g= 0.0001*X(jp,kp);
* Internal vapor/liquid flow with distillate/bottom stream flow rate
col1(j,k)$colm(j,k).. sum(i,D(i,j,k)) =e= V1(j,k) - L1(j,k);
col2(j,k)$colm(j,k).. sum(i,B(i,j,k)) =e= L2(j,k) - V2(j,k);
* Redundant constraint
col3(j,k)$colm(j,k).. V1(j,k) - L1(j,k) + L2(j,k) - V2(j,k) =e= sum(i,F_column(i,j,k));

*Complete recovery of LLK/HHK
sharp1(i,j,k)$colm(j,k).. D(i,j,k) =l= (1- sum(ip$(ord(ip) lt ord(i)),YHK(j,k,ip)))*bigMf;
sharp2(i,j,k)$colm(j,k).. B(i,j,k) =l=  sum(ip$(ord(ip) lE ord(i)),YLK(j,k,ip))*bigMf;

*=============================================================================================
* Reformulation of Underwood Equation
*=============================================================================================
Equations
undw11
undw22
undw33
undw222
undw333
undw444;

positive variable
V1min(j,k)
V2min(j,k)
S(ip,j,k)
;

undw11(j,k)$(colm(j,k))..V1(j,k) =g= V1min(j,k);
undw22(j,k)$(colm(j,k))..V2(j,k) =g= V2min(j,k);
undw33(j,k)$colm(j,k).. V1(j,k) - V2(j,k) =e= V1min(j,k) - V2min(j,k);


undw2(j,k,ip)$(ord(j) le ord(k) and ord(k) lt card(k) and ra(ip,j,k)).. sum(i,uf(i,ip,j,k)) =E= V1(j,k)-V2(j,k);
undw1(j,k,ip)$(ord(j) le ord(k) and ord(k) lt card(k) and ra(ip,j,k)).. sum(i,Ud(i,ip,j,k)) =g= V1min(j,k) - S(ip,j,k);
undw3(j,k,ip)$(ord(j) le ord(k) and ord(k) lt card(k) and ra(ip,j,k)).. sum(i,Ub(i,ip,j,k)) =g= V2min(j,k) - S(ip,j,k);
undw222(j,k,ip)$(ord(j) le ord(k) and ord(k) lt card(k) and ra(ip,j,k)).. sum(i,Ud(i,ip,j,k)) =l= V1min(j,k)  ;
undw333(j,k,ip)$(ord(j) le ord(k) and ord(k) lt card(k) and ra(ip,j,k)).. sum(i,Ub(i,ip,j,k)) =l= V2min(j,k)  ;
undw444(j,k,ip)$(colm(j,k) and ra(ip,j,k)).. S(ip,j,k) =l= bigMF * (1-Z(j,k,ip)+YHK(j,k,ip+1));

undw2a(i,j,k,ip)$(colm(j,k) and ra(ip,j,k)).. alfa(i)*F_column(i,j,k) =e= (alfa(i)-phi(j,k,ip))*Uf(i,ip,j,k);
undw3a(i,j,k,ip)$(colm(j,k) and ra(ip,j,k)).. alfa(i)*B(i,j,k) =e= -(alfa(i)-phi(j,k,ip))*UB(i,ip,j,k);
undw1a(i,j,k,ip)$(colm(j,k) and ra(ip,j,k)).. alfa(i)*D(i,j,k) =e= (alfa(i)-phi(j,k,ip))*Ud(i,ip,j,k);
* Deactivation of variables when node is not active
uNdw2b(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Uf(i,ip,j,k) =l= bigM*x(j,k);
undw2c(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Uf(i,ip,j,k) =g= -bigM*x(j,k);

uNdw1b(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k) ).. Ud(i,ip,j,k) =l= bigM*x(j,k);
undw1c(i,j,k,ip)$(colm(j,k) and ra(ip,j,k) ).. Ud(i,ip,j,k) =g= -bigM*x(j,k);

uNdw3b(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Ub(i,ip,j,k) =l= bigM*x(j,k);
undw3c(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Ub(i,ip,j,k)=g= -bigM*x(j,k);

* These constraints need to be modified when multiple inlets are present
uNdw2d(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Uf(i,ip,j,k) =l= bigM*y(i);
undw2e(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Uf(i,ip,j,k) =g= -bigM*y(i);
uNdw1d(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Ud(i,ip,j,k) =l= bigM*y(i);
undw1e(i,j,k,ip)$(colm(j,k)  and ra(ip,j,k)).. Ud(i,ip,j,k) =g= -bigM*y(i);

undw5(j,k)$(colm(j,k)).. V1(j,k) =l= xc(j,k)*bigMF;
undw7(j,k)$(colm(j,k)).. V1(j,k) =g= xc(j,k)*0.0001;

undw4(ip,i,j,k)$(colm(j,k) and ra(ip,j,k)).. (phi(j,k,ip)-alfa(i))*(phi(j,k,ip)-alfa(i)) =g= delta;
phicut1(i,j,k)$(rad(i,j,k) and colm(j,k)).. phi(j,k,i) =g= alfa(i+1) * xr(j,k,i);
phicut2(i,j,k)$(rad(i,j,k) and colm(j,k)).. phi(j,k,i+1) =l= alfa(i+1) * xr(j,k,i) + alfa('A') * ( 1 - xr(j,k,i));

key1(j,k,i)$(rad(i,j,k) and colm(j,k)).. -phi(j,k,i+1) + phi(j,k,i) =e= dp(j,k,i);
key2(j,k,i)$(rad(i,j,k) and colm(j,k))..dp(j,k,i) =g= dpp(j,k,i)*xr(j,k,i);
key3(j,k,i)$(rad(i,j,k) and colm(j,k))..dp(j,k,i) =l= dp.up(j,k,i)*xr(j,k,i);

* Location of Underwood roots
key4(j,k,ip)$(colm(j,k) and ra(ip,j,k))..phi(j,k,ip)  =g= sum(i,(alfa(i)+delta+0.001)*yhk(j,k,i));
key5(j,k,ip)$(colm(j,k) and ra(ip,j,k) )..phi(j,k,ip) =l= sum(i,(alfa(i) - delta-0.001)*yLK(j,k,i)) ;

*decide xr
key6(j,k,i)$(colm(j,k) and rad(i,j,k))..xr(j,k,i) =e= z(j,k,i+1)*y(i+1);
*deactivation of xr
key7(j,k,i)$(colm(j,k) and rad(i,j,k))..xr(j,k,i) =l= xc(j,k);

*lk must exist
key8(j,k,i)$(colm(j,k))..Ylk(j,k,i) =l= y(i);
*hk must exist
key9(j,k,i)$(colm(j,k))..Yhk(j,k,i) =l= y(i);
*=============================================================================================
*                Separation Task Key Selection
*=============================================================================================
*If column exists, there should be LK/HK
key10(j,k)$(colm(j,k))..sum(i,YLK(j,k,i)) =e= xc(j,k);
key11(j,k)$(colm(j,k))..sum(i,YhK(j,k,i)) =e= xc(j,k);
* Distributed components
key12(j,k,i)$(colm(j,k))..Z(j,k,i) =e= sum(ip$(ord(ip) le ord(i)-1),YLK(j,k,ip)) - sum(ip$(ord(ip) le ord(i)),YHK(j,k,ip));
*=============================================================================================
*                                Set Operation Part
*=============================================================================================
*determine light key - lk must be one of the lk candidates
key13(j,k ,jp,kp,jpp,kpp)$(arct(j,k,jp,kp) and arcb(j,k,jpp,kpp) and colm(j,k))..
sum(i$ilk1(jp,kp,jpp,kpp,i),YLK(j,k,i)) =g= yb(j,k,jpp,kpp) + yt(j,k,jp,kp)-xc(j,k);
*determine heavy key - hk must be one of the hk candidates
key14(j,k ,jp,kp,jpp,kpp)$(arct(j,k,jp,kp) and arcb(j,k,jpp,kpp) and colm(j,k))..
sum(i$ihk1(jp,kp,jpp,kpp,i),YHK(j,k,i)) =g= yb(j,k,jpp,kpp) + yt(j,k,jp,kp)-xc(j,k);
*i is a lk if it is the heaviest existiing lk candidate
key15(i,ip,j,k ,jp,kp,jpp,kpp)$(arct(j,k,jp,kp) and arcb(j,k,jpp,kpp) and ilk1(jp,kp,jpp,kpp,i) and ilk1(jp,kp,jpp,kpp,ip) and ord(ip) gt ord(i))..
Ylk(j,k,i) =l= 1-Y(ip)+(2-yb(j,k,jpp,kpp) - yt(j,k,jp,kp));
*i is a hk if it is the lightest existiing hk candidate
key16(i,ip,j,k ,jp,kp,jpp,kpp)$(arct(j,k,jp,kp) and arcb(j,k,jpp,kpp) and ihk1(jp,kp,jpp,kpp,i) and ihk1(jp,kp,jpp,kpp,ip) and ord(ip) lt ord(i))..
Yhk(j,k,i) =l= 1-Y(ip)+(2-yb(j,k,jpp,kpp) - yt(j,k,jp,kp));

*=============================================================================================
*                                Vapor/Liquid Flow Rate Balance
*=============================================================================================
Equations
heat_exchanger1
heat_exchanger2

heat_exchanger3a
heat_exchanger4a
heat_exchanger5a
heat_exchanger6a

heat_exchanger_act1
heat_exchanger_act2
heat_exchanger_act3
heat_exchanger_act4
heat_exchanger_act5
heat_exchanger_act6
heat_exchanger_act7
heat_exchanger_act8
heat_exchanger_act9
heat_exchanger_act10
heat_exchanger_act11
heat_exchanger_act12
heat_exchanger_existence
;

*===============================================================================
* If upstream column has heat exchanger, feed is liquid, so no change in vapor flow
heat_exchanger1(j,k)$(colm(j,k)).. V1(j,k) - V2(j,k) =l= bigMF*(Winv(j,k));
heat_exchanger2(j,k)$(colm(j,k)).. V1(j,k) - V2(j,k) =g= -bigMF*(Winv(j,k));
* If thermally coupled, vapor/liquid flow rate balance hold
heat_exchanger3a(j,k)$(colm(j,k)).. V1(j,k) - V2(j,k) - sum((jp,kp)$arct(jp,kp,j,k), V1_T(jp,kp,j,k)) + sum((jp,kp)$arcb(jp,kp,j,k), V2_B(jp,kp,j,k)) =l= bigMF*(1-Winv(j,k)) ;
heat_exchanger4a(j,k)$(colm(j,k)).. V1(j,k) - V2(j,k) - sum((jp,kp)$arct(jp,kp,j,k), V1_T(jp,kp,j,k)) + sum((jp,kp)$arcb(jp,kp,j,k), V2_B(jp,kp,j,k)) =g= -bigMF*(1-Winv(j,k)) ;
heat_exchanger5a(j,k)$(colm(j,k)).. L1(j,k) - L2(j,k) - sum((jp,kp)$arct(jp,kp,j,k), L1_T(jp,kp,j,k)) + sum((jp,kp)$arcb(jp,kp,j,k), L2_B(jp,kp,j,k)) =l= bigMF*(1-Winv(j,k));
heat_exchanger6a(j,k)$(colm(j,k)).. L1(j,k) - L2(j,k) - sum((jp,kp)$arct(jp,kp,j,k), L1_T(jp,kp,j,k)) + sum((jp,kp)$arcb(jp,kp,j,k), L2_B(jp,kp,j,k)) =g= -bigMF*(1-Winv(j,k));

* Disaggregated vapor/liquid flow rate calculation
heat_exchanger_act1(jp,kp,j,k)$(snk(j,k) and arct(jp,kp,j,k)).. V1_T(jp,kp,j,k) - V1(jp,kp) =l= 0;
heat_exchanger_act2(jp,kp,j,k)$(snk(j,k) and arct(jp,kp,j,k)).. V1_T(jp,kp,j,k) - V1(jp,kp) =g= -bigMF * (1 - yt(jp,kp,j,k) + sum(p, yp(j,k,p))$(%relax% = 1));
heat_exchanger_act3(jp,kp,j,k)$(snk(j,k) and arcb(jp,kp,j,k)).. V2_B(jp,kp,j,k) - V2(jp,kp) =l= 0;
heat_exchanger_act4(jp,kp,j,k)$(snk(j,k) and arcb(jp,kp,j,k)).. V2_B(jp,kp,j,k) - V2(jp,kp) =g= -bigMF * (1 - yb(jp,kp,j,k));

heat_exchanger_act5(jp,kp,j,k)$(snk(j,k) and arct(jp,kp,j,k)).. L1_T(jp,kp,j,k) - L1(jp,kp) =l= 0;
heat_exchanger_act6(jp,kp,j,k)$(snk(j,k) and arct(jp,kp,j,k)).. L1_T(jp,kp,j,k) - L1(jp,kp) =g= -bigMF * (1 - yt(jp,kp,j,k));
heat_exchanger_act7(jp,kp,j,k)$(snk(j,k) and arcb(jp,kp,j,k)).. L2_B(jp,kp,j,k) - L2(jp,kp) =l= 0;
heat_exchanger_act8(jp,kp,j,k)$(snk(j,k) and arcB(jp,kp,j,k)).. L2_B(jp,kp,j,k) - L2(jp,kp) =g= -bigMF * (1 - yb(jp,kp,j,k) + sum(p, yp(j,k,p))$(%relax% = 1));
* Deactivation of disaggregated vapor/liquid flow rate
heat_exchanger_act9(jp,kp,j,k)$(snk(j,k) and arct(jp,kp,j,k))..  V1_T(jp,kp,j,k) =l= bigMF * (yt(jp,kp,j,k));
heat_exchanger_act10(jp,kp,j,k)$(snk(j,k) and arcb(jp,kp,j,k)).. V2_B(jp,kp,j,k) =l= bigMF * (yb(jp,kp,j,k));
heat_exchanger_act11(jp,kp,j,k)$(snk(j,k) and arct(jp,kp,j,k)).. L1_T(jp,kp,j,k) =l= bigMF * (yt(jp,kp,j,k));
heat_exchanger_act12(jp,kp,j,k)$(snk(j,k) and arcb(jp,kp,j,k)).. L2_B(jp,kp,j,k) =l= bigMF * (yb(jp,kp,j,k));
* If node is not active, it is thermally coupled
heat_exchanger_existence(j,k)$(snk(j,k)).. 1-Winv(j,k)  =l=  X(j,k) ;
*=============================================================================================
Equations
* When two active connections
two_mb_MODIFIED1
two_mb_MODIFIED2
two_mb_MODIFIED3
;
* If there are two active connections, upstream columns are stacked, and only liquid
* side stream is recovered. Thus, vapor flow rates are equal
two_mb_MODIFIED1(j,k)$(snk(j,k))..
sum((jp,kp)$(arct(jp,kp,j,k)), V1_T(jp,kp,j,k)) - sum((jp,kp)$(arcb(jp,kp,j,k)), V2_B(jp,kp,j,k))
=l= bigMF*(2 - sum((jp,kp)$(arcb(jp,kp,j,k)), yb(jp,kp,j,k))
                 - sum((jp,kp)$(arct(jp,kp,j,k)), yt(jp,kp,j,k))
         );
two_mb_MODIFIED2(j,k)$(snk(j,k))..
sum((jp,kp)$(arct(jp,kp,j,k)), V1_T(jp,kp,j,k)) - sum((jp,kp)$(arcb(jp,kp,j,k)), V2_B(jp,kp,j,k))
=g= -bigMF*(2 - sum((jp,kp)$(arcb(jp,kp,j,k)), yb(jp,kp,j,k))
                 - sum((jp,kp)$(arct(jp,kp,j,k)), yt(jp,kp,j,k))
         );
*If there are two active arcs, it should be thermally coupled
two_mb_MODIFIED3(j,k)$(snk(j,k))..
2 - sum((jp,kp)$(arcb(jp,kp,j,k)), yb(jp,kp,j,k))
  - sum((jp,kp)$(arct(jp,kp,j,k)), yt(jp,kp,j,k)) =g= 1-Winv(j,k);

*===============================================================================
*        Conditional Constraints required when bypass stream relaxation is considered
*===============================================================================
$if %relax% == 1 Equation correction1, correction2, correction3, correction4;
*Even though bypass stream relax the lower bound, it should stay the same for nodes with heat exchanger (Win(jk)=0)
$if %relax% == 1 correction1(jp,kp,j,k)$(colm(jp,kp) and snk(j,k) and arct(jp,kp,j,k))..V1_T(jp,kp,j,k) =g= V1(jp,kp) - bigMF * (1 - YT(jp,kp,j,k) +  Winv(j,k)  );
$if %relax% == 1 correction2(jp,kp,j,k)$(colm(jp,kp) and snk(j,k) and arcb(jp,kp,j,k))..L2_B(jp,kp,j,k) =g= L2(jp,kp) - bigMF * (1 - YB(jp,kp,j,k) +  Winv(j,k)  );
*===============================================================================
* If two connections are active, thermal coupling is enforced, then lower bound on V1_T is relaxed.
* Then, V1_T can be anything If liquid side stream product is required
$if %relax% == 1 correction3(jp,kp,j,k)$(colm(jp,kp) and snk(j,k) and arct(jp,kp,j,k))..V1_T(jp,kp,j,k) =g= V1(jp,kp)- bigMF * [2 - YT(jp,kp,j,k) - sum((jpp,kpp)$(colm(jpp,kpp) and arcb(jpp,kpp,j,k)), yb(jpp,kpp,j,k))] ;
$if %relax% == 1 correction4(jp,kp,j,k)$(colm(jp,kp) and snk(j,k) and arcb(jp,kp,j,k))..V2_B(jp,kp,j,k) =g= V2(jp,kp)- bigMF * [2 - YB(jp,kp,j,k) - sum((jpp,kpp)$(colm(jpp,kpp) and arct(jpp,kpp,j,k)), yt(jpp,kpp,j,k))] ;
*===============================================================================
*        Heat duty calculation
*===============================================================================
equations heat1,heat2,heat3,heat4,heat6
$if %vaporbypass% == 0 heat5
$if %vaporbypass% == 1 heat5a, heat5b
;

heat1(jpp,kpp)$(snk(jpp,kpp) )..QB(jpp,kpp) - sum((jp,kp)$(arcb(jp,kp,jpp,kpp)), V2_B(jp,kp,jpp,kpp) ) =l= bigMF * (Winv(jpp,kpp));
heat2(jpp,kpp)$(snk(jpp,kpp) )..QB(jpp,kpp) - sum((jp,kp)$(arcb(jp,kp,jpp,kpp)), V2_B(jp,kp,jpp,kpp) ) =g= -bigMF * (Winv(jpp,kpp));
heat3(jpp,kpp)$(snk(jpp,kpp) )..QT(jpp,kpp) - sum((jp,kp)$(arct(jp,kp,jpp,kpp)), V1_T(jp,kp,jpp,kpp) ) =l= bigMF * (Winv(jpp,kpp));
heat4(jpp,kpp)$(snk(jpp,kpp))..QT(jpp,kpp) - sum((jp,kp)$(arct(jp,kp,jpp,kpp)), V1_T(jp,kp,jpp,kpp) ) =g= -bigMF * (Winv(jpp,kpp));
heat6(j,k)$snk(j,k)..sum((jp,kp)$(arcb(jp,kp,j,k)), Qb(j,k))  =l= bigMF * (1-Winv(j,k));
$if %vaporbypass% == 0 heat5(j,k)$snk(j,k)..sum((jp,kp)$(arct(jp,kp,j,k)), QT(j,k))  =l= bigMF * (1-Winv(j,k));
$if %vaporbypass% == 1 heat5a(jp,kp,j,k)$(arct(jp,kp,j,k) and snk(j,k))..QT(j,k) =g= V1(jp,kp) - V1_T(jp,kp,j,k) - bigMF * (2 - Winv(j,k) - YT(jp,kp,j,k));
$if %vaporbypass% == 1 heat5b(jp,kp,j,k)$(arct(jp,kp,j,k) and snk(j,k))..QT(j,k) =l= V1(jp,kp) - V1_T(jp,kp,j,k) + bigMF * (2 - Winv(j,k) - YT(jp,kp,j,k));

*====================================================================================
*        Product Specifications/Objective Function
*====================================================================================
obj.. cost=e= sum((i,j,k,p)$(snk(j,k) and ord(k)=card(k) and (ord(p)>1)), price(i)*F_bypass(i,j,k,p))
         - ReactorCost('1')*sum(i, FF1(i))
         - ReactorCost('2')*sum(i, FF2(i))
         - ReactorCost('3')*sum(i, FF3(i))
         - %Energycost%*sum((j,k)$snk(j,k), QT(j,k)+QB(j,k))
         - %Capitalcost%*sum((j,k)$colm(j,k),[V1(j,k) + V2(j,k)]);
;
* Purity constraint
purup(i,p)$(IPP(i,p)).. sum((j,k)$(snk(j,k)), F_bypass(i,j,k,p)) =g= pp(i,p) * sum((ip,j,k)$(snk(j,k)), F_bypass(ip,j,k,p));
purlo(i,p)$(IPNP(i,p)).. sum((j,k)$(snk(j,k)), F_bypass(i,j,k,p)) =l= 0 ;
recovery_const(i,p)$(IRP(i,p))..sum((j,k)$(snk(j,k)), F_bypass(i,j,k,p)) =g= rr(i,p) * f0(i);

* Pure components into pure component nodes
pur12(ip,j,k)$(ord(k)=card(k) and ord(ip)=ord(j)).. F(ip,j,k) =g= 1*sum(i,F(i,j,k));
* Optional Minimum recovery of keys
recy2(i,j,k)$(colm(j,k)).. D(i,j,k) =g= rlk(i)*F_column(i,j,k) - bigMF*(1-YLK(j,k,i));
recy3(i,j,k)$(colm(j,k)).. B(i,j,k) =g= (1-rhk(i))*F_column(i,j,k) - bigMF*(1-YHK(j,k,i));

*===============================================================================
*        Terminal Node Constraints
*===============================================================================
* If there is one active incoming stream while no column inside, there should be heat exchanger
Equation Final_product_constraint1, Final_product_constraint2;
Final_product_constraint1(j,k)$(snk(j,k)).. 1 =g= -sum((jp,kp)$(arct(jp,kp,j,k)), YT(jp,kp,j,k)) + sum((jp,kp)$(arcb(jp,kp,j,k)), YB(jp,kp,j,k))
                                            -(XC(j,k))$(colm(j,k)) + Winv(j,k) ;
Final_product_constraint2(j,k)$(snk(j,k)).. 1 =g= sum((jp,kp)$(arct(jp,kp,j,k)), YT(jp,kp,j,k)) - sum((jp,kp)$(arcb(jp,kp,j,k)), YB(jp,kp,j,k))
                                            -(XC(j,k))$(colm(j,k)) + Winv(j,k);


*===============================================================================
*        Root flowing
*===============================================================================
Equation
RF1
RF2
;
RF1(ip,j,k,jp,kp)$(arct(j,k,jp,kp) and ra(ip,j,k) and ra(ip,jp,kp))..phi(j,k,ip) =l= phi(jp,kp,ip) + (alfa('A')-1)*[(2 - XC(j,k) - XC(jp,kp)) + (1-Winv(jp,kp)) + (sum((jpp,kpp)$arct(jpp,kpp,jp,kp), YT(jpp,kpp,jp,kp))+sum((jpp,kpp)$arcb(jpp,kpp,jp,kp),YB(jpp,kpp,jp,kp))-1)$(ord(jp) ne ord(kp) and ord(jp) ne 1) ];
RF2(ip,j,k,jp,kp)$(arcb(j,k,jp,kp) and ra(ip,j,k) and ra(ip,jp,kp))..phi(j,k,ip) =g= phi(jp,kp,ip) - (alfa('A')-1)*[(2 - XC(j,k) - XC(jp,kp)) + (1-Winv(jp,kp)) + (sum((jpp,kpp)$arct(jpp,kpp,jp,kp), YT(jpp,kpp,jp,kp))+sum((jpp,kpp)$arcb(jpp,kpp,jp,kp),YB(jpp,kpp,jp,kp))-1)$(ord(jp) ne ord(kp) and ord(jp) ne 1) ];


*===============================================================================
*        Variable Bounds
*===============================================================================
F.up(i,j,k)=10;
D.up(i,j,k)=F.up(i,j,k);
B.up(i,j,k)=F.up(i,j,k);
V1.up(j,k) = 50;
V2.up(j,k)=50;
L1.up(j,k) = 50;
L2.up(j,k)=50;
Bs.up(i,j,k,jp,kp)=10;
Ds.up(i,j,k,jp,kp)=10;

Ud.lo(i,ip,j,k)=-50;
Ud.up(i,ip,j,k)=50;
Uf.lo(i,ip,j,k)=-50;
Uf.up(i,ip,j,k)=50;
Ub.lo(i,ip,j,k)=-50;
Ub.up(i,ip,j,k)=50;

Phi.up(j,k,i)=smax(ip,alfa(ip));
loop(i,
phi.up(j,k,i)$(ord(j)=ord(i))=alfa(i)-delta;)
;
dp.up(j,k,i)$(colm(j,k))=phi.up(j,k,i)-phi.lo(j,k,i);

display phi.up,phi.lo,phi.l;

*==============================================================================
*                Reactor Modeling Part
*==============================================================================

Equation
reactor1
reactor2
reactor3
reactor4
reactor5
reactor6
reactor7
*reactor8
reactor9

reactor11
reactor12
reactor13
reactor14
reactor15
reactor16
reactor17
reactor18

;
reactor1..F1('A') =e= FF1('A')  - E1;
reactor2..F1('B') =e= FF1('B')  + E1;
reactor3..F1('C') =e= FF1('C')  + E1;
F1.fx('D') = 0;
reactor4..F2('A') =e= FF2('A')  - E2;
reactor5..F2('C') =e= FF2('C')  + E2;
F2.fx('B') = 0;
reactor6..F2('D') =e= FF2('D')  + E2;

reactor7..F3('A') =e= FF3('A')  - E3;
*reactor8..F3('B') =e= FF3('B')  + E3;
reactor9..F3('D') =e= FF3('D')  + 2*E3;

F3.fx('C') = 0;
F3.fx('B') = 0;

reactor11.. E1 =e= x1 * (FF1('A'))  ;
reactor12.. E2 =e= x2 * (FF2('A')) ;
reactor13.. E3 =e= x3 * (FF3('A')) ;

reactor14.. y1+y2+y3 =e= 1;
reactor15(i).. F1(i) + F2(i) + F3(i) =e= f0(i);
reactor16(i)$(sameas(i,'A')).. FF1(i)  =e= 1*y1;
reactor17(i)$(sameas(i,'A')).. FF2(i)  =e= 1*y2;
reactor18(i)$(sameas(i,'A')).. FF3(i)  =e= 1*y3;

FF1.fx(i)$(sameas(i,'B') or sameas(i,'C') or sameas(i,'D')) = 0;
FF2.fx(i)$(sameas(i,'B') or sameas(i,'C') or sameas(i,'D')) = 0;
FF3.fx(i)$(sameas(i,'B') or sameas(i,'C') or sameas(i,'D')) = 0;

*QB.lo('3','3') =0.1;

*equation
*add1;
*add1(i).. F_column(i,'3','3') =e= 0;
*==============================================================================
*                Model Definition
*==============================================================================
model cut /
all
-Rf1
-Rf2
*-heat1
*-heat6
*-heat_exchanger_act1
*-heat_exchanger_act2
*-heat_exchanger_act3
*-heat_exchanger_act4
*-heat_exchanger_act5
*-heat_exchanger_act6
*-heat_exchanger_act7
*-heat_exchanger_act8
/
;

*==============================================================================
*                Solver Options
*==============================================================================
*cut.OptFile = 1;
option mip=cplex
option nlp =conopt4;
option minlp=baron;
option reslim = 1000;
option rminlp=baron;
option optcr=0.0000;

*==============================================================================
*                Final Additional Constraints and Solve statement
*==============================================================================
*y2.fx=1;
solve cut using minlp max cost;

*==============================================================================
*                Data Preparation
*==============================================================================


parameter
product_flow(i,p)
FL(i,j,k)
FBL(i,j,k,p)
W_inv_par(j,k)
W_par(j,k)
Revenue
Reactor_cost
Energy_cost
Capital_cost
data_save(j,k,*)
;

product_flow(i,p) = sum((j,k), F_bypass.l(i,j,k,p));
FL(i,j,k) = F_column.l(i,j,k);
FBL(i,j,k,p) = F_bypass.l(i,j,k,p);
W_inv_par(j,k)$(y0.l(j,k) ne 1) = Winv.l(j,k);
W_par(j,k)$(snk(j,k) and y0.l(j,k) ne 1 and X.l(j,k)=1)= 1 - Winv.l(j,k);
Revenue = sum((i,j,k,p)$(snk(j,k) and ord(k)=card(k) and (ord(p)>1)), price(i)*F_bypass.l(i,j,k,p));
Reactor_cost = ReactorCost('1')*sum(i, FF1.l(i)) + ReactorCost('2')*sum(i, FF2.l(i)) + ReactorCost('3')*sum(i, FF3.l(i)) ;
Energy_cost =  %Energycost%*sum((j,k)$snk(j,k), QT.l(j,k)+QB.l(j,k));
Capital_cost =  %Capitalcost%*sum((j,k)$colm(j,k),[V1.l(j,k) + V2.l(j,k)]);
data_save(j,k,'V1')$(colm(j,k)) = V1.l(j,k);
data_save(j,k,'V2')$(colm(j,k)) = V2.l(j,k);
data_save(j,k,'L1')$(colm(j,k)) = L1.l(j,k);
data_save(j,k,'L2')$(colm(j,k)) = L2.l(j,k);


option FL:3:1:2, FBL:3:3:1;

execute_unload 'Ex1_result.gdx';

display
x.l
phi.l
ylk.l
yhk.l
W_par
W_inv_par
*Y_H.l
*Y_L.l
Qt.l
Qb.l
V1.l
V2.l
L1.l
L2.l
V1_T.l
V2_B.l
L1_T.l
L2_B.l
f0.l
F1.l
F2.l
F3.l
F.l
FL
FBL
F_bypass.l
D.l
B.l
yp.l
yt.l
yb.l
product_flow
Revenue
Reactor_cost
Energy_cost
Capital_cost
data_save
  ;

