option limrow = 0;
option limcol = 0;
option solprint = off;

set T ordered /0*200 /;
set R ordered /1*40/;

scalar scale / 1.0e9 /;



parameter x0(R);
$call csv2gdx X0Data.csv id=x0  index=1 values=2..lastCol useHeader=y trace=0
$gdxIn X0Data.gdx
$load  x0
$gdxIn 
display x0;


parameter z0(R);
$call csv2gdx Z0Data.csv id=z0  index=1 values=2..lastCol useHeader=y trace=0
$gdxIn Z0Data.gdx
$load  z0
$gdxIn 
display z0;

parameter zbar(R);
$call csv2gdx ZbarData.csv id=zbar  index=1 values=2..lastCol useHeader=y trace=0
$gdxIn ZbarData.gdx
$load  zbar
$gdxIn 
display zbar;


parameter gamma(R);

$call csv2gdx GammaData.csv id=gamma  index=1 values=2..lastCol useHeader=y trace=0
$gdxIn GammaData.gdx
$load  gamma 
$gdxIn

display gamma;


parameter theta(R);
$call csv2gdx ThetaData.csv id=theta  index=1 values=2..lastCol useHeader=y trace=0
$gdxIn ThetaData.gdx
$load  theta 
$gdxIn
display theta;


parameter delta / 0.02 /;
parameter p_e / 25 /;
parameter p_a / 41.945/;

parameter alpha / 0.045007414 /;
parameter kappa / 2.094215255 /;
parameter zeta / 1.66e-04 /;
parameter dt /1 /; 


parameter x_agg(T);
parameter w_agg(T);

parameter z_agg(T);
parameter u_agg(T);
parameter v_agg(T);
parameter c_agg(T);

positive variable z(T,R);
positive variable u(T,R);
positive variable v(T,R);
         variable x(T,R);
         variable w(T);
         variable obj;

equation zdot_def(T,R);
equation xdot_def(T,R);
equation w_def(T);
equation obj_def;

zdot_def(T,R)$(ord(T) < card(T))..
 ( z(T+1,R) - z(T,R))/dt =e= (u(T,R) - v(T,R));

xdot_def(T,R)$(ord(T) < card(T))..
  (x(T+1,R) - x(T,R))/dt =e= (-gamma(R)*u(T,R) - alpha*x(T,R) + alpha*gamma(R)*(zbar(R)/scale - z(T,R)));

w_def(T)$(ord(T) < card(T))..
  w(T) =e= sum(R, u(T,R) + v(T,R));

obj_def..
  obj =e= sum(T$(ord(T) < card(T)), exp(-delta*(ord(T)*dt-dt))*(-p_e*sum(R, kappa*z(T,R) - (x(T+1,R) - x(T,R))/dt) + p_a*sum(R, theta(R)*z(T,R)) - scale*zeta/2*sqr(w(T)))*dt);
    
model amazon / all /;


file results_x / "amazon_data_x.dat" /;
file results_w / "amazon_data_w.dat" /;

file regionresults_z / "amazon_data_z.dat" /;
file regionresults_u / "amazon_data_u.dat" /;
file regionresults_v / "amazon_data_v.dat" /;

regionresults_u.pw = 163840000;
regionresults_v.pw = 163840000;
regionresults_z.pw = 163840000;
x.fx(T,R)$(ord(T) = 1) = x0(R) / scale; 
z.fx(T,R)$(ord(T) = 1) = z0(R) / scale;
z.up(T,R)$(ord(T) > 1) = zbar(R) / scale ;
u.fx(T,R)$(ord(T) = card(T)) = 0;
v.fx(T,R)$(ord(T) = card(T)) = 0;
w.fx(T)$(ord(T) = card(T)) = 0;

option qcp = cplex;
solve amazon using qcp maximizing  obj;

x_agg(T) = sum(R, x.l(T,R));
w_agg(T) = w.l(T);
z_agg(T) = sum(R, z.l(T,R));
u_agg(T) = sum(R, u.l(T,R));
v_agg(T) = sum(R, v.l(T,R));
c_agg(T) = sum(R, u.l(T,R)*v.l(T,R));

put results_w;
put 'T':4 system.tab 'w_agg':20/;
loop(T, put T.tl:4 system.tab w_agg(T):20:16 /);
putclose;

put results_x;
put 'T':4 system.tab 'x_agg':20/;
loop(T, put T.tl:4 system.tab x_agg(T):20:16 /);
putclose;

put regionresults_z;
put 'T/R':4;
loop(R, put system.tab R.tl:10);
put /;
loop(T,
  put T.tl:4;
  loop(R, put system.tab (z.l(T,R)):16:10);

  put /;
);
putclose;

put regionresults_u;
put 'T/R':4;
loop(R, put system.tab R.tl:10);
put /;
loop(T,
  put T.tl:4;
  loop(R, put system.tab (u.l(T,R)):16:10);

  put /;
);
putclose;

put regionresults_v;
put 'T/R':4;
loop(R, put system.tab R.tl:10);
put /;
loop(T,
  put T.tl:4;
  loop(R, put system.tab (v.l(T,R)):16:10);

  put /;
);
putclose;

