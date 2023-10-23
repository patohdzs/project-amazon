option limrow = 0;
option limcol = 0;
option solprint = off;

set T ordered /0*200 /;
set R ordered /1*24/;

scalar scale / 1.0e9 /;



parameter x0(R);
$call csv2gdx X0Data.csv id=x0  index=1 values=2..lastCol useHeader=y trace=0
$gdxIn X0Data.gdx
$load  x0
$gdxIn 
display x0;


parameter z0(R) /
1 5589.185181
2 787674.0828
3 9646.861745
4 128634.5031
5 16248.65369
6 77127.05958
7 413615.3759
8 1415637.311
9 2975490.403
10 7911410.246
11 311951.5001
12 491859.8041
13 548889.5383
14 1193672.493
15 1416617.977
16 6055143.679
17 6746442.015
18 84410.49025
19 2448200.526
20 7385903.416
21 8456499.239
22 5225548.432
23 860690.5049
24 2311937.326



/;

parameter zbar(R) /
1 5923957.508
2 16269303.58
3 8317106.102
4 11492722.45
5 3641993.119
6 27941871.29
7 26031526.11
8 26494123.96
9 25883433.49
10 18981310.08
11 986453.6577
12 23983732.12
13 28379164.51
14 27505450.39
15 28080696.69
16 27530403.94
17 9828455.437
18 4377913.526
19 8346530.632
20 21218696.15
21 20770966.56
22 16507489.38
23 1577160.408
24 3458325.658





/;

 

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
parameter p_a / 44.75/;

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


file results_x / "results/amazon_data_x.dat" /;
file results_w / "results/amazon_data_w.dat" /;

file regionresults_z / "results/amazon_data_z.dat" /;
file regionresults_u / "results/amazon_data_u.dat" /;
file regionresults_v / "results/amazon_data_v.dat" /;

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

