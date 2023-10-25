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


parameter z0(R) /
1 250.80078773591563
2 3534.568207268975
3 564439.1039092198
4 650.6406752988532
5 968.820082036218
6 27535.746897144163
7 4716.298642953472
8 22055.50292252441
9 300211.1750508643
10 99294.63938700294
11 641899.4880127095
12 250060.1707892678
13 1721312.1577844007
14 2009.1043394730968
15 54520.0607255271
16 49367.59583538092
17 352615.37500308617
18 687987.2643486774
19 1556095.6496074467
20 4628754.300016958
21 7056730.818600473
22 300993.4347421814
23 376492.1242193523
24 501972.5719503983
25 597713.3877686
26 705355.6161430755
27 371144.01942469337
28 1554126.330330825
29 7503364.662832906
30 559098.2865366042
31 7471.360220528475
32 1520541.845735727
33 2627528.4312725803
34 4986665.422304678
35 5394735.127975366
36 4499513.552096191
37 2483766.692381097
38 1636997.8774710726
39 2771996.557329408
40 830023.379807841




/;

parameter zbar(R) /
1 668418.5913065937
2 2779482.9439592776
3 7456860.319674995
4 630703.6475810524
5 3338123.2486141897
6 3247777.0276685245
7 12789009.089098873
8 13615218.935454436
9 14176077.614825662
10 15433469.722759722
11 15345830.215443855
12 9283388.110591304
13 3969407.115423103
14 3408237.8629657542
15 14216603.335133215
16 15582609.418698223
17 15008209.424406173
18 14781262.755443327
19 15108658.049888996
20 15458704.737096386
21 13904836.584494572
22 903477.2922439108
23 13618843.667688578
24 16050793.19262713
25 15655699.743307205
26 15478777.16134089
27 15752105.83800765
28 15841353.955941873
29 13707565.814692184
30 805065.3054766154
31 1903753.478378784
32 6521153.62854427
33 8856108.22030659
34 14717726.870318938
35 14049906.508759096
36 15382032.80416708
37 4573449.39131777
38 3521558.523282204
39 4781681.083117344
40 1915466.8266412937





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

