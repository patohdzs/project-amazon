option limrow = 0;
option limcol = 0;
option solprint = off;

set T ordered /0*200 /;
set R ordered /1*78/;

scalar scale / 1.0e9 /;



parameter x0(R);
$call csv2gdx X0Data.csv id=x0  index=1 values=2..lastCol useHeader=y trace=0
$gdxIn X0Data.gdx
$load  x0
$gdxIn 
display x0;


parameter z0(R) /
1 5713.628934761566
2 147689.3937465045
3 9344.426290023148
4 577.0929400978703
5 1350.812047242954
6 3758.080485774173
7 8963.86542448523
8 625307.1947267982
9 7103.935333340407
10 2542.926411305254
11 3581.564117650438
12 115672.09184089689
13 852.1239077835021
14 2187.0705355650616
15 8018.458580218519
16 26785.50539495797
17 63523.93014856061
18 9171.070457588336
19 357997.8998547132
20 362379.7711735249
21 98531.7254632279
22 713183.2333392987
23 1062176.598537366
24 15396.529777370208
25 30235.674703547968
26 36685.85576082633
27 29872.3815064708
28 293433.5588127973
29 286460.17902218533
30 762008.1618069287
31 1048023.9231448268
32 1466554.9830308636
33 2505504.298582405
34 3630546.1151278983
35 300993.4347421814
36 1677.500344874387
37 22228.959840250005
38 13781.180081558976
39 9385.082155257762
40 33392.67025352242
41 94994.45112957398
42 67688.62727610428
43 429667.5655363259
44 543244.4131700981
45 2574056.481091863
46 3907754.71387018
47 1033201.850958964
48 230923.08910011404
49 237030.25479343388
50 383720.6725936014
51 142002.6034327703
52 642573.0423107159
53 422712.32938820776
54 143328.89572640735
55 775932.8880314889
56 112998.60999727668
57 2824844.17491085
58 1805485.450413916
59 1933.989422476656
60 82476.5008257816
61 1443602.715707997
62 968643.7866233932
63 3398858.5621535047
64 1226145.4334957227
65 2640896.539930064
66 3210661.266963465
67 1466631.4985449116
68 1992156.500752608
69 35954.02398489008
70 1003991.536193392
71 1756907.883689507
72 662782.7461972998
73 1942158.6855006733
74 1263981.622725933
75 502778.80999226007
76 860690.5049079219
77 2149764.6691846712
78 159391.58816178655




/;

parameter zbar(R) /
1 3205198.859834969
2 1368030.8055876978
3 1215945.072311783
4 827842.6730102487
5 3933940.238911268
6 1812371.995813367
7 6085749.816471196
8 5610324.102973549
9 3419085.232796117
10 4898020.869117723
11 5746261.74412789
12 4427045.0974597465
13 1615504.4326193312
14 7080140.335864618
15 6980452.768150298
16 6077666.371600857
17 6710199.81977433
18 7189723.396478963
19 6843905.459975344
20 6547659.59833465
21 5635336.488127723
22 2951317.8562700404
23 2270927.73951364
24 2026488.6864479545
25 7016105.196144057
26 6865172.989221514
27 6769856.859815905
28 6473803.057185335
29 6272576.532592798
30 6187918.567517067
31 6783239.11990195
32 6917198.278746616
33 6847818.215978191
34 6911246.265614724
35 903477.2922439108
36 3923651.056674422
37 7153919.837935575
38 7137828.745054413
39 7128421.011116844
40 6963686.102870945
41 6934925.333682166
42 7072091.702277897
43 7131869.105470084
44 7088549.145154162
45 7152090.871228486
46 5545534.502206791
47 1930475.1565786377
48 5771525.0021374775
49 7134636.21968513
50 7130323.314335648
51 6982591.444101469
52 6749146.300699924
53 6857692.649058391
54 6951490.588027891
55 6925245.291118768
56 7070400.547221225
57 6219363.373477953
58 2352445.7779283584
59 761650.1796125177
60 3616263.3468388943
61 4046993.580471642
62 3862731.19142241
63 6899631.38409216
64 7068235.401118808
65 7014755.113686711
66 6973697.21585846
67 6957779.728332206
68 3624090.7578613865
69 436805.8600466905
70 2849125.8629063475
71 4401703.506077482
72 2630069.582148917
73 4152444.650218677
74 4952058.725588936
75 973560.1652094162
76 1577160.4082438345
77 3177990.268652912
78 270291.9381979981






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

