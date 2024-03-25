option limrow = 0;
option limcol = 0;
option solprint = off;

set T /0*200/;
set R /1*78/;
set Y /1*200/;
set S /1*16/;
alias(S,S1);
 
set dep(T,S,S1) /
    0.2.1   0.3.1   0.4.1   0.5.1   0.6.1   0.7.1   0.8.1   0.9.1   0.10.1  0.11.1  0.12.1  0.13.1  0.14.1  0.15.1  0.16.1 
    1.2.1   1.3.1   1.4.1   1.5.1   1.6.1   1.7.1   1.8.1   1.10.9  1.11.9  1.12.9  1.13.9  1.14.9  1.15.9  1.16.9
    2.2.1   2.3.1   2.4.1   2.6.5   2.7.5   2.8.5   2.10.9  2.11.9  2.12.9  2.14.13 2.15.13 2.16.13 
    3.2.1   3.4.3   3.6.5   3.8.7   3.10.9  3.12.11 3.14.13 3.16.15    
/;

scalar scale / 1.0e11 /;

scalar p_a_low / 35.71 /;
scalar p_a_high / 44.25 /;

parameter scenario(Y) ;
$call csv2gdx mc_1.csv id=scenario index=1 values=2..lastCol useHeader=y trace=0
$gdxIn mc_1.gdx
$load  scenario
$gdxIn 
display scenario;

parameter prob(S);


parameter prob_low(S) /
1 0.2490479683700778
2 0.1034951587545163
3 0.025092369985659297
4 0.12141156316570699
5 0.0250923699856593
6 0.010427464364350469
7 0.02943619682459999
8 0.14242953822897736
9 0.0250923699856593
10 0.010427464364350469
11 0.0025281355861598875
12 0.012232598737622108
13 0.029436196824599994
14 0.01223259873762211
15 0.034531998531497965
16 0.16708600755294128
/;

parameter prob_high(S) /
1 0.060381604721641786
2 0.0250923699856593
3 0.006083637525409759
4 0.02943619682459999
5 0.006083637525409759
6 0.0025281355861598875
7 0.007136797030724112
8 0.03453199853149797
9 0.07083447287710407
10 0.029436196824599994
11 0.007136797030724112
12 0.03453199853149797
13 0.08309687314386374
14 0.034531998531497965
15 0.09748205987595236
16 0.4716752254536568
/;

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
parameter p_e / 6.9 /;
table p_a(S,T) 
	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200
1	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
2	35.71	35.71	35.71	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25
3	35.71	35.71	35.71	44.25	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
4	35.71	35.71	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25
5	35.71	35.71	44.25	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
6	35.71	35.71	44.25	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25
7	35.71	35.71	44.25	44.25	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
8	35.71	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25
9	35.71	44.25	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
10	35.71	44.25	35.71	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25
11	35.71	44.25	35.71	44.25	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
12	35.71	44.25	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25
13	35.71	44.25	44.25	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
14	35.71	44.25	44.25	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25
15	35.71	44.25	44.25	44.25	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71	35.71
16	35.71	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25	44.25;

parameter alpha / 0.045007414 /;
parameter kappa / 2.094215255 /;
parameter zeta / 1.66e-04 /;

parameter x_agg(S,T);
parameter z_agg(S,T);
parameter u_agg(S,T);
parameter v_agg(S,T);
parameter c_agg(S,T);
parameter objval(S);

positive variable z(S,T,R);
positive variable u(S,T,R);
positive variable v(S,T,R);
         variable x(S,T,R);
         variable w(S,T);
         variable obj;

equation zdot_def(S,T,R);
equation xdot_def(S,T,R);
equation w_def(S,T);
equation non_x_def(T,S,S1,R);
equation non_z_def(T,S,S1,R);
equation non_u_def(T,S,S1,R);
equation non_v_def(T,S,S1,R);
equation non_w_def(T,S,S1);
equation obj_def;

zdot_def(S,T,R)$(ord(T) < card(T))..
  z(S,T+1,R) - z(S,T,R) =e= u(S,T,R) - v(S,T,R);

xdot_def(S,T,R)$(ord(T) < card(T))..
  x(S,T+1,R) - x(S,T,R) =e= -gamma(R)*u(S,T,R) - alpha*x(S,T,R) + alpha*gamma(R)*(zbar(R)/scale - z(S,T,R));

w_def(S,T)$(ord(T) < card(T))..
  w(S,T) =e= sum(R, u(S,T,R) + v(S,T,R));

non_x_def(dep(T,S,S1),R)..
  x(S,T,R) =e= x(S1,T,R);

non_z_def(dep(T,S,S1),R)..
  z(S,T,R) =e= z(S1,T,R);

non_u_def(dep(T,S,S1),R)..
  u(S,T,R) =e= u(S1,T,R);

non_v_def(dep(T,S,S1),R)..
  v(S,T,R) =e= v(S1,T,R);

non_w_def(dep(T,S,S1))..
  w(S,T) =e= w(S1,T);

obj_def..
  obj =e= sum(S, prob(S)*sum(T$(ord(T) < card(T)), exp(-delta*(ord(T)-1))*(-p_e*sum(R, kappa*z(S,T+1,R) - (x(S,T+1,R) - x(S,T,R))) + p_a(S,T)*sum(R, theta(R)*z(S,T+1,R)) - scale*zeta/2*sqr(w(S,T)))));

model amazon / all /;

file results / 25sites_output_dynamic_200.dat /;
file regionresults / 25sites_output_regions_dynamic_200.dat /;
file regionresults_x / "amazon_data_x.dat" /;
file regionresults_w / "amazon_data_w.dat" /;
file regionresults_z / "amazon_data_z.dat" /;
file regionresults_u / "amazon_data_u.dat" /;
file regionresults_v / "amazon_data_v.dat" /;

results.nr = 2;
regionresults.nr = 2;
regionresults.pw = 16384;

x.l(S,T,R)$(ord(T) = 1) = x0(R) / scale;
z.l(S,T,R)$(ord(T) = 1) = z0(R) / scale;
z.up(S,T,R)$(ord(T) > 1) = zbar(R) / scale;




loop(Y,
* fix boundary values

  x.fx(S,T,R)$(ord(T) = 1) = x.l(S,T,R);
  z.fx(S,T,R)$(ord(T) = 1) = z.l(S,T,R);
  u.fx(S,T,R)$(ord(T) = card(T)) = 0;
  v.fx(S,T,R)$(ord(T) = card(T)) = 0;
  w.fx(S,T)$(ord(T) = card(T)) = 0;

* update the probabilities based on the actual scenario
* update the initial prices based on the actual scenario

  prob(S) = prob_low(S)$(scenario(Y) = 1) + prob_high(S)$(scenario(Y) = 2);
  p_a(S,'0') = p_a_low$(scenario(Y) = 1) + p_a_high$(scenario(Y) = 2);
  
  option qcp = cplex;
  solve amazon using qcp maximizing obj;

* output the results
  x_agg(S,T) = scale*sum(R, x.l(S,T,R));
  z_agg(S,T) = scale*sum(R, z.l(S,T,R));
  u_agg(S,T) = scale*sum(R, u.l(S,T,R));
  v_agg(S,T) = scale*sum(R, v.l(S,T,R));
  c_agg(S,T) = scale*scale*sum(R, u.l(S,T,R)*v.l(S,T,R));
  objval(S) = sum(T$(ord(T) < card(T)), exp(-delta*(ord(T)-1))*(-p_e*sum(R, kappa*z.l(S,T+1,R) - (x.l(S,T+1,R) - x.l(S,T,R))) + p_a(S,T)*sum(R, theta(R)*z.l(S,T+1,R)) - scale*zeta/2*sqr(w.l(S,T))));

  put results;
  put // 'Y = ' Y.tl:4 '  p_a = ' p_a('1','0') //;
  loop(S,
    put // 'S = ' S.tl:4 '  obj = ' objval(S) //;
    put 'T':4 system.tab 'x_agg':10 system.tab 'z_agg':10 system.tab 'u_agg':10 system.tab 'v_agg':10 system.tab 'c_agg':10 /;
    loop(T, 
      put T.tl:4 system.tab x_agg(S,T):10:6 system.tab z_agg(S,T):10:6 system.tab u_agg(S,T):10:6 system.tab v_agg(S,T):10:6 system.tab c_agg(S,T):10:6 /
    );
  );

* output region, x0, z0, u0, v0
  put regionresults_z;
  put // 'Y = ' Y.tl:4 '  p_a = ' p_a('1','0') //;
  loop(R, 
    put system.tab (z.l('1','0',R)*scale):10:6 /
  );
  
  put regionresults_x;
  put // 'Y = ' Y.tl:4 '  p_a = ' p_a('1','0') //;
  loop(R, 
    put system.tab (x.l('1','0',R)*scale):10:6 /
  );
  
  put regionresults_u;
  put // 'Y = ' Y.tl:4 '  p_a = ' p_a('1','0') //;
  loop(R, 
    put system.tab (u.l('1','0',R)*scale):10:6 /
  );  


  put regionresults_v;
  put // 'Y = ' Y.tl:4 '  p_a = ' p_a('1','0') //;
  loop(R, 
    put system.tab (v.l('1','0',R)*scale):10:6 /
  );

  put regionresults_w;
  put // 'Y = ' Y.tl:4 '  p_a = ' p_a('1','0') //;
  loop(R, 
    put system.tab (w.l('1','0')*scale):10:6 /
  );

$ontext
  loop(S,
    put // 'S = ' S.tl:4 '  obj = ' objval(S) //;
    put 'T/R':4;
    loop(R, 
      put system.tab R.tl:10
    );
    put /;
    loop(T,
      put T.tl:4;
      loop(R, 
        put system.tab (scale*z.l(S,T,R)):10:6);
      put /;
    );
  );
$offtext

* update based on the dynamics
* the updated x and z have the correct values from the nonanticipativity constraints
* the updated u and v would not satisfy the nonanticipativity constraints for the next round
  x.l(S,T,R)$(ord(T)<card(T)) = x.l(S,T+1,R);
  z.l(S,T,R)$(ord(T)<card(T)) = z.l(S,T+1,R);
* u.l(S,T,R)$(ord(T)<card(T)) = u.l(S,T+1,R);
* v.l(S,T,R)$(ord(T)<card(T)) = v.l(S,T+1,R);
* w.l(S,T)$(ord(T)<card(T)) = w.l(S,T+1);
);

putclose results;
putclose regionresults_z;
putclose regionresults_x;
putclose regionresults_u;
putclose regionresults_v;
putclose regionresults_w;