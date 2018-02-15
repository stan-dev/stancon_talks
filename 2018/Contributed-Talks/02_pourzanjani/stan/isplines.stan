functions {
  real I1(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -16) return(ret + -0.9999999999 * x^4/256 + -47.99999999 * x^3/192 + -767.9999999 * x^2/128 + -4095.9999998976 * x/64 - 254.999999926333);
	 else ret = ret + 1.00000006033338;

	 if(x >= -16) return(1);

	 return(ret);
}

real I2(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -16) return(ret + 448 * x^4/131072 + 22272 * x^3/98304 + 365568 * x^2/65536 + 1986560 * x/32768 - -246.875);
	 else ret = ret + 0.875;

	 if(x >= -16 && x < -12) return(ret + -0.9999999999 * x^4/2048 + -36 * x^3/1536 + -431.9999999 * x^2/1024 + -1727.9999999568 * x/512 - 10.00000002685);
	 else ret = ret + 0.124999987212497;

	 if(x >= -12) return(1);

	 return(ret);
}

real I3(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= -20 && x < -16) return(ret + -1408 * x^4/1179648 + -75264 * x^3/884736 + -1320960 * x^2/589824 + -7577600 * x/294912 - 107.638888888889);
	 else ret = ret + 0.361111111111086;

	 if(x >= -16 && x < -12) return(ret + 224 * x^4/294912 + 8832 * x^3/221184 + 112128 * x^2/147456 + 464896 * x/73728 - -20);
	 else ret = ret + 0.583333333333329;

	 if(x >= -12 && x < -8) return(ret + -0.9999999999 * x^4/4608 + -24 * x^3/3456 + -192 * x^2/2304 + -511.9999999872 * x/1152 - 0.83333333365);
	 else ret = ret + 0.0555555552388891;

	 if(x >= -8) return(1);

	 return(ret);
}

real I4(real x) { 
	 real ret = 0;
	 if(x < -20) return(0);

	 if(x >= -20 && x < -16) return(ret + 1 * x^4/6144 + 60 * x^3/4608 + 1200 * x^2/3072 + 8000 * x/1536 - -26.0416666666667);
	 else ret = ret + 0.0416666666666714;

	 if(x >= -16 && x < -12) return(ret + -3 * x^4/6144 + -132 * x^3/4608 + -1872 * x^2/3072 + -8383.9999999999 * x/1536 - 16.6666666666656);
	 else ret = ret + 0.458333333333599;

	 if(x >= -12 && x < -8) return(ret + 3 * x^4/6144 + 84 * x^3/4608 + 720 * x^2/3072 + 1984 * x/1536 - -3.125);
	 else ret = ret + 0.458333333333332;

	 if(x >= -8 && x < -4) return(ret + -0.9999999999 * x^4/6144 + -12 * x^3/4608 + -47.99999999 * x^2/3072 + -63.9999999984 * x/1536 - 2.66666633219614e-10);
	 else ret = ret + 0.0416666664520834;

	 if(x >= -4) return(1);

	 return(ret);
}

real I5(real x) { 
	 real ret = 0;
	 if(x < -16) return(0);

	 if(x >= -16 && x < -12) return(ret + 1 * x^4/6144 + 48 * x^3/4608 + 768 * x^2/3072 + 4096 * x/1536 - -10.6666666666667);
	 else ret = ret + 0.0416666666666643;

	 if(x >= -12 && x < -8) return(ret + -3 * x^4/6144 + -96 * x^3/4608 + -960 * x^2/3072 + -2816 * x/1536 - 2.875);
	 else ret = ret + 0.458333333333332;

	 if(x >= -8 && x < -4) return(ret + 3 * x^4/6144 + 48 * x^3/4608 + 192 * x^2/3072 + 256 * x/1536 - -0.666666666666666);
	 else ret = ret + 0.458333333333333;

	 if(x >= -4 && x < 0) return(ret + -0.9999999999 * x^4/6144 - -0.0416666666625);
	 else ret = ret + 0.0416666666625;

	 if(x >= 0) return(1);

	 return(ret);
}

real I6(real x) { 
	 real ret = 0;
	 if(x < -12) return(0);

	 if(x >= -12 && x < -8) return(ret + 1 * x^4/6144 + 36 * x^3/4608 + 432 * x^2/3072 + 1728 * x/1536 - -3.375);
	 else ret = ret + 0.0416666666666661;

	 if(x >= -8 && x < -4) return(ret + -3 * x^4/6144 + -60 * x^3/4608 + -336 * x^2/3072 + -30720 * x/147456 - -0.666666666666666);
	 else ret = ret + 0.458333333333333;

	 if(x >= -4 && x < 0) return(ret + 3 * x^4/6144 + 12 * x^3/4608 + -48 * x^2/3072 + 64 * x/1536 - -0.458333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 0 && x < 4) return(ret + -1 * x^4/6144 + 12 * x^3/4608 + -48 * x^2/3072 + 64 * x/1536 - 0);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 4) return(1);

	 return(ret);
}

real I7(real x) { 
	 real ret = 0;
	 if(x < -8) return(0);

	 if(x >= -8 && x < -4) return(ret + 1 * x^4/6144 + 24 * x^3/4608 + 192 * x^2/3072 + 512 * x/1536 - -0.666666666666667);
	 else ret = ret + 0.0416666666666666;

	 if(x >= -4 && x < 0) return(ret + -3 * x^4/6144 + -24 * x^3/4608 + 256 * x/1536 - -0.458333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 0 && x < 4) return(ret + 3 * x^4/6144 + -24 * x^3/4608 + 256 * x/1536 - 0);
	 else ret = ret + 0.458333333333333;

	 if(x >= 4 && x < 8) return(ret + -1 * x^4/6144 + 24 * x^3/4608 + -192 * x^2/3072 + 512 * x/1536 - 0.625);
	 else ret = ret + 0.0416666666666666;

	 if(x >= 8) return(1);

	 return(ret);
}

real I8(real x) { 
	 real ret = 0;
	 if(x < -4) return(0);

	 if(x >= -4 && x < 0) return(ret + 1 * x^4/6144 + 12 * x^3/4608 + 48 * x^2/3072 + 64 * x/1536 - -0.0416666666666667);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 0 && x < 4) return(ret + -3 * x^4/6144 + 12 * x^3/4608 + 48 * x^2/3072 + 64 * x/1536 - 0);
	 else ret = ret + 0.458333333333333;

	 if(x >= 4 && x < 8) return(ret + 3 * x^4/6144 + -60 * x^3/4608 + 336 * x^2/3072 + -320 * x/1536 - 0.208333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 8 && x < 12) return(ret + -1 * x^4/6144 + 36 * x^3/4608 + -432 * x^2/3072 + 1728 * x/1536 - 3.33333333333333);
	 else ret = ret + 0.0416666666666661;

	 if(x >= 12) return(1);

	 return(ret);
}

real I9(real x) { 
	 real ret = 0;
	 if(x < 0) return(0);

	 if(x >= 0 && x < 4) return(ret + 1 * x^4/6144 - 0);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 4 && x < 8) return(ret + -3 * x^4/6144 + 48 * x^3/4608 + -192 * x^2/3072 + 256 * x/1536 - 0.208333333333333);
	 else ret = ret + 0.458333333333333;

	 if(x >= 8 && x < 12) return(ret + 3 * x^4/6144 + -96 * x^3/4608 + 960 * x^2/3072 + -2816 * x/1536 - -3.33333333333333);
	 else ret = ret + 0.458333333333332;

	 if(x >= 12 && x < 16) return(ret + -1 * x^4/6144 + 48 * x^3/4608 + -768 * x^2/3072 + 4096 * x/1536 - 10.625);
	 else ret = ret + 0.0416666666666643;

	 if(x >= 16) return(1);

	 return(ret);
}

real I10(real x) { 
	 real ret = 0;
	 if(x < 4) return(0);

	 if(x >= 4 && x < 8) return(ret + 1 * x^4/6144 + -12 * x^3/4608 + 48 * x^2/3072 + -64 * x/1536 - -0.0416666666666667);
	 else ret = ret + 0.0416666666666667;

	 if(x >= 8 && x < 12) return(ret + -3 * x^4/6144 + 84 * x^3/4608 + -720 * x^2/3072 + 1984 * x/1536 - 2.66666666666667);
	 else ret = ret + 0.458333333333332;

	 if(x >= 12 && x < 16) return(ret + 3 * x^4/6144 + -132 * x^3/4608 + 1872 * x^2/3072 + -8383.9999999999 * x/1536 - -17.1249999999992);
	 else ret = ret + 0.458333333333599;

	 if(x >= 16 && x < 20) return(ret + -1 * x^4/6144 + 60 * x^3/4608 + -1200 * x^2/3072 + 8000 * x/1536 - 26);
	 else ret = ret + 0.0416666666666714;

	 if(x >= 20) return(1);

	 return(ret);
}

real I11(real x) { 
	 real ret = 0;
	 if(x < 8) return(0);

	 if(x >= 8 && x < 12) return(ret + 1 * x^4/4608 + -24 * x^3/3456 + 192 * x^2/2304 + -512 * x/1152 - -0.888888888888889);
	 else ret = ret + 0.0555555555555558;

	 if(x >= 12 && x < 16) return(ret + -224 * x^4/294912 + 8832 * x^3/221184 + -112128 * x^2/147456 + 464896 * x/73728 - 19.4166666666667);
	 else ret = ret + 0.583333333333329;

	 if(x >= 16 && x < 20) return(ret + 1408 * x^4/1179648 + -75264 * x^3/884736 + 1320960 * x^2/589824 + -7577600 * x/294912 - -108);
	 else ret = ret + 0.361111111111086;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20) return(1);

	 return(ret);
}

real I12(real x) { 
	 real ret = 0;
	 if(x < 12) return(0);

	 if(x >= 12 && x < 16) return(ret + 1 * x^4/2048 + -36 * x^3/1536 + 432 * x^2/1024 + -1728 * x/512 - -10.125);
	 else ret = ret + 0.125;

	 if(x >= 16 && x < 20) return(ret + -448 * x^4/131072 + 22272 * x^3/98304 + -365568 * x^2/65536 + 1986560 * x/32768 - 246);
	 else ret = ret + 0.875;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20) return(1);

	 return(ret);
}

real I13(real x) { 
	 real ret = 0;
	 if(x < 16) return(0);

	 if(x >= 16 && x < 20) return(ret + 1 * x^4/256 + -48 * x^3/192 + 768 * x^2/128 + -4095.9999999999 * x/64 - -255.999999999975);
	 else ret = ret + 1.00000000000637;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20 && x < 20) return(ret + 0 - 0);
	 else ret = ret + 0;

	 if(x >= 20) return(1);

	 return(ret);
}

real[] ispline(vector x, matrix c) {

	 int N;
	 real ret[dims(x)[1]];
	 N = dims(x)[1];

	 for(n in 1:N) ret[n] = c[n][1]*I1(x[n]) + c[n][2]*I2(x[n]) + c[n][3]*I3(x[n]) + c[n][4]*I4(x[n]) + c[n][5]*I5(x[n]) + c[n][6]*I6(x[n]) + c[n][7]*I7(x[n]) + c[n][8]*I8(x[n]) + c[n][9]*I9(x[n]) + c[n][10]*I10(x[n]) + c[n][11]*I11(x[n]) + c[n][12]*I12(x[n]) + c[n][13]*I13(x[n]);

	 return(ret);
}
}
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  simplex[13] c;
}
model {
  y ~ normal(ispline(x, rep_matrix(to_row_vector(c), N)), 0.05);
}
generated quantities {
  vector[81] seq;
  real fhat[81];
  for(n in 1:81) seq[n] = (n-1)*0.5 - 20;
  fhat = ispline(seq, rep_matrix(to_row_vector(c), 101));
}
