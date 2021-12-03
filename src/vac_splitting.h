// kt2 * dP/(dkt2 dx)
// Conventions: 
//   P(x) = R(x) / [(1-x)]_+ + A(x) + D(x) delta(1-x)
//   Pij : i->j
double Rqq(double x, double kt2);
double Dqq(double kt2);
double Aqg(double x, double kt2);
double Rgg(double x, double kt2);
double Agg(double x, double kt2);
double Dgg(double kt2, double nf);
double Agq(double x, double kt2);
