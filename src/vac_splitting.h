// kt2 * dP/(dkt2 dx)
// Conventions: 
//   P(x) = R(x) / [(1-x)]_+ + A(x) + D(x) delta(1-x)
//   Pij : i->j
double Rqq(double x, double kt2);
double Dqq(double kt2);
double Aqg(double x, double kt2);
double Rgg(double x, double kt2);
double Agg(double x, double kt2);
double Dgg(double kt2);
double Agq(double x, double kt2);
double Dgq(double kt2);
// (kt2+1mx2m2) * dP/(d(kt2+1nx2m2) dx)
double RQQ(double x, double M2overQ2);
double AQQ(double x, double M2overQ2);
double DQQ(double M2overQ2);
// (kt2+m2) * dP/(d(kt2+m2) dx)
double AgQ(double x, double M2overQ2);
double DgQ(double M2overQ2);
// (kt2+x2m2) * dP/(d(kt2+x2m2) dx)
double AQg(double x, double M2overQ2);
