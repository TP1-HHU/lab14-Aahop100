#include "EM.hxx"
#include <cmath>
#include <iostream>
//----------------------------------
using namespace std;
//----------------------------------
class particle{
public:
	particle(){x=0; y=0; z=0;  // location at t = 0
		         px=0; py=0; pz=0;  // momentum at t = -dt/2
						 gam=0;}  // gamma at t= -dt/2
	// Print particle position, momentum u and gamma-factor
	void print_info(const double t) const;
	void boris_push(const double Ex, const double Ey, const double Ez,
			            const double Bx, const double By, const double Bz,
	                const double dt);
  double get_x() const{return x;}
private:
	double x,y,z;
	double px,py,pz;
  double gam;
	// Calculate gamma-factor from current momentum
  void calc_gam();
};
//----------------------------------
//----------------------------------
int main(){
  const double FWHM = 10;
  const double a0 = 1.44;
  const double T0 = 200;
  const double Tend = 800;
  const double dt = 0.001;
  const int N = int(Tend/dt + 0.5);

  double Ex,Ey,Ez, Bx, By, Bz;
  double t=0;

  EM em(a0, FWHM,T0);
  particle p;

  for(int i=0; i<N; i++){
	  em.fields(Ex,Ey,Ez,Bx,By,Bz, t, p.get_x());
	  p.boris_push(Ex, Ey, 0, Bx, 0, Bz, dt);
	  if (i%5 == 0 ) p.print_info(t);
	  t += dt;
  }
  return 0;
}
//----------------------------------
void particle::calc_gam(){
	gam = sqrt(1.0 + px*px + py*py + pz*pz);
}
//----------------------------------
void particle::boris_push(const double Ex, const double Ey, const double Ez,
		                      const double Bx, const double By, const double Bz,
													const double dt)
{
	// p Update nach Boris
	// p(t-dt/2) -> p(t+dt/2)
  const double dth = 0.5 * dt;
	px +=  Ex * dth;
	py +=  Ey * dth;
	pz +=  Ez * dth;

  calc_gam();

	double tx = Bx * dth / gam;
	double ty = By * dth / gam;
	double tz = Bz * dth / gam;

	double psx = py * tz - pz * ty;
	double psy = pz * tx - px * tz;
	double psz = px * ty - py * tx;

	psx += px;
	psy += py;
	psz += pz;

	double n = 1 + tx*tx + ty*ty + tz*tz;
	double sx = 2*tx /n;
	double sy = 2*ty /n;
	double sz = 2*tz /n;

  double ppx = psy * sz - psz * sy;
  double ppy = psz * sx - psx * sz;
  double ppz = psx * sy - psy * sx;

  ppx += px;
  ppy += py;
  ppz += pz;

  px = ppx + Ex * dth;
  py = ppy + Ey * dth;
  pz = ppz + Ez * dth;

  // x Update nach Leap-Frog
  calc_gam();

  x += px / gam * dt;
  y += py / gam * dt;
  z += pz / gam * dt;
}
//--------------------------------------
void particle::print_info(const double t) const
{
	cout << t << "\t" << x << "\t" << y << "\t" << z
			 << "\t" << px << "\t" << py << "\t" << pz << "\t" << gam << endl;
}
