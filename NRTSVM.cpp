/*C++ CODE - MANGEAT MATTHIEU - 2026*/
/*NON RECIPROCAL TWO SPECIES VICSEK MODEL*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

using namespace std;

//Personal libraries
#include "lib/random_OMP.cpp"
#include "lib/special_functions.cpp"

template <typename T>
T modulo(const T &x, const T &LX);

///////////////////////////////
///// OPENMP VECTOR CLASS /////
///////////////////////////////

template <typename T>
class omp_vector : public vector<T>
{
	private:
	
	omp_lock_t lock;
	
	public:
	
	omp_vector();
	void atomic_push_back(const T &p);
};

template <typename T>
omp_vector<T>::omp_vector()
{
	omp_init_lock(&lock);
}

template <typename T>
void omp_vector<T>::atomic_push_back(const T &p)
{
	omp_set_lock(&lock);
	vector<T>::push_back(p);
	omp_unset_lock(&lock);
}

//////////////////////////
///// PARTICLE CLASS /////
//////////////////////////

class particle
{
	public:
	
	double x,y; //Position of the particle.
	double dx, dy; //Displacement of the particle.
	double theta, ctheta, stheta; //Orientation of the particle.
	double theta_avg; //Average orientation in the neighborhood.
	double omega; //angular velocity.
	int species; //Species = 0 or 1

	particle(const int &LX, const int &LY, const int &species0, const int &init); //Initial position and orientation.
	int band_formation(const double &x0, const double &alpha, const double &kappa, const int &LX);
	void update(const double &v0, const double &eta, const int &LX, const int &LY); //Update position and orientation.
};

//Creation of bands.
int particle::band_formation(const double &x0, const double &alpha, const double &kappa, const int &LX)
{
	const double r=ran();
	double XX=0.;
	if (r<kappa)
	{
		XX=modulo(x0+alpha*r/kappa,1.);
	}
	else
	{
		XX=modulo(x0+(1-alpha)*(r-kappa)/(1-kappa)+alpha,1.);
		theta=-M_PI+2*M_PI*ran(); //Gaseous phase -> change the spin to random.
	}
	return int(XX*LX);
}


//Creation of particles.
particle::particle(const int &LX, const int &LY, const int &species0, const int &init)
{
	static const double alpha=0.125, kappa=0.8;
	
	//Spin of the particle.
	species=species0;
	omega=0;
	
	//Displacement.
	dx=0.;
	dy=0.;
	
	if (init==0) //APF state.
	{
		theta=M_PI*species;
		//x=0.5*LX*(1-species)+LX*ran()/8;
		x=band_formation(0.5*species,alpha,kappa,LX);
		y=LY*ran();
	}
	else if (init==1) //PF state.
	{
		theta=0;
		//x=0.5*LX*(1-species)+LX*ran()/8;
		x=band_formation(0.5*species,alpha,kappa,LX);
		y=LY*ran();
	}
	else if (init==2) //Chiral orientation.
	{
		theta=M_PI/2.*species;
		x=LX*ran();
		y=LY*ran();
	}
	else if (init==3) //Disordered configuration.
	{
		theta=-M_PI+2*M_PI*ran();
		x=LX*ran();
		y=LY*ran();
	}
	else
	{
		cerr << "BAD INIT VALUE: " << init << endl;
		abort();
	}
	
	ctheta=cos(theta);
	stheta=sin(theta);
}

void particle::update(const double &v0, const double &eta, const int &LX, const int &LY)
{
	//Update angle.
	omega=theta;
	
	theta=theta_avg + 2*M_PI*eta*(ran()-0.5);
	if (theta<-M_PI) theta+=2*M_PI;
	else if (theta>=M_PI) theta-=2*M_PI;
	
	ctheta=cos(theta);
	stheta=sin(theta);
	
	omega=theta-omega;
	if (omega<-M_PI) omega+=2*M_PI;
	else if (omega>=M_PI) omega-=2*M_PI;
	
	//Update the displacement.
	dx += v0*ctheta;
	dy += v0*stheta;
	
	//Update the position with the new angle (+periodic condition).
	x += v0*ctheta;
	y += v0*stheta;	
	
	if (x>=LX) x-=LX;
	else if (x<0) x+=LX;
		
	if (y>=LY) y-=LY;
	else if (y<0) y+=LY;
}

////////////////////////////////
///// AVERAGES IN SUBBOXES /////
////////////////////////////////



///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Modulo function.
template <typename T>
T modulo(const T &x, const T &L)
{
	if (x<0)
	{
		return x+L;
	}
	else if (x>=L)
	{
		return x-L;
	}
	else
	{
		return x;
	}
}

//Distance between two particles in periodic space.
double distance(const particle &part1, const particle &part2, const int &LX, const int &LY)
{
	const double DX=fabs(part1.x-part2.x);
	const double DY=fabs(part1.y-part2.y);	
	return square(min(DX,LX-DX)) + square(min(DY,LY-DY));
}

//Order parameter VS.
vector<double> VS(const int &Npart, const vector<particle> &PART)
{
	vector<double> vs(2,0.);
	for (int i=0; i<Npart; i++)
	{
		vs[0]+=PART[i].ctheta;
		vs[1]+=PART[i].stheta;
	}
	vs[0]/=Npart;
	vs[1]/=Npart;
	return vs;	
}

//Order parameter VA.
vector<double> VA(const int &Npart, const vector<particle> &PART)
{
	vector<double> va(2,0.);
	for (int i=0; i<Npart; i++)
	{
		const int spin=1-2*PART[i].species;
		va[0]+=spin*PART[i].ctheta;
		va[1]+=spin*PART[i].stheta;
	}
	va[0]/=Npart;
	va[1]/=Npart;
	return va;	
}

//mean orientation of S species
double THETAS(const int &Npart, const vector<particle> &PART, const int &s)
{
	vector<double> mag(2,0.);
	for (int i=0; i<Npart; i++)
	{
		if (PART[i].species==s)
		{
			mag[0]+=PART[i].ctheta;
			mag[1]+=PART[i].stheta;
		}
	}
	return atan2(mag[1],mag[0]);
}

//turning activity of S species
double OMEGAS(const int &Npart, const vector<particle> &PART, const int &s)
{
	double omega=0;
	for (int i=0; i<Npart; i++)
	{
		if (PART[i].species==s)
		{
			omega+=sin(PART[i].omega);
		}
	}
	return 2*omega/Npart;
}



//Mean-square displacement.
double MSD(const int &Npart, const vector<particle> &PART)
{
	double DX2=0., DY2=0.;
	for (int i=0; i<Npart; i++)
	{
		DX2+=PART[i].dx*PART[i].dx/Npart;
		DY2+=PART[i].dy*PART[i].dy/Npart;
	}	
	return DX2+DY2;
}

//Average angle near each particles at time t.
void THETA_AVG(const vector< vector<double> > &Jcoup, const int &LX, const int &LY, const int &Npart, vector<particle> &PART)
{
	//Create a matrix with particle locations, box[i][j] regroups the particle indices with i<x<i+1 and j<y<j+1.
	vector< vector< omp_vector<int> > > box(LX);
	for (int x = 0; x < LX; x++)
	{
		box[x].resize(LY);
	}
	
	#pragma omp parallel for default(shared) schedule(static)
	for (int i=0; i<Npart; i++)
	{
		box[int(PART[i].x)][int(PART[i].y)].atomic_push_back(i);
	}
	
	#pragma omp parallel for collapse(2) default(shared) schedule(static)
	for (int x = 0; x < LX; x++)
	{
		for (int y = 0; y < LY; y++)
		{
			sort(box[x][y].begin(), box[x][y].end());
		}
	}
	
	//Calculate the average angle of neighbours.
	#pragma omp parallel for default(shared) schedule(static)
	for (int i=0; i<Npart; i++)
	{
		const int X0=int(PART[i].x), Y0=int(PART[i].y);
		
		//Take the average magnetisation for particles in neighbour boxes and with a distance smaller than 1.	
		double MX=0., MY=0.;		
		for (int XN=X0-1;XN<=X0+1;XN++)
		{
			for (int YN=Y0-1;YN<=Y0+1;YN++)
			{
				const omp_vector<int> &neighbours=box[modulo(XN,LX)][modulo(YN,LY)];
				for (int l=0; l<neighbours.size(); l++)
				{
					const int j=neighbours[l];
					if (distance(PART[i],PART[j],LX,LY)<1)
					{
						double J=Jcoup[PART[i].species][PART[j].species];
						MX+=J*PART[j].ctheta;
						MY+=J*PART[j].stheta;
					}
				}
			}
		}
		PART[i].theta_avg=atan2(MY,MX);
	}
}

//Update the number fluctuations at time t.
void updateFluctuations(vector<double> &n, vector<double> &n2, const int &L0, const int &rx, const int &ry, const int &LX, const int &LY, const int &Npart, const vector<particle> &PART)
{
	//Density of each sector at time t
	vector<int> RHO=vector<int>(LX*LY,0);
	for (int i=0; i<Npart; i++)
	{
		RHO[int(PART[i].x)*LY+int(PART[i].y)]++;
	}
	
	static const int Nboxes=10;
	//Averages on all sub-boxes
	#pragma omp parallel for default(shared) schedule(static)
	for (int l=1; l<L0; l++)
	{
		//Select Nboxes different sub-boxes of size l
		for (int nbox=0; nbox<Nboxes; nbox++)
		{
			//Random position for the bottom left corner
			const int x0=int((LX-rx*l)*ran2());
			const int y0=int((LY-ry*l)*ran2());
			
			//Density and magnetization in this sub-box
			long unsigned int rho=0;
			
			for (int x=x0; x<x0+rx*l; x++)
			{
				for (int y=y0; y<y0+ry*l; y++)
				{
					rho+=RHO[x*LY+y];
				}
			}
			
			//Add to n,n2
			n[l]+=double(rho)/Nboxes;
			n2[l]+=double(rho*rho)/Nboxes;
		}
	}
}

//Export the number fluctuations in a file.
void exportFluctuations(const vector<double> &n, const vector<double> &n2, const int &L0, const int &Nav, const double &v0, const double &eta, const double &rho0, const double &JAB, const double &JBA, const int &LX, const int &LY, const int &init, const int &RAN)
{
	int returnSystem=system("mkdir -p data_TSVM_fluctuations/");
	stringstream ss;	
	ss << "./data_TSVM_fluctuations/TSVM_fluctuations_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";	
	string nameAV = ss.str();			
	ofstream fileAV(nameAV.c_str(),ios::trunc);	
	fileAV.precision(6);
	
	fileAV << "#The number of averages are: " << Nav << endl;

	for (int l=1; l<L0; l++)
	{
		double N=n[l]/Nav;
		double delN=n2[l]/Nav-N*N;
		fileAV << l << "\t" << N << "\t" << delN << endl;
	}
	
	fileAV.close();
}

//Update proba of omega.
void updateProba(vector<double> &PA, vector<double> &PB, const int &Nomega, const double &domega, const int &Npart, const vector<particle> &PART)
{
	const int nthreads = omp_get_max_threads();
	vector<vector<int>> pA(nthreads, vector<int>(Nomega, 0)), pB(nthreads, vector<int>(Nomega, 0));

	#pragma omp parallel for default(shared) schedule(static)
	for (int i=0; i<Npart; i++)
	{
		const int tid = omp_get_thread_num();
		const int k=int((M_PI+PART[i].omega)/domega);
		if (PART[i].species==0)
		{
			pA[tid][k]++;
		}
		else
		{
			pB[tid][k]++;
		}
	}
	
	#pragma omp parallel for default(shared) schedule(static)
	for (int k=0; k<Nomega; k++)
	{
		for (int tid=0; tid<nthreads; tid++)
		{
			PA[k]+=2*double(pA[tid][k])/Npart;
			PB[k]+=2*double(pB[tid][k])/Npart;
		}
	}
}

//Export proba of omega in a file.
void exportProba(const vector<double> &PA, const vector<double> &PB, const int &Nomega, const double &domega, const int &Nav, const double &v0, const double &eta, const double &rho0, const double &JAB, const double &JBA, const int &LX, const int &LY, const int &init, const int &RAN)
{
	int returnSystem=system("mkdir -p data_TSVM_proba_omega/");
	stringstream ss;	
	ss << "./data_TSVM_proba_omega/TSVM_proba_omega_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";	
	string nameAV = ss.str();			
	ofstream fileAV(nameAV.c_str(),ios::trunc);	
	fileAV.precision(6);
	
	fileAV << "#The number of averages are: " << Nav << endl;

	for (int k=0; k<Nomega; k++)
	{
		fileAV << -M_PI+(k+0.5)*domega << "\t" << PA[k]/(Nav*domega) << "\t" << PB[k]/(Nav*domega) << endl;
	}
	
	fileAV.close();
}

//Update correlation function.
void updateCorr(vector<double> &CORR, const int &Rmax, const int &LX, const int &LY, const int &Npart, const vector<particle> &PART)
{
	const int NX=2*LX, NY=2*LY;
	const int Ngrid=NX*NY;
	vector<double> MAx(Ngrid,0), MAy(Ngrid,0), MBx(Ngrid,0), MBy(Ngrid,0);
	vector<int> hasA(Ngrid, 0), hasB(Ngrid, 0);
	
	for (int i=0; i<Npart; i++)
	{
		const int k=int(2*PART[i].x)*NY+int(2*PART[i].y);
		if (PART[i].species==0)
		{
			MAx[k]+=PART[i].ctheta;
			MAy[k]+=PART[i].stheta;
			hasA[k]=1;
		}
		else
		{
			MBx[k]+=PART[i].ctheta;
			MBy[k]+=PART[i].stheta;
			hasB[k]=1;
		}
	}
	
	vector<double> PHI(Ngrid,0);
	vector<int> valid(Ngrid, 0);
	
	for (int k=0; k<Ngrid; k++)
	{
		if (hasA[k] and hasB[k])
		{
			PHI[k]=atan2(MAy[k],MAx[k])-atan2(MBy[k],MBx[k]); //thetaA-thetaB without periodicity.
			valid[k]=1;
		}
	}
	
	#pragma omp parallel for schedule(static)
	for (int r=0; r<Rmax+1; r++) //delta r.
	{
		double corr=0.;
		int count=0;
		for (int k=0; k<Ngrid; k++) //init position (averaged on it)
		{
			if (valid[k])
			{
				const int x=k/NY, y=k%NY;
				const int kxp=k+r*NY - (x>=NX-r)*Ngrid;
				const int kxm=k-r*NY + (x<r)*Ngrid;
				const int kyp=k+r - (y>=NY-r)*NY;
				const int kym=k-r + (y<r)*NY;
				
				if (valid[kxp])
				{
					corr+=cos(PHI[kxp]-PHI[k]);
					count++;
				}
				if (valid[kxm])
				{
					corr+=cos(PHI[kxm]-PHI[k]);
					count++;
				}
				if (valid[kyp])
				{
					corr+=cos(PHI[kyp]-PHI[k]);
					count++;
				}
				if (valid[kym])
				{
					corr+=cos(PHI[kym]-PHI[k]);
					count++;
				}
			}
		}
		if (count>0)
		{
			CORR[r]+=corr/count;
		}
	}
}

//Export correlation function in a file.
void exportCorr(const vector<double> &CORR, const int &Rmax, const int &Nav, const double &v0, const double &eta, const double &rho0, const double &JAB, const double &JBA, const int &LX, const int &LY, const int &init, const int &RAN)
{
	int returnSystem=system("mkdir -p data_TSVM_correlation/");
	stringstream ss;	
	ss << "./data_TSVM_correlation/TSVM_correlation_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";	
	string nameAV = ss.str();			
	ofstream fileAV(nameAV.c_str(),ios::trunc);	
	fileAV.precision(6);
	
	fileAV << "#The number of averages are: " << Nav << endl;

	for (int r=0; r<Rmax+1; r++)
	{
		fileAV << double(r)/Rmax << "\t" << CORR[r]/Nav << endl; //plot C(2r/L) in real coordinates.
	}
	
	fileAV.close();
}

//Export densities and mean orientations in a file.
void exportDensity(const double &v0, const double &eta, const double &rho0, const double &JAB, const double &JBA, const int &LX, const int &LY, const int &t, const int &init, const int &RAN, const int &Npart, const vector<particle> &PART)
{
	const int Nsites=LX*LY;
	
	//Local density of each species.
	vector<int> rhoA(Nsites,0), rhoB(Nsites,0);
	vector<double> mAx(Nsites,0.), mAy(Nsites,0.), mBx(Nsites,0.), mBy(Nsites,0.);
	
	for (int i=0; i<Npart; i++)
	{
		if (PART[i].species==0)
		{
			rhoA[int(PART[i].y)*LX+int(PART[i].x)]++;
			mAx[int(PART[i].y)*LX+int(PART[i].x)]+=PART[i].ctheta;
			mAy[int(PART[i].y)*LX+int(PART[i].x)]+=PART[i].stheta;
		}
		else
		{
			rhoB[int(PART[i].y)*LX+int(PART[i].x)]++;
			mBx[int(PART[i].y)*LX+int(PART[i].x)]+=PART[i].ctheta;
			mBy[int(PART[i].y)*LX+int(PART[i].x)]+=PART[i].stheta;
		}
	}
	
	
	vector<int16_t> rho(Nsites,0);
	vector<float> thetaA(Nsites,0), thetaB(Nsites,0);
	
	for (int i=0;i<Nsites;i++)
	{
		if (rhoA[i]>rhoB[i] or (rhoA[i]==rhoB[i] and ran_hash()<0.5))
		{
			rho[i]=rhoA[i];
		}
		else
		{
			rho[i]=-rhoB[i];
		}
		if (rhoA[i]==0)
		{
			thetaA[i]=NAN;
		}
		else
		{
			thetaA[i]=atan2(mAy[i],mAx[i]);
		}
		if (rhoB[i]==0)
		{
			thetaB[i]=NAN;
		}
		else
		{
			thetaB[i]=atan2(mBy[i],mBx[i]);
		}
	}
	
	//Creation of the file.
	static const int dossier=system("mkdir -p ./data_TSVM_dynamics/");
	
	stringstream ssDensity;
	ssDensity << "./data_TSVM_dynamics/TSVM_rho_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".bin";
	string nameDensity = ssDensity.str();
	ofstream fileDensity(nameDensity.c_str(),ios::binary);
	fileDensity.write(reinterpret_cast<const char*>(rho.data()), Nsites*sizeof(int16_t));
	fileDensity.close();
	
	stringstream ssThetaA;
	ssThetaA << "./data_TSVM_dynamics/TSVM_thetaA_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".bin";
	string nameThetaA = ssThetaA.str();
	ofstream fileThetaA(nameThetaA.c_str(),ios::binary);
	fileThetaA.write(reinterpret_cast<const char*>(thetaA.data()), Nsites*sizeof(float));
	fileThetaA.close();
	
	stringstream ssThetaB;
	ssThetaB << "./data_TSVM_dynamics/TSVM_thetaB_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".bin";
	string nameThetaB = ssThetaB.str();
	ofstream fileThetaB(nameThetaB.c_str(),ios::binary);
	fileThetaB.write(reinterpret_cast<const char*>(thetaB.data()), Nsites*sizeof(float));
	fileThetaB.close();
}

//Export particle positions and orientations in a file.
void exportParticles(const double &v0, const double &eta, const double &rho0, const double &JAB, const double &JBA, const int &LX, const int &LY, const int &t, const int &init, const int &RAN, const int &Npart, const vector<particle> &PART)
{
	//Creation of the file.
	static const int dossier=system("mkdir -p ./data_TSVM_particles/");
	stringstream ssPart;
	ssPart << "./data_TSVM_particles/TSVM_particles_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".bin";
	string namePart = ssPart.str();
	ofstream filePart(namePart.c_str(),ios::binary);
	
	//Write in the file.
	for (int i=0;i<Npart;i++)
	{
		float part[3]={static_cast<float>(PART[i].x),static_cast<float>(PART[i].y),static_cast<float>(PART[i].theta)};
		filePart.write(reinterpret_cast<const char*>(part), 3*sizeof(float));
	}
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &v0, double &eta, double &rho0, double &JAB, double &JBA, int &LX, int &LY, int &tmax, int &init, int &RAN, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		
		if (strstr(argv[i], "-v0=" ))
		{
			v0=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-eta=" ))
		{
			eta=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-JAB=" ))
		{
			JAB=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-JBA=" ))
		{
			JBA=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
		}
		else if (strstr(argv[i], "-threads=" ))
		{
			THREAD_NUM=atoi(argv[i]+9);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: v0=self-propulsion velocity, eta=noise parameter, rho0=average density, JAB=coupling of A feeling B, JBA=coupling of B feeling A, LX*LY=size of the box.
	double v0=0.0015, eta=0.02, rho0=64., JAB=0.25, JBA=-0.25;
	int LX=12, LY=12;
	
	//Numerical parameters: init=initial condition, tmax=maximal time, RAN=index of RNG, THREAD_NUM=number of threads.
	int init=3, tmax=70000, RAN=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,v0,eta,rho0,JAB,JBA,LX,LY,tmax,init,RAN,THREAD_NUM);
	
	//OpenMP.
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;
	
	//Number of particles.
	const int Npart=int(rho0*LX*LY);
	
	//Coupling constant.
	const vector< vector<double>> Jcoup={{1.,JAB},{JBA,1.}};

	//Start the random number generator.
	init_gsl_ran();
	cout << "GSL index (initial time) = " << RAN << "\n";
	for (int k=0; k<THREAD_NUM; k++)
	{
		gsl_rng_set(GSL_r[k],THREAD_NUM*RAN+k);
		gsl_rng_set(GSL_r2[k],THREAD_NUM*(RAN+1234567)+k);
	}
	//set_hash_seed(12345u * RAN);
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_TSVM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_TSVM_averages/TSVM_averages_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	//Initial state.
	vector<particle> PART;
	for (int k=0; k<Npart/2; k++)
	{
		particle partplus(LX,LY,0,init);
		PART.push_back(partplus);
		particle partminus(LX,LY,1,init);
		PART.push_back(partminus);
	}
	
	const double teq=10000, texp=tmax-2000;
	
	//Averages.
	vector<double> DTHETA(0), OMEGA(0), VSNORM(0), VANORM(0), MSD2(0);
	
	//Number fluctuations.
	int L0=min(LX,LY), rx=max(1,LX/LY), ry=max(1,LY/LX);
	vector<double> n(L0,0.), n2(L0,0.);
	
	//Probability of omega.
	const int Nomega=1001;
	const double domega=2*M_PI/Nomega;
	vector<double> PA(Nomega,0.), PB(Nomega,0.);
	
	//Correlation function.
	const int Rmax=min(LX,LY);
	vector<double> CORR(Rmax+1,0.);
	
	int Nav=0;
	bool breakdown=false;
	
	//Time evolution.
	for(int t=0; t<=tmax; t++)
	{
		//Export data.
		const vector<double> vs=VS(Npart,PART), va=VA(Npart,PART); //Order parameter VS and VA (sometime useful).
		const double vs_norm=sqrt(square(vs[0])+square(vs[1])), va_norm=sqrt(square(va[0])+square(va[1])); //norm of VS and VA.
		
		const double msd=MSD(Npart,PART); //Mean-square displacement.
		
		double dtheta=THETAS(Npart,PART,0)-THETAS(Npart,PART,1); //Averaged angles.
		if (dtheta<-M_PI) dtheta+=2*M_PI;
		else if (dtheta>=M_PI) dtheta-=2*M_PI;
	
		const double omega=0.5*(OMEGAS(Npart,PART,0)+OMEGAS(Npart,PART,1)); //Turning activity;
		
		fileAverages << t << " " << dtheta << " " << omega << " " <<  vs_norm << " " << va_norm << " " << msd << endl;
		
		DTHETA.push_back(dtheta);
		OMEGA.push_back(omega);
		VSNORM.push_back(vs_norm);
		VANORM.push_back(va_norm);
		MSD2.push_back(msd);
		
		if (t>=texp)
		{
			exportParticles(v0,eta,rho0,JAB,JBA,LX,LY,t,init,RAN,Npart,PART);
			exportDensity(v0,eta,rho0,JAB,JBA,LX,LY,t,init,RAN,Npart,PART);
		}
		
		if (t>teq)
		{
			updateFluctuations(n,n2,L0,rx,ry,LX,LY,Npart,PART);
			updateProba(PA,PB,Nomega,domega,Npart,PART);
			updateCorr(CORR,Rmax,LX,LY,Npart,PART);
			Nav++;
		}
		
		if (t%(tmax/1000)==0 or t==tmax)
		{
			cout << double(100*t)/double(tmax) << "% steps performed -t=" << t << " -dtheta=" << dtheta << " -Omega=" << omega << " -vs=" << vs_norm << " -va=" << va_norm << " -MSD=" << msd << running_time.TimeRun(" ") << endl;
			if (Nav>0)
			{
				exportFluctuations(n,n2,L0,Nav,v0,eta,rho0,JAB,JBA,LX,LY,init,RAN);
				exportProba(PA,PB,Nomega,domega,Nav,v0,eta,rho0,JAB,JBA,LX,LY,init,RAN);
				exportCorr(CORR,Rmax,Nav,v0,eta,rho0,JAB,JBA,LX,LY,init,RAN);
			}
		}
		
		//UNCOMMENT TO CHARACTERIZE THE BREAKING TIME EFFICIENTLY.
		/*if (fabs(cos(dtheta))>0.9 and not breakdown)
		{
			breakdown=true;
			//tmax=t+1000;
			tmax=1000*(int((t+500)/1000)+1);
		}*/
		
		//Calculate the average orientation in the neighborhood of each particles.
		THETA_AVG(Jcoup,LX,LY,Npart,PART);
		
		//Update orientations and positions sequentially.
		#pragma omp parallel for default(shared) schedule(static)
		for (int i=0; i<Npart; i++)
		{
			PART[i].update(v0,eta,LX,LY);			
		}
	}
	
	stringstream ssDTHETA;
	ssDTHETA << "./data_TSVM_averages/TSVM_dtheta_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".bin";
	string nameDTHETA = ssDTHETA.str();
	ofstream fileDTHETA(nameDTHETA.c_str(),ios::binary);
	fileDTHETA.write(reinterpret_cast<const char*>(DTHETA.data()), DTHETA.size()*sizeof(double));
	fileDTHETA.close();
	
	stringstream ssOMEGA;
	ssOMEGA << "./data_TSVM_averages/TSVM_omega_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".bin";
	string nameOMEGA = ssOMEGA.str();
	ofstream fileOMEGA(nameOMEGA.c_str(),ios::binary);
	fileOMEGA.write(reinterpret_cast<const char*>(OMEGA.data()), OMEGA.size()*sizeof(double));
	fileOMEGA.close();

	stringstream ssVS;
	ssVS << "./data_TSVM_averages/TSVM_vs_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".bin";
	string nameVS = ssVS.str();
	ofstream fileVS(nameVS.c_str(),ios::binary);
	fileVS.write(reinterpret_cast<const char*>(VSNORM.data()), VSNORM.size()*sizeof(double));
	fileVS.close();

	stringstream ssVA;
	ssVA << "./data_TSVM_averages/TSVM_va_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".bin";
	string nameVA = ssVA.str();
	ofstream fileVA(nameVA.c_str(),ios::binary);
	fileVA.write(reinterpret_cast<const char*>(VANORM.data()), VANORM.size()*sizeof(double));
	fileVA.close();

	stringstream ssMSD;
	ssMSD << "./data_TSVM_averages/TSVM_msd_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".bin";
	string nameMSD = ssMSD.str();
	ofstream fileMSD(nameMSD.c_str(),ios::binary);
	fileMSD.write(reinterpret_cast<const char*>(MSD2.data()), MSD2.size()*sizeof(double));
	fileMSD.close();
	
	return 0;
}
