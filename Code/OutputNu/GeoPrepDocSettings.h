#pragma once

#include <stdlib.h>
#include <direct.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

struct GeoPrepDocSettings
{

	bool auto_mode;
	string model_file_path;
	string model_file_name;
	

	double rect_bak;
	double rect_hmin;
	double rect_coeff;

	double rectsp_bak;
	double rectsp_hmin;
	double rectsp_coeff;

	double rectloop_bak;
	double rectloop_rcoeff;
	double rectloop_coeff;

	double impulse_tcoeff;
	double impulse_t;

	int line_mrad_coeff;
	int line_infc;
	int line_infm;

	int time_gaps_for_loop;

	vector<double> logscale;

	int polygon_mrad_coeff;

	bool enable_recalculations;

	double slae3d_eps;

	bool enable_2drecalculations;

	bool vfem_direct_output;

	bool meshed;

	bool emf_in_line;

	bool enable_pariterscheme;
	double epsgfp;

	bool un3dmg;

	bool optXZ;
	bool optAZ;
	bool optUZ;

	double Cxyz[3];
	double C2xyz[3];
	int Cplmax;
	double Sxyz[3];

	bool objframe;

	int N_Periods;
	int N_StepsPerPeriod;

	int N_Vect;
	int MaxIter;

	bool enable_block_relaxation;
	bool flagRunGMRES;
	bool enable_smooth_emf2d_on_z0;

	vector<int> mesh3d_h2;

	double smooth_alpha;

	bool smooth1;

	int line_coeff;

	int FindLLtPrecond;

	bool GPMessages;

	bool OptBaseMesh2d;

	bool b_in_line;

	bool line_zones;

	double slae2d_eps;

	bool ced_cr;

	bool calcE;
	int nLayerE;

	double OmegaMax;
	double FreqMed;

	int MaxIterBlockRelax;

	bool Resultant;

	bool vX; double cX;
	bool vY; double cY;
	bool vZ; double cZ;

	bool ObjectTitles;

	int MaxIterFit;
	int MaxIterReg;

	bool ForInv;

	char PathTo2d[256];

private:

	string pws;
	char pwk[4];
	char GetCharFromPws(unsigned int i)
	{
		int c;
		char ch[2];
		i-=1;
		if (i>=26)
			i=i%26;
		ch[0]=pws[i<<1];
		ch[1]=pws[(i<<1)+1];
		sscanf_s(ch, "%d", &c);
		if (c<0x5A)
			c+=0x64;
		return char(c);
	}

public:

	GeoPrepDocSettings()
	{
		auto_mode=false;

		rect_bak	=100000;
		rect_hmin	=0.1;
		rect_coeff	=1.05;

		rectsp_bak	=150000;
		rectsp_hmin	=1;
		rectsp_coeff=1.1;

		rectloop_bak	=1000000;
		rectloop_rcoeff	=25;
		rectloop_coeff	=1.1;

		impulse_tcoeff	=1.05;
		impulse_t		=100;

		line_mrad_coeff	=20;
		line_infc		=450;
		line_infm		=1000;

		time_gaps_for_loop=0;

		logscale.resize(1);
		logscale[0]=5;

		polygon_mrad_coeff=10;

		enable_recalculations=false;
		enable_2drecalculations=false;

		slae3d_eps=1e-6;

		vfem_direct_output=false;

		meshed=false;

		emf_in_line=false;

		enable_pariterscheme=true;
		epsgfp=1e-5;

		un3dmg=true;

		optXZ=true;
		optAZ=true;
		optUZ=true;

		Cxyz[0]=Cxyz[1]=Cxyz[2]=1.3;
		C2xyz[0]=C2xyz[1]=C2xyz[2]=1.5;
		Cplmax=1;
		Sxyz[0]=Sxyz[1]=Sxyz[2]=2.0;

		objframe=true;

		N_Periods=1;
		N_StepsPerPeriod=8;

		N_Vect=10;
		MaxIter=10000;

		enable_block_relaxation=false;
		flagRunGMRES=true;

		enable_smooth_emf2d_on_z0=false;

		mesh3d_h2.resize(0);

		smooth_alpha = 0.01;

		smooth1 = false;

		line_coeff = 1;

		FindLLtPrecond=false;

		GPMessages=false;

		OptBaseMesh2d=true;

		b_in_line=false;

		line_zones=false;

		slae2d_eps=1e-20;

		ced_cr=false;

		calcE = false;
		nLayerE = 6;

		OmegaMax = 0.01;
		FreqMed=10.0;

		MaxIterBlockRelax = 100;

		Resultant=false;

		vX=false; cX=0;
		vY=false; cY=0;
		vZ=false; cZ=0;

		ObjectTitles=false;

		MaxIterFit=1000;
		MaxIterReg=1000;

		ForInv=false;

		PathTo2d[0]='\0';
	}

	int Read(const char* fname)
	{
			int n;
			ifstream inf;
			
			inf.open(fname);

			if (!inf)
				return 0;

			inf.exceptions(ios_base::badbit|ios_base::failbit);

			inf>>rect_bak;							inf.ignore(1000, '\n');
			inf>>rect_hmin;							inf.ignore(1000, '\n');
			inf>>rect_coeff;						inf.ignore(1000, '\n');

			inf>>rectsp_bak;						inf.ignore(1000, '\n');
			inf>>rectsp_hmin;						inf.ignore(1000, '\n');
			inf>>rectsp_coeff;						inf.ignore(1000, '\n');

			inf>>rectloop_bak;						inf.ignore(1000, '\n');
			inf>>rectloop_rcoeff;					inf.ignore(1000, '\n');
			inf>>rectloop_coeff;					inf.ignore(1000, '\n');

			inf>>impulse_tcoeff;					inf.ignore(1000, '\n');
			inf>>impulse_t;							inf.ignore(1000, '\n');

			inf>>line_mrad_coeff;					inf.ignore(1000, '\n');
			inf>>line_infc;							inf.ignore(1000, '\n');
			inf>>line_infm;							inf.ignore(1000, '\n');

			inf>>time_gaps_for_loop;				inf.ignore(1000, '\n');


			inf>>n;									inf.ignore(1000, '\n');
			logscale.resize(n);
			for (int i=0; i<n; i++)
				inf>>logscale[i]; 
													inf.ignore(1000, '\n');
			
			inf>>polygon_mrad_coeff;				inf.ignore(1000, '\n');
			
			inf>>enable_recalculations;				inf.ignore(1000, '\n');
			
			inf>>slae3d_eps;						inf.ignore(1000, '\n');

			inf>>enable_2drecalculations;			inf.ignore(1000, '\n');

			inf>>vfem_direct_output;				inf.ignore(1000, '\n');

			inf>>meshed;							inf.ignore(1000, '\n');

			inf>>emf_in_line;						inf.ignore(1000, '\n');

			inf>>enable_pariterscheme;				inf.ignore(1000, '\n');
			inf>>epsgfp;							inf.ignore(1000, '\n');
			inf>>un3dmg;							inf.ignore(1000, '\n');

			inf>>optXZ;								inf.ignore(1000, '\n');
			inf>>optAZ;								inf.ignore(1000, '\n');
			inf>>optUZ;								inf.ignore(1000, '\n');
			
			inf>>Cxyz[0];							inf.ignore(1000, '\n');
			inf>>Cxyz[1];							inf.ignore(1000, '\n');
			inf>>Cxyz[2];							inf.ignore(1000, '\n');

			inf>>C2xyz[0];							inf.ignore(1000, '\n');
			inf>>C2xyz[1];							inf.ignore(1000, '\n');
			inf>>C2xyz[2];							inf.ignore(1000, '\n');

			inf>>Cplmax;							inf.ignore(1000, '\n');

			inf>>Sxyz[0];							inf.ignore(1000, '\n');
			inf>>Sxyz[1];							inf.ignore(1000, '\n');
			inf>>Sxyz[2];							inf.ignore(1000, '\n');

			inf>>objframe;							inf.ignore(1000, '\n');
			
			inf>>N_Periods;							inf.ignore(1000, '\n');
			inf>>N_StepsPerPeriod;					inf.ignore(1000, '\n');

			inf>>N_Vect;							inf.ignore(1000, '\n');
			inf>>MaxIter;							inf.ignore(1000, '\n');

			inf>>enable_block_relaxation;			inf.ignore(1000, '\n'); 

			inf>>enable_smooth_emf2d_on_z0;			inf.ignore(1000, '\n'); 

			inf>>n;									inf.ignore(1000, '\n'); 
			mesh3d_h2.resize(n);
			for (int i=0; i<n; i++)
				inf>>mesh3d_h2[i]; 
													inf.ignore(1000, '\n'); 

			inf>>smooth_alpha;						inf.ignore(1000, '\n'); 

			inf>>smooth1;							inf.ignore(1000, '\n'); 

			inf>>line_coeff;						inf.ignore(1000, '\n'); 

			inf>>FindLLtPrecond;					inf.ignore(1000, '\n');

			inf>>GPMessages;						inf.ignore(1000, '\n');

			inf>>OptBaseMesh2d;						inf.ignore(1000, '\n');

			inf>>b_in_line;							inf.ignore(1000, '\n');
			
			inf>>line_zones;						inf.ignore(1000, '\n');

			inf>>slae2d_eps;						inf.ignore(1000, '\n');

			inf>>ced_cr;							inf.ignore(1000, '\n');

			inf>>calcE;								inf.ignore(1000, '\n'); 
			inf>>nLayerE;							inf.ignore(1000, '\n'); 

			inf>>OmegaMax;							inf.ignore(1000, '\n');
			inf>>FreqMed;							inf.ignore(1000, '\n');

			inf>>MaxIterBlockRelax;					inf.ignore(1000, '\n'); 
			inf>>Resultant;							inf.ignore(1000, '\n'); 

			inf>>MaxIterFit;						inf.ignore(1000, '\n'); 
			inf>>MaxIterReg;						inf.ignore(1000, '\n'); 

			inf.close();

			return 0;
	}

	int Write(const char* fname)
	{
			int n;
			ofstream inf;
			inf.exceptions(ios_base::badbit|ios_base::failbit);

			inf.open(fname);			

			inf<<rect_bak;					inf<<"\t";	inf<<"// Rectmesh, Bak\n";
			inf<<rect_hmin;					inf<<"\t";	inf<<"//           Hmin\n";
			inf<<rect_coeff;				inf<<"\t";	inf<<"//           Coeff\n";

			inf<<rectsp_bak;				inf<<"\t";	inf<<"// Rectmesh SP, Bak\n";
			inf<<rectsp_hmin;				inf<<"\t";	inf<<"//              Hmin\n";
			inf<<rectsp_coeff;				inf<<"\t";	inf<<"//              Coeff\n";

			inf<<rectloop_bak;				inf<<"\t";	inf<<"// Rectmesh for loop, Bak\n";
			inf<<rectloop_rcoeff;			inf<<"\t";	inf<<"//                    R-coeff\n";
			inf<<rectloop_coeff;			inf<<"\t";	inf<<"//                    Coeff\n";

			inf<<impulse_tcoeff;			inf<<"\t";	inf<<"// Time, Coeff\n";
			inf<<impulse_t;					inf<<"\t";	inf<<"//       Divisor\n";

			inf<<line_mrad_coeff;			inf<<"\t";	inf<<"// Line, Loops for line\n";
			inf<<line_infc;					inf<<"\t";	inf<<"//       Infinity coeff for loops\n";
			inf<<line_infm;					inf<<"\t";	inf<<"//       Infinity coeff for mesh\n";

			inf<<time_gaps_for_loop;		inf<<"\t";	inf<<"// Time, Gaps\n";

			n=(int)logscale.size();
			inf<<n;							inf<<"\t";	inf<<"// Log scale size\n";
			for (int i=0; i<n; i++)
				inf<<logscale[i]<<" ";
											inf<<"\t";	inf<<"// Log scale\n";
			
			inf<<polygon_mrad_coeff;		inf<<"\t";	inf<<"// Loops for square loop\n";
			
			inf<<enable_recalculations;		inf<<"\t";	inf<<"// Enable recalculations 3D\n";
			
			inf<<slae3d_eps;				inf<<"\t";	inf<<"// SLAE3D eps\n";

			inf<<enable_2drecalculations;	inf<<"\t";	inf<<"// Enable recalculations 2D\n";

			inf<<vfem_direct_output;		inf<<"\t";	inf<<"// Direct output in VFEM\n";

			inf<<meshed;					inf<<"\t";	inf<<"// Meshed\n";

			inf<<emf_in_line;				inf<<"\t";	inf<<"// EMF2d for line\n";

			inf<<enable_pariterscheme;		inf<<"\t";	inf<<"// Enable parabolic iteration scheme\n";
			inf<<epsgfp;					inf<<"\t";	inf<<"// Residual for parabolic iteration scheme\n";
			inf<<un3dmg;					inf<<"\t";	inf<<"// Use new 3d mesh generator\n";

			inf<<optXZ;						inf<<"\t";	inf<<"// Use xz optimization\n";
			inf<<optAZ;						inf<<"\t";	inf<<"// Use above z optimization\n";
			inf<<optUZ;						inf<<"\t";	inf<<"// Use under z optimization\n";
			
			inf<<Cxyz[0];					inf<<"\t";	inf<<"//   xy-contour\n";
			inf<<Cxyz[1];					inf<<"\t";	inf<<"//   xz-contour\n";
			inf<<Cxyz[2];					inf<<"\t";	inf<<"//   yz-contour\n";

			inf<<C2xyz[0];					inf<<"\t";	inf<<"//   cube in xy\n";
			inf<<C2xyz[1];					inf<<"\t";	inf<<"//   cube in xz\n";
			inf<<C2xyz[2];					inf<<"\t";	inf<<"//   cube in yz\n";

			inf<<Cplmax;					inf<<"\t";	inf<<"//   maximum number of t-nodes\n";

			inf<<Sxyz[0];					inf<<"\t";	inf<<"//   correct x-edges in hex area\n";
			inf<<Sxyz[1];					inf<<"\t";	inf<<"//   correct y-edges in hex area\n";
			inf<<Sxyz[2];					inf<<"\t";	inf<<"//   correct z-edges in hex area\n";

			inf<<objframe;					inf<<"\t";	inf<<"// Object frame\n";
			
			inf<<N_Periods;					inf<<"\t";	inf<<"// Number of periods\n";
			inf<<N_StepsPerPeriod;			inf<<"\t";	inf<<"// Steps for period\n";

			inf<<N_Vect;					inf<<"\t";	inf<<"// Vectors for GMRes\n";
			inf<<MaxIter;					inf<<"\t";	inf<<"// Max iter\n";

			inf<<enable_block_relaxation;	inf<<"\t";	inf<<"// Use block relaxation\n";

			inf<<enable_smooth_emf2d_on_z0;	inf<<"\t";	inf<<"// Smooth EMF2d\n";

			n=(int)mesh3d_h2.size();
			inf<<n;							inf<<"\t";	inf<<"// H2 for 3D mesh\n";
			mesh3d_h2.resize(n);
			for (int i=0; i<n; i++)
				inf<<mesh3d_h2[i]<<" ";
											inf<<"\t";	inf<<"// H2 types\n";

			inf<<smooth_alpha;				inf<<"\t";	inf<<"// Alfa\n";

			inf<<smooth1;					inf<<"\t";	inf<<"// Smooth with 1 point\n";

			inf<<line_coeff;				inf<<"\t";	inf<<"// Line coeff\n";

			inf<<FindLLtPrecond;			inf<<"\t";	inf<<"// Find LLt\n";

			inf<<GPMessages;				inf<<"\t";	inf<<"// Show messages\n";

			inf<<OptBaseMesh2d;				inf<<"\t";	inf<<"// Use base mesh optimization\n";

			inf<<b_in_line;					inf<<"\t";	inf<<"// Calculate EMF for line\n";
			
			inf<<line_zones;				inf<<"\t";	inf<<"// Use line zones\n";

			inf<<slae2d_eps;				inf<<"\t";	inf<<"// SLAE2D eps\n";

			inf<<ced_cr;					inf<<"\t";	inf<<"// Use CED resultants\n";

			inf<<calcE;						inf<<"\t";	inf<<"// Use RotE\n";
			inf<<nLayerE;					inf<<"\t";	inf<<"// Start time number for RotE\n";

			inf<<OmegaMax;					inf<<"\t";	inf<<"// Omega Max\n";
			inf<<FreqMed;					inf<<"\t";	inf<<"// Frequency Medium\n";

			inf<<MaxIterBlockRelax;			inf<<"\t";	inf<<"// Max iter block relax\n";
			inf<<Resultant;					inf<<"\t";	inf<<"// Resultant\n";

			inf<<MaxIterFit;				inf<<"\t";	inf<<"// Max iter fit\n";
			inf<<MaxIterReg;				inf<<"\t";	inf<<"// Max iter reg\n";

			inf.close();

			return 0;
	}
}; 

