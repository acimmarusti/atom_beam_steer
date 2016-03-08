/* Propagation of atom beam through 2D optical molasses and 2D laser lens and into optical cavity
 * 
 * Language: C++ using ROOT: Data analysis framework
 *
 * Author:		Andres Cimmarusti
 *			Joint Quantum Institute
 *			University of Maryland
 *			College Park, MD 20742
 *
 *Contact:		candres@umd.edu
 *
 *License:		GNU GPLv3
 *			
 *Dedication:		To Maya Kabkab
 *
 *Acknowledgements:
 *			David Norris
 *			Steve Rolston
 *			Luis Orozco
 *			Howard Carmichael
 *
 * Last update: 2011-02-25
 */

/* TODO:

   - Parallelizing the code has become urgent as it takes too long to run after the addition of the laser lens beams.

   - Code needs clean-up. Could be made more modular.

   - From several tests I've run it is quite clear that the molasses and lens stage creates the most elements (dt shortest). Improvements in other stages are of little help, in terms of performance.
*/

/*****C++ STL*****/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <ctime>
#include <algorithm>


/*****ROOT libraries*****/

#include "TROOT.h"
#include "TApplication.h"

/*****ROOT math/physics libraries*****/

#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TRotation.h"
#include "TVector3.h"


/*****ROOT histogram/graph libraries*****/

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"


/*****ROOT geometry package libraries*****/

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TGeoBBox.h"
#include "TGeoHype.h"
#include "TGeoCompositeShape.h"


/*****Other header files*****/

#include "read_inputs.h"


const double pi = TMath::Pi();		// Pi
const double h = TMath::H();		// Planck's constant
const double kb = TMath::K();		// Boltzmann constant
const double mproton = 1.67262158e-27;	// Mass of a proton (Kg)
const double mneutron = 1.67492729e-27;	// Mass of a neutron (Kg)


/*****Function prototypes*****/

double rayleigh_length(double wavelength, double waist);
double spot_size(double wavelength, double waist, double z);
void rand_state_atom(int &atomic_state, double proba);
void spontaneous_decay(int &state_atom, int kick);

/*****Main*****/

int lvis_mol_lens(int arg = 0, string region = "", string par_type = "", int ite = 0, double step = 0) {

  /*****Input variable declaration*****/


  /*****Atomic parameters*****/

  int natoms;			// Number of atoms
  int nkeep;			// Number of atom trajectories to keep
  int nplot;			// Number of atom trajectories to plot
  int atomic_num;		// Z -> Atomic number
  int mass_num;			// A -> Mass number (isotope)
  double gam;			// Linewidth in Hz

  /*****Position and velocity of the atoms to keep*****/

  vector<double> time_atom;
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> vx;
  vector<double> vy;
  vector<double> vz;

  /*****Vacuum chamber geometry*****/

  double molcham_rmin;	        // Molasses chamber internal radius
  double molcham_rmax;	        // Molasses chamber external radius
  double winrad;		// Radius of a chamber window
  double motchamsize;		// Magneto-Optical trap chamber size
  double hole;			// Diameter of hole in mirror

  /*****LVIS parameters*****/

  double theta0;		// Initial polar angle
  double phi0;			// Initial azimuthal angle
  double v0;			// Initial speed
  double dv0;			// Initial longitudinal speed fwhm
  int sampling;			// Sampling outside laser regions and cavity
  double grav;			// Gravitional acceleration
  double lam;			// Wavelength

  /*****Optical cavity parameters*****/

  double cavspa;		// Mirror spacing
  double mirdiam;		// Mirror diameter
  double mirwidth;		// Mirror width
  double rcavmod;		// Mode radius
  double zcav;			// Mode position

  /*****2D molasses beams*****/

  int nmolbeams;		// Number molasses beams
  double posw0mol;		// Waist position
  double waistmol;		// Waist
  double satmol;		// On-resonance saturation parameter at waist
  double deltamol;		// Detuning
  double centermol[3];	        // Position molasses center

  /*****2D laser lens beams*****/

  int nlensbeams;		// Number of laser lens beams
  double posw0lens;		// Waist position
  double waistlens;		// Waist
  double satlens;		// On-resonance saturation parameter at waist
  double deltalens;		// Detuning
  double centerlens[3];	        // Position laser lens center
  int lenslock;			// Lock spot size at center?
  double lenssize;		// Fixed spot size at center
  int atom6;			// Use 6-level atom approximation

  /*****Tools and flags to read from configuration file*****/

  int config_exist;
  string input = "";
  string useful_in = "";
  string useless_in = "";
  int natoms_flag = 0;
  int nkeep_flag = 0;
  int nplot_flag = 0;
  int atomic_num_flag = 0;
  int mass_num_flag = 0;
  int gam_flag = 0;
  int molcham_rmin_flag = 0;
  int molcham_rmax_flag = 0;
  int winrad_flag = 0;
  int motcham_flag = 0;
  int theta_flag = 0;
  int phi_flag = 0;
  int v0_flag = 0;
  int dv0_flag = 0;
  int sampling_flag = 0;
  int g_flag = 0;
  int lam_flag = 0;
  int cavspa_flag = 0;
  int mirdiam_flag = 0;
  int mirwidth_flag = 0;
  int rcavmod_flag = 0;
  int hole_flag = 0;
  int zcav_flag = 0;
  int nmol_flag = 0;
  int w0molpos_flag = 0;
  int w0mol_flag = 0;
  int s0mol_flag = 0;
  int detunmol_flag = 0;
  int cmolx_flag = 0;
  int cmoly_flag = 0;
  int cmolz_flag = 0;
  int nlens_flag = 0;
  int w0lenspos_flag = 0;
  int w0lens_flag = 0;
  int s0lens_flag = 0;
  int detunlens_flag = 0;
  int clensx_flag = 0;
  int clensy_flag = 0;
  int clensz_flag = 0;
  int lenslock_flag = 0;
  int lenssize_flag = 0;
  int atom6_flag = 0;

  /*****Reading inputs from configuration file*****/

  fstream config_file;

  config_file.open ("lvis_mol_lens.conf", ios::in);

  if (config_file.is_open()) {

    while (getline (config_file, input)) {
      
      size_t ini = input.find_first_not_of(" \t\n\v");

      size_t eq_pos = input.find_last_of("=");

      if (ini != string::npos && input[ini] == '#') continue;

      useful_in = input.substr(eq_pos + 1);

      useless_in = input.substr(0, eq_pos);

      read_input(natoms_flag, "natoms", useless_in, useful_in, natoms);
      read_input(nkeep_flag, "nkeep", useless_in, useful_in, nkeep);
      read_input(nplot_flag, "nplot", useless_in, useful_in, nplot);
      read_input(atomic_num_flag, "atomic_num", useless_in, useful_in, atomic_num);
      read_input(mass_num_flag, "mass_num", useless_in, useful_in, mass_num);
      read_input(gam_flag, "gam", useless_in, useful_in, gam);
      read_input(molcham_rmin_flag, "molcham_rmin", useless_in, useful_in, molcham_rmin);
      read_input(molcham_rmax_flag, "molcham_rmax", useless_in, useful_in, molcham_rmax);
      read_input(winrad_flag, "winrad", useless_in, useful_in, winrad);
      read_input(motcham_flag, "motchamsize", useless_in, useful_in, motchamsize);
      read_input(hole_flag, "hole", useless_in, useful_in, hole);
      read_input(theta_flag, "theta0", useless_in, useful_in, theta0);
      read_input(phi_flag, "phi0", useless_in, useful_in, phi0);
      read_input(v0_flag, "v0", useless_in, useful_in, v0);
      read_input(dv0_flag, "dv0", useless_in, useful_in, dv0);
      read_input(sampling_flag, "sampling", useless_in, useful_in, sampling);
      read_input(g_flag, "g", useless_in, useful_in, grav);
      read_input(lam_flag, "lam", useless_in, useful_in, lam);
      read_input(cavspa_flag, "cavspa", useless_in, useful_in, cavspa);
      read_input(mirdiam_flag, "mirdiam", useless_in, useful_in, mirdiam);
      read_input(mirwidth_flag, "mirwidth", useless_in, useful_in, mirwidth);
      read_input(rcavmod_flag, "rcavmod", useless_in, useful_in, rcavmod);
      read_input(zcav_flag, "zcav", useless_in, useful_in, zcav);
      read_input(nmol_flag, "nmolbeams", useless_in, useful_in, nmolbeams);
      read_input(w0molpos_flag, "w0molpos", useless_in, useful_in, posw0mol);
      read_input(w0mol_flag, "w0mol", useless_in, useful_in, waistmol);
      read_input(s0mol_flag, "s0mol", useless_in, useful_in, satmol);
      read_input(detunmol_flag, "detunmol", useless_in, useful_in, deltamol);
      read_input(cmolx_flag, "cmolx", useless_in, useful_in, centermol[0]);
      read_input(cmoly_flag, "cmoly", useless_in, useful_in, centermol[1]);
      read_input(cmolz_flag, "cmolz", useless_in, useful_in, centermol[2]);
      read_input(nlens_flag, "nlensbeams", useless_in, useful_in, nlensbeams);
      read_input(w0lenspos_flag, "w0lenspos", useless_in, useful_in, posw0lens);
      read_input(w0lens_flag, "w0lens", useless_in, useful_in, waistlens);
      read_input(s0lens_flag, "s0lens", useless_in, useful_in, satlens);
      read_input(detunlens_flag, "detunlens", useless_in, useful_in, deltalens);
      read_input(clensx_flag, "clensx", useless_in, useful_in, centerlens[0]);
      read_input(clensy_flag, "clensy", useless_in, useful_in, centerlens[1]);
      read_input(clensz_flag, "clensz", useless_in, useful_in, centerlens[2]);
      read_input(lenslock_flag, "lenslock", useless_in, useful_in, lenslock);
      read_input(lenssize_flag, "lenssize", useless_in, useful_in, lenssize);
      read_input(atom6_flag, "atom6", useless_in, useful_in, atom6);
    }
		
    config_exist = 1;

  } else {
  
    config_exist = 0;
		
    cout << "WARNING: Configuration file not present, or could not be read" /*<< endl << "the program will create one"*/ << endl;

    return 1;
  }

  config_file.close ();

  if (natoms < nkeep || nplot > nkeep)

    return 2;


  /*****Parsing command line inputs*****/

  int lenslockin = 0;

  if (arg == 0)

    cout << "using laser parameters from configuration file" << endl << endl;

  else {
    if (region == "molasses" || region == "mol") {
      if (par_type == "detun" || par_type == "detuning" || par_type == "delta")

	deltamol = ite * step;

      if (par_type == "waist" || par_type == "w0" || par_type == "focus")

	waistmol = (ite + 1) * step;
    }

    if (region == "lens" || region == "laserlens") {
      if (par_type == "detun" || par_type == "detuning" || par_type == "delta")

	deltalens = ite * step;

      if (par_type == "waist" || par_type == "w0" || par_type == "focus") {
	waistlens = (ite + 1) * step;

	lenslockin = 1;
      }
    }
  }


  /*****Preliminary calculations and conversions*****/

  double matom = atomic_num * mproton + (mass_num - atomic_num) * mneutron;	// Mass of an atom

  theta0 = theta0 * pi / 180;					// Conversion from degrees to radians

  phi0 = phi0 * pi / 180;

  gam = gam * 2 * pi;							// Linewidth


  /*****Molasses beams characteristics*****/

  double *wmol = new double [nmolbeams];			// Beam spot size

  double *w0mol = new double [nmolbeams];			// Beam waist

  double *beamdivermol = new double [nmolbeams];	// Stereo beam divergence angle

  double *detunmol = new double [nmolbeams];		// Detunings

  double *gs0mol = new double [nmolbeams];		// Gaussian beam on-resonance saturation parameter

  double *s0mol = new double [nmolbeams];			// On-resonance saturation parameter (I / Isat) at waist

  double *smol = new double [nmolbeams];			// Saturation parameter

  double *srmol = new double [nmolbeams];			// Scattering rates

  TVector3 kmol[4];							// Wavevectors molasses

  TVector3 w0molpos[4];						// Molasses beam waist position from center molasses reference frame

  TRotation molrotz[4];						// Rotations to transform from center molasses reference frame to beam frame

  TRotation molrotx[4];

  double molwidth = 0;						// Average beam width at molasses center

  int i;

  for (i = 0 ; i < nmolbeams ; i ++) {
    kmol[i].SetX(1);

    kmol[i].SetMag(2 * pi / lam);

    kmol[i].SetTheta(pi / 2);

    if (i < 2)

      kmol[i].SetPhi(i * pi);

    else

      kmol[i].SetPhi((2 * i - 3) * pi / 2);

    w0molpos[i] = - kmol[i].Unit();

    w0molpos[i].SetMag(posw0mol);

    /* Preserve order, first rotate around z then around x */
    molrotz[i].TRotation::RotateZ(kmol[i].Phi() + pi / 2);

    molrotx[i].TRotation::RotateX(kmol[i].Theta());

    w0mol[i] = waistmol;

    s0mol[i] = satmol;

    detunmol[i] = deltamol * gam;

    molwidth += spot_size(lam, w0mol[i], w0molpos[i].Mag()) / nmolbeams;

    beamdivermol[i] = atan(lam / (pi * w0mol[i])) * 180 / pi;
  }


  /*****Laser lens beams characteristics*****/

  double *wlens = new double [nlensbeams];		// Beam spot size

  double *w0lens = new double [nlensbeams];		// Beam waist

  double *beamdiverlens = new double [nlensbeams];	// Stereo beam divergence angle

  double *detunlens = new double [nlensbeams];		// Detunings

  double *gs0lens = new double [nlensbeams];		// Gaussian beam on-resonance saturation parameter

  double *s0lens = new double [nlensbeams];		// On-resonance saturation parameter (I / Isat) at waist

  double *slens = new double [nlensbeams];		// Saturation parameter

  double *srlens = new double [nlensbeams];		// Scattering rates

  double *clebsch = new double [nlensbeams];		// Clebsch-Gordan coefficients

  TVector3 klens[4];							// Wavevectors laser lens

  TVector3 w0lenspos[4];						// Lens beam waist position from center reference frame

  TRotation lensrotz[4];						// Rotations to transform from center lens reference frame to beam frame

  TRotation lensrotx[4];

  double lenswidth = 0;	 					// Average beam width at laser lens center

  for (i = 0 ; i < nlensbeams ; i ++) {

    klens[i].SetX(1);

    klens[i].SetMag(2 * pi / lam);

    klens[i].SetTheta(pi / 2);

    if (i < 2)

      klens[i].SetPhi(i * pi);

    else

      klens[i].SetPhi((2 * i - 3) * pi / 2);

    w0lenspos[i] = - klens[i].Unit();

    w0lens[i] = waistlens;

    if (lenslock == 1 &&  lenslockin == 1) {
      posw0lens = pi * waistlens * sqrt(pow(lenssize, 2) - pow(waistlens, 2)) / lam;

      /* cannot be greater than the geometry boundary */
      if (posw0lens > 0.5)

	return 3;
    }

    w0lenspos[i].SetMag(posw0lens);
      
    /* Preserve order, first rotate around z then around x */
    lensrotz[i].TRotation::RotateZ(klens[i].Phi() + pi / 2);

    lensrotx[i].TRotation::RotateX(klens[i].Theta());

    s0lens[i] = satlens / pow(w0lens[i] / spot_size(lam, w0lens[i], w0lenspos[i].Mag()), 2);

    detunlens[i] = (deltalens) * gam;

    lenswidth += spot_size(lam, w0lens[i], w0lenspos[i].Mag()) / nlensbeams;

    beamdiverlens[i] = atan(lam / (pi * w0lens[i])) * 180 / pi;

    clebsch[i] = 1;
  }


  /*****Beam geometry (all lengths must be in centimeters)*****/


  /*****Definition of top volume corresponding to a cube circumscribing the internal region of the molasses vacuum chamber*****/

  TGeoManager *mgr = new TGeoManager("apparatus", "Simple geometry");

  TGeoMedium *medium = 0;

  TGeoVolume *apparatus = mgr->MakeBox("A", medium, 50, 50, 50);

  //	TGeoVolume *apparatus = mgr->MakeBox("A", medium, 100 * molcham_rmin, 100 * molcham_rmin, 100 * molcham_rmin);

  mgr->SetTopVolume(apparatus);


  /*****Molasses vacuum chamber*****/
  /*
    TGeoVolume *chamber = mgr->MakeSph("C", medium, molcham_rmin * 100, molcham_rmax * 100, 50);

    chamber->SetLineColor(5);

    apparatus->AddNode(chamber, 1);

  */
  /*****Rotations and translations to be performed on hyperboloid base*****/


  /*****Rotations are in GEANT3 style, 3 pairs of angles, corresponding to angular position of the new axes*****/

  TGeoRotation *rot1 = new TGeoRotation("rot1", 90., 90., 0., 0., 90., 0.);

  TGeoRotation *rot2 = new TGeoRotation("rot2", 90., 270., 180., 0., 90., 0.);

  TGeoRotation *rot3 = new TGeoRotation("rot3", 90., 0., 180., 0., 90., 90.);

  TGeoRotation *rot4 = new TGeoRotation("rot4", 90., 0., 0., 0., 90., 270.); // TODO: make it so that the angles are put in terms of k


  /*****Beam translations in centimeters*****/

  double transmol[4][3];

  double translens[4][3];

  for (int ii = 0 ; ii < 4 ; ii++)

    for (int jj = 0 ; jj < 3 ; jj++) {

      transmol[ii][jj] = (w0molpos[ii](jj) + centermol[jj]) * 100;

      translens[ii][jj] = (w0lenspos[ii](jj) + centerlens[jj]) * 100;
    }


  /*****Combination of translations and rotations (molasses)*****/

  TGeoCombiTrans *combi1 = new TGeoCombiTrans("combi1", transmol[0][0], transmol[0][1], transmol[0][2], rot1);

  TGeoCombiTrans *combi2 = new TGeoCombiTrans("combi2", transmol[1][0], transmol[1][1], transmol[1][2], rot2);

  TGeoCombiTrans *combi3 = new TGeoCombiTrans("combi3", transmol[2][0], transmol[2][1], transmol[2][2], rot3);

  TGeoCombiTrans *combi4 = new TGeoCombiTrans("combi4", transmol[3][0], transmol[3][1], transmol[3][2], rot4);

  combi1->RegisterYourself();

  combi2->RegisterYourself();

  combi3->RegisterYourself();

  combi4->RegisterYourself();


  /*****Combination of translations and rotations (laser lens)*****/

  TGeoCombiTrans *combi11 = new TGeoCombiTrans("combi11", translens[0][0], translens[0][1], translens[0][2], rot1);

  TGeoCombiTrans *combi21 = new TGeoCombiTrans("combi21", translens[1][0], translens[1][1], translens[1][2], rot2);

  TGeoCombiTrans *combi31 = new TGeoCombiTrans("combi31", translens[2][0], translens[2][1], translens[2][2], rot3);

  TGeoCombiTrans *combi41 = new TGeoCombiTrans("combi41", translens[3][0], translens[3][1], translens[3][2], rot4);

  combi11->RegisterYourself();

  combi21->RegisterYourself();

  combi31->RegisterYourself();

  combi41->RegisterYourself();


  /*****Center conversions to centimeters*****/

  double cmolcm[3];

  double clenscm[3];

  for (int kk = 0 ; kk < 3 ; kk++) {

    cmolcm[kk] = centermol[kk] * 100;

    clenscm[kk] = centerlens[kk] * 100;
  }

  /*****Basic shapes: a hyperboloid for each beam and a box as constraint (molasses)*****/

  TGeoBBox *molbox = new TGeoBBox("molbox", (w0molpos[0].Mag() + w0molpos[1].Mag()) * 100 / 2, (w0molpos[2].Mag() + w0molpos[3].Mag()) * 100 / 2, molwidth * 100, cmolcm);

  TGeoHype *beam1 = new TGeoHype("beam1", 0, 0, w0mol[0] * 100, beamdivermol[0], w0molpos[0].Mag() * 100);

  TGeoHype *beam2 = new TGeoHype("beam2", 0, 0, w0mol[1] * 100, beamdivermol[1], w0molpos[1].Mag() * 100);

  TGeoHype *beam3 = new TGeoHype("beam3", 0, 0, w0mol[2] * 100, beamdivermol[2], w0molpos[2].Mag() * 100);

  TGeoHype *beam4 = new TGeoHype("beam4", 0, 0, w0mol[3] * 100, beamdivermol[3], w0molpos[3].Mag() * 100);


  /*****Basic shapes: a hyperboloid for each beam and a box as constraint (laser lens)*****/

  TGeoBBox *lensbox = new TGeoBBox("lensbox", (w0lenspos[0].Mag() + w0lenspos[1].Mag()) * 100 / 2, (w0lenspos[2].Mag() + w0lenspos[3].Mag()) * 100 / 2, lenswidth * 100, clenscm);

  TGeoHype *beam11 = new TGeoHype("beam11", 0, 0, w0lens[0] * 100, beamdiverlens[0], w0lenspos[0].Mag() * 100);

  TGeoHype *beam21 = new TGeoHype("beam21", 0, 0, w0lens[1] * 100, beamdiverlens[1], w0lenspos[1].Mag() * 100);

  TGeoHype *beam31 = new TGeoHype("beam31", 0, 0, w0lens[2] * 100, beamdiverlens[2], w0lenspos[2].Mag() * 100);

  TGeoHype *beam41 = new TGeoHype("beam41", 0, 0, w0lens[3] * 100, beamdiverlens[3], w0lenspos[3].Mag() * 100);


  /*****Composite shape: translation + rotation of each hyperboloid intersected with box (molasses)*****/

  TGeoCompositeShape *beams = new TGeoCompositeShape("beams", "(beam1:combi1 + beam2:combi2 + beam3:combi3 + beam4:combi4) * molbox");

  TGeoVolume *molasses = new TGeoVolume("M", beams);

  molasses->SetLineColor(kRed);

  apparatus->AddNode(molasses, 1);


  /*****Composite shape: translation + rotation of each hyperboloid intersected with box (laser lens)*****/

  TGeoCompositeShape *lensbeams = new TGeoCompositeShape("lensbeams", "(beam11:combi11 + beam21:combi21 + beam31:combi31 + beam41:combi41) * lensbox");

  TGeoVolume *laserlens = new TGeoVolume("L", lensbeams);

  laserlens->SetLineColor(kRed);

  apparatus->AddNode(laserlens, 2);

  mgr->CloseGeometry();


  TVector3 cmol(1); 	// Position molasses center

  TVector3 clens(1); 	// Position lens center

  TVector3 rhole(1);	// Position of LVIS hole center in molasses reference frame

  double dmothole;	// Distance from MOT to hole in mirror

  double tube;		// In vacuo tube length

  double thetam;	// Minimum polar angle to clear cavity


  if (nmolbeams != 0) {

    if (nlensbeams != 0) {
      //mgr->SetNsegments(120);

      //apparatus->Raytrace();

      //apparatus->Draw();

      clens.SetXYZ(centerlens[0], centerlens[1], centerlens[2]);
    }

    cmol.SetXYZ(centermol[0], centermol[1], centermol[2]);


    /*****Position of LVIS hole center in molasses reference frame*****/

    double d0 = sqrt(3.) * molwidth;		// Distance to LVIS hole center from molasses center

    rhole.SetTheta(pi - theta0);

    rhole.SetPhi(pi + phi0);

    rhole.SetMag(d0);


    /*****Position of the intersection between the molasses chamber and the in vacuo tube*****/


    /*****Change center molasses vector to tube rotation-only reference frame*****/

    TVector3 cmoltube(cmol);			// Position molasses center in tube rotation-only reference frame

    cmoltube.RotateZ(phi0);

    cmoltube.RotateX(theta0);


    /*****Set parameters in this frame for position vector for the intersection between chamber and tube*****/

    TVector3 tubepos;				// Position intersection between molasses chamber and LVIS tube

    tubepos.SetY(cmoltube.Y());

    tubepos.SetX(cmoltube.X());

    tubepos.SetZ(sqrt(pow(molcham_rmax, 2) - pow(cmoltube.X(), 2) - pow(cmoltube.Y(), 2)));		// Keep magnitude equal to chamber radius


    /*****Change back to center chamber reference frame*****/

    tubepos.RotateX(- theta0);

    tubepos.RotateZ(- phi0);


    /*****Find tube's length*****/

    TVector3 tubeposmol;				// Position intersection between molasses chamber and LVIS tube in molasses reference frame

    tubeposmol = tubepos - cmol;

    tube = tubeposmol.Mag() - rhole.Mag();

    dmothole = tube + motchamsize / 2;

    thetam = atan((cavspa / 2 + mirwidth) / (TMath::Abs(2 * cmol.Z()) + zcav));

	
    /*****Sanity check*****/

    TVector3 sanity; 

    sanity = tubeposmol.Cross(rhole);

    cout << "sanity check : (" << sanity.X() << ", " << sanity.Y() << ", " << sanity.Z() << ")" << endl;
  } else {
    cout << "No molasses beams" << endl;

    cmol.SetXYZ(0, 0, 0);

    rhole.SetXYZ(0, 0, 0);

    dmothole = 0.075 - zcav;
  }


  /*****Standard deviation LVIS initial speed*****/

  double v0dev = dv0 / (2 * sqrt(2 * log(2)));


  /*****LVIS initial transverse speed width*****/

  double dvtmax = hole / dmothole * v0;


  /*****Initialization of atomic beam averages*****/

  double *height = new double [sampling];			// Z-position of atomic beam

  double *rtrms = new double [sampling];			// Normalized root mean square position average of atomic beam

  double *xavg = new double [sampling];			// Average y position of atomic beam

  double *yavg = new double [sampling];			// Average x position of atomic beam

  double *ravg = new double [sampling];			// Average position magnitude of atomic beam

  double *vzavg = new double [sampling];			// Average z velocity of atomic beam

  double *vavg = new double [sampling];			// Average speed of atomic beam

  int *entries = new int [sampling];			// Number of entries per bin

  double hstep = (zcav - centermol[2] + molwidth / 2) / sampling;		// Bin height size

  int zbin;

  for (zbin = 0 ; zbin < sampling ; zbin++) {
    height[zbin] = (zbin + 0.5) * hstep - molwidth / 2 + centermol[2];

    rtrms[zbin] = 0;

    xavg[zbin] = 0;

    yavg[zbin] = 0;

    ravg[zbin] = 0;

    vzavg[zbin] = 0;

    vavg[zbin] = 0;

    entries[zbin] = 0;
  }


  int nbins;						// Number of bins for histograms

  if (natoms < 1e4)

    nbins = natoms / 10;

  else

    if (natoms < 1e5)

      nbins = natoms / 100;

    else

      nbins = natoms / 1000;


  /*****histograms characterization atom beam*****/

  //	TH2D *rinitest = new TH2D("rinitest", "Initial x-y position distribution", nbins, - hole / 2, hole / 2, nbins, - hole / 2, hole / 2);

  //	TH1D *xrinitest = new TH1D("xrinitest", "Initial x position histogram", nbins, - hole / 2, hole / 2);

  //	TH1D *yrinitest = new TH1D("yrinitest", "Initial y position histogram", nbins, - hole / 2, hole / 2);

  TH2D *xypos = new TH2D("xypos", "x-y position distribution;x (m);y (m)", nbins, - 4 * hole, 4 * hole, nbins, - 4 * hole, 4 * hole);

  TH2D *xyposlens = new TH2D("xyposlens", "x-y position distribution after laser lens;x (m);y (m)", nbins, - 4 * hole, 4 * hole, nbins, - 4 * hole, 4 * hole);

  //	TH1D *vxmot = new TH1D("vxmot", "X velocity distribution leaving MOT", nbins, - dvtmax, dvtmax);

  //	TH1D *vymot = new TH1D("vymot", "Y velocity distribution leaving MOT", nbins, - dvtmax, dvtmax);

  //	TH1I *numit_mol = new TH1I("numit_mol", "Number of iterations molasses", natoms / 10, 0, 2.5e6);

  //	TH1I *numit_lens = new TH1I("numit_lens", "Number of iterations laser lens", natoms / 10, 0, 2.5e6);

  //	TH1I *numit_cav = new TH1I("numit_cav", "Number of iterations propagation cavity", natoms / 10, 0, 1e3);

  //	TH1D *xmol = new TH1D("xmol", "x distribution leaving molasses", nbins, - hole, hole);

  //	TH1D *ymol = new TH1D("ymol", "y distribution leaving molasses", nbins, - hole, hole);

  TH1D *vxmol = new TH1D("vxmol", "x velocity distribution leaving molasses", nbins, - 0.5, 0.5);

  //	TH1D *vymol = new TH1D("vymol", "y velocity distribution leaving molasses", nbins, - 0.5, 0.5);

  TH1D *vzmol = new TH1D("vzmol", "z velocity distribution leaving molasses", nbins, (v0 - dv0) * cos(theta0), (v0 - dv0) * cos(theta0));

  TH1D *vxmoldiv = new TH1D("vxmoldiv", "x diverging velocity distribution leaving molasses", nbins * 100, (v0 - dv0) * sin(theta0) * cos(phi0), (v0 - dv0) * sin(theta0) * cos(phi0));

  TH1D *vzmoldiv = new TH1D("vzmoldiv", "z diverging velocity distribution leaving molasses", nbins * 100, (v0 - dv0) * cos(theta0), (v0 - dv0) * cos(theta0));

  TH1D *xmoldiv = new TH1D("xmoldiv", "x diverging distribution leaving molasses", nbins * 100, - hole, hole);

  //	TMultiGraph *x_vs_t = new TMultiGraph();

  //	TMultiGraph *y_vs_t = new TMultiGraph();

  //	TMultiGraph *z_vs_t = new TMultiGraph();

  TMultiGraph *z_vs_x = new TMultiGraph();

  //	TMultiGraph *z_vs_y = new TMultiGraph();

  //	TMultiGraph *vx_vs_t = new TMultiGraph();

  //	TMultiGraph *vy_vs_t = new TMultiGraph();

  //	TMultiGraph *vz_vs_t = new TMultiGraph();

  TMultiGraph *vx_vs_z = new TMultiGraph();

  //	TMultiGraph *vy_vs_z = new TMultiGraph();


  /*****Random numbers*****/

  TRandom3 randnum;

  TRandom3 randnum1(0);

  TRandom3 randnum2(0);

  TRandom3 randnum3(0);

  TRandom3 randnum4(0);

  TRandom3 randnum5(0);

  TRandom3 randnum6(0);

  double prob[2];				// Storage for random numbers


  /*****Atom counters and/or indicators*****/

  int natoms_mod = 0;				// Number of atoms that go through the cavity mode

  int benchmark = 0;				// Number of atoms randomly generated in geometric projection of the cavity mode

  int natoms_mol = 0;				// Number of atoms that go through molasses region

  int natoms_lens = 0;			// Number of atoms that go through laser lens region

  double vzmodsum = 0; 			// Sum z-velocity of atoms going through the cavity mode

  double vt2avgmol = 0;			// Average transverse speed squared after molasses

  int doppler_count = 0;			// Counter to estimate Doppler cooling limit

  int sr_count = 0;				// Counter to estimate average scattering rate at molasses center

  double avg_sr = 0;				// Average scattering rate at molasses center

  double phot_scatt = 0;			// Average photon scattering


  /*****Recoil velocity parameters*****/

  TVector3 vrec_abs(1);				// Absorption recoil velocity

  TVector3 vrec_spon;					// Spontaneous emission recoil velocity

  TVector3 spon_dir;					// Spontaneous emission direction vector

  double vrec = h / (matom * lam);		// Recoil velocity magnitude

  vrec_abs.SetTheta(pi / 2);

  vrec_abs.SetMag(vrec);


  /*****Individual trajectory per atom*****/

  time_t start;					// Start time marker

  time(& start);

  for (int j = 0 ; j < natoms ; j ++) {
    double t = 0;

    double randcont = randnum.Rndm();


    /*****Initial position spread in tube reference frame*****/

    TVector3 rini(1);			// Initial position atom in tube reference frame

    double radini; 				// Initial radial position of the atom in the hole

    do {
      rini.SetX(randnum1.Uniform(- hole / 2, hole / 2));

      rini.SetY(randnum2.Uniform(- hole / 2, hole / 2));

      radini = sqrt(pow(rini.X(), 2) + pow(rini.Y(), 2));
    }
    while (radini > hole / 2);

    rini.SetZ(0);

    //		rinitest->Fill(rini.X(), rini.Y());

    //		xrinitest->Fill(rini.X());

    //		yrinitest->Fill(rini.Y());


    /*****Initial velocity in tube reference frame*****/

    TVector3 vini(1);				// Initial velocity atom in tube reference frame


    /*****Position of edges in reference frame of the atom*****/

    double xedgep = sqrt(pow(hole / 2, 2) - pow(rini.Y(), 2));

    double xedgem = - sqrt(pow(hole / 2, 2) - pow(rini.Y(), 2));

    double yedgep = sqrt(pow(hole / 2, 2) - pow(rini.X(), 2));

    double yedgem = - sqrt(pow(hole / 2, 2) - pow(rini.X(), 2));

    double xp = - rini.X() + xedgep;

    double yp = - rini.Y() + yedgep;

    double xm = - rini.X() + xedgem;

    double ym = - rini.Y() + yedgem;


    /*****Position dependent transverse velocity spreads*****/

    double dvxp = xp / dmothole * v0;

    double dvyp = yp / dmothole * v0;

    double dvxm = xm / dmothole * v0;

    double dvym = ym / dmothole * v0;


    vini.SetX(randnum3.Uniform(- dvxp, - dvxm));

    vini.SetY(randnum4.Uniform(- dvyp, - dvym));

    vini.SetZ(randnum5.Gaus(v0, v0dev));

    //		vini.SetXYZ(0, 0, v0);		// All atoms have the same velocity with no spread (for testing purposes only)

    //		vxmot->Fill(vini.X());

    //		vymot->Fill(vini.Y());


    /*****Transformation to molasses reference frame*****/

    if (nmolbeams != 0) {

      vini.RotateX(- theta0);

      vini.RotateZ(- phi0);

      rini.RotateX(- theta0);

      rini.RotateZ(- phi0);
    }

    TVector3 v(vini);				// Velocity atom

    TVector3 rmol; 				// Position atom from molasses reference frame

    rmol = rhole + rini;

    TVector3 r; 					// Position atom

    r = rmol + cmol;				// Change to center chamber frame of reference

    TVector3 g(0, 0, -grav);		// Gravitational acceleration


    /*****Settings for no molasses beams*****/

    int atom_out; 					// Indicator atom out of molasses range	

    if (nmolbeams == 0) {

      atom_out = 1;


      /*****Benchmark: percentage of atoms generated in cavity mode projection area*****/

      double yconst;		// Geometric Y constraint for the cavity mode

      if (hole > cavspa)

	yconst = cavspa / 2;

      else

	yconst = hole / 2;

      if ((rini.Y() > - yconst && rini.Y() < yconst) && (rini.X() > - rcavmod && rini.X() < rcavmod))

	benchmark ++;
    }
    else

      atom_out = 0;


    double pos[3];			// Position atom

    double dir[3];			// Unit velocity atom

    int q;

    for (q = 0 ; q < 3 ; q ++) {

      pos[q] = r(q) * 100;

      dir[q] = v(q) / v.Mag();
    }

    string molcham = "/A_1";			// Geometry string for molasses chamber


    /*****Update geometry navigator, find distance to next region and propagate*****/

    TGeoNode *ini_node = gGeoManager->InitTrack(pos, dir);

    TGeoNode *fin_node = gGeoManager->FindNextBoundaryAndStep(2 * sqrt(3.) * molcham_rmin * 100);

    double dist_atom_mol = gGeoManager->GetStep();		// Distance of atom to molasses beams

    string location = gGeoManager->GetPath();		// Location of the atom within the geometry hierarchy

    const double *steppoint = gGeoManager->GetCurrentPoint();

    if (j < 10) {

      cout << "- atom " << j << endl;

      cout << "Molasses" << endl;		

      cout << "initial position : (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << endl;

      cout << "direction : (" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << endl;

      cout << "position after step : (" << steppoint[0] << ", " << steppoint[1] << ", " << steppoint[2] << ")" << endl;

      cout << "path after step : " << location << endl;
    }

    if (location == molcham || gGeoManager->IsOutside() || nmolbeams == 0)

      atom_out = 1;

    double dtp[2]; 			// Varying precision time step

    dtp[0] = (dist_atom_mol / 100) / v.Mag() / sampling;


    /*****Undo update geometry navigator and probe location of atom*****/

    ini_node = gGeoManager->InitTrack(pos, dir);

    location = gGeoManager->GetPath();

    if (j < 10) {

      cout << "path before step : " << location << endl;

      cout << "Is atom out ? " << atom_out << endl;

      cout << "Distance to molasses : " << dist_atom_mol << endl;
    }


    /*****LVIS propagation prior to molasses beams*****/

    int tracker = 0;

    double dtp2over2;			// Squared precision time step divided by 2

    while (atom_out == 0 && location == molcham) {

      tracker++;

      if (j < nkeep) {

	time_atom.push_back(t);

	x.push_back(r(0));

	y.push_back(r(1));

	z.push_back(r(2));

	vx.push_back(v(0));

	vy.push_back(v(1));

	vz.push_back(v(2));
      }

      dtp2over2 = pow(dtp[0], 2) / 2; 

      t += dtp[0]; 

      v += g * dtp[0];

      r += v * dtp[0];

      r += g * dtp2over2;

      rmol = r - cmol;

      for (q = 0 ; q < 3 ; q ++) {

	pos[q] = r(q) * 100;

	dir[q] = v(q) / v.Mag();
      }


      /*****Update geometry navigator and probe location of atom*****/

      ini_node = gGeoManager->InitTrack(pos, dir);

      location = gGeoManager->GetPath();
    }


    /*****LVIS propagation through molasses beams*****/

    string mol_region = "/A_1/M_1";		// Geometry string for molasses chamber

    TVector3 rbeam; 			// Position atom in beam reference frame

    double dt;				// Photon scattering region time step

    double dt2over2;			// Squared time step divided by 2

    double tot_s;				// Total saturation parameter

    double tot_sr;				// Total scattering rate

    double time_start_mol = t;	// Time stamp when entering 2-D molasses

    tracker = 0;

    while (atom_out == 0 && location == mol_region) {

      tot_sr = 0;

      tot_s = 0;

      if (tracker == 0)

	natoms_mol++;

      tracker++;


      /*****Random numbers: 1 next scattering time, 1 absorption probability, 2 reemission direction*****/

      randnum6.RndmArray(2, prob);

      for (i = 0 ; i < nmolbeams ; i ++) {

	/*****Gaussian beam correction to on-resonance saturation parameter*****/

	rbeam = rmol - w0molpos[i];

	rbeam.Transform(molrotz[i]);

	rbeam.Transform(molrotx[i]);

	wmol[i] = spot_size(lam, w0mol[i], rbeam.Z());

	gs0mol[i] = s0mol[i] * pow(w0mol[i] / wmol[i], 2) * exp(- 2 * pow(rbeam.Perp() / wmol[i], 2));


	/*****Saturation parameter per beam*****/

	smol[i] = gs0mol[i] / (1 + pow(2 * (detunmol[i] - kmol[i].Dot(v)) / gam , 2));

	tot_s += smol[i];
      }


      /*****Scattering rate per beam*****/

      for (i = 0 ; i < nmolbeams ; i ++) {
	srmol[i] = (gam / 2) * smol[i] / (1 + tot_s);

	tot_sr += srmol[i];
      }

      avg_sr += tot_sr;

      sr_count++;


      /*****Time it takes for the atom to decay to ground state and scatter a new photon*****/

      dt = - log(prob[0]) / tot_sr;

      t += dt;


      /*****From which molasses beam did the atom get a kick?*****/

      if (nmolbeams == 2) {
	if (prob[1] < srmol[0] / tot_sr)

	  vrec_abs.SetPhi(0);

	else

	  vrec_abs.SetPhi(pi);
      } else {
	if (prob[1] < srmol[0] / tot_sr)

	  vrec_abs.SetPhi(0);

	else {
	  if (prob[1] < (srmol[0] + srmol[1]) / tot_sr)

	    vrec_abs.SetPhi(pi);

	  else {
	    if (prob[1] < (srmol[0] + srmol[1] + srmol[2]) / tot_sr)

	      vrec_abs.SetPhi(pi / 2);

	    else

	      vrec_abs.SetPhi(3 * pi / 2);
	  }
	}
      }

      /*****Velocity changes due to absorption and spontaneous emission*****/

      spon_dir.SetX(randnum1.Uniform(- 1, 1));

      spon_dir.SetY(randnum2.Uniform(- 1, 1));

      spon_dir.SetZ(randnum3.Uniform(- 1, 1));

      vrec_spon = vrec * spon_dir.Unit();

      v += vrec_abs;

      v += vrec_spon;

      v += g * dt;


      /*****Position changes due to change in velocity*****/

      dt2over2 = pow(dt, 2) / 2;

      r += v * dt;

      r += g * dt2over2;

      rmol = r - cmol;


      /*****Update position and velocity*****/

      for (q = 0 ; q < 3 ; q ++) {

	pos[q] = r(q) * 100;

	dir[q] = v(q) / v.Mag();
      }

      if (j < nkeep) {

	time_atom.push_back(t);

	x.push_back(r(0));

	y.push_back(r(1));

	z.push_back(r(2));

	vx.push_back(v(0));

	vy.push_back(v(1));

	vz.push_back(v(2));
      }


      /*****Update geometry navigator and probe location of atom*****/

      ini_node = gGeoManager->InitTrack(pos, dir);

      location = gGeoManager->GetPath();
    }

    //		numit_mol->Fill(tracker);

    xypos->Fill(r.X(), r.Y());

    for (q = 0 ; q < 3 ; q ++) {

      pos[q] = r(q) * 100;

      dir[q] = v(q) / v.Mag();
    }


    /*****Final check for spread anomaly*****/

    double vxvzratio = TMath::Abs(v.X()) / v.Z();

    if (vxvzratio > tan(pi / 90)) {

      vxmoldiv->Fill(v.X());

      vzmoldiv->Fill(v.Z());

      xmoldiv->Fill(r.X());
    }

    /*****Doppler temperature: Sampling velocity in molasses*****/

    double vt2 = pow(v.X(), 2) + pow(v.Y(), 2);

    if (atom_out == 0) {

      //			xmol->Fill(r.X());

      //			ymol->Fill(r.Y());

      vxmol->Fill(v.X());

      //			vymol->Fill(v.Y());

      vzmol->Fill(v.Z());

      phot_scatt += tracker / (t - time_start_mol);

      if (sqrt(vt2) < 1.) {

	vt2avgmol += vt2;

	doppler_count++;
      }
    }


    /*****Update geometry navigator, find distance to next region and propagate*****/

    ini_node = gGeoManager->InitTrack(pos, dir);

    fin_node = gGeoManager->FindNextBoundaryAndStep(2 * sqrt(3.) * molcham_rmin * 100);

    double dist_atom_lens = gGeoManager->GetStep();	// Distance of atom to lens beams

    location = gGeoManager->GetPath();

    steppoint = gGeoManager->GetCurrentPoint();

    if (j < 10) {

      cout << "Laser lens" << endl;		

      cout << "initial position : (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << endl;

      cout << "direction : (" << dir[0] << ", " << dir[1] << ", " << dir[2] << ")" << endl;

      cout << "position after step : (" << steppoint[0] << ", " << steppoint[1] << ", " << steppoint[2] << ")" << endl;

      cout << "path after step : " << location << endl;
    }

    if (location == molcham || gGeoManager->IsOutside() || nlensbeams == 0)

      atom_out = 1;

    dtp[0] = (dist_atom_lens / 100) / v.Mag() / sampling;


    /*****Undo update geometry navigator and probe location of atom*****/

    ini_node = gGeoManager->InitTrack(pos, dir);

    location = gGeoManager->GetPath();

    if (j < 10) {

      cout << "path before step : " << location << endl;

      cout << "Is atom out ? " << atom_out << endl;

      cout << "Distance to laser lens : " << dist_atom_lens << endl << endl;
    }

    /*****LVIS propagation prior to laser lens*****/

    tracker = 0;

    TVector3 rlens; 			// Position atom from lens reference frame

    while (atom_out == 0 && location == molcham) {
      tracker++;

      if (j < nkeep) {

	time_atom.push_back(t);

	x.push_back(r(0));

	y.push_back(r(1));

	z.push_back(r(2));

	vx.push_back(v(0));

	vy.push_back(v(1));

	vz.push_back(v(2));
      }

      dtp2over2 = pow(dtp[0], 2) / 2; 

      t += dtp[0]; 

      v += g * dtp[0];

      r += v * dtp[0];

      r += g * dtp2over2;

      rlens = r - clens;

      for (q = 0 ; q < 3 ; q ++) {

	pos[q] = r(q) * 100;

	dir[q] = v(q) / v.Mag();
      }


      /*****Update geometry navigator and probe location of atom*****/

      ini_node = gGeoManager->InitTrack(pos, dir);

      location = gGeoManager->GetPath();
    }


    /*****LVIS propagation through laser lens*****/


    /*****Randomly choose state of the atom (x and y quantization axes)*****/

    tracker = 0;

    int state_atomx;			// Ground state of the atom (x quantization axis)

    int state_atomy;			// Ground state of the atom (y quantization axis)

    rand_state_atom(state_atomx, 0.5);

    rand_state_atom(state_atomy, 0.5);

    int kick_atom;				// Polarization for beam from which atom absorbs photon (1 -> sigma plus, 0 -> sigma minus)

    int quant_axis = 0;			// Quantization axis to use (0 -> x, 1 -> y)
	
    string lens_region = "/A_1/L_2";	// Geometry string for molasses region

    while (atom_out == 0 && location == lens_region) {

      tot_sr = 0;

      tot_s = 0;

      if (tracker == 0)

	natoms_lens++;

      tracker++;


      /*****Clebsch-Gordan coefficients for the possible transitions*****/

      if (atom6 == 1) {

	if (state_atomx == 0) {

	  clebsch[0] = 1. / 3.;

	  clebsch[1] = 1.;
	} else {

	  clebsch[0] = 1.;

	  clebsch[1] = 1. / 3.;
	}

	if (state_atomy == 0) {
	  clebsch[2] = 1. / 3.;

	  clebsch[3] = 1.;
	} else {

	  clebsch[2] = 1.;

	  clebsch[3] = 1. / 3.;
	}
      }


      /*****Random numbers: 1 next scattering time, 1 absorption probability, 2 reemission direction*****/

      randnum6.RndmArray(2, prob);

      for (i = 0 ; i < nlensbeams ; i ++) {

	/*****Gaussian beam correction to on-resonance saturation parameter*****/

	rbeam = rlens - w0lenspos[i];

	rbeam.Transform(lensrotz[i]);

	rbeam.Transform(lensrotx[i]);

	wlens[i] = spot_size(lam, w0lens[i], rbeam.Z());

	gs0lens[i] = s0lens[i] * pow(w0lens[i] / wlens[i], 2) * exp(- 2 * pow(rbeam.Perp() / wlens[i], 2));


	/*****Saturation parameter per beam*****/

	slens[i] = clebsch[i] * gs0lens[i] / (1 + pow(2 * (detunlens[i] - klens[i].Dot(v)) / gam , 2));

	tot_s += slens[i];
      }


      /*****Scattering rate per beam*****/

      for (i = 0 ; i < nlensbeams ; i ++) {
	srlens[i] = (gam / 2) * slens[i] / (1 + tot_s);

	tot_sr += srlens[i];
      }


      /*****Time to next scattering: Waiting time distribution*****/

      dt = - log(prob[0]) / tot_sr;

      t += dt;


      /*****From which molasses beam did the atom get a kick?*****/

      if (nlensbeams == 2) {
	if (prob[1] < srlens[0] / tot_sr) {

	  vrec_abs.SetPhi(0);

	  kick_atom = 1;

	} else {

	  vrec_abs.SetPhi(pi);

	  kick_atom = 0;
	}
      } else {

	if (prob[1] < srlens[0] / tot_sr) {

	  vrec_abs.SetPhi(0);

	  kick_atom = 1;

	  quant_axis = 0;
	} else {

	  if (prob[1] < (srlens[0] + srlens[1]) / tot_sr) {

	    vrec_abs.SetPhi(pi);

	    kick_atom = 0;

	    quant_axis = 0;

	  } else {

	    quant_axis = 1;

	    if (prob[1] < (srlens[0] + srlens[1] + srlens[2]) / tot_sr) {

	      vrec_abs.SetPhi(pi / 2);

	      kick_atom = 1;

	    } else {

	      vrec_abs.SetPhi(3 * pi / 2);

	      kick_atom = 0;
	    }
	  }
	}
      }


      /*****Velocity changes due to absorption and spontaneous emission*****/

      spon_dir.SetX(randnum4.Uniform(- 1, 1));

      spon_dir.SetY(randnum5.Uniform(- 1, 1));

      spon_dir.SetZ(randnum6.Uniform(- 1, 1));

      vrec_spon = vrec * spon_dir.Unit();

      v += vrec_abs;

      v += vrec_spon;

      v += g * dt;


      /*****Position changes due to change in velocity*****/

      dt2over2 = pow(dt, 2) / 2;

      r += v * dt;

      r += g * dt2over2;

      rlens = r - clens;


      /*****Atom goes back to ground state*****/

      if (quant_axis == 0) {

	spontaneous_decay(state_atomx, kick_atom);

	rand_state_atom(state_atomy, 0.5);
      } else {

	spontaneous_decay(state_atomy, kick_atom);

	rand_state_atom(state_atomx, 0.5);
      }


      /*****Update position and velocity*****/

      for (q = 0 ; q < 3 ; q ++) {

	pos[q] = r(q) * 100;

	dir[q] = v(q) / v.Mag();
      }

      if (j < nkeep) {

	time_atom.push_back(t);

	x.push_back(r(0));

	y.push_back(r(1));

	z.push_back(r(2));

	vx.push_back(v(0));

	vy.push_back(v(1));

	vz.push_back(v(2));
      }


      /*****Update geometry navigator and probe location of atom*****/

      ini_node = gGeoManager->InitTrack(pos, dir);

      location = gGeoManager->GetPath();
    }

    //		numit_lens->Fill(tracker);

    xyposlens->Fill(r.X(), r.Y());


    /*****LVIS propagation to cavity*****/

    tracker = 0;

    dtp[0] = TMath::Abs((zcav - r.Z()) / v.Z()) / sampling;

    dtp[1] = 2 * rcavmod / TMath::Abs(v.Z()) / 4;

    int prec = 0;				// Precision chooser

    double radxz;				// Atom radial position from cavity center reference frame

    while (r.Z() < zcav + mirdiam / 2 && r.Mag() < molcham_rmin) {

      tracker++;

      if (j < nkeep) {

	time_atom.push_back(t);

	x.push_back(r(0));

	y.push_back(r(1));

	z.push_back(r(2));

	vx.push_back(v(0));

	vy.push_back(v(1));

	vz.push_back(v(2));
      }

      dtp2over2 = pow(dtp[prec], 2) / 2;

      t += dtp[prec];

      v += g * dtp[prec];

      r += v * dtp[prec];

      r += g * dtp2over2;

      if (r.Z() > zcav - mirdiam / 2) {

	radxz = sqrt(pow(r.Z() - zcav, 2) + pow(r.X(), 2));

	prec = 1;

	if ((r.Y() > - cavspa / 2 && r.Y() < cavspa / 2) && radxz < rcavmod) {

	  vzmodsum += v.Z();

	  natoms_mod ++;

	  break;
	}
      }

    }

    //		numit_cav->Fill(tracker);

    if (j < nplot) {

      /*****Graph position vs time per atom*****/

      //			TGraph *xvst = new TGraph(time_atom.size(), &(time_atom[0]), &(x[0]));

      //			TGraph *yvst = new TGraph(time_atom.size(), &(time_atom[0]), &(y[0]));

      //			TGraph *zvst = new TGraph(time_atom.size(), &(time_atom[0]), &(z[0]));

      TGraph *zvsx = new TGraph(time_atom.size(), &(x[0]), &(z[0]));

      //			TGraph *zvsy = new TGraph(time_atom.size(), &(y[0]), &(z[0]));

      //			xvst->SetLineColor(static_cast<int>(randcont * 50));

      //			yvst->SetLineColor(static_cast<int>(randcont * 50));

      //			zvst->SetLineColor(static_cast<int>(randcont * 50));

      zvsx->SetLineColor(static_cast<int>(randcont * 50));

      //			zvsy->SetLineColor(static_cast<int>(randcont * 50));

      //			x_vs_t->Add(xvst);

      //			y_vs_t->Add(yvst);

      //			z_vs_t->Add(zvst);

      z_vs_x->Add(zvsx);

      //			z_vs_y->Add(zvsy);


      /*****Graph velocity vs time per atom*****/

      //			TGraph *vxvst = new TGraph(time_atom.size(), &(time_atom[0]), &(vx[0]));

      //			TGraph *vyvst = new TGraph(time_atom.size(), &(time_atom[0]), &(vy[0]));

      //			TGraph *vzvst = new TGraph(time_atom.size(), &(time_atom[0]), &(vz[0]));

      //			TGraph *vxvsz = new TGraph(time_atom.size(), &(z[0]), &(vx[0]));

      //			TGraph *vyvsz = new TGraph(time_atom.size(), &(z[0]), &(vy[0]));

      //			vxvst->SetLineColor(static_cast<int>(randcont * 50));

      //			vyvst->SetLineColor(static_cast<int>(randcont * 50));

      //			vzvst->SetLineColor(static_cast<int>(randcont * 50));

      //			vxvsz->SetLineColor(static_cast<int>(randcont * 50));

      //			vyvsz->SetLineColor(static_cast<int>(randcont * 50));

      //			vx_vs_t->Add(vxvst);

      //			vy_vs_t->Add(vyvst);

      //			vz_vs_t->Add(vzvst);

      //			vx_vs_z->Add(vxvsz);

      //			vy_vs_z->Add(vyvsz);
    }


    if (j < nkeep) {

      size_t stopper = 0;		// Array positioning

      /*****Calculating atomic beam averages as a function of z-position*****/

      for (zbin = 0 ; zbin < sampling ; zbin++) {

	while (stopper < z.size()) {

	  if (z[stopper] > height[zbin] - hstep / 2 && z[stopper] <= height[zbin] + hstep / 2) {

	    rtrms[zbin] += pow(x[stopper], 2) + pow(y[stopper], 2);

	    ravg[zbin] += sqrt(pow(x[stopper], 2) + pow(y[stopper], 2) + pow(z[stopper], 2));

	    xavg[zbin] += x[stopper];

	    yavg[zbin] += y[stopper];

	    vzavg[zbin] += vz[stopper];

	    vavg[zbin] += sqrt(pow(vx[stopper], 2) + pow(vy[stopper], 2) + pow(vz[stopper], 2));

	    entries[zbin]++;

	    break;
	  }

	  stopper++;
	}
      }


      /*****Clear vector containers for next atom*****/

      time_atom.clear();

      x.clear();

      y.clear();

      z.clear();

      vx.clear();

      vy.clear();

      vz.clear();
    }

  }

  /*****Calculation and graphing of atomic beam position averages as function of height*****/

  int stopmol = 0;

  int stoplens = 0;

  int outmol = 0;

  int outlens = 0;

  for (zbin = 0 ; zbin < sampling ; zbin++) {

    rtrms[zbin] = sqrt(rtrms[zbin] / entries[zbin]);

    ravg[zbin] = ravg[zbin] / entries[zbin];

    xavg[zbin] = xavg[zbin] / entries[zbin];

    yavg[zbin] = yavg[zbin] / entries[zbin];

    vzavg[zbin] = vzavg[zbin] / entries[zbin];

    vavg[zbin] = vavg[zbin] / entries[zbin];

    if (height[zbin] > molwidth + cmol.Z() && stopmol == 0) {

      outmol = zbin;

      stopmol = 1;
    }

    if (height[zbin] > lenswidth + clens.Z() && stoplens == 0) {

      outlens = zbin;

      stoplens = 1;
    }
  }

  TGraph *atom_trms = new TGraph(sampling, height, rtrms);

  double themol_avg = acos(vzavg[outmol] / vavg[outmol]); 		// Average polar angle of atoms after molasses

  double thelens_avg = acos(vzavg[outlens] / vavg[outlens]);		// Average polar angle of atoms after lens


  /*****Timing*****/

  time_t end;			// End time marker

  time (& end);

  double timing = difftime (end, start);

  cout.precision(5);


  if (nmolbeams != 0) {

    /*****Doppler cooling limit*****/

    vt2avgmol = vt2avgmol / doppler_count;

    double Tdoppler = matom * vt2avgmol / (2 * kb);


    /*****2D molasses outputs*****/

    cout << "Minimum polar angle to clear cavity : " << thetam * 180 / pi << endl;

    cout << "In vacuo tube length required : " << tube << endl;

    cout << "Average molasses beam spot size : " << molwidth << endl;

    cout << "Atoms that enter molasses : " << (double) natoms_mol / natoms * 100 << "%" << endl;

    cout << "Average transverse position of atoms after molasses : (" << xavg[outmol] << ", " << yavg[outmol] << ")" << endl;

    cout << "Average polar angle of atoms after molasses : " << themol_avg * 180 / pi << endl;

    cout << "Doppler speed : " << sqrt(vt2avgmol) * 100 << " cm/s" << endl;

    cout << "Doppler temperature : " << Tdoppler * 1e6 << " uK" << endl;

    cout << "Average scattering rate : " << avg_sr / sr_count << endl;

    cout << "Average rate of photons scattered : " << phot_scatt / natoms_mol << endl;


    /*****File output: average trajectory*****/

    ofstream avgtraj;

    stringstream avgtraj_filename;

    avgtraj_filename << "avgtraj_" << region << "_" << par_type << ite << ".txt"; 

    avgtraj.open(avgtraj_filename.str().c_str(), ios::out | ios::trunc);

    for (zbin = 0 ; zbin < sampling ; zbin++) {

      avgtraj << height[zbin] << "  " << xavg[zbin] << "  " << yavg[zbin] << "  " << rtrms[zbin] << "\n";
    }

    avgtraj.close();

    if (nlensbeams != 0) {

      /*****Find atomic beam focus*****/

      vector<double> focus_range(rtrms + outmol, rtrms + (sampling - 1));

      vector<double>::iterator min_pos;

      min_pos = min_element(focus_range.begin(), focus_range.end());

      int actual_min_pos = distance(focus_range.begin(), min_pos) + outmol;


      /*****Laser lens outputs*****/

      cout << "Average lens beam spot size : " << lenswidth << endl;

      cout << "On-resonance saturation parameter at waist for laser lens : " << s0lens[0] << endl;

      cout << "Atoms that enter laser lens : " << (double) natoms_lens / natoms * 100 << "%" << endl;

      cout << "Average transverse position of atoms after lens : (" << xavg[outlens] << ", " << yavg[outlens] << ")" << endl;

      cout << "Average polar angle of atoms after lens : " << thelens_avg * 180 / pi << endl;

      cout << "Position of the minimum : " << height[actual_min_pos] << endl;

      cout << "Minimum standard deviation : " << *min_pos << endl;


      /*****File output: beam parameters*****/

      ofstream par_expl;

      stringstream par_filename;

      par_filename << "par_" << region << "_" << par_type << ".txt";

      par_expl.open(par_filename.str().c_str(), ios::out | ios::app);

      par_expl << w0mol[0] << "  " << w0molpos[0].Mag() << "  " << s0mol[0] << "  " << (double) detunmol[0] / gam << "  " << themol_avg << "  " << w0lens[0] << "  " << w0lenspos[0].Mag() << "  " << s0lens[0] << "  " << (double) detunlens[0] / gam << "  " << thelens_avg << "  " << rtrms[outlens] << "  " << (double) natoms_mod / natoms * 100 << "\n";

      par_expl.close();
    }
  } else

    cout << "Geometrical bound on # of atoms : " << (double) benchmark / natoms * 100 << "%" << endl;

  cout << "Atoms that go through : " << (double) natoms_mod / natoms * 100 << "%" << endl;

  if (natoms_mod != 0)

    cout << "Average Vz of atoms going through : " << (double) (vzmodsum / natoms_mod) << " m/s" << endl;

  int hours = (int) (timing / 3600);

  int minutes = (int) ((timing - hours * 3600) / 60);

  int remain = (int) (timing - hours * 3600 - minutes * 60);

  cout << "Processing time " << hours << ":" << minutes << ":" << remain << endl;

  cout << "Or simply " << timing << " seconds" << endl;


  /*****Drawing*****/

  TCanvas *zofx = new TCanvas("zofx", "z vs x", 1);

  zofx->cd();

  zofx->SetFillColor(10);

  z_vs_x->Draw("AC");

  //	TCanvas *zofy = new TCanvas("zofy", "z vs y", 1);

  //	zofy->cd();

  //	zofy->SetFillColor(10);

  //	z_vs_y->Draw("AC");


  //	TCanvas *vxofz = new TCanvas("vxofz", "vx vs z", 1);

  //	vxofz->cd();

  //	vxofz->SetFillColor(10);

  //	vx_vs_z->Draw("AC");


  //	TCanvas *vyofz = new TCanvas("vyofz", "vy vs z", 1);

  //	vyofz->cd();

  //	vyofz->SetFillColor(10);

  //	vy_vs_z->Draw("AC");


  TCanvas *atom_rtrms = new TCanvas("atom_rtrms", "Transverse rms position vs z", 1);

  atom_rtrms->cd();

  atom_rtrms->SetFillColor(10);

  atom_trms->Draw("AC");


  TCanvas *xy2d = new TCanvas("xy2d", "x-y position histogram", 1);

  xy2d->cd();

  xy2d->SetFillColor(10);

  xypos->Draw();


  TCanvas *xy2dlens = new TCanvas("xy2dlens", "x-y position histogram (lens)", 1);

  xy2dlens->cd();

  xy2dlens->SetFillColor(10);

  xyposlens->Draw();


  //	TCanvas *vxmothist = new TCanvas("vxmothist", "x velocity histogram after MOT", 1);

  //	vxmothist->cd();

  //	vxmothist->SetFillColor(10);

  //	vxmot->Draw();

  //	TCanvas *vymothist = new TCanvas("vymothist", "y velocity histogram after MOT", 1);

  //	vymothist->cd();

  //	vymothist->SetFillColor(10);

  //	vymot->Draw();


  //	TCanvas *xmolhist = new TCanvas("xmolhist", "x histogram after molasses", 1);

  //	xmolhist->cd();

  //	xmolhist->SetFillColor(10);

  //	xmol->Draw();

  //	TCanvas *ymolhist = new TCanvas("ymolhist", "y histogram after molasses", 1);

  //	ymolhist->cd();

  //	ymolhist->SetFillColor(10);

  //	ymol->Draw();


  TCanvas *vxmolhist = new TCanvas("vxmolhist", "x velocity histogram after molasses", 1);

  vxmolhist->cd();

  vxmolhist->SetFillColor(10);

  vxmol->Draw();

  //	TCanvas *vymolhist = new TCanvas("vymolhist", "y velocity histogram after molasses", 1);

  //	vymolhist->cd();

  //	vymolhist->SetFillColor(10);

  //	vymol->Draw();

  TCanvas *vzmolhist = new TCanvas("vzmolhist", "z velocity histogram after molasses", 1);

  vzmolhist->cd();

  vzmolhist->SetFillColor(10);

  vzmol->Draw();

  TCanvas *vzmoldivhist = new TCanvas("vzmoldivhist", "z velocity histogram after molasses", 1);

  vzmoldivhist->cd();

  vzmoldivhist->SetFillColor(10);

  vzmoldiv->Draw();

  TCanvas *vxmoldivhist = new TCanvas("vxmoldivhist", "x velocity histogram after molasses", 1);

  vxmoldivhist->cd();

  vxmoldivhist->SetFillColor(10);

  vxmoldiv->Draw();

  TCanvas *xmoldivhist = new TCanvas("xmoldivhist", "x histogram after molasses", 1);

  xmoldivhist->cd();

  xmoldivhist->SetFillColor(10);

  xmoldiv->Draw();


  //	TCanvas *nit_mol = new TCanvas("nit_mol", "Number of iterations molasses", 1);

  //	nit_mol->cd();

  //	nit_mol->->SetFillColor(10);

  //	numit_mol->Draw();

  //	TCanvas *nit_lens = new TCanvas("nit_lens", "Number of iterations laser lens", 1);

  //	nit_lens->cd();

  //	nit_lens->SetFillColor(10);

  //	numit_lens->Draw();

  //	TCanvas *nit_cav = new TCanvas("nit_cav", "Number of iterations propagation cavity", 1);

  //	nit_cav->cd();

  //	nit_cav->SetFillColor(10);

  //	numit_cav->Draw();

  gStyle->SetOptStat("em");

  /*****Delete all dynamically allocated arrays*****/

  delete[] detunmol;

  delete[] s0mol;

  delete[] gs0mol;

  delete[] smol;

  delete[] srmol;

  delete[] wmol;

  delete[] w0mol;

  delete[] beamdivermol;

  delete[] detunlens;

  delete[] s0lens;

  delete[] gs0lens;

  delete[] slens;

  delete[] srlens;

  delete[] wlens;

  delete[] w0lens;

  delete[] beamdiverlens;

  delete[] ravg;

  delete[] xavg;

  delete[] yavg;

  delete[] vzavg;

  delete[] vavg;

  delete[] rtrms;

  delete[] entries;

  delete[] height;

  return 0;
}


/*****Personal functions*****/

double rayleigh_length(double wavelength, double waist)
{
  return pi * pow (waist, 2) / wavelength;
}

double spot_size(double wavelength, double waist, double z)
{
  return waist * sqrt(1 + pow(z / rayleigh_length(wavelength, waist), 2));
}

void rand_state_atom(int &atomic_state, double proba)
{
  TRandom3 rnum(0);

  if (rnum.Rndm() < proba)

    atomic_state = 0;

  else

    atomic_state = 1;
}

void spontaneous_decay(int &state_atom, int kick)
{
  if (state_atom == 0 && kick == 0)

    state_atom = 0;

  else {
    if (state_atom == 0 && kick == 1)

      rand_state_atom(state_atom, 1 / 3.);

    else {
      if (state_atom == 1 && kick == 0)

	rand_state_atom(state_atom, 2 / 3.);

      else
	state_atom = 1;
    }
  }
}
