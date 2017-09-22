/* Copyright (c) Wolfgang Brehm @ DESY 2017                                  */
/* This file is part of a crystallographic refinement program                */
/* All rights reserved                                                       */

#ifndef RECIPROCALSPACE_TARGET_FUNCTION_H
#define RECIPROCALSPACE_TARGET_FUNCTION_H

//#include "discrete_gaussian.hpp"
#include "electron_density_table.hpp"
#include "atom.hpp"

#include <Eigen/Dense>

#include <openbabel/atom.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>

#include "wmath.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <iostream>
#include <limits>
#include <random>
#include <string>

namespace{
  constexpr double pi          = 3.14159265358979323846;
  constexpr double LOG_DBL_MAX = 709.78271289338402993962517939506;
  constexpr double R           = 8.3144598;
  using std::max_element;
  using std::floor;
  using std::abs;
  using std::vector;
  using std::complex;
  using std::cout;
  using std::cerr;
  using std::cin;
  using std::endl;
  using std::numeric_limits;
  using OpenBabel::OBConversion;
  using OpenBabel::OBMol;
  using OpenBabel::OBForceField;
  using OpenBabel::OBConformerData;
  using OpenBabel::OBConversion;
  using OpenBabel::vector3;
  using std::fill;
}

constexpr double structure_factor(
    const uint8_t& z,
    const uint8_t& e,
    const  double& x){
  const size_t i = proton_electron_table_index(z,e);
  const double a0 = coeff_a[i][0];
  const double a1 = coeff_a[i][1];
  const double a2 = coeff_a[i][2];
  const double a3 = coeff_a[i][3];
  const double b0 = coeff_bp[i][0]; // b²/2
  const double b1 = coeff_bp[i][1];
  const double b2 = coeff_bp[i][2];
  const double b3 = coeff_bp[i][3];
  return a0*exp(-b0*x*x) // maybe its b0/(2pi)
        +a1*exp(-b1*x*x) // TODO find out
        +a2*exp(-b2*x*x) // too tired...
        +a3*exp(-b3*x*x);// probaly correct
}

constexpr double diff_structure_factor(
    const uint8_t& z,
    const uint8_t& e,
    const  double& x){
  const size_t i = proton_electron_table_index(z,e);
  const double a0 = coeff_a[i][0];
  const double a1 = coeff_a[i][1];
  const double a2 = coeff_a[i][2];
  const double a3 = coeff_a[i][3];
  const double b0 = coeff_bp[i][0]; // b²/2
  const double b1 = coeff_bp[i][1];
  const double b2 = coeff_bp[i][2];
  const double b3 = coeff_bp[i][3];
  return 2*a0*b0*x*exp(-b0*x*x)
        +2*a1*b1*x*exp(-b1*x*x)
        +2*a2*b1*x*exp(-b2*x*x)
        +2*a3*b1*x*exp(-b3*x*x);
}

constexpr double partial_diff_structure_factor(
    const uint8_t& z,
    const uint8_t& e,
    const  double& x){
  const size_t i = proton_electron_table_index(z,e);
  const double a0 = coeff_a[i][0];
  const double a1 = coeff_a[i][1];
  const double a2 = coeff_a[i][2];
  const double a3 = coeff_a[i][3];
  const double b0 = coeff_bp[i][0]; // b²/2
  const double b1 = coeff_bp[i][1];
  const double b2 = coeff_bp[i][2];
  const double b3 = coeff_bp[i][3];
  return 2*a0*b0*exp(-b0*x*x)
        +2*a1*b1*exp(-b1*x*x)
        +2*a2*b1*exp(-b2*x*x)
        +2*a3*b1*exp(-b3*x*x);
}

class reciprocalspace_target_function{
  private:
  const double& spacing;
  const size_t px;
  const size_t py;
  const size_t pz;
  const uint16_t* data;
  OpenBabel::OBMol molecule;
  OpenBabel::OBForceField* pFF;
  std::mt19937_64 mt;
  double min_score=numeric_limits<double>::max();
  size_t num_calls=0;
  public:
  reciprocalspace_target_function(const OpenBabel::OBMol& mol,
                                  const double& spacing,
                                  const size_t& px,
                                  const size_t& py,
                                  const size_t& pz,
                                  const uint16_t* data)
                                  :molecule(mol),spacing(spacing),
                                   px(px),py(py),pz(pz),data(data)
  {
    //cerr << "setting up target function" << endl;
    pFF = OpenBabel::OBForceField::FindForceField("GAFF");
    pFF->SetLogFile(&std::cerr);
    pFF->SetLogLevel(OBFF_LOGLVL_NONE);
    pFF->Setup(molecule);
    //min_score = numeric_limits<double>::max();
    //cerr << "target function set up" << endl;
  }
  double const operator()(const Eigen::VectorXd& x,Eigen::VectorXd& grad){
    cerr << "reciprocalspace_target_function::operator()" << endl;
    ++num_calls;
    grad = Eigen::VectorXd::Zero(x.size());
    {
      size_t i=0;
      auto atom=molecule.BeginAtoms();
      while (atom!=molecule.EndAtoms()){
        (**atom).SetVector(x(i),x(i+1),x(i+2));
        ++atom;
        i+=3;
      }
    }
    //std::mt19937_64 mt;
    double score;
    pFF->Setup(molecule);
    //cerr << "molecule energy= " << pFF->Energy() << endl;
    pFF->GetCoordinates(molecule);
    //cerr << "pFF->GetCoordinates(molecule)" << endl;
    OBConformerData *cd =
    dynamic_cast<OBConformerData*>
    (molecule.GetData(OpenBabel::OBGenericDataType::ConformerData));
    //cerr << "survived dynamic cast" << endl;
    if (cd) {
      vector<double> energies = cd->GetEnergies();
      //cerr << "cd->GetEnergies();" << endl;
      //cerr << energies.size() << endl;
      score = pFF->Energy(true);
      //cerr << "energy = " << energies[0] << endl;
      // do something with forces (= negative gradients)
    }
    {
      size_t i=0;
      //vector<vector<vector3> > confForces = cd->GetForces();
      //cerr << confForces.size() << endl;
      vector<vector3> forces = cd->GetForces()[0];
      while (i!=molecule.NumAtoms()){
        grad(3*i+0)-=forces[i][0]*1000/(R*298);
        grad(3*i+1)-=forces[i][1]*1000/(R*298);
        grad(3*i+2)-=forces[i][2]*1000/(R*298);
        ++i;
      }
    }
    //score=0;
    score*=1000/(R*298);
    cerr << "molecule score=" << score << endl;
    //return score;
    complex<double>* dx = new complex<double>[molecule.NumAtoms()];
    complex<double>* dy = new complex<double>[molecule.NumAtoms()];
    complex<double>* dz = new complex<double>[molecule.NumAtoms()];
    for (size_t z=0;z!=pz;++z){
    //cerr << z << endl;
    for (size_t y=0;y!=py;++y){
    //cerr << z << " " << y << endl;
    for (size_t x=0;x!=px;++x){
      const size_t n= x+y*px+z*px*py;
      complex<double> F(0.0,0.0);
      fill(dx,dx+molecule.NumAtoms(),complex<double>(0.0,0.0));
      fill(dy,dy+molecule.NumAtoms(),complex<double>(0.0,0.0));
      fill(dz,dz+molecule.NumAtoms(),complex<double>(0.0,0.0));
      {
      size_t i=0;
      for (auto atom=molecule.BeginAtoms(); atom!=molecule.EndAtoms(); ++atom){
        const double sx   = 2.0*pi/(px*spacing); //2π is an arbitrary constant
        const double sy   = 2.0*pi/(py*spacing);
        const double sz   = 2.0*pi/(pz*spacing);
        const double G[3] = {(double(x)-px/2)*sx,
                             (double(y)-py/2)*sy,
                             (double(z)-pz/2)*sz};
        const double r[3] = {(**atom).GetX(),
                             (**atom).GetY(),
                             (**atom).GetZ()}; 
        const double abs_G = sqrt(G[0]*G[0]+G[1]*G[1]+G[2]*G[2]);
        const double G_r   = G[0]*r[0]+G[1]*r[1]+G[2]*r[2];
        const double f_G   = structure_factor(
                             (**atom).GetAtomicNum(),
                             (**atom).GetAtomicNum()-(**atom).GetFormalCharge(),
                             abs_G);
        F+=complex<double>(f_G*cos(G_r),f_G*sin(G_r));
        dx[i]+=complex<double>(-f_G*G[0]*sin(G_r),f_G*G[0]*cos(G_r));
        dy[i]+=complex<double>(-f_G*G[1]*sin(G_r),f_G*G[1]*cos(G_r));
        dz[i]+=complex<double>(-f_G*G[2]*sin(G_r),f_G*G[2]*cos(G_r));
        ++i;
      }
      }
      const double k = data[n];
      const double c = 1.0;
      F*=c;
      const double l = norm(F);//=F.real()*F.real()+F.imag()*F.imag()
      //const size_t p = std::poisson_distribution<size_t>(l)(mt);
      //cout.write(reinterpret_cast<const char*>(&p),2);
      //cout << p << " ";
      score += lgamma(k+1)+l-(k>0?k*log(l):0);
      for (size_t i=0;i!=molecule.NumAtoms();++i){
        dx[i]*=c;
        dy[i]*=c;
        dz[i]*=c;
        const double dlx = 2*(F.real()*dx[i].real()+F.imag()*dx[i].imag());
        const double dly = 2*(F.real()*dy[i].real()+F.imag()*dy[i].imag());
        const double dlz = 2*(F.real()*dz[i].real()+F.imag()*dz[i].imag());
        grad(3*i+0)+=dlx-(k>0?k*dlx/l:0); // this derivative is wrong :(
        grad(3*i+1)+=dly-(k>0?k*dly/l:0);
        grad(3*i+2)+=dlz-(k>0?k*dlz/l:0);
      }
    }
    //cout << endl;
    }
    }
    cerr << "score= " << score << endl;
    return score;
    if (score<min_score){
      min_score=score;
      std::ofstream stepfile(std::to_string(num_calls)+".mol");
      OBConversion conv(&cin,&stepfile);
      conv.SetOutFormat("MOL");
      conv.Write(&molecule);
    }
    return score;
  }
  //~reciprocalspace_target_function(){}
};

#endif // RECIPROCALSPACE_TARGET_FUNCTION_H
