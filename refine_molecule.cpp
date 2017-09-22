/* Copyright (c) Wolfgang Brehm @ DESY 2017                                  */
/* This file is part of a crystallographic refinement program                */
/* All rights reserved                                                       */

// asu is needed, get it at https://github.com/1ykos/asu
// wmath is needed, get it at https://github.com/1ykos/wmath
// LBFGSpp is needed, get it at https://github.com/yixuan/LBFGSpp
// openbabel is needed,
// Eigen headers are needed,
// C++14 and C++17 are supported.

//#define BOOST_THREAD_PROVIDES_FUTURE

#include "wmath.hpp"
#include <asu> // symmetry operators and helper functions
#include <encode>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <exception>
#include <future>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <streambuf>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <LBFGS.h>

/*#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/thread/future.hpp>
#include <boost/thread/thread.hpp>
*/

#include <openbabel/atom.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

//#include "atomic_electron_density.hpp"
//#include "discrete_gaussian_lut.hpp"

#include "reciprocalspace_target_function.hpp"

using std::array;
using std::async;
using std::bitset;
using std::cerr;
using std::cin;
using std::cout;
using std::enable_if;
using std::endl;
using std::exception;
using std::fill;
using std::fixed;
using std::future;
using std::get;
using std::ifstream;
using std::is_floating_point;
using std::isnan;
using std::launch;
using std::make_shared;
using std::max;
using std::min;
using std::numeric_limits;
using std::ref;
using std::scientific;
using std::setprecision;
using std::setw;
using std::shared_ptr;
using std::sprintf;
using std::stod;
using std::stoi;
using std::streambuf;
using std::streamsize;
using std::string;
using std::stringstream;
using std::tanh;
using std::tuple;
using std::unordered_map;
using std::vector;

using SYMMETRY::asu;
using SYMMETRY::asu_hashed;
using SYMMETRY::MillerIndex;
using SYMMETRY::ReciprocalCell;
using SYMMETRY::DirectCell;

//boost::asio::io_service ioService;
//boost::thread_group threadpool;
//boost::asio::io_service::work work(ioService);

using OpenBabel::OBConversion;
using OpenBabel::OBMol;
using OpenBabel::OBForceField;

using Eigen::VectorXd;

using LBFGSpp::LBFGSSolver;
using LBFGSpp::LBFGSParam;

int main(int argc, char *argv[])
{
  std::ios::sync_with_stdio(false);
  const double spacing = 0.5; // Angstroem
  const double min_padding = 8; // Angstroem
  const size_t num_photons = 65535;
  const size_t px = 128;
  const size_t py = 128;
  const size_t pz = 128;
  OBConversion conv(&cin,&cout);
  //cerr << "test" << endl;
  int return_code = conv.SetInFormat("PDB");
  if(!return_code) return return_code;
  //cerr << "so far so good" << endl;
  OBMol molecule;
  return_code = conv.Read(&molecule);
  if(!return_code) return return_code;
  uint16_t* data = new uint16_t[px*py*pz]();
  ifstream datafile(argv[1]);
  datafile.read(reinterpret_cast<char*>(data),px*py*pz*2);
  reciprocalspace_target_function target(molecule,spacing,px,py,pz,data);
  VectorXd x(molecule.NumAtoms()*3);
  VectorXd grad(molecule.NumAtoms()*3);
  {
    size_t i=0;
    for (auto atom=molecule.BeginAtoms(); atom!=molecule.EndAtoms(); ++atom){
      x(i++)=(**atom).GetX();
      x(i++)=(**atom).GetY();
      x(i++)=(**atom).GetZ();
    }
  }
  //target(x,grad);return 0;
  //cerr << "so far so good" << endl;
  /*const double temp = target(x,grad); 
  cerr << grad(1) << endl;
  x(1)+=0.0001;
  cerr << (target(x,grad)-temp)*10000 << endl;
  return 0;*/
  /*const double temp = target(x,grad);
  VectorXd temp_grad(molecule.NumAtoms()*3);
  target(x,temp_grad);
  //cout << temp_grad(0) << " " << temp_grad(1) << " " << temp_grad(2) << endl;
  for (size_t i=0;i!=x.size();++i){
    VectorXd temp_x = x;
    temp_x(i)+=1e-4;
    cout << temp_grad(i) << " " << (target(temp_x,grad)-temp)*1e4 << endl;
  }
  return 0;*/
  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-1;
  param.max_iterations = 1024;

  // Create solver and function object
  LBFGSSolver<double> solver(param);
  double fx;
  solver.minimize(target,x,fx);
  delete[] data;
  return 0;
}
