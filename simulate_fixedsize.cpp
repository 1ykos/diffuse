/* Copyright (c) Wolfgang Brehm 2016                                         */
/* This file is part of a scaling and merging progarm for XFEL data          */
/* All rights reserved                                                       */

// the asu library is needed, get it at https://github.com/1ykos/asu
#define BOOST_THREAD_PROVIDES_FUTURE

#include <asu> // symmetry operators and helper functions
#include <encode>

// dlib is needed, get it at http://dlib.net/
#include <dlib/optimization.h>
#include <dlib/matrix.h>

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

#include <fftw3.h>

/*#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/thread/future.hpp>
#include <boost/thread/thread.hpp>
*/

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>

//#include "atomic_electron_density.hpp"
//#include "discrete_gaussian_lut.hpp"

namespace{
  constexpr double pi = 3.14159265358979323846;
  constexpr double LOG_DBL_MAX = 709.78271289338402993962517939506;
}

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
using std::scientific;
using std::is_floating_point;
using std::isnan;
using std::launch;
using std::make_shared;
using std::max;
using std::min;
using std::numeric_limits;
using std::ref;
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

int main(int argc, char *argv[])
{
  std::ios::sync_with_stdio(false);
  const double spacing = 0.5; // Angstroem
  const double min_padding = 8; // Angstroem
  const size_t num_photons = 65535;
  const size_t dim = 128;
  const size_t dim2= dim*dim;
  const size_t dim3= dim*dim2;
  OBConversion conv(&cin,&cout);
  cerr << "test" << endl;
  int return_code = conv.SetInFormat("PDB");
  if(!return_code) return return_code;
  cerr << "so far so good" << endl;
  OBMol molecule;
  return_code = conv.Read(&molecule);
  if(!return_code) return return_code;
  double min_x = numeric_limits<double>::max();
  double min_y = numeric_limits<double>::max();
  double min_z = numeric_limits<double>::max();
  double max_x = numeric_limits<double>::lowest();
  double max_y = numeric_limits<double>::lowest();
  double max_z = numeric_limits<double>::lowest();
  for (auto atom=molecule.BeginAtoms(); atom!=molecule.EndAtoms(); ++atom){
    const double x=(**atom).GetX();
    const double y=(**atom).GetY();
    const double z=(**atom).GetZ();
    min_x=x<min_x?x:min_x;min_y=x<min_y?y:min_y;min_z=x<min_z?z:min_z;
    max_x=x>max_x?x:max_x;max_y=y>max_y?y:max_y;max_z=z>max_z?z:max_z;
  }
  cerr << "Molecule has " << molecule.NumHvyAtoms() << " heavy atoms" << endl;
  cerr << "Fourier transform of size " << dim3 << " incoming ... " << endl;
  fftw_complex* in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*dim3);
  fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*dim3);   
  fftw_plan plan=fftw_plan_dft_3d(dim,dim,dim,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fill( in[0], in[0]+dim3,0.0);
  fill(out[0],out[0]+dim3,0.0);
  for (auto atom=molecule.BeginAtoms(); atom!=molecule.EndAtoms(); ++atom){
    const double x=(**atom).GetX();
    const double y=(**atom).GetY();
    const double z=(**atom).GetZ();
    const int32_t px=round((x-min_x+min_padding)/spacing);
    const int32_t py=round((y-min_y+min_padding)/spacing);
    const int32_t pz=round((z-min_z+min_padding)/spacing);
    cerr << px << " " << py << " " << pz << endl;
    for (int dz = -12;dz!=13;++dz){
      for (int dy = -12;dy!=13;++dy){
        for (int dx = -12;dx!=13;++dx){
          const size_t n=(px+dx)
                        +(py+dy)*dim
                        +(pz+dz)*dim2;
          in[n][0]+=exp(-0.25*dx*dx-0.25*dy*dy-0.25*dz*dz)
                   *(**atom).GetAtomicNum()-(**atom).GetFormalCharge();
        }
      }
    }
  }
  //cerr << "first element in input array = " << in[0][0]*in[0][0]+in[0][1]*in[0][1] << endl;
  //in[0][0]=0.0;
  //in[0][1]=0.0;
  fftw_execute(plan); 
  //cout << scientific;
  //cout << setprecision(0);
  double total=0.0;
  for (size_t n=0;n!=dim3;++n){
    total+=out[n][0]*out[n][0]+out[n][1]*out[n][1];
  }
  std::mt19937_64 mt(std::random_device{}());
  for (size_t z = 0; z!=dim; ++z){
    //for (size_t z = 0; z!=1; ++z){
    for (size_t y = 0; y!=dim; ++y){
      for (size_t x = 0; x!=dim; ++x){
        const size_t n = x+y*dim+z*dim2;
        const double l = (out[n][0]*out[n][0]+out[n][1]*out[n][1])
                         /total*num_photons;
        std::poisson_distribution<uint16_t> photons(l);
        uint16_t p = photons(mt);
        cout.write(reinterpret_cast<char*>(&p),2);
//        cerr << p << endl;
        //cout << photons(mt) << " ";
      }
      //cout << endl;
    }
  }
  fftw_destroy_plan(plan);
  fftw_free(in);fftw_free(out);
  return 0;
}
