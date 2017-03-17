/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(body,PairBody)

#else

#ifndef LMP_PAIR_BODY_H
#define LMP_PAIR_BODY_H

#include "pair.h"

#include "waterbox/qm_interpolation.h"
#include "waterbox/read_energy_force_new.h"

#include <array>
#include <string>
#include <vector>
#include <python2.7/Python.h>

namespace LAMMPS_NS {

class PairBody : public Pair {
 public:
  PairBody(class LAMMPS *);
  ~PairBody();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  double cut_global;
  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4;

  class AtomVecBody *avec;
  class BodyNparticle *bptr;

  double **discrete;            // list of all sub-particles for all bodies
  int ndiscrete;                // number of discretes in list
  int dmax;                     // allocated size of discrete list
  int *dnum;                    // number of discretes per line, 0 if uninit
  int *dfirst;                  // index of first discrete per each line
  int nmax;                     // allocated size of dnum,dfirst vectors

  void allocate();
  void body2space(int);

  class WaterAtom {
   public:
     std::vector<std::array<double,3>> atoms;
  };

 private:
  std::vector<double> computeForceAndTorque(const WaterAtom& lhs,
                             const WaterAtom& rhs);
  // invoke the python
  class PythonWater {
  public:
    PythonWater() {
      Py_Initialize();
      if (!Py_IsInitialized()) {
        fprintf(stderr, "can't initialize python environment.\n");
        exit(1);
      }

      PyRun_SimpleString("import sys");
      PyRun_SimpleString("sys.path.insert(0, '/home/zengping/work/code/waterbox/QMCPP/QM/python')");
      PyRun_SimpleString("sys.path.insert(0, '/usr/lib64/python2.7/lib-dynload/')");
      PyRun_SimpleString("sys.path.insert(0, '/usr/lib64/python2.7/site-packages/')");
        PyRun_SimpleString("sys.path.insert(0, '/usr/lib64/python2.7/')");
      PyRun_SimpleString("print sys.path");

      pModule_ = PyImport_ImportModule("interp");
      if(!pModule_) {
        PyErr_Print();
        exit(1);
      }
      pDict_ = PyModule_GetDict(pModule_);
      if(!pDict_) {
        fprintf(stderr, "can't find dictionary.\n");
        exit(1);
      }

      pClassQM_ = PyDict_GetItemString(pDict_, "QMInterpolation");
      if(!pClassQM_) {
        fprintf(stderr, "can't find QMInterpolation class.\n");
        exit(1);
      }

      // Construct QM Instance
      PyObject* tmp = PyString_FromString("wtr");
      pInstanceArgs_ = PyTuple_Pack(1, tmp);
      if (!pInstanceArgs_) {
        PyErr_Print();
        exit(1);
      }
      Py_DECREF(tmp);

      pInstanceQM_ = PyInstance_New(pClassQM_, pInstanceArgs_, NULL);
      if(!pInstanceQM_)
      {
        PyErr_Print();
        fprintf(stderr, "can't create QM instance.\n");
        exit(1);
      }
    }

    std::vector<double> Calculate(const WaterAtom& lhs,
                                  const WaterAtom& rhs) {
      assert(lhs.atoms.size() == 3);
      assert(rhs.atoms.size() == 3);
      PyObject* lhs_p = Py_BuildValue(
        "[(s, [dddd]), (s, [dddd]), (s, [dddd])]",
        "O", 0.0, lhs.atoms[0][0], lhs.atoms[0][1], lhs.atoms[0][2],
        "H", 0.0, lhs.atoms[1][0], lhs.atoms[1][1], lhs.atoms[1][2],
        "H", 0.0, lhs.atoms[2][0], lhs.atoms[2][1], lhs.atoms[2][2]);
      if(!lhs_p) {
        PyErr_Print();
        exit(1);
      }

      PyObject* rhs_p = Py_BuildValue(
        "[(s, [dddd]), (s, [dddd]), (s, [dddd])]",
        "O", 0.0, rhs.atoms[0][0], rhs.atoms[0][1], rhs.atoms[0][2],
        "H", 0.0, rhs.atoms[1][0], rhs.atoms[1][1], rhs.atoms[1][2],
        "H", 0.0, rhs.atoms[2][0], rhs.atoms[2][1], rhs.atoms[2][2]);
      if(!rhs_p) {
        PyErr_Print();
        exit(1);
      }
      PyObject *pReturn = NULL;
      pReturn = PyObject_CallMethod(pInstanceQM_,
                                    (char*)"calculate",
                                    (char*)"(OO)", lhs_p, rhs_p);

      if (!pReturn) {
        PyErr_Print();
        exit(1);
      }
      std::vector<double> values(6);
      if (!PyArg_ParseTuple(pReturn, "dddddd", &values[0],
                            &values[1], &values[2], &values[3],
                            &values[4], &values[5])) {
        fprintf(stderr, "Failed to parse the return val\n");
        exit(1);
      }

      Py_DECREF(lhs_p);
      Py_DECREF(rhs_p);
      Py_DECREF(pReturn);
      return values;

    }

    ~PythonWater() {
      Py_DECREF(pInstanceArgs_);
      Py_DECREF(pInstanceQM_);
      Py_DECREF(pClassQM_);
      Py_DECREF(pDict_);
      Py_DECREF(pModule_);
      Py_Finalize();
    }

  private:
    PyObject* pModule_;
    PyObject* pDict_;
    PyObject* pClassQM_;
    PyObject* pInstanceQM_;
    PyObject* pInstanceArgs_;
  };

  class WaterDb {
   public:
    WaterDb() {
      std::cout << "initialize water db" << std::endl;
      std::string model = "wtr";
      std::string database_file = "/home/fuqiuyu/work/QEM_MD/body/sethbrin/QM/version2/Dimer_deltaEForceTorque_data_mp2_wtr_wtr.txt.gz";
      energe_ = new database::EnergeForceDatabase(database_file.c_str(), model.c_str());
      energe_->read_file();
      interpolation_ = new QM::QMInterpolation(model, *energe_);
      std::cout << "initialize water db finished" << std::endl;
    }

    std::vector<double> Calculate(const WaterAtom& lhs,
                                  const WaterAtom& rhs) {
      std::vector<std::vector<double>> lhs_point = {
        {lhs.atoms[0][0], lhs.atoms[0][1], lhs.atoms[0][2]},
        {lhs.atoms[1][0], lhs.atoms[1][1], lhs.atoms[1][2]},
       {lhs.atoms[2][0], lhs.atoms[2][1], lhs.atoms[2][2]}};
      std::vector<std::vector<double>> rhs_point = {
          {rhs.atoms[0][0], rhs.atoms[0][1], rhs.atoms[0][2]},
           {rhs.atoms[1][0], rhs.atoms[1][1], rhs.atoms[1][2]},
           {rhs.atoms[2][0], rhs.atoms[2][1], rhs.atoms[2][2]}};

      interpolation_->calculate(lhs_point, rhs_point);
      std::vector<double> res;
      res.insert(res.end(), interpolation_->m_force.begin(), interpolation_->m_force.end());
      res.insert(res.end(), interpolation_->m_torque.begin(), interpolation_->m_torque.end());

      return res;
    }

    ~WaterDb() {
      delete interpolation_;
      delete energe_;
    }

   private:
    database::EnergeForceDatabase* energe_;
    QM::QMInterpolation* interpolation_;
  };


 private:
  //PythonWater python_water_db_;
    WaterDb water_db_;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair body requires atom style body

Self-explanatory.

E: Pair body requires body style nparticle

This pair style is specific to the nparticle body style.

*/
