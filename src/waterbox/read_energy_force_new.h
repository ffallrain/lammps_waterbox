/*
 * read energe force data from database
 *
 * Author: Ping Zeng(zengping10212044@gmail.com)
 */
#ifndef _READ_ENERGY_FORCE_NEW_H_
#define _READ_ENERGY_FORCE_NEW_H_

#include <cmath>
#include <map>
#include <vector>
#include <cassert>

#include "grid_structures.h"
#include "gzstream/gzstream.h"
#include "common.h"

namespace database
{
  class Atom
  {
  public:
    explicit Atom(std::string& line);

    std::string m_name;
    std::vector<double> m_coord;
  };


  class Config
  {
  public:
    explicit Config();
    explicit Config(int n1, int n2);

    double get_prop(std::string name, double w);
    void setEnergy(double energe);
    void setForce(std::array<double, 3> ff);
    void setTorque(std::array<double, 3> tq);
    void setName(std::string name);
    int m_n1;
    int m_n2;
    double m_energy;
    std::string m_name;
    std::vector<Atom> m_fmole1;
    std::vector<Atom> m_fmole2;
    std::vector<Atom> m_xmole2;
    double m_Fx, m_Fy, m_Fz;
    double m_Tx, m_Ty, m_Tz;
  };

  class Molecule
  {
  public:
    explicit Molecule(){}
    void addAtom(std::string& line);

    std::vector<Atom> m_atoms;
    double m_energy; //energe
  };

  class DataBase
  {
  };

  class EnergeForceDatabase
  {
  public:
    explicit EnergeForceDatabase(const char* filename, const char* fragtype);
    ~EnergeForceDatabase();
    Config*& at_all(int dim1, int dim2, int dim3, int dim4);
    void read_file();

    double get_prop(int i, int j, int si, int ni, std::string name, double w, double eheigh);
    const std::vector<Config*>& getEnergy(int i, int j, int k) {
        return m_all_config[i][j][k];
    }

  private:

    // read_mole from the first two line
    void read_mole(int ind, igzstream& datafile);

    std::vector<std::vector<std::vector<std::vector<Config*>>>> m_all_config;
    std::string m_filename;
    std::string m_fragtype;
    DataStructure* m_pDS;
    int m_natoms1;
    int m_natoms2;

    Molecule m_mole1;
    Molecule m_mole2;

  };

  /**
   * Return the weights of the two configurations at one bisector direction.
   *   Linear interpolation is used.
   *   NOTE: The two vectors of the two configurations are calculated here
   *       but in future, they should be retrieve from a pre-calculated
   *       constant data list or dictionary of fixed directions.
   */
  std::pair<double, double> weights_for_2_configs(const vector<double>& vec, const vector<database::Atom> config1, const vector<database::Atom> config2, double cut=0.0000001);

}

#endif
