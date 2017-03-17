/*
 * define some constant of QM calc
 *
 * Author: Ping Zeng(zengping10212044@gmail.com)
 */
#ifndef _GRID_STRUCTURES_H_
#define _GRID_STRUCTURES_H_

#include <cmath>
#include <map>
#include <vector>
#include <cassert>

namespace database
{
  typedef union dataStructureParam
  {
    double dim1;
    double dim2[2];
  } DataStructureParam;

  const double D2R = 3.14159265358979/180.0;
  const std::map<char, double> Mass_table = {{'O', 15.999}, {'H', 1.008}};

  class DataStructure;
  typedef double (DataStructure::*Fun_lj)(double, double, double);

  class WaterStructure;

  class DataStructure
  {
  public:
    explicit DataStructure(const char* fftype);

    void initialize();

    void set_grid_data(const char* structure_type);

    std::map<int, std::vector<int> > m_grid_data;
    std::vector<std::string> m_symface;
    int m_n1;
    int m_n2;
    std::vector<double> m_R_NDX;
    std::vector<double> m_DR;
    std::vector<double> m_PHI_angles;
    std::map<int, std::vector<double> > m_THETA_angles;
    std::vector<int> m_nConf;
    std::vector<int> m_nNorm;
    std::vector<int> m_NTheta;
    int m_nGrid;
    const char* m_fftype;

    virtual void set_symmetry();
    virtual void set_num_of_atoms();
    virtual void set_R();
    virtual void set_phi();
    virtual void set_theta();
    virtual void set_nConf();
    virtual void degree2radius();

    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec1(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);
    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec2(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);
    virtual std::vector<double> calt_dvec(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);

  private:
  };

  class PrpStructure: public DataStructure
  {
  public:
    explicit PrpStructure(const char* fftype):DataStructure(fftype){}
    void set_theta();
    void set_symmetry();
    void set_num_of_atoms();

    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec1(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);
    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec2(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);
    virtual std::vector<double> calt_dvec(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);
  };

  class WaterStructure : public PrpStructure
  {
  public:
    explicit WaterStructure(const char* fftype):PrpStructure(fftype){}
    void set_num_of_atoms();
    void set_theta();

    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec1(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);
    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec2(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);
    virtual std::vector<double> calt_dvec(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c);

    //virtual void set_num_of_atoms();
  private:
  };


  DataStructure* getDataStructure(const char* fragtype);

  class OrientDataStructure
  {
  public:
    explicit OrientDataStructure(const char* fftype);

    void initialize();

    std::vector<double> m_PHI_angles;
    std::map<int, std::vector<double> > m_THETA_angles;
    std::vector<int> m_NTheta;
    int m_nGrid;
    const char* m_fftype;

    virtual void set_phi();
    virtual void set_theta();
    virtual void degree2radius();

  private:
  };

  class OrientWaterStructure2 : public OrientDataStructure {
      public:
          explicit OrientWaterStructure2(const char* fftype): OrientDataStructure(fftype){}
          void set_theta();
          void set_phi();
  };

  class OrientWaterStructure3 : public OrientDataStructure {
      public:
          explicit OrientWaterStructure3(const char* fftype): OrientDataStructure(fftype){}
          void set_theta();
          void set_phi();
  };

  OrientDataStructure* getOrientDataStructure2(const char* fragtype);
  OrientDataStructure* getOrientDataStructure3(const char* fragtype);
}
#endif
