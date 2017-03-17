/*
 * The QM interpolation
 *
 * Author: Ping Zeng(zengping10212044@gmail.com)
 */
#ifndef _QM_INTERPOLATION_H_
#define _QM_INTERPOLATION_H_

#include <vector>
#include <string>
#include <map>
#include <array>

#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>

#include "common.h"
#include "grid_structures.h"
#include "read_energy_force_new.h"

namespace QM
{
  class Atom
  {
  public:
    std::vector<double> m_x;

    virtual ~Atom(){}

    virtual Atom* copy() { return NULL;}
  };

  // TODO complete the struct
  class GjfAtom : public Atom
  {
  public:
    explicit GjfAtom(std::string line);
    ~GjfAtom(){}

    Atom* copy();

    std::string m_name;
  };

  // TODO complete the struct
  class PdbAtom : public Atom
  {
  public:
    explicit PdbAtom(std::string line);
    ~PdbAtom(){}

    Atom* copy();

    int m_i_atm;
    std::string m_a_name;
    std::string m_a_res;
    std::string m_a_chn;
    int m_i_res;
  };

  class InpAtom : public Atom
  {
  public:
    explicit InpAtom(std::string line);
    InpAtom(const std::vector<double>& points, std::string name);

    ~InpAtom(){}
    Atom* copy();

    std::string m_name;
  };

  class Coordinates
  {
  public:
    explicit Coordinates(int n1, int n2, std::string& fragtype, std::string& name);

    void addAtom(std::string& line, const std::string& ftype);

    void addAtom(std::vector<double>& points, std::string name, const std::string& ftype);

    void ReorientToOrigin(double cut);

    void MirrorAll();

    void MirrorBackProperty();
    void ReorientToOldVec();

    void indexing_auto3();
    void indexing_orient_auto3(int ri);
    void calt_conf_energy(database::EnergeForceDatabase& allconfig, bool isForce, double ehigh);

    double get_interp_energy();

    std::vector<double> get_interp_force();
    std::vector<double> get_interp_torque();

    void reverse_force_toque();

    std::vector<int> m_r_ndxs;
    std::vector<double> m_vbis;
    std::vector<double> m_vnrm;
    std::vector<double>  m_rel_x;

    std::map<int, std::vector<int>> m_dgrid_ndx_layer;
    std::map<int, std::vector<int>> m_dtheta_ndx_layer;

    std::vector<double> m_orientVec;
    double m_orient_ang1;
    double m_orient_ang2;

    std::map<int, database::OrientDataStructure*> m_orient_DS;
    std::map<int, double> m_orient_tr;
    std::map<int, std::map<int, std::vector<int>>> m_dgrid_sub_ndx;
    std::map<int, std::map<int, std::vector<int>>> m_dtheta_ndx;

    std::vector<double> m_force;
    std::vector<double> m_torque;

    std::map<std::string, double> m_properties;
    int m_n1;
    int m_n2;
    std::string m_fragtype;
    std::string m_name;
    database::DataStructure* m_pDS;
    bool m_is_oriented;
    std::map<std::string, int> m_facendx;
    std::vector<int> m_symm;
    int m_center;
    double m_r;
    double m_ang1;
    double m_ang2;
    std::vector<double> m_center2;
    std::vector<double> m_origin_com;
    std::vector<boost::shared_ptr<Atom> > m_atoms;
    std::vector<int> m_operateNdx;
    std::vector<std::pair<std::vector<double>, double>> m_operations;

    std::vector<double> m_center_coord_normalized;
    bool m_exit_before = false;

  private:
    void spherical_x();
    void spherical_orient();

    std::vector<double> get_com(const std::vector<boost::shared_ptr<Atom>>& atoms) {
        double totalM = 0;
        std::vector<double> x = {0, 0, 0};
        for (int i = 0; i < atoms.size(); i++) {
            for (int k = 0; k < 3; k++) {
                x[k] += atoms[i]->m_x[k] * TMASS[i];
            }
            totalM += TMASS[i];
        }
        for (int k = 0; k < 3; k++) {
            x[k] /= totalM;
        }
        return x;
    }


    std::vector<double> norm_prob(const std::vector<database::Atom>& config,
            std::vector<int> ndx, std::string prob="wtr") {
        std::vector<double> vec;
        if (prob == "wtr") {
            vec = get_normal_unit(subtraction(config[ndx[1]].m_coord, config[ndx[0]].m_coord),
                    subtraction(config[ndx[2]].m_coord, config[ndx[0]].m_coord));
        }
        return vec;
    }

    void insert_map_element_2_map(const std::map<int, std::vector<int>>& element,
            const int key,
            std::map<int, std::map<int, std::vector<int>>>& m) {
        m[key] = {};
        for (const auto& item : element) {
            m[key][item.first] = item.second;
        }
    }

    vector<int>  get_map_keys(const std::map<int, std::vector<int>>& m) {
        vector<int> res;
        for (const auto& item : m) {
            res.push_back(item.first);
        }
        return res;
    }

    // get index of a struct
    int get_index(double r, const vector<double>& vec);

    int index(int r, const vector<int>& vec) {
        int res = -1;
        for (auto item : vec) {
            res ++;
            if (r == item) {
                return res;
            }
        }
        assert(res != vec.size());
        return res;
    }

    constexpr static double R2D = 180.0 / 3.14159265358979;
    constexpr static double PI4 = 0.78539816339744817;
    const std::vector<double> TMASS = {15.999, 1.008, 1.008};
  };

  class QMInterpolation
  {
      public:
          explicit QMInterpolation(std::string fftype, database::EnergeForceDatabase& allconfig);

          std::vector<std::string> process(std::string filename);

          void calculate(const std::vector<std::vector<double>>& lhs,
                         const std::vector<std::vector<double>>& rhs);

          std::vector<double> m_force;
          std::vector<double> m_torque;
          double m_energy;
      private:
          std::string m_fftype;
          database::EnergeForceDatabase& m_allconfig;

          const static std::vector<int> m_aa_ndx;
          const static std::vector<int> m_prob_ndx;

  };
}
#endif
