#pragma once

#include "../flowed-hmc-generalized.h"
// #include "staples.h"

namespace qlat{
  namespace Generalized{

    // =========================
    // for observations
    namespace Obs{

      using VectFunc2=std::function<const std::vector<int>(const int, const int)>;
      using VectFunc3=std::function<const std::vector<int>(const int, const int, const int)>;
      using StapleFunction = std::function<const ColorMatrix(const GaugeField&, const Coordinate&, const int)>;

      // ---------------------------

      double frobenius_sq( const qlat::ColorMatrix& m ){
        const qlat::ColorMatrix mmdagger = m * matrix_adjoint(m);
        return matrix_trace(mmdagger).real();
      }

      template <class T>
      double force_norm(const T& gm_force){
        double squared_norm = 0.0;
        //
        const qlat::Geometry geo = geo_reform(gm_force.geo());
        qfor(index, geo.local_volume(), {
            const qlat::Coordinate xl = geo.coordinate_from_index(index);
            const qlat::Vector<qlat::ColorMatrix> gm_force_v = gm_force.get_elems_const(xl);
            for (int mu = 0; mu < 4; ++mu) {
              squared_norm += 2.0 * frobenius_sq(gm_force_v[mu]) / ( 4.0 * geo.total_volume() );
            }
            //
          }); // end qfor
        qlat::glb_sum(squared_norm);
        return std::sqrt(squared_norm);
      }

      std::vector<int> convert2staple( const VectFunc2& loop, const int mu, const int nu){
        std::vector<int> tmp = loop(mu,nu);
        tmp.erase(tmp.begin());
        return tmp;
      }
      std::vector<int> convert2staple( const VectFunc3& loop, const int mu, const int nu, const int rho){
        std::vector<int> tmp = loop(mu,nu,rho);
        tmp.erase(tmp.begin());
        return tmp;
      }


      // ---------------------------

      double calculate_W_two_dir_no_comm(const GaugeField& gf, const VectFunc2& w,
                                         const bool is_mu_nu_dir_identical = false);
      double calculate_W_three_dir_no_comm(const GaugeField& gf, const VectFunc3& w);
      double calculate_W3r_no_comm(const GaugeField& gf);
      double calculate_W3c_no_comm(const GaugeField& gf);
      double calculate_W4r_no_comm(const GaugeField& gf);
      double calculate_W4c_no_comm(const GaugeField& gf);
      double calculate_W6_no_comm(const GaugeField& gf);
      double calculate_W7_no_comm(const GaugeField& gf);
      //
      ColorMatrix calculate_staple_two_dir_no_comm(const GaugeField& gf,
                                                   const VectFunc2& w,
                                                   const Coordinate& xl,
                                                   const int mu);
      ColorMatrix calculate_staple_three_dir_no_comm(const GaugeField& gf,
                                                     const VectFunc3& loop,
                                                     const Coordinate& xl,
                                                     const int mu);
      ColorMatrix staple_G0_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G1r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G1c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G1_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G2r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G2c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G2_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G3r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G3c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G3_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G4r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G4c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G4_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G5_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G6_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);
      ColorMatrix staple_G7_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0);

      inline void set_flow_all_kernel_mask_mu_no_comm( FieldM<ColorMatrix, 4>& kernel, const GaugeField& gf_ext, const std::vector<FlowType>& flow_types );
      inline void set_flow_all_kernel_coeff_mask_mu_no_comm( FieldM<array<double, 8>, DIMN>& coeff, const GaugeField& gf_ext, const std::vector<FlowType>& flow_types, const ColorMatrixConstants& cmcs );

      std::vector<FlowType> w_type(const int i, const int rect_chair=-1/*rect=0, chair=1*/);

      // ---------------------------

      void observe_W0_W7( std::ofstream& ofs, const GaugeField& gf, const FlowInfo2& fi ){

        // const int max_size = fi.get_max_flow_size();
        const int max_size = 2;
        const Coordinate expand_left(max_size,max_size,max_size,max_size);
        const Coordinate expand_right(max_size,max_size,max_size,max_size);
        const Geometry geo_ext = geo_reform(gf.geo(), 1, expand_left, expand_right);
        GaugeField gf_ext;
        gf_ext.init(geo_ext);
        gf_ext = gf;
        refresh_expanded(gf_ext);

        const VectFunc2 plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1};
        };
        const VectFunc2 rect = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, nu, -mu-1, -nu-1, -nu-1};
        };
        const VectFunc3 chair = [](const int mu, const int nu, const int rho){
          return std::vector<int>{mu, nu, rho, -mu-1, -rho-1, -nu-1};
        };
        const VectFunc2 twisted_rect = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, nu, mu, -nu-1, -mu-1, -nu-1};
        };
        const VectFunc3 twisted_chair = [](const int mu, const int nu, const int rho){
          return std::vector<int>{mu, nu, -mu-1, rho, mu, -rho-1, -mu-1, -nu-1};
        };
        const VectFunc2 twice_plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1, mu, nu, -mu-1, -nu-1};
        };

        // w0
        const double w0 = calculate_W_two_dir_no_comm(gf_ext, plaq, true);
        // w1
        const double w1r = calculate_W_two_dir_no_comm(gf_ext, rect);
        const double w1c = calculate_W_three_dir_no_comm(gf_ext, chair);
        // w2
        const double w2r = calculate_W_two_dir_no_comm(gf_ext, twisted_rect);
        const double w2c = calculate_W_three_dir_no_comm(gf_ext, twisted_chair);
        // w3
        const double w3r = calculate_W3r_no_comm(gf_ext);
        const double w3c = calculate_W3c_no_comm(gf_ext);
        // w4
        const double w4r = calculate_W4r_no_comm(gf_ext);
        const double w4c = calculate_W4c_no_comm(gf_ext);
        // w5
        const double w5 = calculate_W_two_dir_no_comm(gf_ext, twice_plaq, true);
        // w6
        const double w6 = calculate_W6_no_comm(gf_ext);
        // w7
        const double w7 = calculate_W7_no_comm(gf_ext);

        // IO::wb(ofs, w0);
        // IO::wb(ofs, w1r+w1c);
        // // IO::wb(ofs, w1c);
        // IO::wb(ofs, w2r+w2c);
        // // IO::wb(ofs, w2c);
        // IO::wb(ofs, w3r+w3c);
        // // IO::wb(ofs, w3c);
        // IO::wb(ofs, w4r+w4c);
        // // IO::wb(ofs, w4c);
        // IO::wb(ofs, w5);
        // IO::wb(ofs, w6);
        // IO::wb(ofs, w7);
        std::stringstream ss;
        ss << w0 << std::endl
           << w1r << std::endl
           << w1c << std::endl
           << w2r << std::endl
           << w2c << std::endl
           << w3r << std::endl
           << w3c << std::endl
           << w4r << std::endl
           << w4c << std::endl
           << w5 << std::endl
           << w6 << std::endl
           << w7 << std::endl;
        IO::print_to(ofs, ss.str() );
      }


      void sd_matrix
      (
       std::ofstream& ofs,
       const GaugeField& gf,
       const int max_size = 2 // NEED TO EXTEND FOR LARGER KERNELS
       ){
        constexpr int nops = 7;
        // constexpr int nops = 1;
        std::array<Generalized::Obs::StapleFunction, nops> staples{{
            Generalized::Obs::staple_G0_no_comm
            , Generalized::Obs::staple_G1_no_comm
            , Generalized::Obs::staple_G2_no_comm
            , Generalized::Obs::staple_G3_no_comm
            , Generalized::Obs::staple_G4_no_comm
            , Generalized::Obs::staple_G5_no_comm
            // , Generalized::Obs::staple_G6_no_comm
            , Generalized::Obs::staple_G7_no_comm
          }};

        double matrix[nops][nops];

        const Coordinate expand_left(max_size,max_size,max_size,max_size);
        const Coordinate expand_right(max_size,max_size,max_size,max_size);
        const Geometry geo_ext = geo_reform(gf.geo(), 1, expand_left, expand_right);
        GaugeField gf_ext;
        gf_ext.init(geo_ext);
        gf_ext = gf;
        refresh_expanded(gf_ext);
        assert(staples.size()==nops);

        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, DIMN*nops*nops> fd1,fd2,fd3;
        fd1.init(geo);
        fd2.init(geo);
        fd3.init(geo);

        const std::function<int(int,int,int)> idx_loc = [](const int i, const int j, const int mu){
          return DIMN*nops*i + DIMN*j + mu;
        };

        qacc_for(index, geo.local_volume(), {
            const Coordinate xl = geo.coordinate_from_index(index);
            Vector<double> fd1x = fd1.get_elems(xl);
            Vector<double> fd2x = fd2.get_elems(xl);
            Vector<double> fd3x = fd3.get_elems(xl);
            for(int mu=0; mu<DIMN; ++mu){
              ColorMatrix staple_xl_mu[nops];
              for(int i=0; i<nops; ++i) staple_xl_mu[i] = staples[i](gf_ext,xl,mu);

              const ColorMatrix U_xl_mu = gf_ext.get_elem(xl,mu);

              ColorMatrix UG[nops];
              for(int i=0; i<nops; ++i) UG[i] = U_xl_mu * staple_xl_mu[i];

              ColorMatrix UGUG[nops][nops];
              for(int i=0; i<nops; ++i) for(int j=i; j<nops; ++j) UGUG[i][j]= UG[j] * UG[i];

              ColorMatrix GGdagger[nops][nops];
              for(int i=0; i<nops; ++i) for(int j=i; j<nops; ++j){
                  GGdagger[i][j] = staple_xl_mu[j] * matrix_adjoint(staple_xl_mu[i]);
                }
              for(int i=0; i<nops; ++i) for(int j=i; j<nops; ++j) {
                  fd1x[idx_loc(i,j,mu)] = matrix_trace(UGUG[i][j]).real();
                  fd2x[idx_loc(i,j,mu)] = matrix_trace(GGdagger[i][j]).real();
                  fd3x[idx_loc(i,j,mu)] = matrix_trace(UG[j]).imag() * matrix_trace(UG[i]).imag();
                }
            }
          });

        for(int i=0; i<nops; ++i) for(int j=i; j<nops; ++j) {
            double sum = 0.0;
            for (long index = 0; index < geo.local_volume(); ++index) {
              const Coordinate xl = geo.coordinate_from_index(index);
              Vector<double> fd1x = fd1.get_elems(xl);
              Vector<double> fd2x = fd2.get_elems(xl);
              Vector<double> fd3x = fd3.get_elems(xl);
              for(int mu=0; mu<DIMN; ++mu){ // @@@@ DEBUG !!!! WHY FACTOR 4? -> removed
                sum += -fd1x[idx_loc(i,j,mu)] / geo.total_volume();
                sum += fd2x[idx_loc(i,j,mu)] / geo.total_volume();
                sum += - (2.0/3.0) * fd3x[idx_loc(i,j,mu)] / geo.total_volume();
              }
            }
            glb_sum(sum);
            matrix[i][j] = sum;
          }

        {
          for(int i=0; i<nops; ++i){
            for(int j=0; j<nops; ++j) {
              double d = (i<=j? matrix[i][j] : matrix[j][i]);
              IO::wb(ofs, d);
            }
          }
        }
      }


      void write_force
      (
       std::ofstream& ofs,
       const GaugeField& gf,
       const StapleFunction& stapleG,
       // const FlowInfo2& fi
       const int max_size
       ){
        //const int max_size = fi.get_max_flow_size();
        const Coordinate expand_left(max_size,max_size,max_size,max_size);
        const Coordinate expand_right(max_size,max_size,max_size,max_size);
        const Geometry geo_ext = geo_reform(gf.geo(), 1, expand_left, expand_right);
        GaugeField gf_ext;
        gf_ext.init(geo_ext);
        gf_ext = gf;
        refresh_expanded(gf_ext);

        const Geometry geo = geo_reform(gf.geo());
        GaugeMomentum fm;
        fm.init(geo);

        qacc_for(index, geo.local_volume(), {
            const Coordinate xl = geo.coordinate_from_index(index);
            Vector<ColorMatrix> fmx = fm.get_elems(xl);
            for(int mu=0; mu<DIMN; ++mu){
              const ColorMatrix staple_xl_mu = stapleG(gf_ext,xl,mu);
              const ColorMatrix U_xl_mu = gf_ext.get_elem(xl,mu);
              fmx[mu] = -make_tr_less_anti_herm_matrix( U_xl_mu * staple_xl_mu );
            }
          });

        IO::print_matrix_field(ofs, fm);
      }


      // [d^b_{y,nu} W_j] [d^b_{y,nu} d^a_{x,mu} W_k] [d^a_{x,mu} W_i]
      void dw_ddw_dw_no_comm
      (
       std::ofstream& ofs,
       const GaugeField& gf_ext,
       const int nops,
       const std::vector<std::vector<FlowType>>& ws,
       const ColorMatrixConstants& cmcs,
       const Geometry& geo_ext
       ){
        const Geometry geo = geo_reform(gf_ext.geo());

        const std::function<int(const int,const int,const int)> tsr_idx
          = [nops](const int i, const int j, const int k){ return (i*nops + j)*nops + k; };

        std::vector<double> dwj_ddwk_dwi_ijk(nops*nops*nops,0);
        std::vector<FieldM<array<double, 8>, DIMN>> dw_adj_vect_ext(nops);
        for(int j=0; j<nops; ++j) {
          FieldM<array<double, 8>, DIMN> dw_adj_vect;
          set_flow_all_kernel_coeff_mask_mu_no_comm(dw_adj_vect, gf_ext, ws[j], cmcs);
          //
          dw_adj_vect_ext[j].init(geo_ext);
          dw_adj_vect_ext[j] = dw_adj_vect;
          refresh_expanded(dw_adj_vect_ext[j]);
        }

        for(int k=0; k<nops; ++k){
          //
          for(const FlowType& w_k_elem : ws[k] ){
            for(int mu=0; mu<DIMN; ++mu){
              for(int mask=0; mask<w_k_elem.n_masks; ++mask){
                Field<AdjointColorMatrix> ddw_k_elem;
                {
                  Field<array<ColorMatrix, 2> > dw_k_elem;
                  set_dw_mat_mask_mu_no_comm(dw_k_elem, gf_ext, mask, mu, w_k_elem);
                  set_n_mat_mask_mu_no_comm(ddw_k_elem, dw_k_elem, gf_ext, mask, w_k_elem, mu);
                }
                // --------
                // for i,j
                // --------
                for(int j=0; j<nops; ++j){
                  for(int i=0; i<nops; ++i){
                    //
                    FieldM<double, DIMN> sum_field;
                    sum_field.init(geo);
                    set_zero(sum_field);
                    //
                    // for y
                    qacc_for(y_qlat_idx, geo.local_volume(), {
                        const Coordinate yl = geo.coordinate_from_index(y_qlat_idx);
                        const Vector<AdjointColorMatrix> ddw_k_elem_vect = ddw_k_elem.get_elems_const(yl);
                        const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,w_k_elem);
                        //
                        // --- for (x,nu) ---
                        for(int m_prime=0; m_prime<w_k_elem.multiplicity_reduced(mask,y_idx); ++m_prime){
                          Coordinate xl;
                          int nu;
                          set_xl_nu_from_mask_mu_yl_m_prime(xl, nu, mask, mu, yl, m_prime, geo, w_k_elem);
                          // ---
                          const AdjointColorMatrix ddw_k_elem_ab = ddw_k_elem_vect[m_prime] * w_k_elem.multiplicative_const;
                          //
                          const array<double, 8> dw_j = dw_adj_vect_ext[j].get_elem(yl,nu);
                          const Eigen::Matrix<double, 8, 1> dw_j_b = Eigen::Map<const Eigen::Matrix<double, 8, 1>>(&dw_j[0], 8, 1);
                          const array<double, 8> dw_i = dw_adj_vect_ext[i].get_elem(xl,mu);
                          const Eigen::Matrix<double, 8, 1> dw_i_a = Eigen::Map<const Eigen::Matrix<double, 8, 1>>(&dw_i[0], 8, 1);

                          sum_field.get_elem(yl,nu) += dw_j_b.dot( ddw_k_elem_ab.em().transpose() * dw_i_a);
                        } // end for m_prime
                      }); // end for y
                    //
                    // ----------
                    // reduction for y
                    // ----------
                    double tsr_cmp = 0.0;
                    const double* ptr_begin = &(sum_field.field[0]);
                    for(long l=0; l<DIMN*geo.local_volume(); ++l) tsr_cmp += *(ptr_begin + l);
                    dwj_ddwk_dwi_ijk[tsr_idx(i,j,k)] += tsr_cmp;
                  }
                } // for i,j
              }}} // end for mask, w_k_elem, mu
        } // end for k

        for(double& elem: dwj_ddwk_dwi_ijk){
          glb_sum(elem);
          elem /= geo.total_volume();
          IO::wb(ofs,elem); // i,j,k
        }
      }


      void set_w_type
      (
       std::vector<std::vector<FlowType>>& ws,
       int& max_flow_size,
       const std::vector<FlowType>& w_type
       ){
        for(const FlowType& flow_type : w_type) max_flow_size = std::max(max_flow_size, flow_type.flow_size);
        ws.push_back(w_type);
      }


      void write_dw_ddw_dw_wrapper
      (
       std::ofstream& ofs,
       const GaugeField& gf
       ){
        const Geometry geo = geo_reform(gf.geo(), 1);
        const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();

        int max_flow_size = 0;

        constexpr int nops = 7;
        // constexpr int nops = 1;
        std::vector<std::vector<FlowType>> ws;
        set_w_type(ws, max_flow_size, Obs::w_type(0));
        set_w_type(ws, max_flow_size, Obs::w_type(1,0));
        // set_w_type(ws, max_flow_size, Obs::w_type(1));
        // set_w_type(ws, max_flow_size, Obs::w_type(1,1));
        set_w_type(ws, max_flow_size, Obs::w_type(2,0));
        // set_w_type(ws, max_flow_size, Obs::w_type(2,1));
        set_w_type(ws, max_flow_size, Obs::w_type(3,0));
        // set_w_type(ws, max_flow_size, Obs::w_type(3,1));
        set_w_type(ws, max_flow_size, Obs::w_type(4,0));
        // set_w_type(ws, max_flow_size, Obs::w_type(4,1));
        set_w_type(ws, max_flow_size, Obs::w_type(5));
        // set_w_type(ws, max_flow_size, Obs::w_type(6));
        set_w_type(ws, max_flow_size, Obs::w_type(7));
        assert(static_cast<int>(ws.size())==nops);

        const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
        const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
        const Geometry geo_ext = geo_reform(geo, 1, expand_left, expand_right);
        GaugeField gf_ext;
        gf_ext.init(geo_ext);
        gf_ext = gf;
        refresh_expanded(gf_ext);

        dw_ddw_dw_no_comm(ofs, gf_ext, nops, ws, cmcs(), geo_ext);
      }



      // --------------------------------
      // --------- DEFINITIONS ----------
      // --------------------------------


      ColorMatrix staple_G0_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const std::function<const std::vector<int>(const int, const int)>
          plaq = [](const int mu, const int nu){
            return std::vector<int>{mu, nu, -mu-1, -nu-1};
          };
        return 2.0 * calculate_staple_two_dir_no_comm(gf_ext,plaq,xl,mu0);
      }


      ColorMatrix staple_G1r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const std::function<const std::vector<int>(const int, const int)>
          lrect = [](const int mu, const int nu){
            return std::vector<int>{mu,nu,nu,-mu-1,-nu-1,-nu-1};
          };
        const std::function<const std::vector<int>(const int, const int)>
          srect1 = [](const int mu, const int nu){
            return std::vector<int>{mu,mu,nu,-mu-1,-mu-1,-nu-1};
          };
        const std::function<const std::vector<int>(const int, const int)>
          srect2 = [](const int mu, const int nu){
            return std::vector<int>{mu,nu,-mu-1,-mu-1,-nu-1,mu};
          };
        ColorMatrix res;
        set_zero(res);
        res += calculate_staple_two_dir_no_comm(gf_ext,lrect,xl,mu0);
        res += calculate_staple_two_dir_no_comm(gf_ext,srect1,xl,mu0);
        res += calculate_staple_two_dir_no_comm(gf_ext,srect2,xl,mu0);
        return 2.0 * res;
      }


      ColorMatrix staple_G1c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const std::function<const std::vector<int>(const int, const int, const int)>
          chair4 = [](const int mu, const int nu, const int rho){
            return std::vector<int>{mu,nu,rho,-mu-1,-rho-1,-nu-1};
          };
        const std::function<const std::vector<int>(const int, const int, const int)>
          chair3 = [](const int mu, const int nu, const int rho){
            return std::vector<int>{mu,nu,-mu-1,rho,-nu-1,-rho-1};
          };
        const std::function<const std::vector<int>(const int, const int, const int)>
          chair5 = [](const int mu, const int nu, const int rho){
            return std::vector<int>{mu,nu,rho,-nu-1,-mu-1,-rho-1};
          };

        ColorMatrix res;
        set_zero(res);
        res += calculate_staple_three_dir_no_comm(gf_ext,chair4,xl,mu0);
        res += calculate_staple_three_dir_no_comm(gf_ext,chair3,xl,mu0);
        res += calculate_staple_three_dir_no_comm(gf_ext,chair5,xl,mu0);
        return 2.0 * res;
      }


      ColorMatrix staple_G1_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        return staple_G1r_no_comm(gf_ext,xl,mu0) + staple_G1c_no_comm(gf_ext,xl,mu0);
      }


      ColorMatrix staple_G2r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const std::function<const std::vector<int>(const int, const int)>
          twisted_lrect = [](const int mu, const int nu){
            return std::vector<int>{mu,nu,-mu-1,nu,mu,-nu-1,-mu-1,-nu-1};
          };
        const std::function<const std::vector<int>(const int, const int)>
          twisted_mrect = [](const int mu, const int nu){
            return std::vector<int>{mu,nu,-mu-1,-nu-1,mu,-nu-1,-mu-1,nu};
          };
        const std::function<const std::vector<int>(const int, const int)>
          twisted_srect1 = [](const int mu, const int nu){
            return std::vector<int>{mu,nu,mu,-nu-1,-mu-1,nu,-mu-1,-nu-1};
          };
        const std::function<const std::vector<int>(const int, const int)>
          twisted_srect2 = [](const int mu, const int nu){
            return std::vector<int>{mu,nu,-mu-1,-nu-1,-mu-1,nu,mu,-nu-1};
          };

        ColorMatrix res;
        set_zero(res);
        res += calculate_staple_two_dir_no_comm(gf_ext,twisted_lrect,xl,mu0);
        res += calculate_staple_two_dir_no_comm(gf_ext,twisted_mrect,xl,mu0);
        res += calculate_staple_two_dir_no_comm(gf_ext,twisted_srect1,xl,mu0);
        res += calculate_staple_two_dir_no_comm(gf_ext,twisted_srect2,xl,mu0);
        return 2.0 * res;
      }


      ColorMatrix staple_G2c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const std::function<const std::vector<int>(const int, const int, const int)>
          twisted_lchair = [](const int mu, const int nu, const int rho){
            return std::vector<int>{mu,nu,-mu-1,rho,mu,-rho-1,-mu-1,-nu-1};
          };
        const std::function<const std::vector<int>(const int, const int, const int)>
          twisted_mchair = [](const int mu, const int nu, const int rho){
            return std::vector<int>{mu,nu,-mu-1,-nu-1,mu,rho,-mu-1,-rho-1};
          };
        const std::function<const std::vector<int>(const int, const int, const int)>
          twisted_schair3 = [](const int mu, const int nu, const int rho){
            return std::vector<int>{mu,nu,-mu-1,-nu-1,rho,nu,-rho-1,-nu-1};
          };
        const std::function<const std::vector<int>(const int, const int, const int)>
          twisted_schair7 = [](const int mu, const int nu, const int rho){
            return std::vector<int>{mu,nu,rho,-nu-1,-rho-1,nu,-mu-1,-nu-1};
          };

        ColorMatrix res;
        set_zero(res);
        res += calculate_staple_three_dir_no_comm(gf_ext,twisted_lchair,xl,mu0);
        res += calculate_staple_three_dir_no_comm(gf_ext,twisted_mchair,xl,mu0);
        res += calculate_staple_three_dir_no_comm(gf_ext,twisted_schair3,xl,mu0);
        res += calculate_staple_three_dir_no_comm(gf_ext,twisted_schair7,xl,mu0);
        return 2.0 * res;
      }


      ColorMatrix staple_G2_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        return staple_G2r_no_comm(gf_ext,xl,mu0) + staple_G2c_no_comm(gf_ext,xl,mu0);
      }


      ColorMatrix staple_G3r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const VectFunc2 plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1};
        };
        const VectFunc2 plaq_below = [](const int mu, const int nu){
          return std::vector<int>{mu, -nu-1, -mu-1, nu};
        };

        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu0);
        const Coordinate xl_minus_mu = coordinate_shifts(xl, -mu0-1);
        for(int nu=0; nu<DIMN; ++nu){
          if(mu0==nu) continue;
          const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
          const Coordinate xl_minus_nu = coordinate_shifts(xl, -nu-1);
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,nu) );
            std::complex<double> factor = 0.0;
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_mu, plaq(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq(mu0,nu) ) );
            tmp *= factor;
            res += tmp;
          }
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq_below,mu0,nu) );
            std::complex<double> factor = 0.0;
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq_below(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq_below(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_mu, plaq_below(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq_below(mu0,nu) ) );
            tmp *= factor;
            res += tmp;
          }
        } // end for nu

        return 2.0 * res;
      }



      ColorMatrix staple_G3c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const VectFunc2 plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1};
        };
        const VectFunc2 plaq_below = [](const int mu, const int nu){
          return std::vector<int>{mu, -nu-1, -mu-1, nu};
        };

        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu0);
        for(int nu=0; nu<DIMN; ++nu){
          if(mu0==nu) continue;
          const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
          const Coordinate xl_minus_nu = coordinate_shifts(xl, -nu-1);
          //
          for(int rho=nu+1; rho<DIMN; ++rho){
            if(mu0==rho) continue;
            const Coordinate xl_plus_rho = coordinate_shifts(xl, rho);
            const Coordinate xl_minus_rho = coordinate_shifts(xl, -rho-1);
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,nu) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq(mu0,rho) ) ); // (1)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq(mu0,-rho-1) ) ); // (2)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(rho,nu) ) ); // (5)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-rho-1,nu) ) ); // (6)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,mu0) ) ); // (4)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,mu0) ) ); // (3)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,rho) ) ); // (8)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,-rho-1) ) ); // (7)
              tmp *= factor;
              res += tmp;
            }
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,-nu-1) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(rho,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-rho-1,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq(mu0,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq(mu0,-rho-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,-rho-1) ) );
              tmp *= factor;
              res += tmp;
            }
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,rho) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(nu,rho) ) ); // (10)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-nu-1,rho) ) ); // (9)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_rho, plaq(mu0,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_rho, plaq(mu0,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,nu) ) ); // (11)
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,-nu-1) ) ); // (12)
              tmp *= factor;
              res += tmp;
            }
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,-rho-1) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_rho, plaq(mu0,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_rho, plaq(mu0,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(nu,-rho-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-nu-1,-rho-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,-nu-1) ) );
              tmp *= factor;
              res += tmp;
            }
          } // end for rho
        } // end for nu

        return 2.0 * res;
      }


      ColorMatrix staple_G3_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        return staple_G3r_no_comm(gf_ext,xl,mu0) + staple_G3c_no_comm(gf_ext,xl,mu0);
      }


      ColorMatrix staple_G4r_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const VectFunc2 plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1};
        };
        const VectFunc2 plaq_below = [](const int mu, const int nu){
          return std::vector<int>{mu, -nu-1, -mu-1, nu};
        };

        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu0);
        const Coordinate xl_minus_mu = coordinate_shifts(xl, -mu0-1);
        for(int nu=0; nu<DIMN; ++nu){
          if(mu0==nu) continue;
          const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
          const Coordinate xl_minus_nu = coordinate_shifts(xl, -nu-1);
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,nu) );
            std::complex<double> factor = 0.0;
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_mu, plaq(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq(mu0,nu) ) );
            tmp *= conj(factor);
            res += tmp;
          }
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq_below,mu0,nu) );
            std::complex<double> factor = 0.0;
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq_below(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq_below(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_mu, plaq_below(mu0,nu) ) );
            factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq_below(mu0,nu) ) );
            tmp *= conj(factor);
            res += tmp;
          }
        } // end for nu
        return 2.0 * res;
      }



      ColorMatrix staple_G4c_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const VectFunc2 plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1};
        };
        const VectFunc2 plaq_below = [](const int mu, const int nu){
          return std::vector<int>{mu, -nu-1, -mu-1, nu};
        };

        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu0);
        for(int nu=0; nu<DIMN; ++nu){
          if(mu0==nu) continue;
          const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
          const Coordinate xl_minus_nu = coordinate_shifts(xl, -nu-1);
          //
          for(int rho=nu+1; rho<DIMN; ++rho){
            if(mu0==rho) continue;
            const Coordinate xl_plus_rho = coordinate_shifts(xl, rho);
            const Coordinate xl_minus_rho = coordinate_shifts(xl, -rho-1);
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,nu) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq(mu0,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_nu, plaq(mu0,-rho-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(rho,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-rho-1,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,-rho-1) ) );
              tmp *= conj(factor);
              res += tmp;
            }
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,-nu-1) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(rho,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-rho-1,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq(mu0,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_nu, plaq(mu0,-rho-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,-rho-1) ) );
              tmp *= conj(factor);
              res += tmp;
            }
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,rho) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(nu,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-nu-1,rho) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_rho, plaq(mu0,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_rho, plaq(mu0,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(rho,-nu-1) ) );
              tmp *= conj(factor);
              res += tmp;
            }
            {
              ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,-rho-1) );
              std::complex<double> factor = 0.0;
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_rho, plaq(mu0,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_minus_rho, plaq(mu0,-nu-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(nu,-rho-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_plus_mu, plaq(-nu-1,-rho-1) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(nu,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-nu-1,mu0) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,nu) ) );
              factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(-rho-1,-nu-1) ) );
              tmp *= conj(factor);
              res += tmp;
            }
          } // end for rho
        } // end for nu

        return 2.0 * res;
      }


      ColorMatrix staple_G4_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        return staple_G4r_no_comm(gf_ext,xl,mu0) + staple_G4c_no_comm(gf_ext,xl,mu0);
      }


      ColorMatrix staple_G5_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const VectFunc2 twice_plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1, mu, nu, -mu-1, -nu-1};
        };

        ColorMatrix res = calculate_staple_two_dir_no_comm(gf_ext,twice_plaq,xl,mu0);
        res *= 4.0;
        return res;
      }


      ColorMatrix staple_G6_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const VectFunc2 plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1};
        };
        const VectFunc2 plaq_below = [](const int mu, const int nu){
          return std::vector<int>{mu, -nu-1, -mu-1, nu};
        };

        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu0);
        for(int nu=0; nu<DIMN; ++nu){
          if(mu0==nu) continue;
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,nu) );
            const std::complex<double> factor = matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(mu0,nu) ) );
            tmp *= factor;
            res += tmp;
          }
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq_below,mu0,nu) );
            const std::complex<double> factor = matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq_below(mu0,nu) ) );
            tmp *= factor;
            res += tmp;
          }
        } // end for nu
        res *= 4.0;

        return res;
      }


      ColorMatrix staple_G7_no_comm(const GaugeField& gf_ext, const Coordinate& xl, const int mu0){
        const VectFunc2 plaq = [](const int mu, const int nu){
          return std::vector<int>{mu, nu, -mu-1, -nu-1};
        };
        const VectFunc2 plaq_below = [](const int mu, const int nu){
          return std::vector<int>{mu, -nu-1, -mu-1, nu};
        };

        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu0);
        for(int nu=0; nu<DIMN; ++nu){
          if(mu0==nu) continue;
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq,mu0,nu) );
            const std::complex<double> factor = matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq(mu0,nu) ) );
            tmp *= conj(factor);
            res += tmp;
          }
          {
            ColorMatrix tmp = gf_wilson_line_no_comm( gf_ext, xl_plus_mu, convert2staple(plaq_below,mu0,nu) );
            const std::complex<double> factor = matrix_trace( gf_wilson_line_no_comm( gf_ext, xl, plaq_below(mu0,nu) ) );
            tmp *= conj(factor);
            res += tmp;
          }
        } // end for nu
        res *= 2.0; // complex conjugate is itself

        return res;
      }




      double calculate_W_two_dir_no_comm
      (
       const GaugeField& gf,
       const std::function<const std::vector<int>(const int, const int)>& w,
       const bool is_mu_nu_dir_identical /*= false; when true, we count c.c. in the for loops */
       ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        qacc_for(index, geo.local_volume(), {
            fd.get_elem(index) = 0;
            const Coordinate xl = geo.coordinate_from_index(index);
            for(int mu=0; mu<DIMN; ++mu){
              for(int nu=0; nu<DIMN; ++nu){
                if(mu==nu) continue;

                const std::complex<double> tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(mu, nu) ));
                fd.get_elem(index) += tmp.real();
              }
            }
            if(!is_mu_nu_dir_identical) fd.get_elem(index) *= 2.0; // complex conjugate; note that we already took .real()
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//

      double calculate_W_three_dir_no_comm
      (
       const GaugeField& gf,
       const std::function<const std::vector<int>(const int, const int, const int)>& w
       ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        qacc_for(index, geo.local_volume(), {
            fd.get_elem(index) = 0;
            const Coordinate xl = geo.coordinate_from_index(index);
            for(int mu=0; mu<DIMN; ++mu){
              for(int nu=mu+1; nu<DIMN; ++nu){
                for(int rho=nu+1; rho<DIMN; ++rho){

                  std::complex<double> tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(mu, nu, rho) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(mu, nu, -rho-1) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(mu, rho, nu) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(mu, -rho-1, nu) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(nu, mu, rho) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(nu, mu, -rho-1) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(nu, rho, mu) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(nu, -rho-1, mu) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(rho, mu, nu) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(rho, mu, -nu-1) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(rho, nu, mu) ));
                  fd.get_elem(index) += tmp.real();
                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, w(rho, -nu-1, mu) ));
                  fd.get_elem(index) += tmp.real();
                }
              }
            }
            fd.get_elem(index) *= 2.0; // complex conjugate; note that we already took .real()
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//

      ColorMatrix calculate_staple_two_dir_no_comm
      (
       const GaugeField& gf,
       const std::function<const std::vector<int>(const int, const int)>& loop,
       const Coordinate& xl,
       const int mu
       ){
        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu);
        for(int nu=0; nu<DIMN; ++nu){
          if(mu==nu) continue;
          res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop,mu,nu) );
          res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop,mu,-nu-1) );
        }
        return res;
      }//


      ColorMatrix calculate_staple_three_dir_no_comm
      (
       const GaugeField& gf,
       const std::function<const std::vector<int>(const int, const int, const int)>& loop,
       const Coordinate& xl,
       const int mu
       ){
        ColorMatrix res;
        set_zero(res);

        const Coordinate xl_plus_mu = coordinate_shifts(xl, mu);
        for(int nu=0; nu<DIMN; ++nu){
          for(int rho=nu+1; rho<DIMN; ++rho){
            if(mu==nu) continue;
            if(mu==rho) continue;
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, nu, rho) );
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, nu, -rho-1) );
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, rho, nu) );
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, rho, -nu-1) );
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, -nu-1, rho) );
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, -nu-1, -rho-1) );
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, -rho-1, nu) );
            res += gf_wilson_line_no_comm( gf, xl_plus_mu, convert2staple(loop, mu, -rho-1, -nu-1) );
          }
        }
        return res;
      }//


      double calculate_W3r_no_comm( const GaugeField& gf ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        const std::function<const std::vector<int>(const int, const int)> plaq
          = [](const int mu, const int nu){
            return std::vector<int>{mu, nu, -mu-1, -nu-1};
          };

        qacc_for(index, geo.local_volume(), {
            const Coordinate xl = geo.coordinate_from_index(index);
            fd.get_elem(index) = 0;
            for(int nu=0; nu<DIMN; ++nu){
              const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
              for(int mu=0; mu<DIMN; ++mu){
                if(mu==nu) continue;

                const std::complex<double> tmp
                  = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ))
                  * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(mu, nu) ));
                fd.get_elem(index) += tmp.real()  ;
              }
            }
            fd.get_elem(index) *= 2.0; // complex conjugate; note that we already took .real()
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//


      double calculate_W3c_no_comm( const GaugeField& gf ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        const std::function<const std::vector<int>(const int, const int)> plaq
          = [](const int mu, const int nu){
            return std::vector<int>{mu, nu, -mu-1, -nu-1};
          };

        qacc_for(index, geo.local_volume(), {
            const Coordinate xl = geo.coordinate_from_index(index);
            fd.get_elem(index) = 0;

            for(int mu=0; mu<DIMN; ++mu){
              const Coordinate xl_plus_mu = coordinate_shifts(xl, mu);
              for(int nu=mu+1; nu<DIMN; ++nu){
                const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
                const Coordinate xl_minus_nu = coordinate_shifts(xl, -nu-1);
                for(int rho=nu+1; rho<DIMN; ++rho){
                  const Coordinate xl_plus_rho = coordinate_shifts(xl, rho);
                  const Coordinate xl_minus_rho = coordinate_shifts(xl, -rho-1);

                  std::complex<double> tmp
                    = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(mu, rho) )); // (1)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(mu, -rho-1) )); // (2)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, rho) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_rho, plaq(mu, nu) )); // (3)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(mu, -rho-1) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_minus_rho, plaq(mu, nu) )); // (4)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(nu, mu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(nu, rho) )); // (5)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(nu, mu) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(nu, -rho-1) )); // (6)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(nu, rho) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_rho, plaq(nu, mu) )); // (7)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(nu, -rho-1) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_minus_rho, plaq(nu, mu) )); // (8)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(rho, mu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(rho, nu) )); // (9)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(rho, mu) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(rho, -nu-1) )); // (10)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(rho, nu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(rho, mu) )); // (11)
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(rho, -nu-1) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_minus_nu, plaq(rho, mu) )); // (12)
                  fd.get_elem(index) += tmp.real();

                }
              }
            }
            fd.get_elem(index) *= 2.0; // complex conjugate; note that we already took .real()
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//


      double calculate_W4r_no_comm( const GaugeField& gf ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        const std::function<const std::vector<int>(const int, const int)> plaq
          = [](const int mu, const int nu){
            return std::vector<int>{mu, nu, -mu-1, -nu-1};
          };

        qacc_for(index, geo.local_volume(), {
            fd.get_elem(index) = 0;
            const Coordinate xl = geo.coordinate_from_index(index);
            // rect shape
            for(int nu=0; nu<DIMN; ++nu){
              const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
              for(int mu=0; mu<DIMN; ++mu){
                if(mu==nu) continue;

                const std::complex<double> tmp
                  = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ))
                  * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(nu, mu) ));
                fd.get_elem(index) += tmp.real();
              }
            }
            fd.get_elem(index) *= 2.0; // complex conjugate; note that we already took .real()
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//

      double calculate_W4c_no_comm( const GaugeField& gf ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        const std::function<const std::vector<int>(const int, const int)> plaq
          = [](const int mu, const int nu){
            return std::vector<int>{mu, nu, -mu-1, -nu-1};
          };

        qacc_for(index, geo.local_volume(), {
            fd.get_elem(index) = 0;
            const Coordinate xl = geo.coordinate_from_index(index);
            // chair shape
            for(int mu=0; mu<DIMN; ++mu){
              const Coordinate xl_plus_mu = coordinate_shifts(xl, mu);
              for(int nu=mu+1; nu<DIMN; ++nu){
                const Coordinate xl_plus_nu = coordinate_shifts(xl, nu);
                const Coordinate xl_minus_nu = coordinate_shifts(xl, -nu-1);
                for(int rho=nu+1; rho<DIMN; ++rho){
                  const Coordinate xl_plus_rho = coordinate_shifts(xl, rho);
                  const Coordinate xl_minus_rho = coordinate_shifts(xl, -rho-1);

                  std::complex<double> tmp
                    = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(rho, mu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(-rho-1, mu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, rho) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_rho, plaq(nu, mu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(mu, -rho-1) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_minus_rho, plaq(nu, mu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(nu, mu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(rho, nu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(nu, mu) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(-rho-1, nu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(nu, rho) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_rho, plaq(mu, nu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(nu, -rho-1) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_minus_rho, plaq(mu, nu) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(rho, mu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(nu, rho) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(rho, mu) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_plus_mu, plaq(-nu-1, rho) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(rho, nu) ))
                    * matrix_trace( gf_wilson_line_no_comm( gf, xl_plus_nu, plaq(mu, rho) ));
                  fd.get_elem(index) += tmp.real();

                  tmp = matrix_trace(gf_wilson_line_no_comm( gf, xl, plaq(rho, -nu-1) ))
                    * matrix_trace(gf_wilson_line_no_comm( gf, xl_minus_nu, plaq(mu, rho) ));
                  fd.get_elem(index) += tmp.real();
                }
              }
            }
            fd.get_elem(index) *= 2.0; // complex conjugate; note that we already took .real()
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//

      double calculate_W6_no_comm( const GaugeField& gf ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        const std::function<const std::vector<int>(const int, const int)> plaq
          = [](const int mu, const int nu){
            return std::vector<int>{mu, nu, -mu-1, -nu-1};
          };

        qacc_for(index, geo.local_volume(), {
            fd.get_elem(index) = 0;
            const Coordinate xl = geo.coordinate_from_index(index);
            for(int mu=0; mu<DIMN; ++mu){
              for(int nu=0; nu<DIMN; ++nu){
                if(mu==nu) continue;

                const std::complex<double> tmp
                  = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ));
                fd.get_elem(index) += (tmp*tmp).real();
              }
            }
            // no x2 multiplication
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//

      double calculate_W7_no_comm( const GaugeField& gf ){
        const Geometry geo = geo_reform(gf.geo());
        FieldM<double, 1> fd;
        fd.init(geo);

        const std::function<const std::vector<int>(const int, const int)> plaq
          = [](const int mu, const int nu){
            return std::vector<int>{mu, nu, -mu-1, -nu-1};
          };

        qacc_for(index, geo.local_volume(), {
            fd.get_elem(index) = 0;
            const Coordinate xl = geo.coordinate_from_index(index);
            for(int mu=0; mu<DIMN; ++mu){
              for(int nu=0; nu<DIMN; ++nu){
                if(mu==nu) continue;

                const std::complex<double> tmp
                  = matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(mu, nu) ))
                  * matrix_trace( gf_wilson_line_no_comm( gf, xl, plaq(nu, mu) ));
                fd.get_elem(index) += tmp.real();
              }
            }
            fd.get_elem(index) /= 2.0; // complex conjugate is itself
          });

        double sum = 0;
        for (long index = 0; index < geo.local_volume(); ++index) sum += fd.get_elem(index);
        glb_sum(sum);
        sum /= geo.total_volume();

        return sum;
      }//

      // -------------------------------
      double topological_charge_link_plus(const GaugeField& gf_ext, const Coordinate& xl, const int mu){
        double res = 0.0;
        int counter = 0;
        for(int nu=-DIMN; nu<DIMN; ++nu){
          for(int rho=-DIMN; rho<DIMN; ++rho){
            for(int sigma=-DIMN; sigma<DIMN; ++sigma){

              if(FlowTypeQPlus2::epsilon(mu,nu,rho,sigma)==1){
                // type 1
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,nu,-mu-1,-nu-1,rho,sigma,-rho-1,-sigma-1}) ));
                ++counter;
                // type 2
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,nu,-mu-1,rho,sigma,-rho-1,-sigma-1,-nu-1}) ));
                ++counter;
                // type 3
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,nu,rho,sigma,-rho-1,-sigma-1,-mu-1,-nu-1}) ));
                ++counter;
                // type 4
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,rho,sigma,-rho-1,-sigma-1,nu,-mu-1,-nu-1}) ));
                ++counter;
              }
            } // end for sigma
          } // end for rho
        } // end for nu
        assert(counter==96);

        return res;
      }


      double topological_charge_link_minus(const GaugeField& gf_ext, const Coordinate& xl, const int mu){
        double res = 0.0;
        int counter = 0;
        for(int nu=-DIMN; nu<DIMN; ++nu){
          for(int rho=-DIMN; rho<DIMN; ++rho){
            for(int sigma=-DIMN; sigma<DIMN; ++sigma){

              if(FlowTypeQPlus2::epsilon(mu,nu,rho,sigma)==-1){
                // type 1
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,nu,-mu-1,-nu-1,rho,sigma,-rho-1,-sigma-1}) ));
                ++counter;
                // type 2
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,nu,-mu-1,rho,sigma,-rho-1,-sigma-1,-nu-1}) ));
                ++counter;
                // type 3
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,nu,rho,sigma,-rho-1,-sigma-1,-mu-1,-nu-1}) ));
                ++counter;
                // type 4
                res += real( matrix_trace( gf_wilson_line_no_comm( gf_ext, xl,
                                                                   std::vector<int>{mu,rho,sigma,-rho-1,-sigma-1,nu,-mu-1,-nu-1}) ));
                ++counter;
              }
            } // end for sigma
          } // end for rho
        } // end for nu
        assert(counter==96);

        return res;
      }


      double topological_charge(const GaugeField& gf_ext){
        double res = 0.0;
        const Geometry& geo = gf_ext.geo();

        qfor(index, geo.local_volume(), {
            const Coordinate xl = geo.coordinate_from_index(index);
            for(int mu=-DIMN; mu<DIMN; ++mu) {
              res += topological_charge_link_plus(gf_ext, xl, mu);
              res -= topological_charge_link_minus(gf_ext, xl, mu);
            }
          });
        return res;
      }


      // --------------------------------
      inline std::vector<double> get_gauge_field_infos_2(const GaugeField& gf)
      // Values:
      // 0: avg(plaq_action)
      // 1: avg(plaq_action^2)
      // 2: tot(topo_density, 5LI)
      // 3: avg(topo_density^2, 5LI)
      // 4: tot(topo_density, CL)
      // 5: avg(topo_density^2, CL)
      {
        TIMER("get_gauge_field_infos_2");
        const int info_vec_size = 6;
        const Geometry geo = geo_reform(gf.geo());
        CloverLeafField clf;
        gf_clover_leaf_field(clf, gf);
        CloverLeafField clf1, clf2, clf3, clf4, clf5;
        gf_clover_leaf_field_5(clf1, clf2, clf3, clf4, clf5, gf);
        FieldM<double, 1> paf;
        clf_plaq_action_field(paf, clf1);
        FieldM<double, 1> topf;
        clf_topology_field_5(topf, clf1, clf2, clf3, clf4, clf5);
        FieldM<double, 1> topf2;
        clf_topology_field(topf2, clf);
        std::vector<double> info_vec(info_vec_size, 0.0);
        for (long index = 0; index < geo.local_volume(); ++index) {
          const double pa = paf.get_elem(index);
          const double top = topf.get_elem(index);
          const double top2 = topf2.get_elem(index);
          info_vec[0] += pa;
          info_vec[1] += sqr(pa);
          info_vec[2] += top;
          info_vec[3] += sqr(top);
          info_vec[4] += top2;
          info_vec[5] += sqr(top2);
        }
        glb_sum(get_data(info_vec));
        info_vec[0] *= 1.0 / (double)geo.total_volume();
        info_vec[1] *= 1.0 / (double)geo.total_volume();
        info_vec[3] *= 1.0 / (double)geo.total_volume();
        info_vec[5] *= 1.0 / (double)geo.total_volume();
        return info_vec;
      }


      // inline std::vector<double> get_gauge_field_infos_2(const GaugeField& gf)
      // // Values:
      // // 0: avg(plaq_action)
      // // 1: avg(plaq_action^2)
      // // 2: tot(topo_density)
      // // 3: tot(topo_density)
      // // 4: avg(topo_density^2)
      // {
      //   TIMER("get_gauge_field_infos_2");
      //   const int info_vec_size = 4;
      //   const Geometry geo = geo_reform(gf.geo());
      //   CloverLeafField clf;
      //   gf_clover_leaf_field(clf, gf);
      //   FieldM<double, 1> paf;
      //   clf_plaq_action_field(paf, clf);
      //   FieldM<double, 1> topf;
      //   clf_topology_field(topf, clf);
      //   std::vector<double> info_vec(info_vec_size, 0.0);
      //   for (long index = 0; index < geo.local_volume(); ++index) {
      //     const double pa = paf.get_elem(index);
      //     const double top = topf.get_elem(index);
      //     info_vec[0] += pa;
      //     info_vec[1] += sqr(pa);
      //     info_vec[2] += top;
      //     info_vec[3] += sqr(top);
      //   }
      //   glb_sum(get_data(info_vec));
      //   info_vec[0] *= 1.0 / (double)geo.total_volume();
      //   info_vec[1] *= 1.0 / (double)geo.total_volume();
      //   info_vec[3] *= 1.0 / (double)geo.total_volume();
      //   return info_vec;
      // }

      inline std::string show_gauge_field_info_line(const int i,
                                                    const std::vector<double>& v)
      {
        qassert(v.size() == 4);
        return ssprintf("%5d %24.17E %24.17E %24.17E %24.17E %24.17E %24.17E", 
                        i, v[0], v[1], v[2], v[3], v[4], v[5]);
      }

      // --------------------------------


      void observe_topological_charge
      (
       const GaugeField& gf,
       std::ofstream& ofs_ene_den,
       std::ofstream& ofs_topo5LI,
       std::ofstream& ofs_topoCL,
       const double flow_time_w_total = 2.0,
       const double epsilon_w = 0.01,
       const int interval_w = 20
       ){
        GaugeField gf1;
        gf1 = gf;
        const int FLOW_STEPS_W_TOTAL = std::floor(flow_time_w_total/epsilon_w); //500
        //
        const double FLOW_TIME_W_EACH = epsilon_w * interval_w;//0.1
        //
        double flow_time_W_tmp = 0.0;
        std::vector<double> ene_den_list;
        std::vector<double> topo5LI_list;
        std::vector<double> topoCL_list;

        {
          ene_den_list.push_back( gf_energy_density(gf1) ); // init
          const std::vector<double> info = get_gauge_field_infos_2(gf1);
          topo5LI_list.push_back( info[2] ); // init
          topoCL_list.push_back( info[4] ); // init
        }

        for (int step_W = 0; step_W < FLOW_STEPS_W_TOTAL; step_W+=interval_w) {
          const std::vector<double> ene_den_list_small = gf_wilson_flow(gf1,
                                                                        flow_time_W_tmp,
                                                                        FLOW_TIME_W_EACH,
                                                                        interval_w);
          flow_time_W_tmp += FLOW_TIME_W_EACH;
          const std::vector<double> info = get_gauge_field_infos_2(gf1);
          ene_den_list.insert(ene_den_list.end(),
                              ene_den_list_small.begin(),
                              ene_den_list_small.end());
          topo5LI_list.push_back(info[2]);
          topoCL_list.push_back(info[4]);
        }
        bool is_first = true;

        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          for(const double& elem : ene_den_list) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_ene_den, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo5LI_list) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo5LI, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topoCL_list) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topoCL, ss.str());
        }
      }

      // void observe_topological_charge
      // (
      //  const GaugeField& gf,
      //  std::ofstream& ofs_ene_den,
      //  std::ofstream& ofs_topo_5li,
      //  std::ofstream& ofs_topo_5li_sq,
      //  std::ofstream& ofs_topo_cl,
      //  std::ofstream& ofs_topo_cl_sq,
      //  const double flow_time_w_total = 2.0,
      //  const double epsilon_w = 0.01,
      //  const int interval_w = 20
      //  ){
      //   GaugeField gf1;
      //   gf1 = gf;
      //   const int FLOW_STEPS_W_TOTAL = std::floor(flow_time_w_total/epsilon_w); //500
      //   //
      //   const double FLOW_TIME_W_EACH = epsilon_w * interval_w;//0.1
      //   //
      //   double flow_time_W_tmp = 0.0;
      //   std::vector<double> ene_den_list;
      //   std::vector<double> topo_list;
      //   std::vector<double> topo_list2;

      //   {
      //     ene_den_list.push_back( gf_energy_density(gf1) ); // init
      //     const std::vector<double> info = get_gauge_field_infos_2(gf1);
      //     topo_list.push_back( info[2] ); // init
      //     topo_list2.push_back( info[4] ); // init
      //   }

      //   for (int step_W = 0; step_W < FLOW_STEPS_W_TOTAL; step_W+=interval_w) {
      //     const std::vector<double> ene_den_list_small = gf_wilson_flow(gf1,
      //                                                                   flow_time_W_tmp,
      //                                                                   FLOW_TIME_W_EACH,
      //                                                                   interval_w);
      //     flow_time_W_tmp += FLOW_TIME_W_EACH;
      //     const std::vector<double> info = get_gauge_field_infos_2(gf1);
      //     /* CONTENTS of info:
      //        avg(plaq_action), avg(plaq_action^2), tot(topo_density), avg(topo_density^2) */
      //     ene_den_list.insert(ene_den_list.end(),
      //                         ene_den_list_small.begin(),
      //                         ene_den_list_small.end());
      //     topo_list.push_back(info[2]);
      //     topo_list2.push_back(info[4]);
      //   }
      //   bool is_first = true;

      //   {
      //     std::stringstream ss;
      //     ss << std::scientific << std::setprecision(15);

      //     for(const double& elem : ene_den_list) {
      //       if(is_first) is_first = false;
      //       else ss << ',';
      //       ss << elem;
      //     }
      //     ss << std::endl;
      //     IO::print_to(ofs_ene_den, ss.str());
      //   }
      //   {
      //     std::stringstream ss;
      //     ss << std::scientific << std::setprecision(15);

      //     is_first = true;
      //     for(const double& elem : topo_list) {
      //       if(is_first) is_first = false;
      //       else ss << ',';
      //       ss << elem;
      //     }
      //     ss << std::endl;
      //     IO::print_to(ofs_topo, ss.str());
      //   }
      //   {
      //     std::stringstream ss;
      //     ss << std::scientific << std::setprecision(15);

      //     is_first = true;
      //     for(const double& elem : topo_list2) {
      //       if(is_first) is_first = false;
      //       else ss << ',';
      //       ss << elem;
      //     }
      //     ss << std::endl;
      //     IO::print_to(ofs_topo2, ss.str());
      //   }
      // }

      void parity_trsf
      (
       GaugeField& Pgf, // initialize outside
       const GaugeField& gf
       ){
        const Geometry geo = geo_reform(gf.geo(), 1);
        const Coordinate total_site = geo.total_site();
        const Geometry geo_ext = geo_reform(geo, 1, total_site, total_site);

        GaugeField gf_ext;
        gf_ext.init(geo_ext);
        gf_ext = gf;
        refresh_expanded(gf_ext);

        qacc_for(index, geo.local_volume(), {
            const Coordinate Pxl = geo.coordinate_from_index(index);
            const Coordinate Pxg = geo.coordinate_g_from_l(Pxl);
            const Coordinate xg(total_site[0]-Pxg[0],Pxg[1],Pxg[2],Pxg[3]);
            const Coordinate xl = geo.coordinate_l_from_g(xg);
            {
              const int Pmu=0;
              ColorMatrix& PU_Pxl_Pmu = Pgf.get_elem(Pxl,Pmu);
              const Coordinate xl_shifted = coordinate_shifts(xl, Pmu-1);
              PU_Pxl_Pmu = matrix_adjoint(gf_ext.get_elem(xl_shifted, Pmu));
            }
            for(int Pmu=1; Pmu<DIMN; ++Pmu){
              ColorMatrix& PU_Pxl_Pmu = Pgf.get_elem(Pxl,Pmu);
              PU_Pxl_Pmu = gf_ext.get_elem(xl,Pmu);
            }
          });
      }

      void calculate_sq(double& S_Q_parity,
                        double& action_part_parity,
                        double& det_part_parity,
                        double& S_Q_time_rev,
                        double& action_part_time_rev,
                        double& det_part_time_rev,
                        const GaugeAction& ga,
                        const GaugeField& gfV,
                        const FlowInfo2& fi_theta ){
        GaugeField PgfV;
        PgfV.init(gfV.geo());
        PgfV = gfV;
        Generalized::Obs::parity_trsf(PgfV, gfV);
        double action_part_V;
        double det_part_V;
        double action_part_PV;
        double det_part_PV;
        Generalized::gf_hamilton_flowed(action_part_V, det_part_V,
                                        gfV, ga, fi_theta); // further flow from gfV to obtain SQ
        Generalized::gf_hamilton_flowed(action_part_PV, det_part_PV,
                                        PgfV, ga, fi_theta);
        FlowInfo2 fi_theta_flipped = fi_theta;
        fi_theta_flipped.flip_sign_of_t();
        double action_part_V_inv;
        double det_part_V_inv;
        Generalized::gf_hamilton_flowed(action_part_V_inv, det_part_V_inv,
                                        gfV, ga, fi_theta_flipped);

        action_part_parity = 0.5*(action_part_V-action_part_PV);
        det_part_parity = -0.5*(det_part_V-det_part_PV);
        S_Q_parity = action_part_parity+det_part_parity;
        action_part_time_rev = 0.5*(action_part_V-action_part_V_inv);
        det_part_time_rev = -0.5*(det_part_V-det_part_V_inv);
        S_Q_time_rev = action_part_time_rev+det_part_time_rev;
      }

      void observe_sq_before_and_after_smear
      (
       const GaugeAction& ga,
       const GaugeField& gfV,
       const FlowInfo2& fi_theta,
       std::ofstream& ofs_ene_den,
       std::ofstream& ofs_topo,
       std::ofstream& ofs_topo2,
       std::ofstream& ofs_topo_sq_parity,
       std::ofstream& ofs_topo_action_part_parity,
       std::ofstream& ofs_topo_det_part_parity,
       std::ofstream& ofs_topo_sq_time_rev,
       std::ofstream& ofs_topo_action_part_time_rev,
       std::ofstream& ofs_topo_det_part_time_rev,
       const double flow_time_w_total = 2.0,
       const double epsilon_w = 0.01,
       const int interval_w = 20
       ){
        GaugeField gfV_flowed;
        gfV_flowed = gfV;
        const int FLOW_STEPS_W_TOTAL = std::floor(flow_time_w_total/epsilon_w); //500
        //
        const double FLOW_TIME_W_EACH = epsilon_w * interval_w;//0.1
        //
        double flow_time_W_tmp = 0.0;
        std::vector<double> ene_den_list;
        std::vector<double> topo_list;
        std::vector<double> topo_list2;
        std::vector<double> topo_list_sq_parity;
        std::vector<double> topo_list_action_part_parity;
        std::vector<double> topo_list_det_part_parity;
        std::vector<double> topo_list_sq_time_rev;
        std::vector<double> topo_list_action_part_time_rev;
        std::vector<double> topo_list_det_part_time_rev;

        {
          ene_den_list.push_back( gf_energy_density(gfV_flowed) ); // init
          const std::vector<double> info = get_gauge_field_infos_2(gfV_flowed);
          topo_list.push_back( info[2] ); // init
          topo_list2.push_back( info[4] ); // init

          double sq_parity;
          double action_part_parity;
          double det_part_parity;
          double sq_time_rev;
          double action_part_time_rev;
          double det_part_time_rev;
          calculate_sq( sq_parity, action_part_parity, det_part_parity,
                        sq_time_rev, action_part_time_rev, det_part_time_rev,
                        ga, gfV_flowed, fi_theta);
          topo_list_sq_parity.push_back( sq_parity );
          topo_list_action_part_parity.push_back( action_part_parity );
          topo_list_det_part_parity.push_back( det_part_parity );
          topo_list_sq_time_rev.push_back( sq_time_rev );
          topo_list_action_part_time_rev.push_back( action_part_time_rev );
          topo_list_det_part_time_rev.push_back( det_part_time_rev );
        }

        for (int step_W = 0; step_W < FLOW_STEPS_W_TOTAL; step_W+=interval_w) {
          const std::vector<double> ene_den_list_small = gf_wilson_flow(gfV_flowed,
                                                                        flow_time_W_tmp,
                                                                        FLOW_TIME_W_EACH,
                                                                        interval_w);
          flow_time_W_tmp += FLOW_TIME_W_EACH;

          const std::vector<double> info = get_gauge_field_infos_2(gfV_flowed);
          /* CONTENTS of info:
             avg(plaq_action), avg(plaq_action^2), tot(topo_density), avg(topo_density^2) */
          ene_den_list.insert(ene_den_list.end(),
                              ene_den_list_small.begin(),
                              ene_den_list_small.end());
          topo_list.push_back(info[2]);
          topo_list2.push_back(info[4]);

          double sq_parity;
          double action_part_parity;
          double det_part_parity;
          double sq_time_rev;
          double action_part_time_rev;
          double det_part_time_rev;
          calculate_sq( sq_parity, action_part_parity, det_part_parity,
                        sq_time_rev, action_part_time_rev, det_part_time_rev,
                        ga, gfV_flowed, fi_theta);
          topo_list_sq_parity.push_back( sq_parity );
          topo_list_action_part_parity.push_back( action_part_parity );
          topo_list_det_part_parity.push_back( det_part_parity );
          topo_list_sq_time_rev.push_back( sq_time_rev );
          topo_list_action_part_time_rev.push_back( action_part_time_rev );
          topo_list_det_part_time_rev.push_back( det_part_time_rev );
        }
        bool is_first = true;

        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : ene_den_list) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_ene_den, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list2) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo2, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list_sq_parity) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo_sq_parity, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list_action_part_parity) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo_action_part_parity, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list_det_part_parity) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo_det_part_parity, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list_sq_time_rev) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo_sq_time_rev, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list_action_part_time_rev) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo_action_part_time_rev, ss.str());
        }
        {
          std::stringstream ss;
          ss << std::scientific << std::setprecision(15);

          is_first = true;
          for(const double& elem : topo_list_det_part_time_rev) {
            if(is_first) is_first = false;
            else ss << ',';
            ss << elem;
          }
          ss << std::endl;
          IO::print_to(ofs_topo_det_part_time_rev, ss.str());
        }
      }

//       void observe_topological_charge
//       (
//        const GaugeField& gf,
//        std::ofstream& ofs_ene_den,
//        std::ofstream& ofs_topo,
//        std::ofstream& ofs_topo2,
//        const double flow_time_w_total = 2.0,
//        const double epsilon_w = 0.01,
//        const int interval_w = 20
//        ){
//         GaugeField gf1;
//         gf1 = gf;
//         const int FLOW_STEPS_W_TOTAL = std::floor(flow_time_w_total/epsilon_w); //500
//         //
//         const double FLOW_TIME_W_EACH = epsilon_w * interval_w;//0.1
//         //
//         double flow_time_W_tmp = 0.0;
//         std::vector<double> ene_den_list;
//         std::vector<double> topo_list;
//         std::vector<double> topo_list;

//         {
//           ene_den_list.push_back( gf_energy_density(gf1) ); // init
//           const std::vector<double> info = get_gauge_field_infos(gf1);
//           topo_list.push_back( info[2] ); // init
//         }

//         for (int step_W = 0; step_W < FLOW_STEPS_W_TOTAL; step_W+=interval_w) {
//           const std::vector<double> ene_den_list_small = gf_wilson_flow(gf1,
//                                                                         flow_time_W_tmp,
//                                                                         FLOW_TIME_W_EACH,
//                                                                         interval_w);
//           flow_time_W_tmp += FLOW_TIME_W_EACH;
//           const std::vector<double> info = get_gauge_field_infos(gf1);
//           /* CONTENTS of info:
//              avg(plaq_action), avg(plaq_action^2), tot(topo_density), avg(topo_density^2) */
//           ene_den_list.insert(ene_den_list.end(),
//                               ene_den_list_small.begin(),
//                               ene_den_list_small.end());
//           topo_list.push_back(info[2]);
// ##
//         }
//         bool is_first = true;

//         {
//           std::stringstream ss;
//           ss << std::scientific << std::setprecision(15);

//           for(const double& elem : ene_den_list) {
//             if(is_first) is_first = false;
//             else ss << ',';
//             ss << elem;
//           }
//           ss << std::endl;
//           IO::print_to(ofs_ene_den, ss.str());
//         }
//         {
//           std::stringstream ss;
//           ss << std::scientific << std::setprecision(15);

//           is_first = true;
//           for(const double& elem : topo_list) {
//             if(is_first) is_first = false;
//             else ss << ',';
//             ss << elem;
//           }
//           ss << std::endl;
//           IO::print_to(ofs_topo, ss.str());
//         }
//       }

      inline void gf_flow_with_obs(const GaugeField& gf0, const FlowInfo2& fi,
                                   std::ofstream& ofs_ene_den,
                                   std::ofstream& ofs_topo,
                                   std::ofstream& ofs_topo2,
                                   const double flow_time_w_total = 2.0,
                                   const double epsilon_w = 0.01,
                                   const int interval_w = 20
                                   )
      // gf0 is the gauge field which we perform HMC evolution
      {
        TIMER("gf_flow_with_obs");
        const int max_flow_size = fi.get_max_flow_size();
        qassert(max_flow_size!=0);
        const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
        const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
        const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
        GaugeField gf_ext;
        gf_ext.init(geo_ext);
        gf_ext = gf0;

        Obs::observe_topological_charge(gf_ext, ofs_ene_den, ofs_topo, ofs_topo2,
                                        flow_time_w_total,
                                        epsilon_w,
                                        interval_w);
        for (int i = 0; i < fi.size(); ++i) {
          const FlowStepInfo2& fsi = fi[i];
          const int mask = fsi.mask;
          const int mu = fsi.mu;
          const double epsilon = fsi.epsilon;
          const FlowType& flow_type = *(fsi.ptr_flow_type);
          refresh_expanded_gf_flow_kernel_mask_mu(gf_ext, mask, mu, flow_type);
          gf_flow_mask_mu_no_comm(gf_ext, gf_ext, mask, mu, epsilon, flow_type);
          if( fsi.is_measurement() ) Obs::observe_topological_charge(gf_ext,
                                                                     ofs_ene_den, ofs_topo, ofs_topo2,
                                                                     flow_time_w_total,
                                                                     epsilon_w,
                                                                     interval_w);
        }
      }

      // inline void gf_flow_with_obs(const GaugeField& gf0, const FlowInfo2& fi,
      //                              std::ofstream& ofs_ene_den,
      //                              std::ofstream& ofs_topo_5li,
      //                              std::ofstream& ofs_topo_5li_sq,
      //                              std::ofstream& ofs_topo_cl,
      //                              std::ofstream& ofs_topo_cl_sq,
      //                              const double flow_time_w_total = 2.0,
      //                              const double epsilon_w = 0.01,
      //                              const int interval_w = 20
      //                              )
      // // gf0 is the gauge field which we perform HMC evolution
      // {
      //   TIMER("gf_flow_with_obs");
      //   const int max_flow_size = fi.get_max_flow_size();
      //   qassert(max_flow_size!=0);
      //   const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      //   const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      //   const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
      //   GaugeField gf_ext;
      //   gf_ext.init(geo_ext);
      //   gf_ext = gf0;

      //   Obs::observe_topological_charge(gf_ext, ofs_ene_den, ofs_topo, ofs_topo2,
      //                                   flow_time_w_total,
      //                                   epsilon_w,
      //                                   interval_w);
      //   for (int i = 0; i < fi.size(); ++i) {
      //     const FlowStepInfo2& fsi = fi[i];
      //     const int mask = fsi.mask;
      //     const int mu = fsi.mu;
      //     const double epsilon = fsi.epsilon;
      //     const FlowType& flow_type = *(fsi.ptr_flow_type);
      //     refresh_expanded_gf_flow_kernel_mask_mu(gf_ext, mask, mu, flow_type);
      //     gf_flow_mask_mu_no_comm(gf_ext, gf_ext, mask, mu, epsilon, flow_type);
      //     if( fsi.is_measurement() ) Obs::observe_topological_charge(gf_ext,
      //                                                                ofs_ene_den,
      //                                                                ofs_topo_5li,
      //                                                                ofs_topo_5li_sq,
      //                                                                ofs_topo_cl,
      //                                                                ofs_topo_cl_sq,
      //                                                                flow_time_w_total,
      //                                                                epsilon_w,
      //                                                                interval_w);
      //   }
      // }

      inline void gf_flow_with_obs_sq(const GaugeAction& ga,
                                      const GaugeField& gf0,
                                      const FlowInfo2& fi,
                                      std::ofstream& ofs_ene_den,
                                      std::ofstream& ofs_topo,
                                      std::ofstream& ofs_topo2,
                                      std::ofstream& ofs_topo_sq_parity,
                                      std::ofstream& ofs_topo_action_part_parity,
                                      std::ofstream& ofs_topo_det_part_parity,
                                      std::ofstream& ofs_topo_sq_time_rev,
                                      std::ofstream& ofs_topo_action_part_time_rev,
                                      std::ofstream& ofs_topo_det_part_time_rev,
                                      const double flow_time_w_total = 2.0,
                                      const double epsilon_w = 0.01,
                                      const int interval_w = 20
                                      )
      // gf0 is the gauge field which we perform HMC evolution
      {
        TIMER("gf_flow_with_obs");
        const int max_flow_size = fi.get_max_flow_size();
        qassert(max_flow_size!=0);
        const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
        const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
        const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
        GaugeField gf_ext;
        gf_ext.init(geo_ext);
        gf_ext = gf0;

        Obs::observe_sq_before_and_after_smear(ga, gf_ext, fi,
                                               ofs_ene_den,
                                               ofs_topo,
                                               ofs_topo2,
                                               ofs_topo_sq_parity,
                                               ofs_topo_action_part_parity,
                                               ofs_topo_det_part_parity,
                                               ofs_topo_sq_time_rev,
                                               ofs_topo_action_part_time_rev,
                                               ofs_topo_det_part_time_rev,
                                               flow_time_w_total,
                                               epsilon_w,
                                               interval_w);

        for (int i = 0; i < fi.size(); ++i) {
          const FlowStepInfo2& fsi = fi[i];
          const int mask = fsi.mask;
          const int mu = fsi.mu;
          const double epsilon = fsi.epsilon;
          const FlowType& flow_type = *(fsi.ptr_flow_type);
          refresh_expanded_gf_flow_kernel_mask_mu(gf_ext, mask, mu, flow_type);
          gf_flow_mask_mu_no_comm(gf_ext, gf_ext, mask, mu, epsilon, flow_type);
        }
        Obs::observe_sq_before_and_after_smear(ga, gf_ext, fi,
                                               ofs_ene_den,
                                               ofs_topo,
                                               ofs_topo2,
                                               ofs_topo_sq_parity,
                                               ofs_topo_action_part_parity,
                                               ofs_topo_det_part_parity,
                                               ofs_topo_sq_time_rev,
                                               ofs_topo_action_part_time_rev,
                                               ofs_topo_det_part_time_rev,
                                               flow_time_w_total,
                                               epsilon_w,
                                               interval_w);
      }

      // // FlowInfo2 will be repeated n_step times
      // template<class Out>
      // inline void gf_flow_with_obs_specific_flow
      // (
      //  const GaugeField& gf0,
      //  const FlowInfo2& fi,
      //  const int n_steps,
      //  const int measurement_interval,
      //  Out& out
      //  )
      // // gf0 is the gauge field which we perform HMC evolution
      // // The flow operation typically smooth the gauge field.
      // {
      //   TIMER("gf_flow_with_obs");
      //   const int max_flow_size = fi.get_max_flow_size();
      //   qassert(max_flow_size!=0);
      //   const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      //   const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      //   const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
      //   GaugeField gf_ext;
      //   gf_ext.init(geo_ext);
      //   gf_ext = gf0;

      //   Obs::observe_topological_charge(gf_ext, "s=0", out);
      //   for (int s = 1; s <= n_steps; ++s) {
      //     for (int i = 0; i < fi.size(); ++i) {
      //       const FlowStepInfo2& fsi = fi[i];
      //       const int mask = fsi.mask;
      //       const int mu = fsi.mu;
      //       const double epsilon = fsi.epsilon;
      //       const FlowType& flow_type = *(fsi.ptr_flow_type);
      //       refresh_expanded_gf_flow_kernel_mask_mu(gf_ext, mask, mu, flow_type);
      //       gf_flow_mask_mu_no_comm(gf_ext, gf_ext, mask, mu, epsilon, flow_type);
      //     } // for i
      //     if(s%measurement_interval==0)Obs::observe_topological_charge(gf_ext,
      //                                                                  "s="+std::to_string(s),
      //                                                                  out);
      //   } // for s
      // }


      inline void set_flow_all_kernel_mask_mu_no_comm
      (
       FieldM<ColorMatrix, 4>& kernel,
       const GaugeField& gf_ext,
       const std::vector<FlowType>& flow_types
       )
      // cf will be initialized
      // gf_ext need proper communication
      {
        TIMER("set_flow_all_staple_mask_mu_no_comm");
        qassert(is_initialized(gf_ext));
        const Geometry geo = geo_reform(gf_ext.geo());
        kernel.init(geo);
        set_zero(kernel); // F'(y,nu) = 0

        qacc_for(index, geo.local_volume(), {
            const Coordinate xl = geo.coordinate_from_index(index);
            for(int mu=0; mu<DIMN; ++mu){
              ColorMatrix tmp;
              set_zero(tmp);
              for(const FlowType& flow_type : flow_types) {
                tmp += gf_flow_staple_no_comm(gf_ext, xl, mu, flow_type) * flow_type.multiplicative_const;
              }
              const ColorMatrix& u0_x_mu = gf_ext.get_elem(xl, mu);
              kernel.get_elem(xl,mu) = u0_x_mu * matrix_adjoint(tmp);
            }); // end for xl
          }
      }

      inline void set_flow_all_kernel_coeff_mask_mu_no_comm
      (
       FieldM<array<double, 8>, DIMN>& coeff,
       const GaugeField& gf_ext,
       const std::vector<FlowType>& flow_types,
       const ColorMatrixConstants& cmcs
       )
      // cf will be initialized
      // gf_ext need proper communication
      {
        TIMER("set_flow_all_staple_mask_mu_no_comm");
        qassert(is_initialized(gf_ext));
        const Geometry geo = geo_reform(gf_ext.geo());
        coeff.init(geo);
        set_zero(coeff); // F'(y,nu) = 0
        //
        const array<ColorMatrix, 8>& ts = cmcs.ts;

        qacc_for(index, geo.local_volume(), {
            const Coordinate xl = geo.coordinate_from_index(index);
            for(int mu=0; mu<DIMN; ++mu){
              ColorMatrix tmp;
              set_zero(tmp);
              for(const FlowType& flow_type : flow_types) {
                tmp += gf_flow_staple_no_comm(gf_ext, xl, mu, flow_type) * flow_type.multiplicative_const;
              }
              const ColorMatrix& u0_x_mu = gf_ext.get_elem(xl, mu);
              const ColorMatrix matrix = u0_x_mu * matrix_adjoint(tmp);
              //
              array<double, 8>& coeff_x_mu = coeff.get_elem(xl,mu);
              for(int a=0; a<8; ++a) coeff_x_mu[a] = 2.0 * real( matrix_trace(ts[a]*matrix) );
            }); // end for xl
          }
      }


      std::vector<FlowType> w_type
      (
       const int i,
       const int rect_chair //-1
       ){
        assert(0<=i&&i<=7);
        std::vector<FlowType> w;

        if(i==0){
          w.push_back(GlobalScopeHash::get_flow_type.at("plaq"));
        }
        else if( i==1 && rect_chair==0 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("srect"));
          w.push_back(GlobalScopeHash::get_flow_type.at("lrect"));
        }
        else if( i==1 && rect_chair==1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_rho3"));
        }
        else if( i==1 && rect_chair==-1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("srect"));
          w.push_back(GlobalScopeHash::get_flow_type.at("lrect"));
          //
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_rho3"));
        }
        else if( i==2 && rect_chair==0 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_lrect"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_mrect"));
        }
        else if( i==2 && rect_chair==1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_outer_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_outer_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_outer_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_inner"));
        }
        else if( i==2 && rect_chair==-1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_lrect"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_mrect"));
          //
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_outer_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_outer_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_outer_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_inner"));
        }
        else if( i==3 && rect_chair==0 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("srect_distinct_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("srect_distinct_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("srect_distinct_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("lrect_distinct"));
          w.push_back(GlobalScopeHash::get_flow_type.at("mrect_distinct"));
        }
        else if( i==3 && rect_chair==1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typea_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typea_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typea_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typeb_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typeb_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typeb_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_inner"));
        }
        else if( i==3 && rect_chair==-1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("srect_distinct_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("srect_distinct_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("srect_distinct_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("lrect_distinct"));
          w.push_back(GlobalScopeHash::get_flow_type.at("mrect_distinct"));
          //
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typea_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typea_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typea_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typeb_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typeb_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_typeb_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("chair_distinct_inner"));
        }
        else if( i==4 && rect_chair==0 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_distinct_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_distinct_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_distinct_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_lrect_distinct"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_mrect_distinct"));
        }
        else if( i==4 && rect_chair==1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typea_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typea_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typea_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typeb_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typeb_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typeb_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_inner"));
        }
        else if( i==4 && rect_chair==-1 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_distinct_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_distinct_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_srect_distinct_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_lrect_distinct"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_mrect_distinct"));
          //
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typea_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typea_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typea_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typeb_rho1"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typeb_rho2"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_typeb_rho3"));
          w.push_back(GlobalScopeHash::get_flow_type.at("twist_chair_distinct_inner"));
        }
        else if( i==5 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("plaq_twice"));
        }
        else if( i==6 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("plaq_squared"));
        }
        else if( i==7 ){
          w.push_back(GlobalScopeHash::get_flow_type.at("plaq_times_reversed"));
        }
        else assert(false);

        return w;
      }

    } // namespace Obs

  }  // namespace Generalized
}  // namespace qlat


