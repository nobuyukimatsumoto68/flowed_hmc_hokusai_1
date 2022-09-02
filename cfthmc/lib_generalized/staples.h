#pragma once

#include "../flowed-hmc-generalized.h"

namespace qlat{
  namespace Generalized{

    // =========================
    // for staples
    namespace Staples{

      // [d^b_{y,nu} W_j] [d^b_{y,nu} d^a_{x,mu} W_k] [d^a_{x,mu} W_i]
      double set_dw_j_ddw_k_dw_i_from_yl_mu_mask
      (
       const GaugeField& gf_ext,
       // const std::vector<FlowType>& w_j,
       const FieldM<ColorMatrix, DIMN>& w_j_cf,
       const FlowType& w_k_elem, // added outside
       // const std::vector<FlowType>& w_i,
       const FieldM<ColorMatrix, DIMN>& w_i_cf,
       const Coordinate& yl,
       const int mu,
       const int mask // of w_k_elem
       )
      {
        double res = 0.0;
        const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
        const Geometry geo = geo_reform(gf_ext.geo());

        const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,w_k_elem);
        const auto beg = w_k_elem.begin() + w_k_elem.begin_mask_yidx(mask,y_idx);
        // essentially m_prime loop [(x,mu),(y,nu) pair loop]
        for(auto itr=beg; itr!=beg+w_k_elem.multiplicity(mask,y_idx); ) { // itr incremented below
          assert(std::distance(itr,beg)<w_k_elem.multiplicity(mask,y_idx));

          const int m_prime = itr->m_prime();
          Coordinate xl;
          int nu;
          set_xl_nu_from_mask_mu_yl_m_prime( xl, nu, mask, mu, yl, m_prime, geo, w_k_elem);

          AdjointColorMatrix ddw_k_elem_ab;
          set_zero(ddw_k_elem_ab);

          while(itr->m_prime()==m_prime){
            const int m = std::distance(beg,itr);
            const array<ColorMatrix, 2> dw_mats = dw_mat_site_no_comm( gf_ext, xl, mu, yl, w_k_elem, m, mask );
            ddw_k_elem_ab += n_mat_site_no_comm( dw_mats, gf_ext, cmcs(), itr, xl, mu, w_k_elem );
            ++itr;
            if(itr==beg+w_k_elem.multiplicity(mask,y_idx)) break;
          }

          const ColorMatrix dw_j = w_j_cf.get_elem(xl,mu);
          const array<double, 8> dw_j_b = basis_projection_anti_hermitian_matrix( matrix_adjoint( dw_j ) ); // @@@@

          const ColorMatrix dw_i = w_i_cf.get_elem(yl,nu);
          const array<double, 8> dw_i_a = basis_projection_anti_hermitian_matrix( matrix_adjoint(dw_i) ); // @@@@

          for(int a=0; a<8; ++a) for(int b=0; b<8; ++b) res += dw_j_b[b] * ddw_k_elem_ab(a,b) * dw_i_a[a];
        }

        return res;
      }

      // factor 2 be careful!! (from real part)

      // [d^b_{y,nu} W_j] [d^b_{y,nu} d^a_{x,mu} W_k] [d^a_{x,mu} W_i]
      double dw_j_ddw_k_dw_i
      (
       const GaugeField& gf_ext,
       const FieldM<ColorMatrix, DIMN>& w_j_cf,
       const std::vector<FlowType>& w_k,
       const FieldM<ColorMatrix, DIMN>& w_i_cf
       ){
        const Geometry geo = geo_reform(gf_ext.geo());
        double res = 0.0;

        FieldM<double, DIMN> tmp;
        tmp.init(geo);
        set_zero(tmp);

        for(int mu=0; mu<DIMN; ++mu){
          // set_n_mat_mask_mu_no_comm
          //  Field<array<ColorMatrix, 2> > dwf;
          // set_dw_mat_mask_mu_no_comm(dwf, gf0_ext, mask, mu, flow_type);
          qacc_for(index, geo.local_volume(), {
              const Coordinate yl = geo.coordinate_from_index(index);
              for(const FlowType& w_k_elem : w_k ){
                for(int mask=0; mask<w_k_elem.n_masks; ++mask){
                  tmp.get_elem(yl,mu) += set_dw_j_ddw_k_dw_i_from_yl_mu_mask( gf_ext, w_j_cf, w_k_elem, w_i_cf, yl, mu, mask );
                } // end for mask
              } // end for w_k_elem
            }); // end for y
        } // end for mu

#ifdef _OPENMP
#pragma omp parallel for reduction(+:res)
#endif
        for(long index=0; index<geo.local_volume(); ++index) {
          const Coordinate yl = geo.coordinate_from_index(index);
          for(int mu=0; mu<DIMN; ++mu){
            res += tmp.get_elem(yl,mu);
          }
        }
        glb_sum(res);

        return res;
      }


      // -----------------

      // [d^b_{y,nu} W_j] [d^b_{y,nu} d^a_{x,mu} W_k] [d^a_{x,mu} W_i]
      // = dw_j_b *  ddw_k_ab * dw_i_a
      double set_dw_j_ddw_k_dw_i_from_yl_mu_mask_v2
      (
       const GaugeField& gf_ext,
       const FieldM<ColorMatrix, DIMN>& w_j_cf,
       const Field<AdjointColorMatrix>& ddw_k_elem_field,
       const FlowType& w_k_elem,
       const FieldM<ColorMatrix, DIMN>& w_i_cf,
       const Coordinate& yl,
       const int mu,
       const int mask // of w_k_elem
       )
      {
        double res = 0.0;
        const Geometry geo = geo_reform(gf_ext.geo());

        const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,w_k_elem);
        const Vector<AdjointColorMatrix> ddw_k_elem_vect = ddw_k_elem_field.get_elems_const(yl);
        // --- for (x,nu) ---
        for(int m_prime=0; m_prime<w_k_elem.multiplicity_reduced(mask,y_idx); ++m_prime){
          Coordinate xl;
          int nu;
          set_xl_nu_from_mask_mu_yl_m_prime(xl, nu, mask, mu, yl, m_prime, geo, w_k_elem);
          // ---

          const AdjointColorMatrix ddw_k_elem_ab = ddw_k_elem_vect[m_prime];

          const ColorMatrix dw_j = w_j_cf.get_elem(xl,mu);
          const array<double, 8> dw_j_b = basis_projection_anti_hermitian_matrix( matrix_adjoint( dw_j ) ); // @@@@ anti-hermitian?

          const ColorMatrix dw_i = w_i_cf.get_elem(yl,nu);
          const array<double, 8> dw_i_a = basis_projection_anti_hermitian_matrix( matrix_adjoint(dw_i) ); // @@@@ anti-hermitian?

          for(int a=0; a<8; ++a) for(int b=0; b<8; ++b) res += dw_j_b[b] * ddw_k_elem_ab(a,b) * dw_i_a[a];
        }

        return res;
      }

      // factor 2 be careful!! (from real part)

      // [d^b_{y,nu} W_j] [d^b_{y,nu} d^a_{x,mu} W_k] [d^a_{x,mu} W_i]
      double dw_j_ddw_k_dw_i_v2
      (
       const GaugeField& gf_ext,
       const FieldM<ColorMatrix, DIMN>& w_j_cf,
       const std::vector<FlowType>& w_k,
       const FieldM<ColorMatrix, DIMN>& w_i_cf
       ){
        const Geometry geo = geo_reform(gf_ext.geo());
        double res = 0.0;

        FieldM<double, DIMN> tmp;
        tmp.init(geo);
        set_zero(tmp);

        for(int mu=0; mu<DIMN; ++mu){
          for(const FlowType& w_k_elem : w_k ){
            for(int mask=0; mask<w_k_elem.n_masks; ++mask){
              Field<AdjointColorMatrix> ddw_k_elem;
              Field<array<ColorMatrix, 2> > dw_k_elem;
              set_dw_mat_mask_mu_no_comm(dw_k_elem, gf_ext, mask, mu, w_k_elem);
              set_n_mat_mask_mu_no_comm(ddw_k_elem, dw_k_elem, gf_ext, mask, w_k_elem, mu); // possible improvement (reduction k_elem)
              //
              qacc_for(y_idx, geo.local_volume(), {
                  const Coordinate yl = geo.coordinate_from_index(y_idx);
                  tmp.get_elem(yl,mu) += set_dw_j_ddw_k_dw_i_from_yl_mu_mask_v2( gf_ext, w_j_cf,
                                                                                 ddw_k_elem,w_k_elem,
                                                                                 w_i_cf, yl, mu, mask );
                }); // end for y
            } // end for w_k_elem
          } // end for mask
        } // end for mu

#ifdef _OPENMP
#pragma omp parallel for reduction(+:res)
#endif
        for(long index=0; index<geo.local_volume(); ++index) {
          const Coordinate yl = geo.coordinate_from_index(index);
          for(int mu=0; mu<DIMN; ++mu){
            res += tmp.get_elem(yl,mu);
          }
        }
        glb_sum(res);

        return res;
      }


      // -----------------

      // [d^b_{y,nu} W_j] [d^b_{y,nu} d^a_{x,mu} W_k] [d^a_{x,mu} W_i]
      // = dw_j_b *  ddw_k_ab * dw_i_a
      double set_dw_j_ddw_k_dw_i_from_yl_mu_mask_v3
      (
       const GaugeField& gf_ext,
       const FieldM<ColorMatrix, DIMN>& w_j_cf,
       const Field<AdjointColorMatrix>& ddw_k_elem_field,
       const FlowType& w_k_elem,
       const FieldM<ColorMatrix, DIMN>& w_i_cf,
       const Coordinate& yl,
       const int mu,
       const int mask // of w_k_elem
       )
      {
        double res = 0.0;
        const Geometry geo = geo_reform(gf_ext.geo());

        const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,w_k_elem);
        const Vector<AdjointColorMatrix> ddw_k_elem_vect = ddw_k_elem_field.get_elems_const(yl);
        // --- for (x,nu) ---
        for(int m_prime=0; m_prime<w_k_elem.multiplicity_reduced(mask,y_idx); ++m_prime){
          Coordinate xl;
          int nu;
          set_xl_nu_from_mask_mu_yl_m_prime(xl, nu, mask, mu, yl, m_prime, geo, w_k_elem);
          // ---

          const AdjointColorMatrix ddw_k_elem_ab = ddw_k_elem_vect[m_prime];

          const ColorMatrix dw_j = w_j_cf.get_elem(xl,mu);
          const array<double, 8> dw_j_b = basis_projection_anti_hermitian_matrix( matrix_adjoint( dw_j ) ); // @@@@ anti-hermitian?

          const ColorMatrix dw_i = w_i_cf.get_elem(yl,nu);
          const array<double, 8> dw_i_a = basis_projection_anti_hermitian_matrix( matrix_adjoint(dw_i) ); // @@@@ anti-hermitian?

          for(int a=0; a<8; ++a) for(int b=0; b<8; ++b) res += dw_j_b[b] * ddw_k_elem_ab(a,b) * dw_i_a[a];
        }

        return res;
      }

      // factor 2 be careful!! (from real part)


      // void reduce_n_mat_field
      // (
      //  Field<AdjointColorMatrix>& ddw_k_tot_field,
      //  const GaugeField& gf_ext
      //  ){
      //   const Geometry geo = geo_reform(gf_ext.geo());
      //   cf.init(geo);
      //   set_zero(ddw_k_tot_field);

      //   for(const FlowType& w_k_elem : w_k ){
      //     for(int mask=0; mask<w_k_elem.n_masks; ++mask){
      //       for(int mu=0; mu<DIMN; ++mu){
      //         Field<AdjointColorMatrix> ddw_k_elem_field;
      //         Field<array<ColorMatrix, 2> > dw_k_elem_field;
      //         set_dw_mat_mask_mu_no_comm(dw_k_elem_field, gf_ext, mask, mu, w_k_elem);
      //         set_n_mat_mask_mu_no_comm(ddw_k_elem_field, dw_k_elem_field, gf_ext, mask, w_k_elem, mu);

      //         qacc_for(y_idx, geo.local_volume(), {
      //             const Coordinate yl = geo.coordinate_from_index(y_idx);
      //             const Vector<AdjointColorMatrix> ddw_k_elem_vect = ddw_k_elem_field.get_elems_const(yl);
      //             // --- for (x,nu) ---
      //             for(int m_prime=0; m_prime<w_k_elem.multiplicity_reduced(mask,y_idx); ++m_prime){
      //               Coordinate xl;
      //               int nu;
      //               set_xl_nu_from_mask_mu_yl_m_prime(xl, nu, mask, mu, yl, m_prime, geo, w_k_elem);

      //               ddw_k_tot_vect.get_elem(xl,mu) += ddw_k_elem_vect[m_prime];
      //             } // end for m_prime
      //           }); // end for y
      //       } // end for mu
      //     } // end for mask
      //   } // end for w_k_elem
      // }



      // ------------------



      std::vector<FlowType> w_type
      (
       const int i,
       const int rect_chair=-1 // rect=0, chair=1
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




    } // namespace Staples

  }  // namespace Generalized
}  // namespace qlat

