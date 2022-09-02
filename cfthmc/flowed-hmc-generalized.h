#pragma once

#include <qlat/hmc.h>
#include <qlat/wilson-flow.h>
#include <qlat/hmc-stats.h>

#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <unordered_map>

#include "lib_generalized/io.h"
#include "lib_generalized/flow_type.hpp"
#include "lib_generalized/flow_type_examples.hpp"


namespace qlat{
  namespace Generalized{

    // ----------
    // LOOKING FOR AN ALTERNATIVE FOR THIS FUNCTIONALITY.
    /* ---
       IMPORTANT NOTE:
       All FlowType objects in the code refer to the objects
       stored in this map though nested references and pointers.
       The memory of these objects, as well as the hash map itself,
       are allocated on each node of MPI (as ordinary).
       --- */

    namespace GlobalScopeHash{
      std::unordered_map<std::string, FlowType> get_flow_type;
      constexpr int n_type = 44; // number of listed flow types
      inline void setup_hash_map(){
        /* ---
           IMPORTANT NOTES:
           - Derived classes are defined in "lib_generatlized/flow_type_examples.hpp".
           - This implementation copies the FlowType objects to the hash map at insert(...),
             which is confirmed by observing that the addresses of the objects in the hash map
             are different from the temporal objects below.
             --- */
        FlowTypePlaq plaq;
        FlowTypeSrect srect;
        FlowTypeLrect lrect;
        FlowTypeChairRho1 chair_rho1;
        FlowTypeChairRho2 chair_rho2;
        FlowTypeChairRho3 chair_rho3;
        FlowTypeQPlus2 q_plus;
        FlowTypeQMinus2 q_minus;
        // FlowTypeQPlus q_plus;
        // FlowTypeQMinus q_minus;
        FlowTypePlaqTwice plaq_twice;
        FlowTypeTwistSrectRho1 twist_srect_rho1;
        // 10
        FlowTypeTwistSrectRho2 twist_srect_rho2;
        FlowTypeTwistSrectRho3 twist_srect_rho3;
        FlowTypeTwistLrect twist_lrect;
        FlowTypeTwistMrect twist_mrect;
        FlowTypeTwistChairOuterRho1 twist_chair_outer_rho1;
        FlowTypeTwistChairOuterRho2 twist_chair_outer_rho2;
        FlowTypeTwistChairOuterRho3 twist_chair_outer_rho3;
        FlowTypeTwistChairInner twist_chair_inner;
        FlowTypePlaqSquared plaq_squared;
        FlowTypePlaqTimesReversed plaq_times_reversed;
        // 20
        FlowTypeSrectDistinctRho1 srect_distinct_rho1;
        FlowTypeSrectDistinctRho2 srect_distinct_rho2;
        FlowTypeSrectDistinctRho3 srect_distinct_rho3;
        FlowTypeLrectDistinct lrect_distinct;
        FlowTypeMrectDistinct mrect_distinct;
        FlowTypeChairDistinctTypeARho1 chair_distinct_typea_rho1;
        FlowTypeChairDistinctTypeARho2 chair_distinct_typea_rho2;
        FlowTypeChairDistinctTypeARho3 chair_distinct_typea_rho3;
        FlowTypeChairDistinctTypeBRho1 chair_distinct_typeb_rho1;
        FlowTypeChairDistinctTypeBRho2 chair_distinct_typeb_rho2;
        // 30
        FlowTypeChairDistinctTypeBRho3 chair_distinct_typeb_rho3;
        FlowTypeChairDistinctInner chair_distinct_inner;
        FlowTypeTwistSrectDistinctRho1 twist_srect_distinct_rho1;
        FlowTypeTwistSrectDistinctRho2 twist_srect_distinct_rho2;
        FlowTypeTwistSrectDistinctRho3 twist_srect_distinct_rho3;
        FlowTypeTwistLrectDistinct twist_lrect_distinct;
        FlowTypeTwistMrectDistinct twist_mrect_distinct;
        FlowTypeTwistChairDistinctTypeARho1 twist_chair_distinct_typea_rho1;
        FlowTypeTwistChairDistinctTypeARho2 twist_chair_distinct_typea_rho2;
        FlowTypeTwistChairDistinctTypeARho3 twist_chair_distinct_typea_rho3;
        // 40
        FlowTypeTwistChairDistinctTypeBRho1 twist_chair_distinct_typeb_rho1;
        FlowTypeTwistChairDistinctTypeBRho2 twist_chair_distinct_typeb_rho2;
        FlowTypeTwistChairDistinctTypeBRho3 twist_chair_distinct_typeb_rho3;
        FlowTypeTwistChairDistinctInner twist_chair_distinct_inner;
        //
        const std::function<void(const FlowType&)> set = [](const FlowType& flow_type){
          get_flow_type.insert({ flow_type.description, flow_type });
        };
        //
        set(static_cast<FlowType>(plaq));
        set(static_cast<FlowType>(srect));
        set(static_cast<FlowType>(lrect));
        set(static_cast<FlowType>(chair_rho1));
        set(static_cast<FlowType>(chair_rho2));
        set(static_cast<FlowType>(chair_rho3));
        set(static_cast<FlowType>(q_plus));
        set(static_cast<FlowType>(q_minus));
        set(static_cast<FlowType>(plaq_twice));
        set(static_cast<FlowType>(twist_srect_rho1));
        set(static_cast<FlowType>(twist_srect_rho2));
        set(static_cast<FlowType>(twist_srect_rho3));
        set(static_cast<FlowType>(twist_lrect));
        set(static_cast<FlowType>(twist_mrect));
        set(static_cast<FlowType>(twist_chair_outer_rho1));
        set(static_cast<FlowType>(twist_chair_outer_rho2));
        set(static_cast<FlowType>(twist_chair_outer_rho3));
        set(static_cast<FlowType>(twist_chair_inner));
        set(static_cast<FlowType>(plaq_squared));
        set(static_cast<FlowType>(plaq_times_reversed));
        set(static_cast<FlowType>(srect_distinct_rho1));
        set(static_cast<FlowType>(srect_distinct_rho2));
        set(static_cast<FlowType>(srect_distinct_rho3));
        set(static_cast<FlowType>(lrect_distinct));
        set(static_cast<FlowType>(mrect_distinct));
        set(static_cast<FlowType>(chair_distinct_typea_rho1));
        set(static_cast<FlowType>(chair_distinct_typea_rho2));
        set(static_cast<FlowType>(chair_distinct_typea_rho3));
        set(static_cast<FlowType>(chair_distinct_typeb_rho1));
        set(static_cast<FlowType>(chair_distinct_typeb_rho2));
        set(static_cast<FlowType>(chair_distinct_typeb_rho3));
        set(static_cast<FlowType>(chair_distinct_inner));
        set(static_cast<FlowType>(twist_srect_distinct_rho1));
        set(static_cast<FlowType>(twist_srect_distinct_rho2));
        set(static_cast<FlowType>(twist_srect_distinct_rho3));
        set(static_cast<FlowType>(twist_lrect_distinct));
        set(static_cast<FlowType>(twist_mrect_distinct));
        set(static_cast<FlowType>(twist_chair_distinct_typea_rho1));
        set(static_cast<FlowType>(twist_chair_distinct_typea_rho2));
        set(static_cast<FlowType>(twist_chair_distinct_typea_rho3));
        set(static_cast<FlowType>(twist_chair_distinct_typeb_rho1));
        set(static_cast<FlowType>(twist_chair_distinct_typeb_rho2));
        set(static_cast<FlowType>(twist_chair_distinct_typeb_rho3));
        set(static_cast<FlowType>(twist_chair_distinct_inner));
      }
    }

    // =========================
    // functions specific to this Generalized version

    inline Coordinate std2Coordinate( const std::array<int, DIMN>& x ){
      return Coordinate(x[0],x[1],x[2],x[3]);
    }

    inline std::array<int,DIMN> Coordinate2std( const Coordinate& x ){
      return std::array<int,DIMN>{{x[0],x[1],x[2],x[3]}};
    }

    inline int get_y_idx( const Coordinate& yg, const int mu, const FlowType& flow_type ){
      const std::array<int, DIMN> y_prime = flow_type.rotate( Coordinate2std(yg), -mu-1 );
      int y_idx;
      std::array<int,DIMN> y_prime_block;
      flow_type.global2unit(y_idx,y_prime_block,y_prime);
      return y_idx;
    }

    inline void get_y_idx_and_block
    (
     int& y_idx,
     std::array<int,DIMN>& y_prime_block,
     const Coordinate& yg,
     const int mu,
     const FlowType& flow_type
     ){
      const std::array<int, DIMN> y_prime = flow_type.rotate( Coordinate2std(yg), -mu-1 );
      flow_type.global2unit(y_idx,y_prime_block,y_prime);
    }

    inline void set_xl_beg_and_shape_tmp_rotated
    (
     Coordinate& xl_beg,
     std::vector<int>& shape_tmp_rotated,
     const Coordinate& xl,
     const int mu,
     const FlowType& flow_type,
     const int which_staple,
     const int which_loop
     ){
      assert(flow_type.number_of_multiplied_loops()>0);
      assert(static_cast<int>(shape_tmp_rotated.size())==flow_type.length_of_loop);
      xl_beg = xl;
      // rotated and shifted
      for(int rho=0; rho<DIMN; ++rho){
        const int shift = flow_type.multiplied_shapes_relative_coordinates(which_staple,which_loop,rho);

        for(int i=0; i<std::abs(shift); ++i){
          if(shift>0) xl_beg = coordinate_shifts(xl_beg, flow_type.rotate(rho,mu));
          else xl_beg = coordinate_shifts(xl_beg, flow_type.rotate(-rho-1,mu));
        }
      }

      // rotated
      for(int i=0; i<flow_type.length_of_loop; ++i) {
        shape_tmp_rotated[i] = flow_type.rotate(flow_type.multiplied_shape(which_staple,which_loop)[i], mu);
      }
    }

    // search for the case x_prime == x by shifting xl_tmp along the loop
    inline void search_link_idx_corresp_to_flowed_link
    (
     int& link_idx,
     bool& is_pos_dir,
     const Coordinate& xl_beg,
     const std::vector<int>& shape_tmp_rotated, // rotated
     const Coordinate& xl,
     const int mu,
     const FlowType& flow_type
     ){
      assert(flow_type.number_of_multiplied_loops()>0);
      assert(static_cast<int>(shape_tmp_rotated.size())==flow_type.length_of_loop);

      Coordinate xl_tmp = xl_beg;
      for( link_idx=0; link_idx<flow_type.length_of_loop; ++link_idx){
        const Coordinate xl_tmp_next = coordinate_shifts(xl_tmp, shape_tmp_rotated[link_idx]);

        if(xl_tmp==xl && shape_tmp_rotated[link_idx]==mu){ // if the link is (x,mu) (positive directed)
          is_pos_dir = true;
          break; // each loop is assumed to be linear
        }
        else if(xl_tmp_next==xl && shape_tmp_rotated[link_idx]==-mu-1){ // if the link is (x,mu) (negative directed)
          is_pos_dir = false;
          break; // each loop is assumed to be linear
        }
        else {
          xl_tmp = xl_tmp_next; // continue for link_idx loop
          continue;
        }
      }
      qassert(link_idx<flow_type.length_of_loop); // must have the flowed link when nonlinear==true
    }


    // =========================
    // "2" is supllied for disambiguation
    class FlowStepInfo2 {
    public:
      const int mask;
      const int mu;
      double epsilon;
      const FlowType* ptr_flow_type;
    private:
      bool is_measurement_after_this_step;
      //
    public:
      FlowStepInfo2()
        : mask(0)
        , mu(0)
        , epsilon(0.0)
        , ptr_flow_type(nullptr)
        , is_measurement_after_this_step(false)
      {}
      FlowStepInfo2(const int mask_, const int mu_, const double epsilon_,
                    const FlowType& flow_type_)
        : mask(mask_), mu(mu_), epsilon(epsilon_)
        , ptr_flow_type(&flow_type_)
        , is_measurement_after_this_step(false)
      {
        assert(typeid(flow_type_)==typeid(Generalized::FlowType));
      }

      void turn_on_measurement_after_this_step() & {
        this->is_measurement_after_this_step = true;
      }

      const bool& is_measurement() const& { return this->is_measurement_after_this_step; }
      void flip_sign_of_t(){ this->epsilon *= -1.0; }
    };

    class FlowInfo2 {
    public:
      std::vector<FlowStepInfo2> v;
      int max_flow_size = 0;

      FlowInfo2( const double epsilon = 0.0,
                 const double epsilon_rect = 0.0,
                 const double epsilon_chair = 0.0,
                 const double epsilon_q = 0.0,
                 const double epsilon_plaq_twice = 0.0,
                 const double epsilon_twist_rect = 0.0,
                 const double epsilon_twist_chair = 0.0,
                 const double epsilon_plaq_squared = 0.0,
                 const double epsilon_plaq_times_reversed = 0.0,
                 const double epsilon_rect_distinct = 0.0,
                 const double epsilon_chair_distinct = 0.0,
                 const double epsilon_twist_rect_distinct = 0.0,
                 const double epsilon_twist_chair_distinct = 0.0
                 );

      FlowInfo2(const bool is_luscher_flow,
                const double beta,
                const double tmax,
                const int nstep,
                const bool is_first_order = true );

      FlowInfo2( const double theta,
                 const int nstep); // theta flow

      FlowInfo2( const double beta,
                 const double tmax,
                 const int nstep );

      FlowInfo2( const double tmax,
                 const int nstep,
                 const std::vector<std::vector<double>>& gamma_list );

      const FlowStepInfo2& operator[](const std::size_t i) const& { return v[i]; }
      const int& get_max_flow_size() const& { return max_flow_size; }
      int size() const& { return static_cast<int>(v.size()); }
      void turn_on_measurement_at_last_step() & {
        (v.end()-1)->turn_on_measurement_after_this_step();
      }
      void flip_sign_of_t(){
        for(FlowStepInfo2& fsi : this->v) fsi.flip_sign_of_t();
        std::vector<FlowStepInfo2> v2;
        for(std::size_t i=0; i<this->v.size(); ++i ) v2.push_back(this->v[this->v.size()-i-1]);
        assert(v2.size() == this->v.size() );
        this->v.clear();
        for(std::size_t i=0; i<v2.size(); ++i ) v.push_back(v2[i]);
      }

      std::vector<double> get_gamma_coeffs( const double beta ) const&;

      std::vector<FlowStepInfo2>::const_iterator begin() const& { return v.begin(); }
      std::vector<FlowStepInfo2>::const_iterator end() const& { return v.end(); }

    private:
      void set_max_flow_size( const std::unordered_map<std::string, FlowType>& hash) &;
      void set_from_hash(const std::unordered_map<std::string, FlowType>& hash,
                         const std::string& description,
                         const double step_size) &; // essentially a push_back function
      void set_from_hash_q(const std::unordered_map<std::string, FlowType>& hash,
                           const double step_size) &; // essentially a push_back function
    };

    // -------------------------------------------------------------------------

    inline std::string show(const FlowInfo2& fi)
    {
      std::ostringstream out;
      for (int i=0; i<fi.size(); ++i) {
        const FlowStepInfo2& fsi = fi[i];
        out << ssprintf("fi[%d]: mask=%d, mu=%d, epsilon=%.15f, description=%s, is_flowed_link_nonlinear=%d, measurement=%d",
                        i, fsi.mask, fsi.mu, fsi.epsilon,
                        fsi.ptr_flow_type->description.c_str(),
                        fsi.ptr_flow_type->is_flowed_link_nonlinear, fsi.is_measurement());
        if (i != fi.size() - 1) out << std::endl;
      }
      return out.str();
    }


    inline int mask_from_coordinate(const Coordinate& xg,
                                    const FlowType& flow_type,
                                    const int mu)
    { return flow_type.mask( get_y_idx(xg,mu,flow_type) ); }


    inline const vector<long>& get_flowed_hmc_indices_mask_flow_size
    (
     const Geometry& geo, const int mask, const FlowType& flow_type, const int mu
     )
    // mu: differentiated link (=flowed link)
    // nu: direction to which rectangles extend
    {
      // static Cache<std::string, vector<long> > cache("flowed_hmc_indices_cache", 8, 2);
      static Cache<std::string, vector<long> > cache("flowed_hmc_indices_cache", 1024, 2);
      const std::string key =
        ssprintf("%s-%d-%s-%d", show(geo.node_site).c_str(), mask,
                 flow_type.description.c_str(), mu);
      vector<long>& vec = cache[key];
      if (vec.size() == 0) {
        long count = 0;
        qfor(index, geo.local_volume(), {
            const Coordinate xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index));
            const int mask_xl = mask_from_coordinate(xg, flow_type, mu);
            if (mask_xl == mask) { count += 1; } // check the mask type
          }); // end qfor
        vec.resize(count);
        count = 0;
        qfor(index, geo.local_volume(), {
            const Coordinate xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index));
            const int mask_xl = mask_from_coordinate(xg, flow_type, mu);
            if (mask_xl == mask) {
              vec[count] = index;
              count += 1;
            }
          }); // end qfor
      } // end if vec.size() == 0
      return cache[key];
    }


    qacc void set_xl_nu_from_mask_mu_yl_m(Coordinate& xl, int& nu,
                                          const int mask, const int mu,
                                          const Coordinate& yl, const int m,
                                          const Geometry& geo, const FlowType& flow_type)
    // mask, mu are the flow parameters
    // xl, mu is the link to be flowed
    // yl, nu is the link to calculate derivative
    // given yl, m is the index of the pair: (xl, mu) (yl, nu)
    {
      const Coordinate yg = geo.coordinate_g_from_l(yl);
      const std::array<int, DIMN> y_prime = flow_type.rotate( Coordinate2std(yg), -mu-1 );
      int y_idx;
      std::array<int,DIMN> y_prime_block;
      flow_type.global2unit(y_idx,y_prime_block,y_prime);

      const auto itr = flow_type.begin() + flow_type.begin_mask_yidx(mask,y_idx) + m;
      const int x_prime_idx = itr->x_idx();
      const int nu_prime = itr->nu();

      std::array<int,DIMN> block_x_prime;
      for(int i=0; i<DIMN; ++i) block_x_prime[i] = y_prime_block[i] - itr->y_block(i);
      const std::array<int,DIMN> x_prime = flow_type.unit2global(x_prime_idx,block_x_prime);
      const Coordinate xg = std2Coordinate( flow_type.rotate( x_prime, mu ) );

      nu = flow_type.rotate( nu_prime, mu );
      xl = geo.coordinate_l_from_g(xg);
    }


    qacc void set_xl_nu_from_mask_mu_yl_m_prime(Coordinate& xl, int& nu,
                                                const int mask, const int mu,
                                                const Coordinate& yl, const int m_prime,
                                                const Geometry& geo, const FlowType& flow_type)
    // mask, mu are the flow parameters
    // xl, mu is the link to be flowed
    // yl, nu is the link to calculate derivative
    // given yl, m_prime is the index of the pair: (xl, mu) (yl, nu)
    {
      const Coordinate yg = geo.coordinate_g_from_l(yl);
      const std::array<int, DIMN> y_prime = flow_type.rotate( Coordinate2std(yg), -mu-1 );
      int y_idx;
      std::array<int,DIMN> y_prime_block;
      flow_type.global2unit(y_idx,y_prime_block,y_prime);

      const auto itr_begin = flow_type.begin() + flow_type.begin_mask_yidx(mask,y_idx);
      const auto itr = std::find_if(itr_begin, itr_begin+flow_type.multiplicity(mask,y_idx),
                                    [&m_prime](const RowOfTable& row)
                                    { return row.m_prime()==m_prime; });

      const int x_prime_idx = itr->x_idx();
      const int nu_prime = itr->nu();
      std::array<int,DIMN> block_x_prime;
      for(int i=0; i<DIMN; ++i) block_x_prime[i] = y_prime_block[i] - itr->y_block(i);
      const std::array<int,DIMN> x_prime = flow_type.unit2global(x_prime_idx,block_x_prime);
      const Coordinate xg = std2Coordinate( flow_type.rotate( x_prime, mu ) );

      nu = flow_type.rotate( nu_prime, mu );
      xl = geo.coordinate_l_from_g(xg);
    }


    // ################## PRIMITVE ########################
    // returns the staple attached to the link U_{x,mu}
    // note that the staple starts from x
    qacc ColorMatrix gf_flow_staple_no_comm
    (
     const GaugeField& gf_ext,
     const Coordinate& xl,
     const int mu,
     const FlowType& flow_type
     )
    {
      ColorMatrix acc;
      set_zero(acc);

      for (int l=0; l<flow_type.n_loops; ++l){
        // main loop
        std::vector<int> shape(flow_type.length_of_loop-1);
        // reversed
        for(int i=1; i<flow_type.length_of_loop; ++i) {
          shape[flow_type.length_of_loop-i-1] = -flow_type.rotate(flow_type.shape(l)[i], mu)-1;
        }

        // sub loops
        std::complex<double> multiplied_factor = 0.0;
        if(flow_type.number_of_multiplied_loops()==0) multiplied_factor = 1.0;
        else{
          for(int k=0; k<flow_type.number_of_multiplied_loops(); ++k){
            Coordinate xl_beg;
            std::vector<int> shape_tmp_rotated(flow_type.length_of_loop);
            set_xl_beg_and_shape_tmp_rotated(xl_beg,shape_tmp_rotated,xl,mu,flow_type,l,k);

            multiplied_factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_beg, shape_tmp_rotated ) );
          } // for k
        }
        // assert(multiplied_factor==1.0); // debug

        acc += std::conj(multiplied_factor) * gf_wilson_line_no_comm(gf_ext, xl, shape);
      } // end of l
      return acc;
    }


    // wrapping gf_flow_staple_no_comm
    inline void set_flow_staple_mask_mu_no_comm
    (
     FieldM<ColorMatrix, 1>& cf,
     const GaugeField& gf_ext,
     const int mask,
     const int mu,
     const FlowType& flow_type
     )
    // cf will be initialized
    // gf_ext need proper communication
    // mask: flow 1:odd / 2:even site
    // mu: flow link direction
    // nu: direction to which rectangles extend
    {
      TIMER("set_flow_staple_mask_mu_no_comm");
      qassert(is_initialized(gf_ext));
      qassert(0 <= mu && mu < 4);
      const Geometry geo = geo_reform(gf_ext.geo());
      cf.init(geo);

      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);
          cf.get_elem(xl) = gf_flow_staple_no_comm(gf_ext, xl, mu, flow_type);
        }); // end qacc_for
    }

    inline void gf_flow_mask_mu_no_comm(GaugeField& gf,
                                        const GaugeField& gf0_ext, const int mask,
                                        const int mu, const double epsilon,
                                        const FlowType& flow_type)
    // mask: masking types
    // mu: flow link direction
    // epsilon: is the flow step size
    // nu: direction to which rectangles extend
    {
      TIMER("gf_flow_mask_mu_no_comm");
      qassert(is_initialized(gf0_ext));
      qassert(0 <= mu && mu < 4);
      const Geometry geo = geo_reform(gf0_ext.geo());
      gf.init(geo);
      gf = gf0_ext;
      FieldM<ColorMatrix, 1> cf;
      set_flow_staple_mask_mu_no_comm(cf, gf0_ext, mask, mu, flow_type);

      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);
          const ColorMatrix& c_x_mu = cf.get_elem(xl); // c_x_mu = Z(U)_{x,mu}
          const ColorMatrix& u0_x_mu = gf0_ext.get_elem(xl, mu);
          const ColorMatrix z_u_x_mu = -make_tr_less_anti_herm_matrix(u0_x_mu * matrix_adjoint(c_x_mu));
          const ColorMatrix e_z_u_x_mu = static_cast<Complex>(epsilon) * z_u_x_mu;
          const ColorMatrix e_u_x_mu = make_matrix_exp(e_z_u_x_mu) * u0_x_mu;
          ColorMatrix& u_x_mu = gf.get_elem(xl, mu);
          u_x_mu = e_u_x_mu;
        }); // end qacc for
    }


    // FOR COMMUNICATION
    inline void set_marks_flow_kernel_mask_mu
    (
     CommMarks& marks,
     const Geometry& geo,
     const std::string& tag
     )
    {
      TIMER_VERBOSE("set_marks_flow_kernel_mask_mu");
      qassert(geo.multiplicity == 4);
      marks.init();
      marks.init(geo);
      set_zero(marks);
      const std::vector<std::string> words = split_line_with_spaces(tag);
      qassert(words.size() == 3);
      const int mask = read_long(words[0]);
      const int mu = read_long(words[1]);
      const std::string description = words[2];
      qassert(description.size()!=0);

      const FlowType& flow_type = GlobalScopeHash::get_flow_type.at(description);

      qassert(0 <= mu && mu < 4);
      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);

          for (int l=0; l<flow_type.n_loops; ++l) {
            std::vector<int> shape(flow_type.length_of_loop-1,0);
            for(int i=1; i<flow_type.length_of_loop; ++i) {
              shape[flow_type.length_of_loop-i-1] = -flow_type.rotate(flow_type.shape(l)[i], mu)-1;
            }
            set_marks_field_path( marks, xl, shape );

            for(int k=0; k<flow_type.number_of_multiplied_loops(); ++k){
              Coordinate xl_beg;
              std::vector<int> shape_tmp_rotated(flow_type.length_of_loop);
              set_xl_beg_and_shape_tmp_rotated(xl_beg,shape_tmp_rotated,xl,mu,flow_type,l,k);

              set_marks_field_path( marks, xl_beg, shape_tmp_rotated );
            }
          }
        }); // end qacc_for
    }


    inline void refresh_expanded_gf_flow_kernel_mask_mu(GaugeField& gf_ext,
                                                        const int mask, const int mu,
                                                        const FlowType& flow_type)
    {
      TIMER("refresh_expanded_gf_flow_kernel_mask_mu");
      const CommPlan& plan = get_comm_plan(set_marks_flow_kernel_mask_mu,
                                           ssprintf("%d %d %s", mask, mu, flow_type.description.c_str()),
                                           gf_ext.geo());
      refresh_expanded(gf_ext, plan);
    }


    inline void gf_flow_inv_mask_mu_no_comm
    (
     GaugeField& gf, const GaugeField& gf1_ext, const int mask, const int mu,
     const double epsilon, const FlowType& flow_type, const int n_iter = 50 // 100 // !!!
     )
    // mu: flow link direction
    // epsilon: is the flow step size
    {
      // constexpr double THRESHOLD_ITER = 1.0e-15;

      TIMER("gf_flow_inv_mask_mu_no_comm");
      qassert(is_initialized(gf1_ext));
      qassert(0 <= mu && mu < 4);
      const Geometry geo = geo_reform(gf1_ext.geo());
      gf.init(geo);
      gf = gf1_ext;

      const Geometry geo_ext = geo_reform(gf1_ext.geo(), 1, gf1_ext.geo().expansion_left, gf1_ext.geo().expansion_right);

      GaugeField gf_tmp; // = gf1_ext except for the link (xl, mu)
      gf_tmp.init(geo_ext);
      gf_tmp = gf1_ext;
      refresh_expanded(gf_tmp);

      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {

          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);

          ColorMatrix& uopt_x_mu = gf_tmp.get_elem(xl, mu); // X=0 at zeroth step
          const ColorMatrix uprime_x_mu = uopt_x_mu;
          ColorMatrix X;
          set_zero(X);
          ColorMatrix fX;
          set_zero(fX);

          for (int n=0; n<n_iter; ++n) {
            const ColorMatrix staple = gf_flow_staple_no_comm(gf_tmp, xl, mu, flow_type);
            const ColorMatrix c_x_mu_dagger = matrix_adjoint(staple);
            X = fX;
            ColorMatrix fX = -make_tr_less_anti_herm_matrix(uopt_x_mu * c_x_mu_dagger);
            // norm = Obs::frobenius_sq(X-fX);
            const ColorMatrix exp_minus_eps_fX = make_matrix_exp(- static_cast<Complex>(epsilon) * fX );
            uopt_x_mu = exp_minus_eps_fX * uprime_x_mu;

            // if( norm/3.0 < THRESHOLD_ITER ) break;
          }
          // if(n==n_iter) {
          //   std::stringstream ss;
          //   ss << "fixed point iteration did not convergved." << std::endl;
          //   Obs::print_to(std::cerr, ss.str());
          //   qassert(false);
          // }
          // const double norm = Obs::frobenius_sq(X-fX);
          // if( norm/std::sqrt(6.0) > THRESHOLD_ITER ) {
          //   std::stringstream ss;
          //   ss << "fixed point iteration did not convergved." << std::endl;
          //   Obs::print_to(std::cerr, ss.str());
          //   qassert(false);
          // }

          ColorMatrix& u_x_mu = gf.get_elem(xl, mu);
          u_x_mu = uopt_x_mu;
        });
    }


    inline void gf_flow
    (
     GaugeField& gf,
     const GaugeField& gf0,
     const FlowInfo2& fi
     )
    // gf0 is the gauge field which we perform HMC evolution
    // gf is supposed to be the desired configuration
    // The flow operation typically smooth the gauge field.
    // Normally call this at the end of Flowed-HMC evolution step.
    {
      TIMER("gf_flow");
      const int max_flow_size = fi.get_max_flow_size();
      qassert(max_flow_size!=0);
      const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Geometry geo = geo_reform(gf0.geo());
      const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
      GaugeField gf_ext;
      gf_ext.init(geo_ext);
      gf_ext = gf0;

      for (int i = 0; i<fi.size(); ++i) {
        const FlowStepInfo2& fsi = fi[i];
        const int mask = fsi.mask;
        const int mu = fsi.mu;
        const double epsilon = fsi.epsilon;
        const FlowType& flow_type = *(fsi.ptr_flow_type);
        refresh_expanded_gf_flow_kernel_mask_mu(gf_ext, mask, mu, flow_type);
        gf_flow_mask_mu_no_comm(gf_ext, gf_ext, mask, mu, epsilon, flow_type);
      }
      gf.init(geo);
      gf = gf_ext;
    }

    inline void gf_flow_inv(GaugeField& gf, const GaugeField& gf1,
                            const FlowInfo2& fi)
    // gf1 is supposed to be the desired configuration
    // gf is the gauge field which we perform HMC evolution
    // The flow operation typically smooth the gauge field
    // Normally call this at the beginning of Flowed-HMC evolution step.
    {
      TIMER("gf_flow_inv");
      const int max_flow_size = fi.get_max_flow_size();
      qassert(max_flow_size!=0);
      const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Geometry geo = geo_reform(gf1.geo());
      const Geometry geo_ext = geo_reform(gf1.geo(), 1, expand_left, expand_right);
      GaugeField gf_ext;
      gf_ext.init(geo_ext);
      gf_ext = gf1;

      for (int i = fi.size()-1; i>=0; --i) {
        const FlowStepInfo2& fsi = fi[i];
        const int mask = fsi.mask;
        const int mu = fsi.mu;
        const double epsilon = fsi.epsilon;
        const FlowType& flow_type = *(fsi.ptr_flow_type);
        refresh_expanded_gf_flow_kernel_mask_mu(gf_ext, mask, mu, flow_type);
        gf_flow_inv_mask_mu_no_comm(gf_ext, gf_ext, mask, mu, epsilon, flow_type);
      }
      gf.init(geo);
      gf = gf_ext;
    }


    // ################## PRIMITVE FUNCTION ########################
    // there is a consistency condition among x,mu,y and m.
    // This is only a part of the formula:
    // d{W(x,mu)}/{ds^b(y,nu)} = sum_{m (d not acting on A)} A * pre_m * T^b * post_m
    //                         + sum_{m (d acting on A)} tr[ pre_m * T^b * post_m ] * U * Cdagger .
    // returns {pre_m, post_m} as a two_component vector of matrices
    array<ColorMatrix, 2> dw_mat_site_no_comm // edited
    (
     const GaugeField& gf_ext,
     const Coordinate& xl,
     const int mu,
     const Coordinate& yl,
     const FlowType& flow_type,
     const int m, // index of the involved pair (x,mu) (y,nu)
     const int mask
     )
    {
      array<ColorMatrix, 2> dw_mats;
      ColorMatrix& dw_pre = dw_mats[0];
      ColorMatrix& dw_post = dw_mats[1];
      const Geometry geo = geo_reform(gf_ext.geo());

      const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,flow_type);
      const auto itr = flow_type.begin() + flow_type.begin_mask_yidx(mask,y_idx) + m;

      const bool is_pos_dir = itr->is_positive();
      const int link_idx = itr->s();

      const int which_loop = itr->which_loop();

      // set shape (and shift x if necessary)
      if(which_loop==-1) { // if differentiated in the main loop
        Coordinate xl_tmp = xl;
        const std::vector<int> shape = flow_type.shape( itr->l() );

        set_unit(dw_pre);
        for(int s=0; s<link_idx + 1-is_pos_dir; ++s) {
          const int dir = flow_type.rotate(shape[s],mu);
          dw_pre *= gf_get_link(gf_ext, xl_tmp, dir);
          xl_tmp = coordinate_shifts(xl_tmp, dir);
        }

        set_unit(dw_post);
        for(std::size_t s=link_idx + 1-is_pos_dir; s<shape.size(); ++s) {
          const int dir = flow_type.rotate(shape[s],mu);
          dw_post *= gf_get_link(gf_ext, xl_tmp, dir);
          xl_tmp = coordinate_shifts(xl_tmp, dir);
        }
      }
      else{ // if differentiated in a sub loop
        Coordinate xl_tmp;
        std::vector<int> shape_tmp_rotated(flow_type.length_of_loop);
        set_xl_beg_and_shape_tmp_rotated(xl_tmp,shape_tmp_rotated,xl,mu,flow_type,itr->l(),which_loop);

        set_unit(dw_pre);
        for(int s=0; s<link_idx + 1-is_pos_dir; ++s) {
          const int dir = shape_tmp_rotated[s];
          dw_pre *= gf_get_link(gf_ext, xl_tmp, dir);
          xl_tmp = coordinate_shifts(xl_tmp, dir);
        }

        set_unit(dw_post);
        for(int s=link_idx + 1-is_pos_dir; s<flow_type.length_of_loop; ++s) {
          const int dir = shape_tmp_rotated[s];
          dw_post *= gf_get_link(gf_ext, xl_tmp, dir);
          xl_tmp = coordinate_shifts(xl_tmp, dir);
        }
      }

      if(!is_pos_dir) dw_pre = -dw_pre;

      return dw_mats;
    }


    // wrapping dw_mat_site_no_comm
    // This is only a part of the formula:
    // d{W(x,mu)}/{ds^b(y,nu)} = sum_{m (d not acting on A)} A * pre_m * T^b * post_m
    //                         + sum_{m (d acting on A)} tr[ pre_m * T^b * post_m ] * U * Cdagger .
    // we set two the two-comonent vectors {pre_m, post_m} for all m
    inline void set_dw_mat_mask_mu_no_comm // edited
    (
     Field<array<ColorMatrix, 2> >& dwf,
     const GaugeField& gf_ext,
     const int mask,
     const int mu,
     const FlowType& flow_type
     )
    // dwf does NOT need to be initialized.
    // It will be initialized with no expansion
    // See set_xl_nu_from_mask_mu_yl_m, n_mat_site_no_comm
    // gf_ext need proper communication
    // mu: flow link direction
    {
      TIMER("set_dw_mat_mask_mu_no_comm");
      qassert(is_initialized(gf_ext));
      qassert(0 <= mu && mu < 4);

      const Geometry geo = geo_reform(gf_ext.geo());
      dwf.init(geo, flow_type.multiplicity_max() );
      qacc_for(index, geo.local_volume(), {
          const Coordinate yl = geo.coordinate_from_index(index);
          Vector<array<ColorMatrix, 2> > dwfv = dwf.get_elems(yl);

          const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,flow_type);
          for (int m = 0; m <flow_type.multiplicity(mask,y_idx); ++m) {
            Coordinate xl;
            int nu = 0; // dummy
            set_xl_nu_from_mask_mu_yl_m(xl, nu, mask, mu, yl, m, geo, flow_type);
            dwfv[m] = dw_mat_site_no_comm(gf_ext, xl, mu, yl, flow_type, m, mask);
          }
        });
    }


    // ################## PRIMITVE FUNCTION ########################
    //  proj^a[ A * pre_m * T^b * post_m ]
    // or
    //  proj^a[ tr[ pre_m * T^b * post_m ] * U * Cdagger ]
    // given pre, post.
    qacc AdjointColorMatrix n_mat_site_no_comm // edited
    (
     const array<ColorMatrix, 2>& dw_mats, // pre, post
     const GaugeField& gf_ext,
     const ColorMatrixConstants& cmcs,
     const std::vector<RowOfTable>::const_iterator& itr,
     const Coordinate& xl,
     const int mu,
     const FlowType& flow_type
     )
    {
      const array<ColorMatrix, 8>& ts = cmcs.ts;
      const ColorMatrix& dw_pre = dw_mats[0];
      const ColorMatrix& dw_post = dw_mats[1];
      AdjointColorMatrix n_mat;

      if(itr->which_loop()==-1){
        std::complex<double> multiplied_factor=0.0;
        if(flow_type.number_of_multiplied_loops()==0) multiplied_factor = 1.0;
        else{ // prepare multiplied factor
          for(int k=0; k<flow_type.number_of_multiplied_loops(); ++k){
            Coordinate xl_beg;
            std::vector<int> shape_tmp_rotated(flow_type.length_of_loop);
            set_xl_beg_and_shape_tmp_rotated(xl_beg,shape_tmp_rotated,xl,mu,flow_type,itr->l(),k);
            multiplied_factor += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_beg, shape_tmp_rotated ) );
          }
        }

        for (int b = 0; b < 8; ++b) {
          const ColorMatrix dw = dw_pre * ts[b] * dw_post;
          const ColorMatrix n_b = make_tr_less_anti_herm_matrix(multiplied_factor * dw);
          const array<double, 8> basis_b = basis_projection_anti_hermitian_matrix(n_b);
          for (int a = 0; a < 8; ++a) n_mat(a, b) = basis_b[a];
        }
      } // end if(which_loop==-1)
      else{
        for (int b = 0; b < 8; ++b) {
          const ColorMatrix dw = dw_pre * ts[b] * dw_post;
          std::complex<double> multiplied_factor = matrix_trace(dw);

          std::vector<int> shape(flow_type.length_of_loop);
          for(int i=0; i<flow_type.length_of_loop; ++i) {
            shape[i] = flow_type.rotate(flow_type.shape(itr->l())[i], mu);
          }

          const ColorMatrix wilson_line = gf_wilson_line_no_comm( gf_ext, xl, shape );

          const ColorMatrix n_b = make_tr_less_anti_herm_matrix(multiplied_factor * wilson_line);
          const array<double, 8> basis_b = basis_projection_anti_hermitian_matrix(n_b);
          for (int a = 0; a < 8; ++a) n_mat(a, b) = basis_b[a];
        }
      } // end if(which_loop!=-1)
      return n_mat;
    }


    // wrapping n_mat_site_no_comm // edited
    inline void set_n_mat_mask_mu_no_comm
    (
     Field<AdjointColorMatrix>& nf, // N(mu) = ( N(x,mu;y,nu) )_{x,y,nu}, to be initialized below with ducf geometry
     const Field<array<ColorMatrix, 2> >& dwf,
     const GaugeField& gf_ext,
     const int mask,
     const FlowType& flow_type,
     const int mu
     )
    {
      TIMER("set_n_mat_mask_mu_no_comm");
      qassert(is_initialized(dwf));
      qassert( dwf.geo().multiplicity == flow_type.multiplicity_max() );
      const Geometry& geo = dwf.geo();
      nf.init(geo, flow_type.multiplicity_reduced_max());
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();

      // --- for y ---
      qacc_for(index, geo.local_volume(), {
          const Coordinate yl = geo.coordinate_from_index(index);
          const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,flow_type);
          // ---

          const Vector<array<ColorMatrix, 2> > dwfv = dwf.get_elems_const(yl);
          Vector<AdjointColorMatrix> nfv = nf.get_elems(yl); // vec{N}(mu,y) = ( N(x,mu;y,nu) )_{x,nu}
          for(int i=0; i<nfv.size(); ++i) set_zero(nfv[i]); // vec{N}(mu,y) = 0

          // --- for x,mu (removing duplicates) ---
          const auto beg = flow_type.begin()+flow_type.begin_mask_yidx(mask,y_idx);
          for(auto itr=beg; itr!=beg+flow_type.multiplicity(mask,y_idx); ++itr) {
            const int m = std::distance(beg,itr);
            Coordinate xl;
            int nu = 0; // dummy
            set_xl_nu_from_mask_mu_yl_m(xl, nu, mask, mu, yl, m, geo, flow_type);

            nfv[itr->m_prime()] += n_mat_site_no_comm(dwfv[m],gf_ext,cmcs(),itr,xl,mu,flow_type);
          }
        });
    }


    // diag part // edited
    inline void set_n_mat_mask_mu_no_comm_diag
    (
     FieldM<AdjointColorMatrix, 1>& nf,
     const GaugeField& gf_ext,
     const int mask,
     const int mu,
     const FlowType& flow_type
     )
    // nf does NOT need to be initialized.
    // It will be initialized with no expansion
    // nf: only xl, mu
    // gf_ext need proper communication
    // mask: flow 1:odd / 2:even site
    // mu: flow link direction
    {
      TIMER("set_n_mat_mask_mu_no_comm(nf)");
      qassert(is_initialized(gf_ext));
      qassert(0 <= mu && mu < 4);
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
      const Geometry geo = geo_reform(gf_ext.geo());
      nf.init(geo);

      Field<array<ColorMatrix, 2> > dwf;
      set_dw_mat_mask_mu_no_comm(dwf, gf_ext, mask, mu, flow_type);

      // --- for x ---
      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);
          const int x_idx = get_y_idx(geo.coordinate_g_from_l(xl),mu,flow_type);
          // ---
          const Vector<array<ColorMatrix, 2> > dwfv = dwf.get_elems_const(xl);
          AdjointColorMatrix& n_mat = nf.get_elem(xl); // N(x,mu)
          set_zero(n_mat);

          // --- for x,mu (removing duplicates) ---
          const auto beg = flow_type.begin()+flow_type.begin_mask_yidx(mask,x_idx);
          for(auto itr=beg; itr!=beg+flow_type.multiplicity(mask,x_idx); ++itr) {
            if( itr->x_idx()==x_idx && itr->nu()==0
                // && itr->is_positive()==1 // CAUTION
                && itr->y_block(0)==0 && itr->y_block(1)==0 && itr->y_block(2)==0 && itr->y_block(3)==0
                ){ // if (xl,mu)==(yl,nu)
              const int m = std::distance(beg,itr);
              n_mat += n_mat_site_no_comm(dwfv[m],gf_ext,cmcs(),itr,xl,mu,flow_type);
            }
          }
        });
    }


    inline void set_ad_x_and_j_n_x_mask_mu_no_comm
    (
     FieldM<AdjointColorMatrix, 2>& f_ad_x_and_j_n_x,
     const FieldM<ColorMatrix, 1>& cf,
     const GaugeField& gf,
     const int mask,
     const int mu,
     const double epsilon,
     const FlowType& flow_type
     )
    {
      TIMER("set_ad_x_and_j_n_x_mask_mu_no_comm");
      qassert(is_initialized(cf));
      qassert(is_initialized(gf));
      qassert(0 <= mu && mu < 4);
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
      const Geometry geo = geo_reform(gf.geo());
      f_ad_x_and_j_n_x.init(geo);

      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);

          const ColorMatrix& u = gf.get_elem(xl, mu);
          const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
          const ColorMatrix x_mat = static_cast<Complex>(-epsilon) * make_tr_less_anti_herm_matrix(u * c_dagger);

          Vector<AdjointColorMatrix> ad_x_and_j_n_x = f_ad_x_and_j_n_x.get_elems(xl);
          AdjointColorMatrix& ad_x_mat = ad_x_and_j_n_x[0];
          AdjointColorMatrix& j_n_x_mat = ad_x_and_j_n_x[1];
          ad_x_mat = make_adjoint_representation(x_mat, cmcs());
          j_n_x_mat = make_diff_exp_map(-ad_x_mat);
        });
    }


    // given x,mu,y,nu,
    // return M(x,mu;y,nu); M(x,mu;y,nu) = delta_{x,mu;y,nu}e^{AdX(x,mu)} - eps * J(-X(x,mu))*N(x,mu;y,nu)
    qacc AdjointColorMatrix m_mat_site_no_comm
    (
     const AdjointColorMatrix& n_mat, // n(x,mu;y,nu)
     const Vector<AdjointColorMatrix> ad_x_and_j_n_x, // AdX(x,mu), J(-X(x,mu))
     const Coordinate& xl,
     const int mu,
     const Coordinate& yl,
     const int nu,
     const double epsilon
     )
    {
      const AdjointColorMatrix& ad_x_mat = ad_x_and_j_n_x[0];
      const AdjointColorMatrix& j_n_x_mat = ad_x_and_j_n_x[1];
      if (mu == nu && xl == yl) { // (x,mu) is one of the flowed links
        return make_matrix_exp(ad_x_mat) - epsilon * j_n_x_mat * n_mat; // eq. (38)
      } else {
        return -epsilon * j_n_x_mat * n_mat; // eq. (39)
      }
    }


    // given mu, for all x,y,nu,
    // set M(x,mu;y,nu); M(x,mu;y,nu) = delta_{x,mu;y,nu}e^{AdX(x,mu)} - eps * J(-X(x,mu))*N(x,mu;y,nu)
    inline void set_m_mat_mask_mu_no_comm
    (
     Field<AdjointColorMatrix>& mf, // M(mu) = ( M(x,mu;y,nu) )_{x,y,nu}
     const Field<AdjointColorMatrix>& nf, // N(mu) = ( M(x,mu;y,nu) )_{x,y,nu}
     const FieldM<AdjointColorMatrix, 2>& f_ad_x_and_j_n_x_ext, // ( {AdX(mu), J(-X(mu))} )_x
     const int mask,
     const int mu, // mu: flow link direction
     const double epsilon,
     const FlowType& flow_type
     )
    // mf and nf have the same structure
    // See set_xl_nu_from_mask_mu_yl_m
    // f_ad_x_and_j_n_x_ext need proper communication
    {
      // --- MISC ---
      TIMER("set_m_mat_mask_mu_no_comm");
      qassert(is_initialized(nf));
      qassert(nf.geo().multiplicity == flow_type.multiplicity_reduced_max());
      qassert(is_initialized(f_ad_x_and_j_n_x_ext));
      qassert(0 <= mu && mu < 4);
      const Geometry& geo = nf.geo();
      mf.init(geo, flow_type.multiplicity_reduced_max());
      // ---

      // --- for y ---
      qacc_for(index, geo.local_volume(), {
          const Coordinate yl = geo.coordinate_from_index(index);
          const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,flow_type);
          // ---

          Vector<AdjointColorMatrix> mfv = mf.get_elems(yl); // vec{M}(mu,y) = ( M(x,mu;y,nu) )_{x,nu}
          const Vector<AdjointColorMatrix> nfv = nf.get_elems_const(yl); // vec{N}(mu,y) = ( n(x,mu;y,nu) )_{x,nu}
          qassert(mfv.size() == nfv.size());

          // --- for (x,nu) ---
          for(int m_prime=0; m_prime<flow_type.multiplicity_reduced(mask,y_idx); ++m_prime){
            Coordinate xl;
            int nu;
            set_xl_nu_from_mask_mu_yl_m_prime(xl, nu, mask, mu, yl, m_prime, geo, flow_type);
            // ---

            const Vector<AdjointColorMatrix> ad_x_and_j_n_x = f_ad_x_and_j_n_x_ext.get_elems_const(xl); // AdX(x,mu), J(-X(x,mu))

            // calculate M(x,mu;y,nu) = delta_{x,mu;y,nu}e^{AdX(x,mu)} - eps * J(-X(x,mu))*N(x,mu;y,nu)
            mfv[m_prime] = m_mat_site_no_comm(nfv[m_prime], ad_x_and_j_n_x, xl, mu, yl, nu, epsilon);
          }
        });
    }

    qacc AdjointColorMatrix mp_mat_site_no_comm
    (
     const AdjointColorMatrix& n_mat, const FieldM<ColorMatrix, 1>& cf,
     const GaugeField& gf, const Coordinate& xl, const int mu,
     const double epsilon, const ColorMatrixConstants& cmcs)
    {
      const ColorMatrix& u = gf.get_elem(xl, mu);
      const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
      const ColorMatrix x_mat = static_cast<Complex>(-epsilon) * make_tr_less_anti_herm_matrix(u * c_dagger);

      const AdjointColorMatrix j_x_mat = make_diff_exp_map(x_mat, cmcs);

      AdjointColorMatrix m_mat;
      set_unit(m_mat);
      m_mat -= epsilon * j_x_mat * n_mat;
      return m_mat;
    }


    inline void set_mp_mat_mask_mu_no_comm
    (
     FieldM<AdjointColorMatrix, 1>& mpf,
     const Field<AdjointColorMatrix>& nf,
     const FieldM<ColorMatrix, 1>& cf,
     const GaugeField& gf, const int mask,
     const int mu, const double epsilon,
     const FlowType& flow_type)
    // only use the first elem in each site of nf
    // mu: flow link direction
    {
      TIMER("set_mp_mat_mask_mu_no_comm(mpf)");
      qassert(is_initialized(nf));
      qassert(is_initialized(cf));
      qassert(is_initialized(gf));
      qassert(0 <= mu && mu < 4);
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
      const Geometry geo = geo_reform(gf.geo());
      mpf.init(geo);

      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);
          const AdjointColorMatrix& n_mat = nf.get_elem(xl, 0);
          AdjointColorMatrix& m_mat = mpf.get_elem(xl);
          m_mat = mp_mat_site_no_comm(n_mat, cf, gf, xl, mu, epsilon, cmcs());
        });
    }


    // given y,mu, for all nu,
    // F'(y,nu) = sum_{x,mu}[ F(x,mu) * M(x,mu;y,nu) ] + (1-delta_{y,nu;x,mu})*F(y,nu)
    qacc void set_gm_force_from_flow_site_no_comm
    (
     GaugeMomentum& gm_force, // F' = ( F'(y,nu) )_{y,nu}; to be initialized
     const GaugeMomentum& gm_force_pre_ext, // F = ( F(y,nu) )_{y,nu}
     const Vector<AdjointColorMatrix> m_mat_vec, // vec{M}(mu;y) = ( M(x,mu;y,nu) )_{x,nu}
     const int mask,
     const int mu,
     const FlowType& flow_type,
     const Coordinate& yl,
     const Geometry& geo
     )
    {
      // --- MISC ---
      const Coordinate yg = geo.coordinate_g_from_l(yl);
      const int mask_yl = mask_from_coordinate(yg, flow_type, mu);
      const int y_idx = get_y_idx(yg,mu,flow_type);
      // ---

      Vector<ColorMatrix> fv = gm_force.get_elems(yl); // vec{F'}(y) = ( F'(y,nu) )_nu
      set_zero(fv); // F'(y,nu) = 0

      // --- for (x, nu) [flowed links] --
      for (int m_prime = 0; m_prime<flow_type.multiplicity_reduced(mask,y_idx); ++m_prime){
        Coordinate xl;
        int nu;
        set_xl_nu_from_mask_mu_yl_m_prime(xl, nu, mask, mu, yl, m_prime, geo, flow_type);
        // ---

        const ColorMatrix& force_x_mu = gm_force_pre_ext.get_elem(xl, mu); // F(x,mu)
        const array<double, 8> basis_b = basis_projection_anti_hermitian_matrix(force_x_mu); // F^b(x,mu)

        array<double, 8> basis_a; // to be added to F(y,nu); Denote f(y,nu) = ( f^a(y,nu) )_a
        set_zero(basis_a); // f(y,nu) = 0
        const AdjointColorMatrix& m_mat = m_mat_vec[m_prime]; // M(x,mu;y,nu)

        for(int a=0; a<8; ++a) for(int b=0; b<8; ++b){
            basis_a[a] += basis_b[b] * m_mat(b, a); // f^a(y,nu) = sum_b[ F^b(x,mu) * M(x,mu;y,nu)^{b,a} ]
          }

        fv[nu] += make_anti_hermitian_matrix(basis_a); // F'(y,nu) = sum_a[ T^a * f^a(y,nu) ]
      }

      const Vector<ColorMatrix> f_pre_v = gm_force_pre_ext.get_elems_const(yl); // vec{F}(y) = ( F(y,nu) )_nu
      for (int nu = 0; nu < 4; ++nu) {
        if (!(mask_yl == mask && mu == nu)) { // if (y,nu) is not one of the flowed links
          fv[nu] += f_pre_v[nu]; // F'(y,nu) += F(y,nu)
        }
      }
    }

    // for all (y,nu),
    // F'(y,nu) = sum_{x,mu}[ F(x,mu) * M(x,mu;y,nu) ] + (1-delta_{y,nu;x,mu})*F(y,nu)
    inline void set_gm_force_from_flow_no_comm
    (
     GaugeMomentum& gm_force, // F' = ( F'(y,nu) )_{y,nu}; to be initialized in **_site_no_comm
     const GaugeMomentum& gm_force_pre_ext, // F = ( F(y,nu) )_{y,nu}
     const Field<AdjointColorMatrix>& mf, // M(mu) = ( M(x,mu;y,nu) )_{x,y,nu}
     const int mask,
     const int mu, // mu: flow link direction
     const FlowType& flow_type
     )
    {
      // --- MISC ---
      TIMER("set_gm_force_from_flow_no_comm");
      qassert(is_initialized(gm_force_pre_ext));
      qassert(is_initialized(mf));
      qassert(mf.geo().multiplicity == flow_type.multiplicity_reduced_max());
      qassert(0 <= mu && mu < 4);
      const Geometry geo = geo_reform(gm_force_pre_ext.geo());
      gm_force.init(geo);
      // ---

      // --- for y [general site] --
      qacc_for(index, geo.local_volume(), {
          const Coordinate yl = geo.coordinate_from_index(index);
          //---
          const Vector<AdjointColorMatrix> m_mat_vec = mf.get_elems_const(yl); // vec{M}(mu;y) = ( M(x,mu;y,nu) )_{x,nu}

          // calculate F'(y,nu) = sum_{x,mu}[ F(x,mu) * M(x,mu;y,nu) ] + (1-delta_{y,nu;x,mu})*F(y,nu)
          set_gm_force_from_flow_site_no_comm(gm_force, gm_force_pre_ext, m_mat_vec, mask, mu, flow_type, yl, geo);
        });
    }


    inline void set_f_det_util_mask_mu
    (
     FieldM<array<double, 8>, 1>& f_e2_dj_x_n_mp_inv,
     FieldM<AdjointColorMatrix, 1>& f_n_e_mp_inv_j_x,
     const FieldM<AdjointColorMatrix, 1>& mpf,
     const Field<AdjointColorMatrix>& nf,
     const FieldM<ColorMatrix, 1>& cf,
     const GaugeField& gf,
     const int mask,
     const int mu,
     const double epsilon,
     const FlowType& flow_type
     )
    // mu: flow link direction
    {
      TIMER("set_f_det_util_mask_mu");
      qassert(is_initialized(mpf));
      qassert(is_initialized(nf));
      qassert(is_initialized(cf));
      qassert(is_initialized(gf));
      qassert(0 <= mu && mu < 4);
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
      const Geometry geo = geo_reform(gf.geo());
      f_n_e_mp_inv_j_x.init(geo);
      f_e2_dj_x_n_mp_inv.init(geo);

      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);

          const ColorMatrix& u = gf.get_elem(xl, mu);
          const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
          const ColorMatrix x_mat = static_cast<Complex>(-epsilon) * make_tr_less_anti_herm_matrix(u * c_dagger);
          const AdjointColorMatrix j_x_mat = make_diff_exp_map(x_mat, cmcs());
          const AdjointColorMatrix& n_mat = nf.get_elem(xl, 0);
          const AdjointColorMatrix mp_inv_mat = matrix_inverse( mpf.get_elem(xl) );
          const AdjointColorMatrix e2_n_mp_inv_mat = sqr(epsilon) * n_mat * mp_inv_mat;

          array<double, 8>& basis = f_e2_dj_x_n_mp_inv.get_elem(xl);
          for (int e = 0; e < 8; ++e) {
            basis[e] = matrix_trace(make_diff_exp_map_diff(x_mat, e, cmcs()), e2_n_mp_inv_mat).real();
          }
          f_n_e_mp_inv_j_x.get_elem(xl) = (-epsilon) * mp_inv_mat * j_x_mat;
        });
    }


    inline void set_gm_force_from_flow_det_no_comm_linear
    (
     GaugeMomentum& gm_force_det,
     const FieldM<array<double, 8>, 1>& f_e2_dj_x_n_mp_inv_ext,
     const FieldM<AdjointColorMatrix, 1>& f_n_e_mp_inv_j_x_ext,
     const Field<AdjointColorMatrix>& nf,
     const Field<array<ColorMatrix, 2> >& ducf, const int mask, const int mu,
     const FlowType& flow_type
     )
    // See set_xl_nu_from_mask_mu_yl_m
    // mask: for flow_size=0,1 : flow 1:odd / 2:even site
    // mu: flow link direction
    //
    // given (x,mu), (y,nu), multiple contributions from kernel loops are summed over in the reduction loop below.
    {
      TIMER("set_gm_force_from_flow_det_no_comm_linear");
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
      qassert(is_initialized(f_e2_dj_x_n_mp_inv_ext));
      qassert(is_initialized(f_n_e_mp_inv_j_x_ext));
      qassert(is_initialized(nf));
      qassert(nf.geo().multiplicity == flow_type.multiplicity_reduced_max());
      qassert(is_initialized(ducf));
      qassert(ducf.geo().multiplicity == flow_type.multiplicity_max());
      qassert(0 <= mu && mu < 4);

      const Geometry geo = geo_reform(nf.geo());
      gm_force_det.init();
      gm_force_det.init(geo);

      qacc_for(index, geo.local_volume(), {
          const Coordinate yl = geo.coordinate_from_index(index);
          const array<ColorMatrix, 8>& ts = cmcs().ts;
          const Vector<AdjointColorMatrix> nfv = nf.get_elems_const(yl);
          const Vector<array<ColorMatrix, 2> > ducfv = ducf.get_elems_const(yl);
          Vector<ColorMatrix> gm_f_v = gm_force_det.get_elems(yl);
          set_zero(gm_f_v);

          const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,flow_type);
          const auto beg = flow_type.begin()+flow_type.begin_mask_yidx(mask,y_idx);
          auto itr = beg; // iterator used in while below
          int m_prime = itr->m_prime();
          int next = itr->next();

          while(std::distance(beg,itr)<flow_type.multiplicity(mask,y_idx)) {
            Coordinate xl;
            int nu;
            set_xl_nu_from_mask_mu_yl_m_prime(xl,nu,mask,mu,yl,m_prime,geo,flow_type);
            const array<double, 8>& e2_dj_x_n_mp_inv = f_e2_dj_x_n_mp_inv_ext.get_elem(xl);
            const AdjointColorMatrix& n_e_mp_inv_j_x_mat = f_n_e_mp_inv_j_x_ext.get_elem(xl);
            const AdjointColorMatrix& n_mat = nfv[m_prime];

            array<double, 8> f_det_basis;
            set_zero(f_det_basis);

            for (int a = 0; a < 8; ++a) {
              for (int e = 0; e < 8; ++e) f_det_basis[a] += n_mat(e, a) * e2_dj_x_n_mp_inv[e];

              ColorMatrix uc;
              set_zero(uc);

              for(int i=0; i<next; ++i){ // reduction for a fixed (y,nu)
                const array<ColorMatrix, 2>& uc_mats = ducfv[std::distance(beg,itr)+i];
                const ColorMatrix& uc_pre = uc_mats[0];
                const ColorMatrix& uc_post = uc_mats[1];
                uc += uc_pre * ts[a] * uc_post;
              }

              for (int c = 0; c < 8; ++c) {
                const ColorMatrix d_n = make_tr_less_anti_herm_matrix(ts[c] * uc); // NEED CHANGE !!!!!
                const array<double, 8> d_n_b = basis_projection_anti_hermitian_matrix(d_n);
                for (int b = 0; b < 8; ++b) f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b];
              }
            } // end for a

            const ColorMatrix f_det = static_cast<Complex>(0.5) * make_anti_hermitian_matrix(f_det_basis);
            gm_f_v[nu] += f_det;
            //
            itr += next;
            m_prime = itr->m_prime();
            next = itr->next();
          } // end while
        }); // end for index
    }


    // given (x,mu), (y,nu), multiple contributions from kernel loops are summed over in the reduction loop below.
    inline void set_gm_force_from_flow_det_no_comm_nonlinear
    (
     GaugeMomentum& gm_force_det,
     const FieldM<array<double, 8>, 1>& f_e2_dj_x_n_mp_inv_ext,
     const FieldM<AdjointColorMatrix, 1>& f_n_e_mp_inv_j_x_ext,
     const Field<AdjointColorMatrix>& nf,
     const GaugeField& gf_ext,
     const int mask,
     const int mu, // mu: flow link direction
     const FlowType& flow_type
     )
    {
      TIMER("set_gm_force_from_flow_det_no_comm_nonlinear");
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
      qassert(is_initialized(f_e2_dj_x_n_mp_inv_ext));
      qassert(is_initialized(f_n_e_mp_inv_j_x_ext));
      qassert(is_initialized(nf));
      qassert(nf.geo().multiplicity == flow_type.multiplicity_reduced_max());
      qassert(is_initialized(gf_ext));
      qassert(0 <= mu && mu < 4);

      const Geometry geo = geo_reform(nf.geo());
      gm_force_det.init();
      gm_force_det.init(geo);

      qacc_for(index, geo.local_volume(), {
          const Coordinate yl = geo.coordinate_from_index(index);
          const array<ColorMatrix, 8>& ts = cmcs().ts;
          const Vector<AdjointColorMatrix> nfv = nf.get_elems_const(yl);
          Vector<ColorMatrix> gm_f_v = gm_force_det.get_elems(yl);
          set_zero(gm_f_v);

          const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,flow_type);
          const auto beg = flow_type.begin()+flow_type.begin_mask_yidx(mask,y_idx);
          auto itr_x = beg; // iterator of the loops for a given y; used in while below
          int m_prime = itr_x->m_prime(); // label of x which the iterator corresponds to
          int next = itr_x->next();

          // run over x, the flowed links (mu is already fixed), for a given y
          while(std::distance(beg,itr_x)<flow_type.multiplicity(mask,y_idx)) {
            Coordinate xl;
            int nu;
            set_xl_nu_from_mask_mu_yl_m_prime(xl,nu,mask,mu,yl,m_prime,geo,flow_type);
            const array<double, 8>& e2_dj_x_n_mp_inv = f_e2_dj_x_n_mp_inv_ext.get_elem(xl);
            const AdjointColorMatrix& n_e_mp_inv_j_x_mat = f_n_e_mp_inv_j_x_ext.get_elem(xl);
            const AdjointColorMatrix& n_mat = nfv[m_prime];

            array<double, 8> f_det_basis;
            set_zero(f_det_basis);

            for (int a = 0; a < 8; ++a) {

              for (int e = 0; e < 8; ++e) f_det_basis[a] += n_mat(e, a) * e2_dj_x_n_mp_inv[e];

              for(auto itr_loop = itr_x; itr_loop!=itr_x+next; ++itr_loop){

                const int x_idx = itr_loop->x_idx();
                assert(x_idx == itr_x->x_idx());
                const bool is_pos_dir_y = itr_loop->is_positive();
                const std::vector<int> shape = flow_type.shape( itr_loop->l() );
                const int link_idx_y = itr_loop->s(); // the link idx s (for a given shape l) at which y is located; T^a insertion

                Coordinate x_prime = xl; // search for the case x_prime == x by shifting x_prime along the loop l
                for(std::size_t s_prime=0; s_prime<shape.size(); ++s_prime){
                  Coordinate x_prime_next = coordinate_shifts(x_prime, flow_type.rotate(shape[s_prime],mu));

                  bool is_pos_dir_another_x;
                  if(x_prime==xl && shape[s_prime]==0){ // if the link is (x,mu) (positive directed)
                    is_pos_dir_another_x = true;
                  }
                  else if(x_prime_next==xl && shape[s_prime]==-1){ // if the link is (x,mu) (negative directed)
                    is_pos_dir_another_x = false;
                  }
                  else{
                    x_prime = x_prime_next;
                    continue; // s_prime loop
                  }

                  const int link_idx_another_x = s_prime; // T^c insertion
                  const bool is_ordering_ca = (link_idx_another_x<link_idx_y
                                               ) || ( (link_idx_another_x==link_idx_y)
                                                      && is_pos_dir_y ); // true if "pre T^c mid T^a post"
                  const int link_m_idx = is_ordering_ca? link_idx_another_x : link_idx_y; // smaller link idx
                  const int link_M_idx = is_ordering_ca? link_idx_y : link_idx_another_x; // larger link idx
                  const bool is_link_m_pos_dir = is_ordering_ca? is_pos_dir_another_x : is_pos_dir_y; // ordering of the m link
                  const bool is_link_M_pos_dir = is_ordering_ca? is_pos_dir_y : is_pos_dir_another_x; // ordering of the M link

                  Coordinate xl_tmp = xl; // shifted in the inner loop

                  ColorMatrix uc_pre;
                  set_unit(uc_pre);
                  for(int s=0; s<link_m_idx + 1-is_link_m_pos_dir; ++s) {
                    const int dir = flow_type.rotate(shape[s],mu);
                    uc_pre *= gf_get_link(gf_ext, xl_tmp, dir);
                    xl_tmp = coordinate_shifts(xl_tmp, dir);
                  }

                  ColorMatrix uc_mid;
                  set_unit(uc_mid);
                  for(int s=link_m_idx + 1-is_link_m_pos_dir; s<link_M_idx + 1-is_link_M_pos_dir; ++s) {
                    const int dir = flow_type.rotate(shape[s],mu);
                    uc_mid *= gf_get_link(gf_ext, xl_tmp, dir);
                    xl_tmp = coordinate_shifts(xl_tmp, dir);
                  }

                  ColorMatrix uc_post;
                  set_unit(uc_post);
                  for(int s=link_M_idx + 1-is_link_M_pos_dir; s<static_cast<int>(shape.size()); ++s) {
                    const int dir = flow_type.rotate(shape[s],mu);
                    uc_post *= gf_get_link(gf_ext, xl_tmp, dir);
                    xl_tmp = coordinate_shifts(xl_tmp, dir);
                  }

                  if(!is_link_m_pos_dir) uc_pre = -uc_pre;
                  if(!is_link_M_pos_dir) uc_mid = -uc_mid;

                  for (int c = 0; c < 8; ++c) {
                    const ColorMatrix uc = is_ordering_ca?
                      uc_pre * ts[c] * uc_mid * ts[a] * uc_post:
                      uc_pre * ts[a] * uc_mid * ts[c] * uc_post;
                    const ColorMatrix d_n = make_tr_less_anti_herm_matrix(uc);
                    const array<double, 8> d_n_b = basis_projection_anti_hermitian_matrix(d_n);
                    for (int b = 0; b < 8; ++b) f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b]; // make this outside?
                  } // end for c

                  x_prime = x_prime_next;
                } // end for s_prime

              } // end for itr_loop
            } // end for a

            const ColorMatrix f_det = static_cast<Complex>(0.5) * make_anti_hermitian_matrix(f_det_basis);
            gm_f_v[nu] += f_det;
            //
            itr_x += next;
            m_prime = itr_x->m_prime();
            next = itr_x->next();
          } // end while
        }); // end for index
    }


    // loops are assumed to be linear in each loop
    inline void set_gm_force_from_flow_det_no_comm_distinct
    (
     GaugeMomentum& gm_force_det,
     const FieldM<array<double, 8>, 1>& f_e2_dj_x_n_mp_inv_ext,
     const FieldM<AdjointColorMatrix, 1>& f_n_e_mp_inv_j_x_ext,
     const Field<AdjointColorMatrix>& nf,
     const GaugeField& gf_ext,
     const Field<array<ColorMatrix, 2> >& ducf,
     const int mask,
     const int mu,
     const FlowType& flow_type
     )
    // See set_xl_nu_from_mask_mu_yl_m
    // mask: for flow_size=0,1 : flow 1:odd / 2:even site
    // mu: flow link direction
    //
    // given (x,mu), (y,nu), multiple contributions from kernel loops are summed over in the reduction loop below.
    {
      TIMER("set_gm_force_from_flow_det_no_comm_distinct");
      const box<ColorMatrixConstants>& cmcs = ColorMatrixConstants::get_instance_box();
      qassert(is_initialized(f_e2_dj_x_n_mp_inv_ext));
      qassert(is_initialized(f_n_e_mp_inv_j_x_ext));
      qassert(is_initialized(nf));
      qassert(nf.geo().multiplicity == flow_type.multiplicity_reduced_max());
      qassert(is_initialized(ducf));
      qassert(ducf.geo().multiplicity == flow_type.multiplicity_max());
      qassert(0 <= mu && mu < 4);
      qassert(flow_type.number_of_multiplied_loops()!=0); // supposed to be called for the multiple-loop cases

      const Geometry geo = geo_reform(nf.geo());
      gm_force_det.init();
      gm_force_det.init(geo);

      qacc_for(index, geo.local_volume(), {
          const Coordinate yl = geo.coordinate_from_index(index);
          const array<ColorMatrix, 8>& ts = cmcs().ts;
          const Vector<AdjointColorMatrix> nfv = nf.get_elems_const(yl);
          const Vector<array<ColorMatrix, 2> > ducfv = ducf.get_elems_const(yl);
          Vector<ColorMatrix> gm_f_v = gm_force_det.get_elems(yl);
          set_zero(gm_f_v);

          const int y_idx = get_y_idx(geo.coordinate_g_from_l(yl),mu,flow_type);
          const auto beg = flow_type.begin()+flow_type.begin_mask_yidx(mask,y_idx);
          auto itr_x = beg; // iterator of the loops for a given y; used in while below
          int m_prime = itr_x->m_prime();
          int next = itr_x->next();

          // run over x, the flowed links (mu is already fixed), for a given y
          while(std::distance(beg,itr_x)<flow_type.multiplicity(mask,y_idx)) {
            Coordinate xl;
            int nu;
            set_xl_nu_from_mask_mu_yl_m_prime(xl,nu,mask,mu,yl,m_prime,geo,flow_type);
            const array<double, 8>& e2_dj_x_n_mp_inv = f_e2_dj_x_n_mp_inv_ext.get_elem(xl);
            const AdjointColorMatrix& n_e_mp_inv_j_x_mat = f_n_e_mp_inv_j_x_ext.get_elem(xl);
            const AdjointColorMatrix& n_mat = nfv[m_prime];

            array<double, 8> f_det_basis;
            set_zero(f_det_basis);

            for (int a = 0; a < 8; ++a) {

              for (int e = 0; e < 8; ++e) f_det_basis[a] += n_mat(e, a) * e2_dj_x_n_mp_inv[e];

              for(auto itr_loop = itr_x; itr_loop!=itr_x+next; ++itr_loop){
                const int which_loop = itr_loop->which_loop();

                if( which_loop==-1){ // d^a_{y,nu} acting on the main loop
                  const array<ColorMatrix, 2>& uc_mats = ducfv[std::distance(beg,itr_loop)];
                  const ColorMatrix& uc_pre = uc_mats[0];
                  const ColorMatrix& uc_post = uc_mats[1];

                  // main loop matrix
                  const ColorMatrix dw_a = uc_pre * ts[a] * uc_post;

                  // sub loop factor for the term where both d^a_{y,nu} and d^c_{x,mu} acting on the main loop
                  std::complex<double> factor_A = 0.0;

                  // dA matrix for the term where d^a_{y,nu} acting on the main loop, but d^c_{x,mu} on a sub loop
                  ColorMatrix dA;
                  set_zero(dA);

                  for(int k=0; k<flow_type.number_of_multiplied_loops(); ++k){

                    Coordinate xl_beg;
                    std::vector<int> shape_tmp_rotated(flow_type.length_of_loop);
                    set_xl_beg_and_shape_tmp_rotated(xl_beg,shape_tmp_rotated,xl,mu,flow_type,itr_loop->l(),k);
                    factor_A += matrix_trace( gf_wilson_line_no_comm( gf_ext, xl_beg, shape_tmp_rotated ) ); // reduction (for k)

                    if(flow_type.is_flowed_link_nonlinear){ // we then have nonzero dA
                      int link_idx;
                      bool is_pos_dir_A;
                      search_link_idx_corresp_to_flowed_link(link_idx,is_pos_dir_A,xl_beg,shape_tmp_rotated,xl,mu,flow_type);

                      // calculate the Wilson lines after the T^c insertion
                      Coordinate xl_tmp = xl_beg; // shifted in the inner loop; reset

                      ColorMatrix uc_pre_A;
                      set_unit(uc_pre_A);
                      for(int s=0; s<link_idx + 1-is_pos_dir_A; ++s) {
                        const int dir = shape_tmp_rotated[s];
                        uc_pre_A *= gf_get_link(gf_ext, xl_tmp, dir);
                        xl_tmp = coordinate_shifts(xl_tmp, dir);
                      }

                      ColorMatrix uc_post_A;
                      set_unit(uc_post_A);
                      for(int s=link_idx + 1-is_pos_dir_A; s<flow_type.length_of_loop; ++s) {
                        const int dir = shape_tmp_rotated[s];
                        uc_post_A *= gf_get_link(gf_ext, xl_tmp, dir);
                        xl_tmp = coordinate_shifts(xl_tmp, dir);
                      }

                      if(!is_pos_dir_A) uc_pre_A = -uc_pre_A;
                      dA += uc_post_A * uc_pre_A; // reduction (for k)

                    } // end if flowed_link_nonlinear==true

                  } // end for k (multiplied loops)

                  for (int c = 0; c < 8; ++c) {
                    const std::complex<double> tr_dA_c = flow_type.is_flowed_link_nonlinear? matrix_trace(ts[c]*dA) : 0.0;
                    const ColorMatrix d_n = make_tr_less_anti_herm_matrix( factor_A * (ts[c]*dw_a) + tr_dA_c * dw_a );
                    const array<double, 8> d_n_b = basis_projection_anti_hermitian_matrix(d_n);
                    for (int b = 0; b < 8; ++b) f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b];
                  }

                } // end if main loop
                else{ // if differentiated on a sub loop
                  // main loop matrix
                  std::vector<int> shape_main(flow_type.length_of_loop);
                  // rotated
                  for(int i=0; i<flow_type.length_of_loop; ++i){
                    shape_main[i] = flow_type.rotate( flow_type.shape(itr_loop->l())[i], mu );
                  }
                  const ColorMatrix mat_w = gf_wilson_line_no_comm( gf_ext, xl, shape_main );

                  // sub loop factors

                  if(!flow_type.is_flowed_link_nonlinear){
                    const array<ColorMatrix, 2>& uc_mats = ducfv[std::distance(beg,itr_loop)];
                    const ColorMatrix& uc_pre = uc_mats[0];
                    const ColorMatrix& uc_post = uc_mats[1];

                    for (int c = 0; c < 8; ++c) {
                      const std::complex<double> tr_dA_a = matrix_trace(uc_pre * ts[a] * uc_post);
                      const ColorMatrix d_n = make_tr_less_anti_herm_matrix( tr_dA_a * (ts[c]*mat_w));
                      const array<double, 8> d_n_b = basis_projection_anti_hermitian_matrix(d_n);
                      for (int b = 0; b < 8; ++b) f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b];
                    }
                  }
                  else{ // we then have nonzero ddA
                    Coordinate xl_beg;
                    std::vector<int> shape_tmp_rotated(flow_type.length_of_loop);
                    set_xl_beg_and_shape_tmp_rotated(xl_beg,shape_tmp_rotated,xl,mu,flow_type,itr_loop->l(),which_loop);

                    // search for the case x_prime == x by shifting xl_tmp along the loop
                    int link_idx_x; // T^c insertion
                    bool is_pos_dir_x = true;
                    search_link_idx_corresp_to_flowed_link(link_idx_x,is_pos_dir_x,xl_beg,shape_tmp_rotated,xl,mu,flow_type);

                    const bool is_pos_dir_y = itr_loop->is_positive();
                    const int link_idx_y = itr_loop->s(); // the link idx at which y is located; T^a insertion

                    const bool is_ordering_ca = (link_idx_x<link_idx_y) || ( (link_idx_x==link_idx_y)
                                                                             && is_pos_dir_y ); // true if "pre T^c mid T^a post"
                    const int link_m_idx = is_ordering_ca? link_idx_x : link_idx_y; // smaller link idx
                    const int link_M_idx = is_ordering_ca? link_idx_y : link_idx_x; // larger link idx
                    const bool is_link_m_pos_dir = is_ordering_ca? is_pos_dir_x : is_pos_dir_y; // ordering of the m link
                    const bool is_link_M_pos_dir = is_ordering_ca? is_pos_dir_y : is_pos_dir_x; // ordering of the M link

                    Coordinate xl_tmp = xl_beg; // shifted in the inner loop

                    ColorMatrix uc_pre;
                    set_unit(uc_pre);
                    for(int s=0; s<link_m_idx + 1-is_link_m_pos_dir; ++s) {
                      const int dir = shape_tmp_rotated[s];
                      uc_pre *= gf_get_link(gf_ext, xl_tmp, dir);
                      xl_tmp = coordinate_shifts(xl_tmp, dir);
                    }

                    ColorMatrix uc_mid;
                    set_unit(uc_mid);
                    for(int s=link_m_idx + 1-is_link_m_pos_dir; s<link_M_idx + 1-is_link_M_pos_dir; ++s) {
                      const int dir = shape_tmp_rotated[s];
                      uc_mid *= gf_get_link(gf_ext, xl_tmp, dir);
                      xl_tmp = coordinate_shifts(xl_tmp, dir);
                    }

                    ColorMatrix uc_post;
                    set_unit(uc_post);
                    for(int s=link_M_idx + 1-is_link_M_pos_dir; s<flow_type.length_of_loop; ++s) {
                      const int dir = shape_tmp_rotated[s];
                      uc_post *= gf_get_link(gf_ext, xl_tmp, dir);
                      xl_tmp = coordinate_shifts(xl_tmp, dir);
                    }

                    for (int c = 0; c < 8; ++c) {
                      if(is_ordering_ca){
                        ColorMatrix ta_post = ts[a] * uc_post;
                        if(!is_link_M_pos_dir) ta_post *= -1.0;
                        const std::complex<double> tr_dA_a = matrix_trace(uc_pre * uc_mid * ta_post);
                        std::complex<double> tr_ddA = matrix_trace(uc_pre * ts[c] * uc_mid * ta_post);
                        if(!is_link_m_pos_dir) tr_ddA *= -1.0;
                        const ColorMatrix d_n = make_tr_less_anti_herm_matrix( tr_dA_a * (ts[c]*mat_w) + tr_ddA * mat_w );
                        const array<double, 8> d_n_b = basis_projection_anti_hermitian_matrix(d_n);
                        for (int b = 0; b < 8; ++b) f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b];
                      }
                      else{
                        ColorMatrix pre_ta = uc_pre * ts[a];
                        if(!is_link_m_pos_dir) pre_ta *= -1.0;
                        const std::complex<double> tr_dA_a = matrix_trace(pre_ta * uc_mid * uc_post);
                        std::complex<double> tr_ddA = matrix_trace(pre_ta * uc_mid * ts[c] * uc_post);
                        if(!is_link_M_pos_dir) tr_ddA *= -1.0;
                        const ColorMatrix d_n = make_tr_less_anti_herm_matrix( tr_dA_a * (ts[c]*mat_w) + tr_ddA * mat_w );
                        const array<double, 8> d_n_b = basis_projection_anti_hermitian_matrix(d_n);
                        for (int b = 0; b < 8; ++b) f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b];
                      } // end else
                    } // end for c
                  } // ene else (flow_type.is_flowed_link_nonlinear == true)

                } // end else (differentiated on a sub loop)

              } // end for itr_loop
            } // end for a

            const ColorMatrix f_det = static_cast<Complex>(0.5) * make_anti_hermitian_matrix(f_det_basis);
            gm_f_v[nu] += f_det;
            //
            itr_x += next;
            m_prime = itr_x->m_prime();
            next = itr_x->next();
          } // end while
        }); // end for index
    }



    inline double mf_ln_det_sum(const Field<AdjointColorMatrix>& mpf, const int mask,
                                const FlowType& flow_type, const int mu)
    {
      TIMER("mf_ln_det_sum");
      const Geometry& geo = mpf.geo();
      FieldM<double, 1> f_ln_det;
      f_ln_det.init(geo);
      set_zero(f_ln_det);

      const vector<long>& flowed_indices = get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_type, mu);
      qacc_for(idx, flowed_indices.size(), {
          const Coordinate xl = geo.coordinate_from_index(flowed_indices[idx]);
          const AdjointColorMatrix& m_mat = mpf.get_elem(xl, 0);
          f_ln_det.get_elem(xl) = std::log(matrix_determinant(m_mat));
        });
      double sum = 0.0;
      for (long index = 0; index < geo.local_volume(); ++index) {
        sum += f_ln_det.get_elem(index);
      }
      return sum;
    }


    inline double gf_flow_and_ln_det_node
    (
     GaugeField& gf,
     const GaugeField& gf0,
     const FlowInfo2& fi
     )
    // Return ln(det(d gf / d gf0)) of the flow (on this node ONLY).
    // And set gf to be the flowed gauge field from gf0.
    // Normally call this function to compute the hamilton for this node.
    // Note the glb_sum is NOT perform for the return value ln_det_node.
    // The ln_det_node returned is to be subtracted from the orignal gauge action.
    {
      TIMER("gf_flow_ln_det_node");
      const int max_flow_size = fi.get_max_flow_size();
      qassert(max_flow_size!=0);
      const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Geometry geo = geo_reform(gf0.geo());
      const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
      GaugeField gf_ext;
      gf_ext.init(geo_ext);
      gf_ext = gf0;

      double ln_det_node = 0.0;
      for(int i=0; i<fi.size(); ++i) {
        const FlowStepInfo2& fsi = fi[i];
        const int mask = fsi.mask;
        const int mu = fsi.mu;
        const double epsilon = fsi.epsilon;
        const FlowType& flow_type = *(fsi.ptr_flow_type);

        refresh_expanded_gf_flow_kernel_mask_mu(gf_ext, mask, mu, flow_type);
        FieldM<ColorMatrix, 1> cf;
        set_flow_staple_mask_mu_no_comm(cf, gf_ext, mask, mu, flow_type);

        FieldM<AdjointColorMatrix, 1> nf;
        set_n_mat_mask_mu_no_comm_diag(nf, gf_ext, mask, mu, flow_type); // diag part

        FieldM<AdjointColorMatrix, 1> mpf;
        set_mp_mat_mask_mu_no_comm(mpf, nf, cf, gf_ext, mask, mu, epsilon, flow_type);
        ln_det_node += mf_ln_det_sum(mpf, mask, flow_type, mu);

        gf_flow_mask_mu_no_comm(gf_ext, gf_ext, mask, mu, epsilon, flow_type);
      }
      gf.init(geo);
      gf = gf_ext;
      return ln_det_node;
    }


    inline double gf_hamilton_flowed_node(const GaugeField& gf0,
                                          const GaugeAction& ga,
                                          const FlowInfo2& fi)
    {
      TIMER("gf_hamilton_flowed_node");
      GaugeField gf;
      const double ln_det_node = gf_flow_and_ln_det_node(gf, gf0, fi);
      return gf_hamilton_node(gf, ga) - ln_det_node;
    }


    inline void gf_hamilton_flowed(double& flowed_action, double& ln_det,
                                   const GaugeField& gf0, const GaugeAction& ga,
                                   const FlowInfo2& fi)
    // effective hamilton for U_0 = flowed_action - ln_det
    {
      TIMER("gf_hamilton_flowed");
      GaugeField gf;
      ln_det = gf_flow_and_ln_det_node(gf, gf0, fi);
      flowed_action = gf_hamilton_node(gf, ga);
      glb_sum(ln_det);
      glb_sum(flowed_action);
    }

    inline void set_gm_force_propagated_from_flow_step
    (
     GaugeMomentum& gm_force_propagated, const GaugeMomentum& gm_force_pre,
     const GaugeField& gf0_ext, const FlowStepInfo2& fsi
     )
    // The gf0_ext must be extended and refreshed and not flowed
    // (for this flow step)
    {
      TIMER("set_gm_force_propagated_from_flow_step");
      qassert(is_initialized(gm_force_pre));
      qassert(is_initialized(gf0_ext));
      const int flow_size = fsi.ptr_flow_type->flow_size;
      const Coordinate expand_left(flow_size,flow_size,flow_size,flow_size);
      const Coordinate expand_right(flow_size,flow_size,flow_size,flow_size);
      const Geometry geo = geo_reform(gf0_ext.geo(), 1);
      const Geometry geo_ext = geo_reform(geo, 1, expand_left, expand_right);
      GaugeMomentum gm_force_pre_ext;
      gm_force_pre_ext.init(geo_ext);
      gm_force_pre_ext = gm_force_pre;
      refresh_expanded(gm_force_pre_ext);
      // clear the output forces
      gm_force_propagated.init();

      const int mask = fsi.mask;
      const int mu = fsi.mu;
      const double epsilon = fsi.epsilon;
      const FlowType& flow_type = *(fsi.ptr_flow_type);

      // To propagate force

      // N
      Field<AdjointColorMatrix> nf;
      Field<array<ColorMatrix, 2> > dwf;
      set_dw_mat_mask_mu_no_comm(dwf, gf0_ext, mask, mu, flow_type);
      set_n_mat_mask_mu_no_comm(nf, dwf, gf0_ext, mask, flow_type, mu);

      //
      FieldM<AdjointColorMatrix, 2> f_ad_x_and_j_n_x_ext;
      f_ad_x_and_j_n_x_ext.init(geo_ext);
      FieldM<ColorMatrix, 1> cf;
      set_flow_staple_mask_mu_no_comm(cf, gf0_ext, mask, mu, flow_type);
      set_ad_x_and_j_n_x_mask_mu_no_comm(f_ad_x_and_j_n_x_ext, cf, gf0_ext, mask, mu, epsilon, flow_type);
      refresh_expanded(f_ad_x_and_j_n_x_ext);

      // M(mu) = ( M(x,mu;y,nu) )_{x,y,nu}
      Field<AdjointColorMatrix> mf;
      set_m_mat_mask_mu_no_comm(mf, nf, f_ad_x_and_j_n_x_ext, mask, mu, epsilon, flow_type);

      // F'(y,nu) = sum_{x,mu}[ F(x,mu) * M(x,mu;y,nu) ] + (1-delta_{y,nu;x,mu})*F(y,nu)
      set_gm_force_from_flow_no_comm(gm_force_propagated, gm_force_pre_ext, mf, mask, mu, flow_type);
    }


    inline void set_gm_force_propagated_and_gm_force_det_from_flow_step
    (
     GaugeMomentum& gm_force_propagated, GaugeMomentum& gm_force_det,
     const GaugeMomentum& gm_force_pre, const GaugeField& gf0_ext,
     const FlowStepInfo2& fsi
     )
    // The gm_force_propagated and gm_force_det need to be summed together.
    // The gf0_ext must be extended and refreshed and not flowed
    // (for this flow step)
    {
      TIMER("set_gm_force_propagated_and_gm_force_det_from_flow_step");
      qassert(is_initialized(gm_force_pre));
      qassert(is_initialized(gf0_ext));
      const int flow_size = fsi.ptr_flow_type->flow_size;
      const Coordinate expand_left(flow_size,flow_size,flow_size,flow_size);
      const Coordinate expand_right(flow_size,flow_size,flow_size,flow_size);
      const Geometry geo = geo_reform(gf0_ext.geo(), 1);
      const Geometry geo_ext = geo_reform(geo, 1, expand_left, expand_right);
      GaugeMomentum gm_force_pre_ext;
      gm_force_pre_ext.init(geo_ext);
      gm_force_pre_ext = gm_force_pre;
      refresh_expanded(gm_force_pre_ext);

      // clear the output forces
      gm_force_det.init();
      gm_force_propagated.init();

      const int mask = fsi.mask;
      const int mu = fsi.mu;
      const double epsilon = fsi.epsilon;
      const FlowType& flow_type = *(fsi.ptr_flow_type);

      // To propagate force
      FieldM<ColorMatrix, 1> cf;
      set_flow_staple_mask_mu_no_comm(cf, gf0_ext, mask, mu, flow_type);
      Field<array<ColorMatrix, 2> > dwf;
      set_dw_mat_mask_mu_no_comm(dwf, gf0_ext, mask, mu, flow_type);
      Field<AdjointColorMatrix> nf;
      set_n_mat_mask_mu_no_comm(nf, dwf, gf0_ext, mask, flow_type, mu);
      FieldM<AdjointColorMatrix, 2> f_ad_x_and_j_n_x_ext;
      f_ad_x_and_j_n_x_ext.init(geo_ext);
      set_ad_x_and_j_n_x_mask_mu_no_comm(f_ad_x_and_j_n_x_ext, cf, gf0_ext,
                                         mask, mu, epsilon, flow_type);

      refresh_expanded(f_ad_x_and_j_n_x_ext);
      Field<AdjointColorMatrix> mf;
      set_m_mat_mask_mu_no_comm(mf, nf, f_ad_x_and_j_n_x_ext, mask, mu, epsilon, flow_type);
      set_gm_force_from_flow_no_comm(gm_force_propagated, gm_force_pre_ext, mf,
                                     mask, mu, flow_type);

      // To include force from determinants
      FieldM<AdjointColorMatrix, 1> mpf;
      set_mp_mat_mask_mu_no_comm(mpf, nf, cf, gf0_ext, mask, mu, epsilon, flow_type);
      FieldM<array<double, 8>, 1> f_e2_dj_x_n_mp_inv_ext;
      f_e2_dj_x_n_mp_inv_ext.init(geo_ext);
      FieldM<AdjointColorMatrix, 1> f_n_e_mp_inv_j_x_ext;
      f_n_e_mp_inv_j_x_ext.init(geo_ext);
      set_f_det_util_mask_mu(f_e2_dj_x_n_mp_inv_ext, f_n_e_mp_inv_j_x_ext, mpf,
                             nf, cf, gf0_ext, mask, mu, epsilon, flow_type);

      refresh_expanded(f_n_e_mp_inv_j_x_ext);
      refresh_expanded(f_e2_dj_x_n_mp_inv_ext);
      if(flow_type.number_of_multiplied_loops()!=0){
        set_gm_force_from_flow_det_no_comm_distinct(gm_force_det, f_e2_dj_x_n_mp_inv_ext,
                                                    f_n_e_mp_inv_j_x_ext, nf, gf0_ext,
                                                    dwf, mask,
                                                    mu, flow_type);
      }
      else if(!flow_type.is_flowed_link_nonlinear){
        set_gm_force_from_flow_det_no_comm_linear(gm_force_det, f_e2_dj_x_n_mp_inv_ext,
                                                  f_n_e_mp_inv_j_x_ext, nf, dwf, mask,
                                                  mu, flow_type);
      }
      else{ // require dducf, which we calculate inside
        set_gm_force_from_flow_det_no_comm_nonlinear(gm_force_det, f_e2_dj_x_n_mp_inv_ext,
                                                     f_n_e_mp_inv_j_x_ext, nf, gf0_ext, mask,
                                                     mu, flow_type);
      }
    }


    inline void set_flowed_gauge_fields(std::vector<GaugeField>& gf_ext_vec,
                                        const GaugeField& gf0,
                                        const FlowInfo2& fi)
    // Fill gf_ext_vec in the flow order (start with gf0)
    // Note gf_ext_vec.size() == fi.size() + 1
    // The gf_ext_vec.back() is the flowed (physical) gauge field.
    // All gauge fields are refreshed except gf_ext_vec.back().
    {
      TIMER("set_flowed_gauge_fields");
      const int n_steps = fi.size();
      clear(gf_ext_vec);
      gf_ext_vec.resize(n_steps + 1);
      const int max_flow_size = fi.get_max_flow_size();
      qassert(max_flow_size!=0);
      const Coordinate expand_left(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Coordinate expand_right(max_flow_size,max_flow_size,max_flow_size,max_flow_size);
      const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
      gf_ext_vec[0].init(geo_ext);
      gf_ext_vec[0] = gf0;
      for(int i=0; i<static_cast<int>(fi.size()); ++i) {
        const FlowStepInfo2& fsi = fi[i];
        refresh_expanded(gf_ext_vec[i]);
        gf_ext_vec[i + 1].init(geo_ext);
        gf_flow_mask_mu_no_comm(gf_ext_vec[i + 1], gf_ext_vec[i], fsi.mask, fsi.mu,
                                fsi.epsilon, *(fsi.ptr_flow_type));
      }
    }


    inline void set_gm_force_propagated_det_from_flow(GaugeMomentum& gm_force,
                                                      const GaugeMomentum& gm_force_pre,
                                                      const std::vector<GaugeField>& gf_ext_vec,
                                                      const FlowInfo2& fi)
    // Propagate the force act on the flowed (physical) gauge field to the unflowed
    // gauge field.
    // Force from the determinant is included.
    // Call set_flowed_gauge_fields(gf_ext_vec, gf0, fi) to obtain gf_ext_vec.
    {
      TIMER("set_gm_force_propagated_det_from_flow");
      const Geometry geo = geo_reform(gm_force_pre.geo());
      gm_force.init(geo);
      gm_force = gm_force_pre;
      for (int i=fi.size()-1; i>=0; --i) {
        GaugeMomentum gm_force_det;
        set_gm_force_propagated_and_gm_force_det_from_flow_step(gm_force, gm_force_det,
                                                                gm_force, gf_ext_vec[i],
                                                                fi[i]);
        gm_force += gm_force_det;
      }
    }


    inline void set_gm_force_propagated_no_det_from_flow
    (
     GaugeMomentum& gm_force,
     const GaugeMomentum& gm_force_pre,
     const std::vector<GaugeField>& gf_ext_vec,
     const FlowInfo2& fi
     )
    // Propagate the force act on the flowed (physical) gauge field to the unflowed gauge field.
    // Force from the determinant is included.
    // Call set_flowed_gauge_fields(gf_ext_vec, gf0, fi) to obtain gf_ext_vec.
    {
      TIMER("set_gm_force_propagated_no_det_from_flow");
      const Geometry geo = geo_reform(gm_force_pre.geo());
      gm_force.init(geo);
      gm_force = gm_force_pre;
      for (int i=fi.size()-1; i>=0; --i) {
        set_gm_force_propagated_from_flow_step(gm_force, gm_force, gf_ext_vec[i],fi[i]);
      }
    }


    inline void set_gm_force_flowed(GaugeMomentum& gm_force,
                                    const GaugeField& gf0,
                                    const GaugeAction& ga,
                                    const FlowInfo2& fi)
    {
      TIMER("set_gm_force_flowed");
      std::vector<GaugeField> gf_ext_vec;
      set_flowed_gauge_fields(gf_ext_vec, gf0, fi);
      gm_force.init();
      set_gm_force(gm_force, gf_ext_vec.back(), ga);
      set_gm_force_propagated_det_from_flow(gm_force, gm_force, gf_ext_vec, fi);
    }


    inline void set_gm_force_flowed_no_det(GaugeMomentum& gm_force,
                                           const GaugeMomentum& gm_force_pre,
                                           const GaugeField& gf0,
                                           const FlowInfo2& fi)
    {
      TIMER("set_gm_force_flowed_no_det");
      std::vector<GaugeField> gf_ext_vec;
      set_flowed_gauge_fields(gf_ext_vec, gf0, fi);
      set_gm_force_propagated_no_det_from_flow(gm_force, gm_force_pre, gf_ext_vec, fi);
    }


    // ------------------------------------------------

    FlowInfo2::FlowInfo2
    (
     const double epsilon /*= 0.0*/,
     const double epsilon_rect /*= 0.0*/,
     const double epsilon_chair /*= 0.0*/,
     const double epsilon_q /*= 0.0*/,
     const double epsilon_plaq_twice /*= 0.0*/,
     const double epsilon_twist_rect /*= 0.0*/,
     const double epsilon_twist_chair /*= 0.0*/,
     const double epsilon_plaq_squared /*=0.0*/,
     const double epsilon_plaq_times_reversed /*=0.0*/,
     const double epsilon_rect_distinct /*=0.0*/,
     const double epsilon_chair_distinct /*=0.0*/,
     const double epsilon_twist_rect_distinct /*=0.0*/,
     const double epsilon_twist_chair_distinct /*=0.0*/
     )
    {
      qassert( GlobalScopeHash::get_flow_type.size()==GlobalScopeHash::n_type );

      if(std::abs(epsilon) > 1.0e-15){
        this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq",epsilon);
      }
      if (std::abs(epsilon_rect) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"srect",epsilon_rect);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect",epsilon_rect);
      }
      if (std::abs(epsilon_chair) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho1",epsilon_chair);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho2",epsilon_chair);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho3",epsilon_chair);
      }
      if (std::abs(epsilon_q) > 1.0e-15) {
        this->set_from_hash_q(GlobalScopeHash::get_flow_type,epsilon_q);
      }
      if (std::abs(epsilon_plaq_twice) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_twice",epsilon_plaq_twice);
      }
      if (std::abs(epsilon_twist_rect) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_srect_rho1",
                            epsilon_twist_rect);
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_srect_rho2",
                            epsilon_twist_rect);
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_srect_rho3",
                            epsilon_twist_rect);
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_lrect",
                            epsilon_twist_rect);
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_mrect",
                            epsilon_twist_rect);
      }
      if (std::abs(epsilon_twist_chair) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_chair_outer_rho1",
                            epsilon_twist_chair);
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_chair_outer_rho2",
                            epsilon_twist_chair);
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_chair_outer_rho3",
                            epsilon_twist_chair);
        this->set_from_hash(GlobalScopeHash::get_flow_type,
                            "twist_chair_inner",
                            epsilon_twist_chair);
      }
      if (std::abs(epsilon_plaq_squared) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_squared",epsilon_plaq_squared);
      }
      if (std::abs(epsilon_plaq_times_reversed) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_times_reversed",epsilon_plaq_times_reversed);
      }
      if (std::abs(epsilon_rect_distinct) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho1",epsilon_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho2",epsilon_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho3",epsilon_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect_distinct",epsilon_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"mrect_distinct",epsilon_rect_distinct);
      }
      if (std::abs(epsilon_chair_distinct) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho1",epsilon_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho2",epsilon_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho3",epsilon_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho1",epsilon_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho2",epsilon_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho3",epsilon_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_inner",epsilon_chair_distinct);
      }
      if (std::abs(epsilon_twist_rect_distinct) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho1",epsilon_twist_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho2",epsilon_twist_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho3",epsilon_twist_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect_distinct",epsilon_twist_rect_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect_distinct",epsilon_twist_rect_distinct);
      }
      if (std::abs(epsilon_twist_chair_distinct) > 1.0e-15) {
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho1",epsilon_twist_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho2",epsilon_twist_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho3",epsilon_twist_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho1",epsilon_twist_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho2",epsilon_twist_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho3",epsilon_twist_chair_distinct);
        this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_inner",epsilon_twist_chair_distinct);
      }
      if(this->size()!=0) this->turn_on_measurement_at_last_step();

      this->set_max_flow_size(GlobalScopeHash::get_flow_type);
    }


    FlowInfo2::FlowInfo2
    (
     const bool is_luscher_flow,
     const double beta,
     const double tmax,
     const int nstep,
     const bool is_first_order /*= true*/
     )
    {
      qassert( GlobalScopeHash::get_flow_type.size()==GlobalScopeHash::n_type );

      const double epsilon_base = tmax/nstep;

      if( is_luscher_flow ){
        for(int i=0; i<nstep; ++i){
          const double epsilon0 = - epsilon_base * (-beta/32.0);
          const double epsilon1 = - epsilon_base * (beta*beta/192.0) * (-4.0/33.0) * (epsilon_base*i);
          const double epsilon2 = - epsilon_base * (beta*beta/192.0) * (12.0/119.0) * (epsilon_base*i);
          const double epsilon3 = - epsilon_base * (beta*beta/192.0) * (1.0/33.0) * (epsilon_base*i);
          const double epsilon4 = - epsilon_base * (beta*beta/192.0) * (-5.0/119.0) * (epsilon_base*i);
          const double epsilon5 = - epsilon_base * (beta*beta/192.0) * (3.0/10.0) * (epsilon_base*i);
          const double epsilon6 = - epsilon_base * (beta*beta/192.0) * (-1.0/5.0) * (epsilon_base*i);
          const double epsilon7 = - epsilon_base * (beta*beta/192.0) * (1.0/9.0) * (epsilon_base*i);
          if (std::abs(epsilon0) > 1.0e-15){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq",epsilon0);
          }
          if (std::abs(epsilon1) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho1",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho2",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho3",epsilon1);
          }
          if (std::abs(epsilon2) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho1",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho2",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho3",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho1",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho2",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho3",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_inner",epsilon2);
          }
          if (std::abs(epsilon3) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho1",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho2",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho3",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect_distinct",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"mrect_distinct",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho1",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho2",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho3",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho1",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho2",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho3",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_inner",epsilon3);
          }
          if (std::abs(epsilon4) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho1",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho2",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho3",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect_distinct",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect_distinct",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho1",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho2",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho3",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho1",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho2",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho3",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_inner",epsilon4);
          }
          if (std::abs(epsilon5) > 1.0e-15 && is_first_order ) this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_twice",epsilon5);
          if (std::abs(epsilon6) > 1.0e-15 && is_first_order ) this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_squared",epsilon6);
          if (std::abs(epsilon7) > 1.0e-15 && is_first_order ) this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_times_reversed",epsilon7);
          this->turn_on_measurement_at_last_step();
        }
      }
      else{ // const double epsilon_base = tmax/nstep;
        for(int i=0; i<nstep; ++i){ // t from 1-tmax to 1 while t = 0,...,tmax
          const double epsilon0 = - epsilon_base * (-beta/32.0);
          const double epsilon1 = - epsilon_base * (beta*beta/192.0) * (-4.0/33.0) * (1.0-tmax+epsilon_base*(i+1));
          const double epsilon2 = - epsilon_base * (beta*beta/192.0) * (12.0/119.0) * (1.0-tmax+epsilon_base*(i+1));
          const double epsilon3 = - epsilon_base * (beta*beta/192.0) * (1.0/33.0) * (1.0-tmax+epsilon_base*(i+1));
          const double epsilon4 = - epsilon_base * (beta*beta/192.0) * (-5.0/119.0) * (1.0-tmax+epsilon_base*(i+1));
          const double epsilon5 = - epsilon_base * (beta*beta/192.0) * (3.0/10.0) * (1.0-tmax+epsilon_base*(i+1));
          const double epsilon6 = - epsilon_base * (beta*beta/192.0) * (-1.0/5.0) * (1.0-tmax+epsilon_base*(i+1));
          const double epsilon7 = - epsilon_base * (beta*beta/192.0) * (1.0/9.0) * (1.0-tmax+epsilon_base*(i+1));
          if (std::abs(epsilon0) > 1.0e-15){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq",epsilon0);
          }
          if (std::abs(epsilon1) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho1",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho2",epsilon1);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho3",epsilon1);
          }
          if (std::abs(epsilon2) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho1",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho2",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho3",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho1",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho2",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho3",epsilon2);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_inner",epsilon2);
          }
          if (std::abs(epsilon3) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho1",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho2",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho3",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect_distinct",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"mrect_distinct",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho1",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho2",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho3",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho1",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho2",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho3",epsilon3);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_inner",epsilon3);
          }
          if (std::abs(epsilon4) > 1.0e-15 && is_first_order ){
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho1",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho2",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho3",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect_distinct",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect_distinct",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho1",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho2",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho3",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho1",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho2",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho3",epsilon4);
            this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_inner",epsilon4);
          }
          if (std::abs(epsilon5) > 1.0e-15 && is_first_order ) this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_twice",epsilon5);
          if (std::abs(epsilon6) > 1.0e-15 && is_first_order ) this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_squared",epsilon6);
          if (std::abs(epsilon7) > 1.0e-15 && is_first_order ) this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_times_reversed",epsilon7);
          this->turn_on_measurement_at_last_step();
        }
      }

      this->set_max_flow_size(GlobalScopeHash::get_flow_type);
    }


    FlowInfo2::FlowInfo2 // theta flow
    (
     const double theta,
     const int nstep
     )
    {
      qassert( GlobalScopeHash::get_flow_type.size()==GlobalScopeHash::n_type );

      const double epsilon_base = 1.0/nstep;

      for(int i=0; i<nstep; ++i){
        const double epsilon_q = - epsilon_base * 2.0*theta/(16.0*32.0*M_PI*M_PI); // 2: c.c.
        if (std::abs(epsilon_q) > 1.0e-15) {
          this->set_from_hash_q(GlobalScopeHash::get_flow_type,epsilon_q);
        }
        this->turn_on_measurement_at_last_step();
      }

      this->set_max_flow_size(GlobalScopeHash::get_flow_type);
    }


    std::vector<double> FlowInfo2::get_gamma_coeffs( const double beta ) const&{
      Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(7,7);
      matrix(0,0) = 16.0/3.0 - 2.0 * beta / 9.0;
      matrix(0,1) = - beta / 6.0;
      matrix(0,2) = beta / 6.0;
      matrix(0,3) = beta / 18.0;;
      matrix(0,4) = - beta / 18.0;;
      matrix(0,5) = 2.0 * beta / 9.0;
      matrix(0,6) = 2.0 * beta / 9.0;
      //
      matrix(1,0) = -5.0 * beta;
      matrix(1,1) = 8.0;
      //
      matrix(2,0) = -25.0 * beta / 3.0;
      matrix(2,2) = 31.0 / 3.0;
      matrix(2,4) = 1.0;
      //
      matrix(3,0) = -20.0 * beta;
      matrix(3,1) = -1.0;
      matrix(3,3) = 11.0;
      //
      matrix(4,0) = -20.0 * beta;
      matrix(4,2) = 1.0;
      matrix(4,4) = 31.0/3.0;
      //
      matrix(5,0) = 8.0 - 2.0 * beta / 3.0;
      matrix(5,5) = 40.0 / 3.0;
      //
      matrix(6,0) = -beta;
      matrix(6,6) = 12.0;
      //

      Eigen::VectorXd vector = Eigen::VectorXd::Zero(7);
      vector(0) = -1.0/6.0;

      Eigen::FullPivLU<Eigen::MatrixXd> lu(matrix);
      Eigen::VectorXd gamma = lu.solve(vector);

      std::vector<double> res(7);
      for(int i=0; i<7; ++i) res[i] = gamma[i];

      return res;
    }


    FlowInfo2::FlowInfo2
    (
     const double beta,
     const double tmax,
     const int nstep
     ){
      qassert( GlobalScopeHash::get_flow_type.size()==GlobalScopeHash::n_type );

      const double epsilon_base = tmax/nstep;

      for(int i=0; i<nstep; ++i){
        const double beta_tmp = beta - i * epsilon_base;
        const std::vector<double> gamma = get_gamma_coeffs(beta_tmp);

        // W0
        double epsilon = -gamma[0]*epsilon_base;
        if (std::abs(epsilon) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq", epsilon);
        }
        // W1
        epsilon = -gamma[1]*epsilon_base;
        if (std::abs(epsilon) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho3", epsilon);
        }
        // W2
        epsilon = -gamma[2]*epsilon_base;
        if (std::abs(epsilon) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect", epsilon);
          //
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_outer_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_inner", epsilon);
        }
        // W3
        epsilon = -gamma[3]*epsilon_base;
        if (std::abs(epsilon) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect_distinct", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"mrect_distinct", epsilon);
          //
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_inner", epsilon);
        }
        // W4
        epsilon = -gamma[4]*epsilon_base;
        if (std::abs(epsilon) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect_distinct", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect_distinct", epsilon);
          //
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho1", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho2", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho3", epsilon);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_inner", epsilon);
        }
        // W5
        epsilon = -gamma[5]*epsilon_base;
        if (std::abs(epsilon) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_twice", epsilon);
        }
        // W7
        epsilon = -gamma[6]*epsilon_base;
        if (std::abs(epsilon) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_times_reversed", epsilon);
        }

        this->turn_on_measurement_at_last_step();
      }

      this->set_max_flow_size(GlobalScopeHash::get_flow_type);
    }



    FlowInfo2::FlowInfo2
    (
     const double tmax,
     const int nstep,
     const std::vector<std::vector<double>>& gamma_list
     )
    {
      qassert( GlobalScopeHash::get_flow_type.size()==GlobalScopeHash::n_type );
      const double epsilon_base = tmax/nstep;
      qassert( static_cast<int>(gamma_list.size())==nstep );
      for(int i=0; i<nstep; ++i ) qassert( gamma_list[i].size()==11 );

      for(int i=0; i<nstep; ++i){
        const double epsilon_plaq = -epsilon_base * gamma_list[i][0];
        const double epsilon_rect = -epsilon_base * gamma_list[i][1];
        const double epsilon_chair = -epsilon_base * gamma_list[i][2];
        const double epsilon_twist_rect = -epsilon_base * gamma_list[i][3]; // bound for double: 1.0 * 1/16
        const double epsilon_twist_chair = -epsilon_base * gamma_list[i][4]; // bound for double: 1.0 * 1/64
        const double epsilon_rect_distinct = -epsilon_base * gamma_list[i][5]; // bound 1.0/(8.0*6.0) ~ 0.021
        const double epsilon_chair_distinct = -epsilon_base * gamma_list[i][6]; // bound 0.5/(8.0*12.0) ~ 0.00521
        const double epsilon_twist_rect_distinct = -epsilon_base * gamma_list[i][7]; // bound 1.0/(8.0*6.0) ~ 0.021
        const double epsilon_twist_chair_distinct = -epsilon_base * gamma_list[i][8]; // bound 0.5/(8.0*12.0) ~ 0.00525
        const double epsilon_plaq_twice = -epsilon_base * gamma_list[i][9]; // bound for double: 0.5 * 1/16 (0.5 comes from the multiclicative factor)
        const double epsilon_plaq_times_reversed = -epsilon_base * gamma_list[i][10]; // bound 1.0/(8.0*6.0) ~ 0.021

        if(std::abs(epsilon_plaq) > 1.0e-15){
          this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq",epsilon_plaq);
        }
        if (std::abs(epsilon_rect) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect",epsilon_rect);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect",epsilon_rect);
        }
        if (std::abs(epsilon_chair) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho1",epsilon_chair);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho2",epsilon_chair);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_rho3",epsilon_chair);
        }
        if (std::abs(epsilon_plaq_twice) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_twice",epsilon_plaq_twice);
        }
        if (std::abs(epsilon_twist_rect) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_srect_rho1",
                              epsilon_twist_rect);
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_srect_rho2",
                              epsilon_twist_rect);
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_srect_rho3",
                              epsilon_twist_rect);
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_lrect",
                              epsilon_twist_rect);
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_mrect",
                              epsilon_twist_rect);
        }
        if (std::abs(epsilon_twist_chair) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_chair_outer_rho1",
                              epsilon_twist_chair);
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_chair_outer_rho2",
                              epsilon_twist_chair);
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_chair_outer_rho3",
                              epsilon_twist_chair);
          this->set_from_hash(GlobalScopeHash::get_flow_type,
                              "twist_chair_inner",
                              epsilon_twist_chair);
        }
        if (std::abs(epsilon_plaq_times_reversed) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"plaq_times_reversed",epsilon_plaq_times_reversed);
        }
        if (std::abs(epsilon_rect_distinct) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho1",epsilon_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho2",epsilon_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"srect_distinct_rho3",epsilon_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"lrect_distinct",epsilon_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"mrect_distinct",epsilon_rect_distinct);
        }
        if (std::abs(epsilon_chair_distinct) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho1",epsilon_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho2",epsilon_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typea_rho3",epsilon_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho1",epsilon_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho2",epsilon_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_typeb_rho3",epsilon_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"chair_distinct_inner",epsilon_chair_distinct);
        }
        if (std::abs(epsilon_twist_rect_distinct) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho1",epsilon_twist_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho2",epsilon_twist_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_srect_distinct_rho3",epsilon_twist_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_lrect_distinct",epsilon_twist_rect_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_mrect_distinct",epsilon_twist_rect_distinct);
        }
        if (std::abs(epsilon_twist_chair_distinct) > 1.0e-15) {
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho1",epsilon_twist_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho2",epsilon_twist_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typea_rho3",epsilon_twist_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho1",epsilon_twist_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho2",epsilon_twist_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_typeb_rho3",epsilon_twist_chair_distinct);
          this->set_from_hash(GlobalScopeHash::get_flow_type,"twist_chair_distinct_inner",epsilon_twist_chair_distinct);
        }
        if(this->size()!=0) this->turn_on_measurement_at_last_step();

        this->set_max_flow_size(GlobalScopeHash::get_flow_type);
      }
    }




    void FlowInfo2::set_max_flow_size
    (
     const std::unordered_map<std::string, FlowType>& hash
     ) &
    {
      for(auto itr = hash.begin(); itr!= hash.end(); ++itr) {
        this->max_flow_size = std::max( itr->second.flow_size, this->max_flow_size );
      }
    }

    void FlowInfo2::set_from_hash
    (
     const std::unordered_map<std::string, FlowType>& hash,
     const std::string& description,
     const double step_size
     ) &
    {
      assert(description!="q_plus" && description !="q_minus");
      const FlowType& flow_type = hash.at(description);
      assert(( flow_type.multiplicative_const*step_size < 1.0 / flow_type.bound_inverse ));
      for(int mu=3; mu>=0; --mu){
        // int rotation_sign = 1;
        for(int mask=0; mask<flow_type.n_masks; ++mask){
          this->v.push_back(FlowStepInfo2(mask, mu,
                                          //rotation_sign*flow_type.multiplicative_const*step_size,
                                          flow_type.multiplicative_const*step_size,
                                          flow_type));
        }
      }
    }

    void FlowInfo2::set_from_hash_q
    (
     const std::unordered_map<std::string, FlowType>& hash,
     const double step_size
     ) &
    {
      const FlowType& q_plus = hash.at("q_plus");
      const FlowType& q_minus = hash.at("q_minus");

      for(int mu=3; mu>=0; --mu){
        const int rotation_sign = 1-(mu%2)*2;
        if(rotation_sign==1){
          for(int mask=0; mask<q_plus.n_masks; ++mask){
            this->v.push_back(FlowStepInfo2(mask, mu,
                                            q_plus.multiplicative_const*step_size,
                                            q_plus));
          }
          for(int mask=0; mask<q_minus.n_masks; ++mask){
            this->v.push_back(FlowStepInfo2(mask, mu,
                                            q_minus.multiplicative_const*step_size,
                                            q_minus));
          }
        }
        else if(rotation_sign==-1){
          for(int mask=0; mask<q_plus.n_masks; ++mask){
            this->v.push_back(FlowStepInfo2(mask, mu,
                                            -q_plus.multiplicative_const*step_size,
                                            q_plus));
          }
          for(int mask=0; mask<q_minus.n_masks; ++mask){
            this->v.push_back(FlowStepInfo2(mask, mu,
                                            -q_minus.multiplicative_const*step_size,
                                            q_minus));
          }
        }
        else assert(false);
      }
    }
  }  // namespace Generalized
}  // namespace qlat

