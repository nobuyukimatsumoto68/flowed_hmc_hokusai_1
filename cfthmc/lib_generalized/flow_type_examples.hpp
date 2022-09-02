// TO DO
// add sign property in the FlowType table
#ifndef FLOW_TYPE_EXAMPLES_HPP
#define FLOW_TYPE_EXAMPLES_HPP

#include <iostream>
#include <iomanip>
#include <array>
#include <string>
#include <iterator>
#include <cassert>
#include <functional>
#include <unordered_map>

#include "unit_info.hpp"
#include "kernel_info.hpp"
#include "mask_info.hpp"
#include "flow_type.hpp"


/*-------------
  IMPORTANT NOTE:
  When an example is added,
  please edit qlat::Generalized::GlobalScopeHash::setup_hash_map at flowed-hmc-generalized.h:31
  and qlat::Generalized::FlowInfo2::FlowInfo2 at flowed-hmc-generalized.h:97.
  Please also be sure of the integer multiplicative_const(=1 if omitted)
  and the boolean is_flowed_link_nonlinear=true/false.
  --------------*/

namespace qlat{
  namespace Generalized{

    class FlowTypePlaq : public FlowType{
    public:
      FlowTypePlaq
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 4,
       const int n_masks_ = 2,
       const std::string description_ = "plaq",
       const int flow_size_ = 1
       );
    };

    class FlowTypeSrect : public FlowType{
    public:
      FlowTypeSrect
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{4,4,4,4}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 6,
       const int n_masks_ = 2,
       const std::string description_ = "srect",
       const int flow_size_ = 2
       );
    };

    class FlowTypeLrect : public FlowType{
    public:
      FlowTypeLrect
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 12,
       const int length_of_loop_ = 6,
       const int n_masks_ = 4,
       const std::string description_ = "lrect",
       const int flow_size_ = 2
       );
    };

    class FlowTypeChairRho1 : public FlowType{
    public:
      FlowTypeChairRho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,2,2}},
       const int n_loops_ = 24,
       const int length_of_loop_ = 6,
       const int n_masks_ = 4,
       const std::string description_ = "chair_rho1",
       const int flow_size_ = 1
       );
    };

    class FlowTypeChairRho2 : public FlowType{
    public:
      FlowTypeChairRho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,1,2}},
       const int n_loops_ = 24,
       const int length_of_loop_ = 6,
       const int n_masks_ = 4,
       const std::string description_ = "chair_rho2",
       const int flow_size_ = 1
       );
    };

    class FlowTypeChairRho3 : public FlowType{
    public:
      FlowTypeChairRho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,1}},
       const int n_loops_ = 24,
       const int length_of_loop_ = 6,
       const int n_masks_ = 4,
       const std::string description_ = "chair_rho3",
       const int flow_size_ = 1
       );
    };


    // class FlowTypeQPlus : public FlowType{
    // public:
    //   FlowTypeQPlus
    //   (
    //    const int mu_=0,
    //    const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
    //    const int n_loops_ = 96, // 96
    //    const int length_of_loop_ = 8,
    //    const int n_masks_ = 2,
    //    const std::string description_ = "q_plus",
    //    const int flow_size_ = 2,
    //    const int multiplicative_const_ = +1
    //    );

    //   static int epsilon(const int mu, const int nu, const int rho, const int sigma );
    //   int sign(const int mu) const& { return epsilon(mu%DIMN,
    //                                                  (mu+1)%DIMN,
    //                                                  (mu+2)%DIMN,
    //                                                  (mu+3)%DIMN); }
    // };


    // class FlowTypeQMinus : public FlowType{
    // public:
    //   FlowTypeQMinus
    //   (
    //    const int mu_=0,
    //    const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
    //    const int n_loops_ = 96, // 96
    //    const int length_of_loop_ = 8,
    //    const int n_masks_ = 2,
    //    const std::string description_ = "q_minus",
    //    const int flow_size_ = 2,
    //    const int multiplicative_const_ = -1
    //    );

    //   const std::function<int(const int, const int, const int, const int)>& epsilon = FlowTypeQPlus::epsilon;
    //   int sign(const int mu) const& { return epsilon(mu%DIMN,
    //                                                  (mu+1)%DIMN,
    //                                                  (mu+2)%DIMN,
    //                                                  (mu+3)%DIMN); }
    // };


    class FlowTypeQPlus2 : public FlowType{
    public:
      FlowTypeQPlus2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 96, // 96
       const int length_of_loop_ = 8,
       const int n_masks_ = 2,
       const std::string description_ = "q_plus",
       const int flow_size_ = 2,
       const int multiplicative_const_ = +1
       );

      static int epsilon(const int mu, const int nu, const int rho, const int sigma );
      int sign(const int mu) const& { return epsilon(mu%DIMN,
                                                     (mu+1)%DIMN,
                                                     (mu+2)%DIMN,
                                                     (mu+3)%DIMN); }
    };


    class FlowTypeQMinus2 : public FlowType{
    public:
      FlowTypeQMinus2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 96, // 96
       const int length_of_loop_ = 8,
       const int n_masks_ = 2,
       const std::string description_ = "q_minus",
       const int flow_size_ = 2,
       const int multiplicative_const_ = -1
       );

      const std::function<int(const int, const int, const int, const int)>& epsilon = FlowTypeQPlus2::epsilon;
      int sign(const int mu) const& { return epsilon(mu%DIMN,
                                                     (mu+1)%DIMN,
                                                     (mu+2)%DIMN,
                                                     (mu+3)%DIMN); }
    };


    class FlowTypePlaqTwice : public FlowType{
    public:
      FlowTypePlaqTwice
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 8,
       const int n_masks_ = 2,
       const std::string description_ = "plaq_twice",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 2,
       const bool is_flowed_link_nonlinear_ = true
       );
    };


    class FlowTypeTwistSrectRho1 : public FlowType{
    public:
      FlowTypeTwistSrectRho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,4,1,1}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 8,
       const int n_masks_ = 4,
       const std::string description_ = "twist_srect_rho1",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false
       );
    };


    class FlowTypeTwistSrectRho2 : public FlowType{
    public:
      FlowTypeTwistSrectRho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,4,1}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 8,
       const int n_masks_ = 4,
       const std::string description_ = "twist_srect_rho2",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false
       );
    };


    class FlowTypeTwistSrectRho3 : public FlowType{
    public:
      FlowTypeTwistSrectRho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,1,4}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 8,
       const int n_masks_ = 4,
       const std::string description_ = "twist_srect_rho3",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false
       );
    };

    class FlowTypeTwistLrect : public FlowType{
    public:
      FlowTypeTwistLrect
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 12,
       const int length_of_loop_ = 8,
       const int n_masks_ = 4,
       const std::string description_ = "twist_lrect",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false
       );
    };

    class FlowTypeTwistMrect : public FlowType{
    public:
      FlowTypeTwistMrect
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 8,
       const int n_masks_ = 2,
       const std::string description_ = "twist_mrect",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1, // DIVIDE TWO CONTRIBUTIONS!!
       const bool is_flowed_link_nonlinear_ = true
       );
    };

    class FlowTypeTwistChairOuterRho1 : public FlowType{
    public:
      FlowTypeTwistChairOuterRho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,2,2}},
       const int n_loops_ = 24,
       const int length_of_loop_ = 8,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_outer_rho1",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false
       );
    };

    class FlowTypeTwistChairOuterRho2 : public FlowType{
    public:
      FlowTypeTwistChairOuterRho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,1,2}},
       const int n_loops_ = 24,
       const int length_of_loop_ = 8,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_outer_rho2",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false
       );
    };

    class FlowTypeTwistChairOuterRho3 : public FlowType{
    public:
      FlowTypeTwistChairOuterRho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,1}},
       const int n_loops_ = 24,
       const int length_of_loop_ = 8,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_outer_rho3",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false
       );
    };

    class FlowTypeTwistChairInner : public FlowType{
    public:
      FlowTypeTwistChairInner
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,2}},
       const int n_loops_ = 24,
       const int length_of_loop_ = 8,
       const int n_masks_ = 2,
       const std::string description_ = "twist_chair_inner",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1, // DIVIDE TWO CONTRIBUTIONS!!
       const bool is_flowed_link_nonlinear_ = true
       );
    };

    class FlowTypePlaqSquared : public FlowType{
    public:
      FlowTypePlaqSquared
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 4,
       const int n_masks_ = 2,
       const std::string description_ = "plaq_squared",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 2,
       const bool is_flowed_link_nonlinear_ = true,
       const int number_of_multiplied_loops_ = 1
       );
    };

    class FlowTypePlaqTimesReversed: public FlowType{
    public:
      FlowTypePlaqTimesReversed
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 4,
       const int n_masks_ = 2,
       const std::string description_ = "plaq_times_reversed",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = true,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeSrectDistinctRho1 : public FlowType{
    public:
      FlowTypeSrectDistinctRho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,4,1,1}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "srect_distinct_rho1",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeSrectDistinctRho2 : public FlowType{
    public:
      FlowTypeSrectDistinctRho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,4,1}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "srect_distinct_rho2",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeSrectDistinctRho3 : public FlowType{
    public:
      FlowTypeSrectDistinctRho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,1,4}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "srect_distinct_rho3",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeLrectDistinct : public FlowType{
    public:
      FlowTypeLrectDistinct
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "lrect_distinct",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };


    class FlowTypeMrectDistinct : public FlowType{
    public:
      FlowTypeMrectDistinct
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 4,
       const int n_masks_ = 2,
       const std::string description_ = "mrect_distinct",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1, // DIVIDE TWO CONTRIBUTIONS!!
       const bool is_flowed_link_nonlinear_ = true,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeChairDistinctTypeARho1 : public FlowType{
    public:
      FlowTypeChairDistinctTypeARho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,2,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "chair_distinct_typea_rho1",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };

    class FlowTypeChairDistinctTypeARho2 : public FlowType{
    public:
      FlowTypeChairDistinctTypeARho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,1,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "chair_distinct_typea_rho2",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };

    class FlowTypeChairDistinctTypeARho3 : public FlowType{
    public:
      FlowTypeChairDistinctTypeARho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,1}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "chair_distinct_typea_rho3",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };

    class FlowTypeChairDistinctTypeBRho1 : public FlowType{
    public:
      FlowTypeChairDistinctTypeBRho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,2,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "chair_distinct_typeb_rho1",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 4
       );
    };

    class FlowTypeChairDistinctTypeBRho2 : public FlowType{
    public:
      FlowTypeChairDistinctTypeBRho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,1,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "chair_distinct_typeb_rho2",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 4
       );
    };

    class FlowTypeChairDistinctTypeBRho3 : public FlowType{
    public:
      FlowTypeChairDistinctTypeBRho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,1}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "chair_distinct_typeb_rho3",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 4
       );
    };

    class FlowTypeChairDistinctInner : public FlowType{
    public:
      FlowTypeChairDistinctInner
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,2}},
       const int n_loops_ = 12,
       const int length_of_loop_ = 4,
       const int n_masks_ = 2,
       const std::string description_ = "chair_distinct_inner",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1, // DIVIDE TWO CONTRIBUTIONS!!
       const bool is_flowed_link_nonlinear_ = true,
       const int number_of_multiplied_loops_ = 2
       );
    };


    class FlowTypeTwistSrectDistinctRho1 : public FlowType{
    public:
      FlowTypeTwistSrectDistinctRho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,4,1,1}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_srect_distinct_rho1",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeTwistSrectDistinctRho2 : public FlowType{
    public:
      FlowTypeTwistSrectDistinctRho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,4,1}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_srect_distinct_rho2",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeTwistSrectDistinctRho3 : public FlowType{
    public:
      FlowTypeTwistSrectDistinctRho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,1,4}},
       const int n_loops_ = 2,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_srect_distinct_rho3",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 1
       );
    };


    class FlowTypeTwistLrectDistinct : public FlowType{
    public:
      FlowTypeTwistLrectDistinct
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{2,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_lrect_distinct",
       const int flow_size_ = 2,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };


    class FlowTypeTwistMrectDistinct : public FlowType{
    public:
      FlowTypeTwistMrectDistinct
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,2}},
       const int n_loops_ = 6,
       const int length_of_loop_ = 4,
       const int n_masks_ = 2,
       const std::string description_ = "twist_mrect_distinct",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1, // DIVIDE TWO CONTRIBUTIONS!!
       const bool is_flowed_link_nonlinear_ = true,
       const int number_of_multiplied_loops_ = 1
       );
    };

    class FlowTypeTwistChairDistinctTypeARho1 : public FlowType{
    public:
      FlowTypeTwistChairDistinctTypeARho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,2,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_distinct_typea_rho1",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };

    class FlowTypeTwistChairDistinctTypeARho2 : public FlowType{
    public:
      FlowTypeTwistChairDistinctTypeARho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,1,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_distinct_typea_rho2",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };

    class FlowTypeTwistChairDistinctTypeARho3 : public FlowType{
    public:
      FlowTypeTwistChairDistinctTypeARho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,1}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_distinct_typea_rho3",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 2
       );
    };

    class FlowTypeTwistChairDistinctTypeBRho1 : public FlowType{
    public:
      FlowTypeTwistChairDistinctTypeBRho1
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,1,2,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_distinct_typeb_rho1",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 4
       );
    };

    class FlowTypeTwistChairDistinctTypeBRho2 : public FlowType{
    public:
      FlowTypeTwistChairDistinctTypeBRho2
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,1,2}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_distinct_typeb_rho2",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 4
       );
    };

    class FlowTypeTwistChairDistinctTypeBRho3 : public FlowType{
    public:
      FlowTypeTwistChairDistinctTypeBRho3
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,1}},
       const int n_loops_ = 4,
       const int length_of_loop_ = 4,
       const int n_masks_ = 4,
       const std::string description_ = "twist_chair_distinct_typeb_rho3",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1,
       const bool is_flowed_link_nonlinear_ = false,
       const int number_of_multiplied_loops_ = 4
       );
    };

    class FlowTypeTwistChairDistinctInner : public FlowType{
    public:
      FlowTypeTwistChairDistinctInner
      (
       const int mu_=0,
       const std::array<int,DIMN> dimensions_ = std::array<int,DIMN>{{1,2,2,2}},
       const int n_loops_ = 12,
       const int length_of_loop_ = 4,
       const int n_masks_ = 2,
       const std::string description_ = "twist_chair_distinct_inner",
       const int flow_size_ = 1,
       const int multiplicative_const_ = 1, // DIVIDE TWO CONTRIBUTIONS!!
       const bool is_flowed_link_nonlinear_ = true,
       const int number_of_multiplied_loops_ = 2
       );
    };



    // -----------------------------------------


    FlowTypePlaq::FlowTypePlaq
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_/*=6*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="plaq"*/,
     const int flow_size_/*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_)
    {
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(nu-1,std::vector<int>{mu_,nu,-mu_-1,-nu-1});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(3+nu-1,std::vector<int>{mu_,-nu-1,-mu_-1,nu});
      const std::function<int(const int, const int, const int, const int)> even_odd
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x1 + x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( even_odd );
      this->initialize();
    }


    FlowTypeSrect::FlowTypeSrect
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*= std::array<int,DIMN>{{4,4,4,4}}*/,
     const int n_loops_/*=6*/,
     const int length_of_loop_/*=6*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="srect"*/,
     const int flow_size_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_)
    {
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(nu-1,std::vector<int>{mu_,nu,nu,-mu_-1,-nu-1,-nu-1});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(3+nu-1,std::vector<int>{mu_,-nu-1,-nu-1,-mu_-1,nu,nu});
      const std::function<int(const int, const int, const int, const int)> even_odd_block2
        = [](const int x1, const int x2, const int x3, const int x4) {
        return ( (x1+1024)/2 + (x2+1024)/2 + (x3+1024)/2 + (x4+1024)/2 )%2;
      };
      this->set_mask( even_odd_block2 );
      this->initialize();
    }


    FlowTypeLrect::FlowTypeLrect
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_/*=12*/,
     const int length_of_loop_/*=6*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="lrect"*/,
     const int flow_size_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_)
    {
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(nu-1,std::vector<int>{mu_,mu_,nu,-mu_-1,-mu_-1,-nu-1});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(3+nu-1,std::vector<int>{mu_,nu,-mu_-1,-mu_-1,-nu-1,mu_});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(6+nu-1,std::vector<int>{mu_,mu_,-nu-1,-mu_-1,-mu_-1,nu});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(9+nu-1,std::vector<int>{mu_,-nu-1,-mu_-1,-mu_-1,nu,mu_});

      const std::function<int(const int, const int, const int, const int)> lrect_mask
        = [](const int x1, const int x2, const int x3, const int x4) {
        return 2*( (x1+1024)%2) + ( x2 + x3 + x4 + 1024)%2;
      };
      this->set_mask( lrect_mask );
      this->initialize();
    }


    FlowTypeChairRho1::FlowTypeChairRho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*= std::array<int,DIMN>{{1,1,2,2}}*/,
     const int n_loops_/*=24*/,
     const int length_of_loop_/*=6*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_rho1"*/,
     const int flow_size_/*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_)
    {
      const int rho = 1;
      int counter = 0;
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,rho,-mu_-1,-rho-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-rho-1,-mu_-1,rho,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,rho,-mu_-1,-rho-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-rho-1,-mu_-1,rho,nu});
        ++counter;
      }
      //
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,nu,-rho-1,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,nu,rho,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,-nu-1,-rho-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,-nu-1,rho,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,rho,-nu-1,-rho-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-rho-1,-nu-1,rho});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,rho,nu,-rho-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,nu,rho});
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
        return 2*( (x3+1024)%2) + ( x4 + 1024 )%2;
      };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeChairRho2::FlowTypeChairRho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,1,2}}*/,
     const int n_loops_/*=24*/,
     const int length_of_loop_/*=6*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_rho2"*/,
     const int flow_size_/*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_)
    {
      const int rho = 2;
      int counter = 0;
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,rho,-mu_-1,-rho-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,rho,-mu_-1,-rho-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-rho-1,-mu_-1,rho,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-rho-1,-mu_-1,rho,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,rho,-mu_-1,-rho-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,rho,-mu_-1,-rho-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-rho-1,-mu_-1,rho,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-rho-1,-mu_-1,rho,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,nu,-rho-1,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,nu,-rho-1,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,nu,rho,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,nu,rho,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,-nu-1,-rho-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,-nu-1,-rho-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,-nu-1,rho,-mu_-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,-nu-1,rho,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,rho,-nu-1,-rho-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,rho,-nu-1,-rho-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-rho-1,-nu-1,rho});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-rho-1,-nu-1,rho});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,rho,nu,-rho-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,rho,nu,-rho-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,nu,rho});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,nu,rho});
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
        return 2*( (x2+1024)%2) + ( x4 + 1024 )%2;
      };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeChairRho3::FlowTypeChairRho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,2,1}}*/,
     const int n_loops_/*=24*/,
     const int length_of_loop_/*=6*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_rho3"*/,
     const int flow_size_/*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_)
    {
      const int rho = 3;
      int counter = 0;
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,rho,-mu_-1,-rho-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-rho-1,-mu_-1,rho,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,rho,-mu_-1,-rho-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-rho-1,-mu_-1,rho,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,nu,-rho-1,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,nu,rho,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,rho,-nu-1,-rho-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-rho-1,-nu-1,rho,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,rho,-nu-1,-rho-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-rho-1,-nu-1,rho});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,rho,nu,-rho-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,nu,rho});
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
        return 2*( (x2+1024)%2) + ( x3 + 1024 )%2;
      };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    // FlowTypeQPlus::FlowTypeQPlus
    // (
    //  const int mu_/*=0*/,
    //  const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
    //  const int n_loops_/*=96*/,
    //  const int length_of_loop_/*=8*/,
    //  const int n_masks_/*=2*/,
    //  const std::string description_/*="q_plus"*/,
    //  const int flow_size_/*=2*/,
    //  const int multiplicative_const_/*=+1*/
    //  )
    //   : UnitInfo(dimensions_)
    //   , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
    //              multiplicative_const_)
    // {
    //   int counter = 0;

    //   // const int nu=1;
    //   // const int rho=2;
    //   // const int sigma=3;
    //   // assert(epsilon(mu_,nu,rho,sigma)==1);
    //   // this->set_shape(counter,
    //   //                 std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,sigma,-rho-1,-sigma-1});
    //   // ++counter;

    //   for(int nu=mu_+1; nu<DIMN; ++nu){
    //     for(int rho=1; rho<DIMN; ++rho){
    //       if(rho==nu) continue;
    //       for(int sigma=rho+1; sigma<DIMN; ++sigma){
    //         if(sigma==nu) continue;
    //         //
    //         if(epsilon(mu_,nu,rho,sigma)==1){ // epsilon_{mu_,nu,rho,sigma} = 1
    //           // type 1
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,sigma,-rho-1,-sigma-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,sigma,rho,-sigma-1,-rho-1});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,sigma,-rho-1,-sigma-1,rho});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-sigma-1,rho,sigma,-rho-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-rho-1,sigma,rho,-sigma-1});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,rho,-sigma-1,-rho-1,sigma});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-rho-1,-sigma-1,rho,sigma});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-sigma-1,-rho-1,sigma,rho});
    //           ++counter;
    //           // type 2
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,rho,sigma,-rho-1,-sigma-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,sigma,rho,-sigma-1,-rho-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,sigma,-rho-1,-sigma-1,rho,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-sigma-1,rho,sigma,-rho-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,sigma,rho,-sigma-1,nu});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,rho,-sigma-1,-rho-1,sigma,nu});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-rho-1,-sigma-1,rho,sigma,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-sigma-1,-rho-1,sigma,rho,nu});
    //           ++counter;
    //           // type 3
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,rho,sigma,-rho-1,-sigma-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,sigma,rho,-sigma-1,-rho-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,sigma,-rho-1,-sigma-1,rho,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-sigma-1,rho,sigma,-rho-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-rho-1,sigma,rho,-sigma-1,-mu_-1,nu});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,rho,-sigma-1,-rho-1,sigma,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-rho-1,-sigma-1,rho,sigma,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-sigma-1,-rho-1,sigma,rho,-mu_-1,nu});
    //           ++counter;
    //           // type 4
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,sigma,-rho-1,-sigma-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,rho,-sigma-1,-rho-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,-rho-1,-sigma-1,rho,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,rho,sigma,-rho-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,sigma,rho,-sigma-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,-sigma-1,-rho-1,sigma,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,-sigma-1,rho,sigma,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,-rho-1,sigma,rho,-nu-1,-mu_-1,nu});
    //           ++counter;
    //         }
    //         else { // epsilon_{mu_,nu,rho,sigma} = -1
    //           // type 1
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,sigma,rho,-sigma-1,-rho-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,rho,sigma,-rho-1,-sigma-1});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-rho-1,sigma,rho,-sigma-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,-sigma-1,-rho-1,sigma});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,sigma,-rho-1,-sigma-1,rho});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-sigma-1,rho,sigma,-rho-1});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-sigma-1,-rho-1,sigma,rho});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-rho-1,-sigma-1,rho,sigma});
    //           ++counter;
    //           // type 2
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,sigma,rho,-sigma-1,-rho-1,-nu-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,rho,sigma,-rho-1,-sigma-1,nu});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-rho-1,sigma,rho,-sigma-1,-nu-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,rho,-sigma-1,-rho-1,sigma,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,sigma,-rho-1,-sigma-1,rho,nu});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-sigma-1,rho,sigma,-rho-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-sigma-1,-rho-1,sigma,rho,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,-sigma-1,rho,sigma,nu});
    //           ++counter;
    //           // type 3
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,sigma,rho,-sigma-1,-rho-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,rho,sigma,-rho-1,-sigma-1,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-rho-1,sigma,rho,-sigma-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,rho,-sigma-1,-rho-1,sigma,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,sigma,-rho-1,-sigma-1,rho,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-sigma-1,rho,sigma,-rho-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-sigma-1,-rho-1,sigma,rho,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-rho-1,-sigma-1,rho,sigma,-mu_-1,nu});
    //           ++counter;
    //           // type 4
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,rho,-sigma-1,-rho-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,sigma,-rho-1,-sigma-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,sigma,rho,-sigma-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,-sigma-1,-rho-1,sigma,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,-rho-1,-sigma-1,rho,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,rho,sigma,-rho-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,-rho-1,sigma,rho,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,-sigma-1,rho,sigma,-nu-1,-mu_-1,nu});
    //           ++counter;
    //         } // end else if epsilon...
    //       } // end for sigma
    //     } // end for rho
    //   } // end for nu
    //   assert(counter==n_loops_);

    //   const std::function<int(const int, const int, const int, const int)> even_odd
    //     = [](const int x1, const int x2, const int x3, const int x4) {
    //     return (x1 + x2 + x3 + x4 + 1024)%2;
    //   };
    //   this->set_mask( even_odd );
    //   this->initialize();
    // }


    // int FlowTypeQPlus::epsilon(const int mu, const int nu, const int rho, const int sigma ){
    //   int sgn = 0;
    //   std::array<int,DIMN> array{{mu,nu,rho,sigma}};
    //   for(auto itr=array.begin(); itr!=array.end(); ++itr){
    //     if(*itr<0){
    //       *itr = -(*itr)-1;
    //       ++sgn;
    //     }
    //   }
    //   assert( (array[0]-array[1])*(array[0]-array[2])*(array[0]-array[3])
    //           *(array[1]-array[2])*(array[1]-array[3])
    //           *(array[2]-array[3])
    //           !=0 );
    //   for(auto itr=array.begin(); itr!=array.end(); ){
    //     const auto smaller = std::find_if(itr+1, array.end(),
    //                                       [&itr](const int alpha){ return (*itr)>alpha; });
    //     if(smaller!=array.end()){
    //       ++sgn;
    //       const int tmp = *itr;
    //       *itr = *smaller;
    //       *smaller = tmp;
    //     }
    //     else ++itr;
    //   }
    //   return 1 - 2*(sgn%2);
    // }


    // FlowTypeQMinus::FlowTypeQMinus
    // (
    //  const int mu_/*=0*/,
    //  const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
    //  const int n_loops_/*=96*/,
    //  const int length_of_loop_/*=8*/,
    //  const int n_masks_/*=2*/,
    //  const std::string description_/*="q_minus"*/,
    //  const int flow_size_/*=2*/,
    //  const int multiplicative_const_/*=-1*/
    //  )
    //   : UnitInfo(dimensions_)
    //   , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
    //              multiplicative_const_)
    // {
    //   int counter = 0;
    //   // const int nu=1;
    //   // const int rho=3;
    //   // const int sigma=2;
    //   // assert(epsilon(mu_,nu,rho,sigma)==-1);
    //   // this->set_shape(counter,
    //   //                 std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,sigma,-rho-1,-sigma-1});
    //   // ++counter;

    //   for(int nu=mu_+1; nu<DIMN; ++nu){
    //     for(int rho=1; rho<DIMN; ++rho){
    //       if(nu==rho) continue;
    //       for(int sigma=rho+1; sigma<DIMN; ++sigma){
    //         if(nu==sigma) continue;
    //         //
    //         if(epsilon(mu_,nu,rho,sigma)==-1){ // epsilon_{mu_,nu,rho,sigma} = -1
    //           // type 1
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,sigma,-rho-1,-sigma-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,sigma,rho,-sigma-1,-rho-1});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,sigma,-rho-1,-sigma-1,rho});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-sigma-1,rho,sigma,-rho-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-rho-1,sigma,rho,-sigma-1});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,rho,-sigma-1,-rho-1,sigma});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-rho-1,-sigma-1,rho,sigma});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-sigma-1,-rho-1,sigma,rho});
    //           ++counter;
    //           // type 2
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,rho,sigma,-rho-1,-sigma-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,sigma,rho,-sigma-1,-rho-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,sigma,-rho-1,-sigma-1,rho,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-sigma-1,rho,sigma,-rho-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,sigma,rho,-sigma-1,nu});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,rho,-sigma-1,-rho-1,sigma,nu});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-rho-1,-sigma-1,rho,sigma,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-sigma-1,-rho-1,sigma,rho,nu});
    //           ++counter;
    //           // type 3
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,rho,sigma,-rho-1,-sigma-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,sigma,rho,-sigma-1,-rho-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,sigma,-rho-1,-sigma-1,rho,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-sigma-1,rho,sigma,-rho-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-rho-1,sigma,rho,-sigma-1,-mu_-1,nu});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,rho,-sigma-1,-rho-1,sigma,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-rho-1,-sigma-1,rho,sigma,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-sigma-1,-rho-1,sigma,rho,-mu_-1,nu});
    //           ++counter;
    //           // type 4
    //           // standard
    //           assert(epsilon(mu_,nu,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,sigma,-rho-1,-sigma-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped
    //           assert(epsilon(mu_,-nu-1,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,rho,-sigma-1,-rho-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped
    //           assert(epsilon(mu_,nu,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,-rho-1,-sigma-1,rho,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,rho,sigma,-rho-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,sigma,rho,-sigma-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,-sigma-1,-rho-1,sigma,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,-sigma-1,rho,sigma,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,-rho-1,sigma,rho,-nu-1,-mu_-1,nu});
    //           ++counter;
    //         }
    //         else { // epsilon_{mu_,nu,rho,sigma} = +1
    //           // type 1
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,sigma,rho,-sigma-1,-rho-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,rho,sigma,-rho-1,-sigma-1});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-rho-1,sigma,rho,-sigma-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,-sigma-1,-rho-1,sigma});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,sigma,-rho-1,-sigma-1,rho});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-sigma-1,rho,sigma,-rho-1});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-nu-1,-sigma-1,-rho-1,sigma,rho});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,nu,-rho-1,-sigma-1,rho,sigma});
    //           ++counter;
    //           // type 2
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,sigma,rho,-sigma-1,-rho-1,-nu-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,rho,sigma,-rho-1,-sigma-1,nu});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-rho-1,sigma,rho,-sigma-1,-nu-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,rho,-sigma-1,-rho-1,sigma,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,sigma,-rho-1,-sigma-1,rho,nu});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-sigma-1,rho,sigma,-rho-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-mu_-1,-sigma-1,-rho-1,sigma,rho,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,-sigma-1,rho,sigma,nu});
    //           ++counter;
    //           // type 3
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,sigma,rho,-sigma-1,-rho-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,rho,sigma,-rho-1,-sigma-1,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-rho-1,sigma,rho,-sigma-1,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,rho,-sigma-1,-rho-1,sigma,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,sigma,-rho-1,-sigma-1,rho,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-sigma-1,rho,sigma,-rho-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,nu,-sigma-1,-rho-1,sigma,rho,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-nu-1,-rho-1,-sigma-1,rho,sigma,-mu_-1,nu});
    //           ++counter;
    //           // type 4
    //           // sigma<->rho permutation
    //           assert(epsilon(mu_,nu,sigma,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,rho,-sigma-1,-rho-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped
    //           assert(epsilon(mu_,-nu-1,rho,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,sigma,-rho-1,-sigma-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // rho flipped
    //           assert(epsilon(mu_,nu,-rho-1,sigma)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,sigma,rho,-sigma-1,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma flipped
    //           assert(epsilon(mu_,nu,rho,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,rho,-sigma-1,-rho-1,sigma,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, rho flipped
    //           assert(epsilon(mu_,-nu-1,sigma,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,sigma,-rho-1,-sigma-1,rho,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, nu flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-sigma-1,rho)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,rho,sigma,-rho-1,-nu-1,-mu_-1,nu});
    //           ++counter;
    //           // sigma<->rho permutation, rho flipped, sigma flipped
    //           assert(epsilon(mu_,nu,-sigma-1,-rho-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-sigma-1,-rho-1,sigma,rho,nu,-mu_-1,-nu-1});
    //           ++counter;
    //           // nu flipped, rho flipped, sigma flipped
    //           assert(epsilon(mu_,-nu-1,-rho-1,-sigma-1)==-1);
    //           this->set_shape(counter,
    //                           std::vector<int>{mu_,-rho-1,-sigma-1,rho,sigma,-nu-1,-mu_-1,nu});
    //           ++counter;
    //         } // end else if epsilon...
    //       } // end for sigma
    //     } // end for rho
    //   } // end for nu
    //   assert(counter==n_loops_);

    //   const std::function<int(const int, const int, const int, const int)> even_odd
    //     = [](const int x1, const int x2, const int x3, const int x4) {
    //     return (x1 + x2 + x3 + x4 + 1024)%2;
    //   };
    //   this->set_mask( even_odd );
    //   this->initialize();
    // }


    FlowTypeQPlus2::FlowTypeQPlus2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_/*=96*/,
     const int length_of_loop_/*=8*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="q_plus"*/,
     const int flow_size_/*=2*/,
     const int multiplicative_const_/*=+1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_)
    {
      int counter = 0;

      for(int nu=-DIMN; nu<DIMN; ++nu){
        for(int rho=-DIMN; rho<DIMN; ++rho){
          for(int sigma=-DIMN; sigma<DIMN; ++sigma){

            if(epsilon(mu_,nu,rho,sigma)==1){
              // type 1
              this->set_shape(counter,
                              std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,sigma,-rho-1,-sigma-1});
              ++counter;
              // type 2
              this->set_shape(counter,
                              std::vector<int>{mu_,nu,-mu_-1,rho,sigma,-rho-1,-sigma-1,-nu-1});
              ++counter;
              // type 3
              this->set_shape(counter,
                              std::vector<int>{mu_,nu,rho,sigma,-rho-1,-sigma-1,-mu_-1,-nu-1});
              ++counter;
              // type 4
              this->set_shape(counter,
                              std::vector<int>{mu_,rho,sigma,-rho-1,-sigma-1,nu,-mu_-1,-nu-1});
              ++counter;
            }

          } // end for sigma
        } // end for rho
      } // end for nu
      assert(counter==n_loops_);

      const std::function<int(const int, const int, const int, const int)> even_odd
        = [](const int x1, const int x2, const int x3, const int x4) {
        return (x1 + x2 + x3 + x4 + 1024)%2;
      };
      this->set_mask( even_odd );
      this->initialize();
    }

    int FlowTypeQPlus2::epsilon(const int mu, const int nu, const int rho, const int sigma ){
      int sgn = 0;
      std::array<int,DIMN> array{{mu,nu,rho,sigma}};
      for(auto itr=array.begin(); itr!=array.end(); ++itr){
        if(*itr<0){
          *itr = -(*itr)-1;
          ++sgn;
        }
      }
      if( (array[0]-array[1])*(array[0]-array[2])*(array[0]-array[3])
          *(array[1]-array[2])*(array[1]-array[3])
          *(array[2]-array[3])
          ==0 ){
        return 0;
      }
      for(auto itr=array.begin(); itr!=array.end(); ){
        const auto smaller = std::find_if(itr+1, array.end(),
                                          [&itr](const int alpha){ return (*itr)>alpha; });
        if(smaller!=array.end()){
          ++sgn;
          const int tmp = *itr;
          *itr = *smaller;
          *smaller = tmp;
        }
        else ++itr;
      }
      return 1 - 2*(sgn%2);
    }


    FlowTypeQMinus2::FlowTypeQMinus2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_/*=96*/,
     const int length_of_loop_/*=8*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="q_minus"*/,
     const int flow_size_/*=2*/,
     const int multiplicative_const_/*=-1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_)
    {
      int counter = 0;
      for(int nu=-DIMN; nu<DIMN; ++nu){
        for(int rho=-DIMN; rho<DIMN; ++rho){
          for(int sigma=-DIMN; sigma<DIMN; ++sigma){
            if(epsilon(mu_,nu,rho,sigma)==-1){
              // type 1
              this->set_shape(counter,
                              std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,sigma,-rho-1,-sigma-1});
              ++counter;
              // type 2
              this->set_shape(counter,
                              std::vector<int>{mu_,nu,-mu_-1,rho,sigma,-rho-1,-sigma-1,-nu-1});
              ++counter;
              // type 3
              this->set_shape(counter,
                              std::vector<int>{mu_,nu,rho,sigma,-rho-1,-sigma-1,-mu_-1,-nu-1});
              ++counter;
              // type 4
              this->set_shape(counter,
                              std::vector<int>{mu_,rho,sigma,-rho-1,-sigma-1,nu,-mu_-1,-nu-1});
              ++counter;
            }
          } // end for sigma
        } // end for rho
      } // end for nu
      assert(counter==n_loops_);

      const std::function<int(const int, const int, const int, const int)> even_odd
        = [](const int x1, const int x2, const int x3, const int x4) {
        return (x1 + x2 + x3 + x4 + 1024)%2;
      };
      this->set_mask( even_odd );
      this->initialize();
    }


    FlowTypePlaqTwice::FlowTypePlaqTwice
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_/*=6*/,
     const int length_of_loop_/*=8*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="plaq_twice"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_,/*=2*/
     const bool is_flowed_link_nonlinear/*true*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear)
    {
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(nu-1,std::vector<int>{mu_,nu,-mu_-1,-nu-1,mu_,nu,-mu_-1,-nu-1});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(3+nu-1,std::vector<int>{mu_,-nu-1,-mu_-1,nu,mu_,-nu-1,-mu_-1,nu});
      const std::function<int(const int, const int, const int, const int)> even_odd
        = [](const int x1, const int x2, const int x3, const int x4) {
        return (x1 + x2 + x3 + x4 + 1024)%2;
      };
      this->set_mask( even_odd );
      this->initialize();
    }


    FlowTypeTwistSrectRho1::FlowTypeTwistSrectRho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,4,1,1}}*/,
     const int n_loops_/*=2*/,
     const int length_of_loop_/*=8*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_srect_rho1"*/,
     const int flow_size_/*=2*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear/*false*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear)
    {
      const int rho = 1;
      this->set_shape(0,std::vector<int>{mu_,rho,-mu_-1,rho,mu_,-rho-1,-mu_-1,-rho-1});
      this->set_shape(1,std::vector<int>{mu_,-rho-1,-mu_-1,-rho-1,mu_,rho,-mu_-1,rho});
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
        return (x2+1024)%4;
      };
      this->set_mask( div4rho );
      this->initialize();
    }


    FlowTypeTwistSrectRho2::FlowTypeTwistSrectRho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*= std::array<int,DIMN>{{1,1,4,1}}*/,
     const int n_loops_ /*= 2*/,
     const int length_of_loop_ /*= 8*/,
     const int n_masks_ /*= 4*/,
     const std::string description_ /*= "twist_srect_rho2"*/,
     const int flow_size_ /*= 2*/,
     const int multiplicative_const_,/*=1*/
     const bool is_flowed_link_nonlinear/*false*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear)
    {
      const int rho = 2;
      this->set_shape(0,std::vector<int>{mu_,rho,-mu_-1,rho,mu_,-rho-1,-mu_-1,-rho-1});
      this->set_shape(1,std::vector<int>{mu_,-rho-1,-mu_-1,-rho-1,mu_,rho,-mu_-1,rho});
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
        return (x3+1024)%4;
      };
      this->set_mask( div4rho );
      this->initialize();
    }


    FlowTypeTwistSrectRho3::FlowTypeTwistSrectRho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*= std::array<int,DIMN>{{1,1,1,4}}*/,
     const int n_loops_ /*= 2*/,
     const int length_of_loop_ /*= 8*/,
     const int n_masks_ /*= 4*/,
     const std::string description_ /*= "twist_srect_rho3"*/,
     const int flow_size_ /*= 2*/,
     const int multiplicative_const_,/*=1*/
     const bool is_flowed_link_nonlinear/*false*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear)
    {
      const int rho = 3;
      this->set_shape(0,std::vector<int>{mu_,rho,-mu_-1,rho,mu_,-rho-1,-mu_-1,-rho-1});
      this->set_shape(1,std::vector<int>{mu_,-rho-1,-mu_-1,-rho-1,mu_,rho,-mu_-1,rho});
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
        return (x4+1024)%4;
      };
      this->set_mask( div4rho );
      this->initialize();
    }


    FlowTypeTwistLrect::FlowTypeTwistLrect
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*= std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_ /*= 12*/,
     const int length_of_loop_ /*= 8*/,
     const int n_masks_ /*= 4*/,
     const std::string description_ /*= "twist_lrect"*/,
     const int flow_size_ /*= 2*/,
     const int multiplicative_const_,/*=1*/
     const bool is_flowed_link_nonlinear/*false*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear)
    {
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(nu-1,std::vector<int>{mu_,nu,mu_,-nu-1,-mu_-1,nu,-mu_-1,-nu-1});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(3+nu-1,std::vector<int>{mu_,nu,-mu_-1,-nu-1,-mu_-1,nu,mu_,-nu-1});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(6+nu-1,std::vector<int>{mu_,-nu-1,mu_,nu,-mu_-1,-nu-1,-mu_-1,nu});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(9+nu-1,std::vector<int>{mu_,-nu-1,-mu_-1,nu,-mu_-1,-nu-1,mu_,nu});

      const std::function<int(const int, const int, const int, const int)> lrect_mask
        = [](const int x1, const int x2, const int x3, const int x4) {
        return 2*( (x1+1024)%2) + ( x2 + x3 + x4 + 1024)%2;
      };
      this->set_mask( lrect_mask );
      this->initialize();
    }


    FlowTypeTwistMrect::FlowTypeTwistMrect
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*= std::array<int,DIMN>{{1,2,2,2}}*/,
     const int n_loops_ /*= 6*/,
     const int length_of_loop_ /*= 8*/,
     const int n_masks_ /*= 2*/,
     const std::string description_ /*= "twist_mrect"*/,
     const int flow_size_ /*= 1*/,
     const int multiplicative_const_ /*= 1*/, // DIVIDE TWO CONTRIBUTIONS!!
     const bool is_flowed_link_nonlinear/*true*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear)
    {
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(nu-1,std::vector<int>{mu_,nu,-mu_-1,-nu-1,mu_,-nu-1,-mu_-1,nu});
      for(int nu=mu_+1; nu<DIMN; ++nu)
        this->set_shape(3+nu-1,std::vector<int>{mu_,-nu-1,-mu_-1,nu,mu_,nu,-mu_-1,-nu-1});
      const std::function<int(const int, const int, const int, const int)> eo_prop_mu
        = [](const int x1, const int x2, const int x3, const int x4) {
        return (x2 + x3 + x4 + 1024)%2;
      };
      this->set_mask( eo_prop_mu );
      this->initialize();
    }

    //

    FlowTypeTwistChairOuterRho1::FlowTypeTwistChairOuterRho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*= std::array<int,DIMN>{{1,1,2,2}}*/,
     const int n_loops_/*=24*/,
     const int length_of_loop_/*=8*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_outer_rho1"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_)
    {
      const int rho = 1;
      int counter = 0;
      for(int nu=rho+1; nu<DIMN; ++nu){
        // standard
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,rho,mu_,-rho-1,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip rho
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-rho-1,mu_,rho,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip nu
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,rho,mu_,-rho-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip rho, mu
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,mu_,rho,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=rho+1; nu<DIMN; ++nu){
        // standard
        this->set_shape(counter,std::vector<int>{mu_,nu,rho,-nu-1,-rho-1,nu,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip rho
        this->set_shape(counter,std::vector<int>{mu_,nu,-rho-1,-nu-1,rho,nu,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip nu
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,rho,nu,-rho-1,-nu-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip rho, nu
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-rho-1,nu,rho,-nu-1,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=rho+1; nu<DIMN; ++nu){
        // standard
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,nu,-rho-1,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip rho
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,-rho-1,nu,rho,-nu-1});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip nu
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,rho,-nu-1,-rho-1,nu});
        ++counter;
      }
      for(int nu=rho+1; nu<DIMN; ++nu){
        // flip rho, nu
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,-rho-1,-nu-1,rho,nu});
        ++counter;
      }
      assert(counter==n_loops_);

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x3+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeTwistChairOuterRho2::FlowTypeTwistChairOuterRho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*= std::array<int,DIMN>{{1,2,1,2}}*/,
     const int n_loops_/*=24*/,
     const int length_of_loop_/*=8*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_outer_rho2"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_)
    {
      const int rho = 2;
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,rho,mu_,-rho-1,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-rho-1,mu_,rho,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,rho,mu_,-rho-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,mu_,rho,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,nu,rho,-nu-1,-rho-1,nu,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,nu,-rho-1,-nu-1,rho,nu,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,rho,nu,-rho-1,-nu-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-rho-1,nu,rho,-nu-1,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,nu,-rho-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,-rho-1,nu,rho,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,rho,-nu-1,-rho-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,-rho-1,-nu-1,rho,nu});
        ++counter;
      }
      assert(counter==n_loops_);

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeTwistChairOuterRho3::FlowTypeTwistChairOuterRho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*= std::array<int,DIMN>{{1,2,2,1}}*/,
     const int n_loops_/*=24*/,
     const int length_of_loop_/*=8*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_outer_rho3"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_)
    {
      const int rho = 3;
      int counter = 0;
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,rho,mu_,-rho-1,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-rho-1,mu_,rho,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,rho,mu_,-rho-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,-rho-1,mu_,rho,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,rho,-nu-1,-rho-1,nu,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-rho-1,-nu-1,rho,nu,-mu_-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,rho,nu,-rho-1,-nu-1,-mu_-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-rho-1,nu,rho,-nu-1,-mu_-1,nu});
        ++counter;
      }
      //
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,rho,nu,-rho-1,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,-rho-1,nu,rho,-nu-1});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,rho,-nu-1,-rho-1,nu});
        ++counter;
      }
      for(int nu=mu_+1; nu<rho; ++nu){
        this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,-rho-1,-nu-1,rho,nu});
        ++counter;
      }
      assert(counter==n_loops_);

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x3 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeTwistChairInner::FlowTypeTwistChairInner
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*= std::array<int,DIMN>{{1,2,2,2}}*/,
     const int n_loops_ /*= 24*/,
     const int length_of_loop_ /*= 8*/,
     const int n_masks_ /*= 2*/,
     const std::string description_ /*= "twist_chair_inner"*/,
     const int flow_size_ /*= 1*/,
     const int multiplicative_const_ /*= 1*/, // DIVIDE TWO CONTRIBUTIONS!!
     const bool is_flowed_link_nonlinear_/*true*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_)
    {
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        for(int rho=nu+1; rho<DIMN; ++rho){
          this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,mu_,rho,-mu_-1,-rho-1});
          ++counter;
          this->set_shape(counter,std::vector<int>{mu_,rho,-mu_-1,-rho-1,mu_,nu,-mu_-1,-nu-1});
          ++counter;
          // flip rho
          this->set_shape(counter,std::vector<int>{mu_,nu,-mu_-1,-nu-1,mu_,-rho-1,-mu_-1,rho});
          ++counter;
          this->set_shape(counter,std::vector<int>{mu_,-rho-1,-mu_-1,rho,mu_,nu,-mu_-1,-nu-1});
          ++counter;
          // flip nu
          this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,mu_,rho,-mu_-1,-rho-1});
          ++counter;
          this->set_shape(counter,std::vector<int>{mu_,rho,-mu_-1,-rho-1,mu_,-nu-1,-mu_-1,nu});
          ++counter;
          // flip rho, nu
          this->set_shape(counter,std::vector<int>{mu_,-nu-1,-mu_-1,nu,mu_,-rho-1,-mu_-1,rho});
          ++counter;
          this->set_shape(counter,std::vector<int>{mu_,-rho-1,-mu_-1,rho,mu_,-nu-1,-mu_-1,nu});
          ++counter;
        }
      }
      assert(counter==n_loops_);

      const std::function<int(const int, const int, const int, const int)> eo_prop_mu
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( eo_prop_mu );
      this->initialize();
    }

    FlowTypePlaqSquared::FlowTypePlaqSquared
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_/*=6*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="plaq_squared"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=2*/,
     const bool is_flowed_link_nonlinear_/*=true*/,
     const int number_of_multiplied_loops_/*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
        this->set_shape(counter,plaq_up);
        const std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
        this->set_shape(counter,plaq_down);
        const std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        add_multiplied_shape(counter,0,plaq_down,relative_coordinates);
        ++counter;
      }
      const std::function<int(const int, const int, const int, const int)> even_odd
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x1 + x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( even_odd );
      this->initialize();
    }


    FlowTypePlaqTimesReversed::FlowTypePlaqTimesReversed
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_/*=6*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="plaq_times_reversed"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=true*/,
     const int number_of_multiplied_loops_/*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
        const std::vector<int> plaq_reversed{nu,mu_,-nu-1,-mu_-1};
        this->set_shape(counter,plaq_up);
        const std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        add_multiplied_shape(counter,0,plaq_reversed,relative_coordinates);
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
        const std::vector<int> plaq_reversed{-nu-1,mu_,nu,-mu_-1};
        this->set_shape(counter,plaq_down);
        const std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        add_multiplied_shape(counter,0,plaq_reversed,relative_coordinates);
        ++counter;
      }
      const std::function<int(const int, const int, const int, const int)> even_odd
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x1 + x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( even_odd );
      this->initialize();
    }



    FlowTypeSrectDistinctRho1::FlowTypeSrectDistinctRho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,4,1,1}}*/,
     const int n_loops_ /*=2*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="srect_distinct_rho1"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 1;
      {
        const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
        this->set_shape(0,plaq_up);
        const std::array<int, DIMN> relative_coordinates{{0,1,0,0}};
        add_multiplied_shape(0,0,plaq_up,relative_coordinates);
      }
      {
        const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
        this->set_shape(1,plaq_down);
        const std::array<int, DIMN> relative_coordinates{{0,-1,0,0}};
        add_multiplied_shape(1,0,plaq_down,relative_coordinates);
      }
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x2+1024)%4;
        };
      this->set_mask( div4rho );
      this->initialize();
    }

    FlowTypeSrectDistinctRho2::FlowTypeSrectDistinctRho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,1,4,1}}*/,
     const int n_loops_ /*=2*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="srect_distinct_rho2"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 2;
      {
        const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
        this->set_shape(0,plaq_up);
        const std::array<int, DIMN> relative_coordinates{{0,0,1,0}};
        add_multiplied_shape(0,0,plaq_up,relative_coordinates);
      }
      {
        const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
        this->set_shape(1,plaq_down);
        const std::array<int, DIMN> relative_coordinates{{0,0,-1,0}};
        add_multiplied_shape(1,0,plaq_down,relative_coordinates);
      }
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x3+1024)%4;
        };
      this->set_mask( div4rho );
      this->initialize();
    }


    FlowTypeSrectDistinctRho3::FlowTypeSrectDistinctRho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,1,1,4}}*/,
     const int n_loops_ /*=2*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="srect_distinct_rho3"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 3;
      {
        const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
        this->set_shape(0,plaq_up);
        const std::array<int, DIMN> relative_coordinates{{0,0,0,1}};
        add_multiplied_shape(0,0,plaq_up,relative_coordinates);
      }
      {
        const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
        this->set_shape(1,plaq_down);
        const std::array<int, DIMN> relative_coordinates{{0,0,0,-1}};
        add_multiplied_shape(1,0,plaq_down,relative_coordinates);
      }
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x4+1024)%4;
        };
      this->set_mask( div4rho );
      this->initialize();
    }


    FlowTypeLrectDistinct::FlowTypeLrectDistinct
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_ /*=6*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="lrect_distinct"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
        this->set_shape(counter,plaq_up);
        std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        relative_coordinates[mu_] = 1;
        add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
        relative_coordinates[mu_] = -1;
        add_multiplied_shape(counter,1,plaq_up,relative_coordinates);
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
        this->set_shape(counter,plaq_down);
        std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        relative_coordinates[mu_] = 1;
        add_multiplied_shape(counter,0,plaq_down,relative_coordinates);
        relative_coordinates[mu_] = -1;
        add_multiplied_shape(counter,1,plaq_down,relative_coordinates);
        ++counter;
      }
      const std::function<int(const int, const int, const int, const int)> lrect_mask
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x1+1024)%2) + ( x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( lrect_mask );
      this->initialize();
    }


    FlowTypeMrectDistinct::FlowTypeMrectDistinct
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,2,2,2}}*/,
     const int n_loops_ /*=6*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=2*/,
     const std::string description_ /*="mrect_distinct"*/,
     const int flow_size_ /*=1*/,
     const int multiplicative_const_ /*=1*/, // DIVIDE TWO CONTRIBUTIONS!!
     const bool is_flowed_link_nonlinear_ /*=true*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
        const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
        std::array<int, DIMN> relative_coordinates{{0,0,0,0}};

        this->set_shape(counter,plaq_up);
        relative_coordinates[nu]=-1;
        add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
        ++counter;

        this->set_shape(counter,plaq_down);
        relative_coordinates[nu]=1;
        add_multiplied_shape(counter,0,plaq_down,relative_coordinates);
        ++counter;
      }
      const std::function<int(const int, const int, const int, const int)> eo_prop_mu
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( eo_prop_mu );
      this->initialize();
    }


    FlowTypeChairDistinctTypeARho1::FlowTypeChairDistinctTypeARho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,1,2,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_distinct_typea_rho1"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 1;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = 1;
          add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down,relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = -1;
          add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down,relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x3+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeChairDistinctTypeARho2::FlowTypeChairDistinctTypeARho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,1,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_distinct_typea_rho2"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 2;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = 1;
          add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down,relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = -1;
          add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down,relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeChairDistinctTypeARho3::FlowTypeChairDistinctTypeARho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,2,1}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_distinct_typea_rho3"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 3;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = 1;
          add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down,relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = -1;
          add_multiplied_shape(counter,0,plaq_up,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down,relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x3 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeChairDistinctTypeBRho1::FlowTypeChairDistinctTypeBRho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,1,2,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_distinct_typeb_rho1"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=4*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 1;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{rho,nu,-rho-1,-nu-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-rho-1,nu,rho,-nu-1},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{nu,rho,-nu-1,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{nu,-rho-1,-nu-1,rho},relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{rho,-nu-1,-rho-1,nu},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-rho-1,-nu-1,rho,nu},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{-nu-1,rho,nu,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-nu-1,-rho-1,nu,rho},relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x3+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeChairDistinctTypeBRho2::FlowTypeChairDistinctTypeBRho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,1,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_distinct_typeb_rho2"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=4*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 2;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{rho,nu,-rho-1,-nu-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-rho-1,nu,rho,-nu-1},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{nu,rho,-nu-1,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{nu,-rho-1,-nu-1,rho},relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{rho,-nu-1,-rho-1,nu},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-rho-1,-nu-1,rho,nu},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{-nu-1,rho,nu,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-nu-1,-rho-1,nu,rho},relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeChairDistinctTypeBRho3::FlowTypeChairDistinctTypeBRho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,2,1}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="chair_distinct_typeb_rho3"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=4*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 3;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{rho,nu,-rho-1,-nu-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-rho-1,nu,rho,-nu-1},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{nu,rho,-nu-1,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{nu,-rho-1,-nu-1,rho},relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{rho,-nu-1,-rho-1,nu},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-rho-1,-nu-1,rho,nu},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{-nu-1,rho,nu,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-nu-1,-rho-1,nu,rho},relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x3 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeChairDistinctInner::FlowTypeChairDistinctInner
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,2,2}}*/,
     const int n_loops_/*=12*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="chair_distinct_inner"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/, // DIVIDE TWO CONTRIBUTIONS!!
     const bool is_flowed_link_nonlinear_/*=true*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;

      const std::array<int, DIMN> relative_coordinates{{0,0,0,0}};

      for(int nu=mu_+1; nu<DIMN; ++nu){
        for(int rho=nu+1; rho<DIMN; ++rho){
          {
            const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
            this->set_shape(counter,plaq_up);
            add_multiplied_shape(counter,0,std::vector<int>{rho,mu_,-rho-1,-mu_-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{-rho-1,mu_,rho,-mu_-1},relative_coordinates);
          }
          ++counter;
          {
            const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
            this->set_shape(counter,plaq_up);
            add_multiplied_shape(counter,0,std::vector<int>{nu,mu_,-nu-1,-mu_-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{-nu-1,mu_,nu,-mu_-1},relative_coordinates);
          }
          ++counter;
          {
            const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
            this->set_shape(counter,plaq_down);
            add_multiplied_shape(counter,0,std::vector<int>{rho,mu_,-rho-1,-mu_-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{-rho-1,mu_,rho,-mu_-1},relative_coordinates);
          }
          ++counter;
          {
            const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
            this->set_shape(counter,plaq_down);
            add_multiplied_shape(counter,0,std::vector<int>{nu,mu_,-nu-1,-mu_-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{-nu-1,mu_,nu,-mu_-1},relative_coordinates);
          }
          ++counter;
        }
      }

      const std::function<int(const int, const int, const int, const int)> eo_prop_mu
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( eo_prop_mu );
      this->initialize();
    }


    // -------------------------------


    FlowTypeTwistSrectDistinctRho1::FlowTypeTwistSrectDistinctRho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,4,1,1}}*/,
     const int n_loops_ /*=2*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="twist_srect_distinct_rho1"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 1;
      {
        const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
        this->set_shape(0,plaq_up);
        const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
        const std::array<int, DIMN> relative_coordinates{{0,1,0,0}};
        add_multiplied_shape(0,0,plaq_up_rev,relative_coordinates);
      }
      {
        const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
        this->set_shape(1,plaq_down);
        const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};
        const std::array<int, DIMN> relative_coordinates{{0,-1,0,0}};
        add_multiplied_shape(1,0,plaq_down_rev,relative_coordinates);
      }
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x2+1024)%4;
        };
      this->set_mask( div4rho );
      this->initialize();
    }

    FlowTypeTwistSrectDistinctRho2::FlowTypeTwistSrectDistinctRho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,1,4,1}}*/,
     const int n_loops_ /*=2*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="twist_srect_distinct_rho2"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 2;
      {
        const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
        this->set_shape(0,plaq_up);
        const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
        const std::array<int, DIMN> relative_coordinates{{0,0,1,0}};
        add_multiplied_shape(0,0,plaq_up_rev,relative_coordinates);
      }
      {
        const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
        this->set_shape(1,plaq_down);
        const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};
        const std::array<int, DIMN> relative_coordinates{{0,0,-1,0}};
        add_multiplied_shape(1,0,plaq_down_rev,relative_coordinates);
      }
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x3+1024)%4;
        };
      this->set_mask( div4rho );
      this->initialize();
    }


    FlowTypeTwistSrectDistinctRho3::FlowTypeTwistSrectDistinctRho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,1,1,4}}*/,
     const int n_loops_ /*=2*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="twist_srect_distinct_rho3"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 3;
      {
        const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
        this->set_shape(0,plaq_up);
        const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
        const std::array<int, DIMN> relative_coordinates{{0,0,0,1}};
        add_multiplied_shape(0,0,plaq_up_rev,relative_coordinates);
      }
      {
        const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
        this->set_shape(1,plaq_down);
        const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};
        const std::array<int, DIMN> relative_coordinates{{0,0,0,-1}};
        add_multiplied_shape(1,0,plaq_down_rev,relative_coordinates);
      }
      const std::function<int(const int, const int, const int, const int)> div4rho
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x4+1024)%4;
        };
      this->set_mask( div4rho );
      this->initialize();
    }


    FlowTypeTwistLrectDistinct::FlowTypeTwistLrectDistinct
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{2,2,2,2}}*/,
     const int n_loops_ /*=6*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=4*/,
     const std::string description_ /*="twist_lrect_distinct"*/,
     const int flow_size_ /*=2*/,
     const int multiplicative_const_ /*=1*/,
     const bool is_flowed_link_nonlinear_ /*=false*/,
     const int number_of_multiplied_loops_ /*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
        this->set_shape(counter,plaq_up);
        const std::vector<int> plaq_up_rev{nu,mu_,-nu-1,-mu_-1};
        std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        relative_coordinates[mu_] = 1;
        add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
        relative_coordinates[mu_] = -1;
        add_multiplied_shape(counter,1,plaq_up_rev,relative_coordinates);
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
        this->set_shape(counter,plaq_down);
        const std::vector<int> plaq_down_rev{-nu-1,mu_,nu,-mu_-1};
        std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
        relative_coordinates[mu_] = 1;
        add_multiplied_shape(counter,0,plaq_down_rev,relative_coordinates);
        relative_coordinates[mu_] = -1;
        add_multiplied_shape(counter,1,plaq_down_rev,relative_coordinates);
        ++counter;
      }
      const std::function<int(const int, const int, const int, const int)> lrect_mask
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x1+1024)%2) + ( x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( lrect_mask );
      this->initialize();
    }


    FlowTypeTwistMrectDistinct::FlowTypeTwistMrectDistinct
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_ /*=std::array<int,DIMN>{{1,2,2,2}}*/,
     const int n_loops_ /*=6*/,
     const int length_of_loop_ /*=4*/,
     const int n_masks_ /*=2*/,
     const std::string description_ /*="twist_mrect_distinct"*/,
     const int flow_size_ /*=1*/,
     const int multiplicative_const_ /*=1*/, // DIVIDE TWO CONTRIBUTIONS!!
     const bool is_flowed_link_nonlinear_ /*=true*/,
     const int number_of_multiplied_loops_ /*=1*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;
      for(int nu=mu_+1; nu<DIMN; ++nu){
        const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
        const std::vector<int> plaq_up_rev{nu,mu_,-nu-1,-mu_-1};
        const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
        const std::vector<int> plaq_down_rev{-nu-1,mu_,nu,-mu_-1};
        std::array<int, DIMN> relative_coordinates{{0,0,0,0}};

        this->set_shape(counter,plaq_up);
        relative_coordinates[nu]=-1;
        add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
        ++counter;

        this->set_shape(counter,plaq_down);
        relative_coordinates[nu]=1;
        add_multiplied_shape(counter,0,plaq_down_rev,relative_coordinates);
        ++counter;
      }
      const std::function<int(const int, const int, const int, const int)> eo_prop_mu
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( eo_prop_mu );
      this->initialize();
    }


    FlowTypeTwistChairDistinctTypeARho1::FlowTypeTwistChairDistinctTypeARho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,1,2,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_distinct_typea_rho1"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 1;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
          const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = 1;
          add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down_rev,relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
          const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = -1;
          add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down_rev,relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x3+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeTwistChairDistinctTypeARho2::FlowTypeTwistChairDistinctTypeARho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,1,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_distinct_typea_rho2"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 2;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
          const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = 1;
          add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down_rev,relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
          const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = -1;
          add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down_rev,relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeTwistChairDistinctTypeARho3::FlowTypeTwistChairDistinctTypeARho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,2,1}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_distinct_typea_rho3"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 3;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
          const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = 1;
          add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down_rev,relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
          const std::vector<int> plaq_up_rev{rho,mu_,-rho-1,-mu_-1};
          const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
          const std::vector<int> plaq_down_rev{-rho-1,mu_,rho,-mu_-1};

          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[nu] = -1;
          add_multiplied_shape(counter,0,plaq_up_rev,relative_coordinates);
          add_multiplied_shape(counter,1,plaq_down_rev,relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x3 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeTwistChairDistinctTypeBRho1::FlowTypeTwistChairDistinctTypeBRho1
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,1,2,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_distinct_typeb_rho1"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=4*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 1;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{nu,rho,-nu-1,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{nu,-rho-1,-nu-1,rho},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{rho,nu,-rho-1,-nu-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-rho-1,nu,rho,-nu-1},relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{-nu-1,rho,nu,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-nu-1,-rho-1,nu,rho},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{rho,-nu-1,-rho-1,nu},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-rho-1,-nu-1,rho,nu},relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x3+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }

    FlowTypeTwistChairDistinctTypeBRho2::FlowTypeTwistChairDistinctTypeBRho2
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,1,2}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_distinct_typeb_rho2"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=4*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 2;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{nu,rho,-nu-1,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{nu,-rho-1,-nu-1,rho},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{rho,nu,-rho-1,-nu-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-rho-1,nu,rho,-nu-1},relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{-nu-1,rho,nu,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-nu-1,-rho-1,nu,rho},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{rho,-nu-1,-rho-1,nu},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-rho-1,-nu-1,rho,nu},relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x4 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeTwistChairDistinctTypeBRho3::FlowTypeTwistChairDistinctTypeBRho3
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,2,1}}*/,
     const int n_loops_/*=4*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=4*/,
     const std::string description_/*="twist_chair_distinct_typeb_rho3"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/,
     const bool is_flowed_link_nonlinear_/*=false*/,
     const int number_of_multiplied_loops_/*=4*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      const int rho = 3;
      int counter = 0;

      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
          this->set_shape(counter,plaq_up);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{nu,rho,-nu-1,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{nu,-rho-1,-nu-1,rho},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{rho,nu,-rho-1,-nu-1},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-rho-1,nu,rho,-nu-1},relative_coordinates);
        }
        ++counter;
      }
      for(int nu=mu_+1; nu<DIMN; ++nu){
        if(nu==rho) continue;
        {
          const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
          this->set_shape(counter,plaq_down);
        }
        {
          std::array<int, DIMN> relative_coordinates{{0,0,0,0}};
          relative_coordinates[mu_] = 1;
          add_multiplied_shape(counter,0,std::vector<int>{-nu-1,rho,nu,-rho-1},relative_coordinates);
          add_multiplied_shape(counter,1,std::vector<int>{-nu-1,-rho-1,nu,rho},relative_coordinates);
          relative_coordinates[mu_] = 0;
          add_multiplied_shape(counter,2,std::vector<int>{rho,-nu-1,-rho-1,nu},relative_coordinates);
          add_multiplied_shape(counter,3,std::vector<int>{-rho-1,-nu-1,rho,nu},relative_coordinates);
        }
        ++counter;
      }

      const std::function<int(const int, const int, const int, const int)> mask_mu_rho_plane
        = [](const int x1, const int x2, const int x3, const int x4) {
          return 2*( (x2+1024)%2) + ( x3 + 1024 )%2;
        };
      this->set_mask( mask_mu_rho_plane );
      this->initialize();
    }


    FlowTypeTwistChairDistinctInner::FlowTypeTwistChairDistinctInner
    (
     const int mu_/*=0*/,
     const std::array<int,DIMN> dimensions_/*=std::array<int,DIMN>{{1,2,2,2}}*/,
     const int n_loops_/*=12*/,
     const int length_of_loop_/*=4*/,
     const int n_masks_/*=2*/,
     const std::string description_/*="twist_chair_distinct_inner"*/,
     const int flow_size_/*=1*/,
     const int multiplicative_const_/*=1*/, // DIVIDE TWO CONTRIBUTIONS!!
     const bool is_flowed_link_nonlinear_/*=true*/,
     const int number_of_multiplied_loops_/*=2*/
     )
      : UnitInfo(dimensions_)
      , FlowType(dimensions_,n_loops_,length_of_loop_,n_masks_,description_,flow_size_,
                 multiplicative_const_,is_flowed_link_nonlinear_,number_of_multiplied_loops_)
    {
      int counter = 0;

      const std::array<int, DIMN> relative_coordinates{{0,0,0,0}};

      for(int nu=mu_+1; nu<DIMN; ++nu){
        for(int rho=nu+1; rho<DIMN; ++rho){
          {
            const std::vector<int> plaq_up{mu_,nu,-mu_-1,-nu-1};
            this->set_shape(counter,plaq_up);
            add_multiplied_shape(counter,0,std::vector<int>{mu_,rho,-mu_-1,-rho-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{mu_,-rho-1,-mu_-1,rho},relative_coordinates);
          }
          ++counter;
          {
            const std::vector<int> plaq_up{mu_,rho,-mu_-1,-rho-1};
            this->set_shape(counter,plaq_up);
            add_multiplied_shape(counter,0,std::vector<int>{mu_,nu,-mu_-1,-nu-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{mu_,-nu-1,-mu_-1,nu},relative_coordinates);
          }
          ++counter;
          {
            const std::vector<int> plaq_down{mu_,-nu-1,-mu_-1,nu};
            this->set_shape(counter,plaq_down);
            add_multiplied_shape(counter,0,std::vector<int>{mu_,rho,-mu_-1,-rho-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{mu_,-rho-1,-mu_-1,rho},relative_coordinates);
          }
          ++counter;
          {
            const std::vector<int> plaq_down{mu_,-rho-1,-mu_-1,rho};
            this->set_shape(counter,plaq_down);
            add_multiplied_shape(counter,0,std::vector<int>{mu_,nu,-mu_-1,-nu-1},relative_coordinates);
            add_multiplied_shape(counter,1,std::vector<int>{mu_,-nu-1,-mu_-1,nu},relative_coordinates);
          }
          ++counter;
        }
      }

      const std::function<int(const int, const int, const int, const int)> eo_prop_mu
        = [](const int x1, const int x2, const int x3, const int x4) {
          return (x2 + x3 + x4 + 1024)%2;
        };
      this->set_mask( eo_prop_mu );
      this->initialize();
    }




  }
}


#endif
