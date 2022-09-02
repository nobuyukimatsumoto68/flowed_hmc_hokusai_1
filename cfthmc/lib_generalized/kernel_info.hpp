#ifndef KERNEL_INFO_HPP
#define KERNEL_INFO_HPP

// INDEX RULES
// loops: l = 0, 1, ..., n_loops; l stands for "loop"
// sites along a loop: s = 0, 1, ..., length_of_loop; s stands for "site"

#include <array>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>

#include "unit_info.hpp"

namespace qlat{
  namespace Generalized{

    class KernelInfo : public virtual UnitInfo {
    public:
      static constexpr int m_mu = 0; // flowed direction; other directions can be mapped with symmetry
    protected:
      const int m_n_loops; // number of loops starting from a flowed link
      const int m_length_of_loop; // length of a loop (including a flowd link)
      bool m_is_initialized; // true after properly initialized
      std::vector<std::vector<int> > m_shapes;
      // ------
      // for separate two loops;
      // there can be multiple loops multiplied to a single staple,
      // but they are assumed to have the length (m_length_of_loop)
      const int m_number_of_multiplied_loops;
      std::vector<std::vector<std::vector<int>>> m_multiplied_shapes;
      std::vector<std::vector<std::array<int, DIMN>>> m_multiplied_shapes_relative_coordinates; // assuming flowed direction = 0
      /* shapes[i] determines a shape of i-th loop;
         e.g., shapes[i] = {+0,+1,-1,-2} for a plaquette. We do not include the reversed loop.
         It is assumed that we assign loops with common shapes at each site. */
    public:
      KernelInfo() = delete;
      KernelInfo
      (
       const std::array<int, DIMN>& dimensions_,
       const int n_loops_,
       const int length_of_loop_,
       const int number_of_multiplied_loops_=0
       )
        : UnitInfo(dimensions_)
        , m_n_loops(n_loops_)
        , m_length_of_loop(length_of_loop_)
        , m_is_initialized(false)
        , m_shapes(m_n_loops, std::vector<int>(m_length_of_loop, 0.0) )
        , m_number_of_multiplied_loops(number_of_multiplied_loops_)
        , m_multiplied_shapes( m_n_loops,
                               std::vector<std::vector<int>>(m_number_of_multiplied_loops,
                                                             std::vector<int>(m_length_of_loop) ) )
        , m_multiplied_shapes_relative_coordinates( m_n_loops,
                                                    std::vector<std::array<int, DIMN>>(m_number_of_multiplied_loops) )
      {}

      // SET
      void set_shape(const int l, const std::vector<int>& shape) & {
        assert(0<=l && l<m_n_loops);
        m_shapes[l] = shape;
      }
      void set_shape(const int l, const std::vector<int>&& shape) & {
        assert(0<=l && l<m_n_loops);
        m_shapes[l] = shape;
      }
      void add_multiplied_shape
      (
       const int l,
       const int i,
       const std::vector<int>& shape,
       const std::array<int, DIMN>& relative_coordinates
       ) & {
        assert(0<=l && l<m_n_loops);
        assert(0<=i && i<m_number_of_multiplied_loops);
        m_multiplied_shapes[l][i] = shape;
        m_multiplied_shapes_relative_coordinates[l][i] = relative_coordinates;
      }
      bool initialize() &;
      // PRINT
      const std::vector<int>& shape(const int l) const& {
        assert(0<=l && l<m_n_loops);
        return m_shapes[l];
      }
      const std::vector<int>& multiplied_shape(const int l, const int i) const& {
        assert(0<=l && l<m_n_loops);
        assert(0<=i && i<m_number_of_multiplied_loops);
        return m_multiplied_shapes[l][i];
      }
      int multiplied_shapes_relative_coordinates(const int l, const int i, const int rho) const& {
        assert(0<=l && l<m_n_loops);
        assert(0<=i && i<m_number_of_multiplied_loops);
        assert(-DIMN<=rho && rho<DIMN);
        return m_multiplied_shapes_relative_coordinates[l][i][rho];
      }
      bool is_initialized() const& { return m_is_initialized; }
      int number_of_multiplied_loops() const& { return m_number_of_multiplied_loops; }

      // CONVERSIONS
      int ls2r(const int l, const int s) const& { return l*m_length_of_loop + s; }
      void ls2r(int& r, const int l, const int s) const& { r = ls2r(l,s); }
      void r2ls(int& l, int& s, const int r) const& {
        this->quotient_remainder(l,s,r,m_length_of_loop);
      }

    };


    bool KernelInfo::initialize() & {
      assert(!m_is_initialized);

      assert(static_cast<int>(m_multiplied_shapes.size())==m_n_loops);
      assert(static_cast<int>(m_multiplied_shapes_relative_coordinates.size())==m_n_loops);

      for(int l=0; l<m_n_loops; ++l) {
        assert(static_cast<int>(m_multiplied_shapes[l].size())==m_number_of_multiplied_loops);
        assert(static_cast<int>(m_multiplied_shapes_relative_coordinates[l].size())==m_number_of_multiplied_loops);
        for(const std::vector<int>& shape : m_multiplied_shapes[l]) assert( static_cast<int>(shape.size())==m_length_of_loop );
        //
        int tmp = 0;
        assert(static_cast<int>(m_shapes[l].size())==m_length_of_loop);
        for(int s=0; s<m_length_of_loop; ++s) {
          int dir = m_shapes[l][s];
          if(dir>=0) dir+=1;
          tmp += dir;
        }
        // if(tmp==0) return m_is_initialized;
        if(tmp!=0) assert(false); // return m_is_initialized;
      }

      m_is_initialized = true;
      return m_is_initialized;
    }
  }
}

#endif
