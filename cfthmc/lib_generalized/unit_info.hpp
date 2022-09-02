#ifndef UNIT_INFO_HPP
#define UNIT_INFO_HPP

#include <array>
#include <iostream>
#include <cassert>

namespace qlat{
  namespace Generalized{

    class UnitInfo {
    public:
      // static constexpr int DIMN = 4;
      static constexpr int DIMN = qlat::DIMN;

      const std::array<int, DIMN> m_dimensions; // dimenstion of the unit; e.g., {{2,2,2,2}}
      const int m_unit_size;
    private:
      std::array<int, DIMN> m_lattice_size;
      std::array<int, DIMN> m_n_unit;
      bool m_is_lattice_size_set;

    public:
      UnitInfo() = delete;
      UnitInfo( const std::array<int, DIMN>& dimensions_ )
        : m_dimensions(dimensions_)
        , m_unit_size( set_unit_size() )
        , m_is_lattice_size_set(false)
      {}

      const int& operator[](const int mu) const& { return m_dimensions[mu]; }
      const int& size() const& { return m_unit_size; }

      inline void shift(std::array<int, DIMN>& x_global, const int mu) const&;

      // ROTATIONS
      int rotate( const int nu, const int mu) const&;
      void rotate( int& nu_prime, const int nu, const int mu ) const& { nu_prime = rotate(nu,mu); }

      void rotate( std::array<int, DIMN>& x_prime, const std::array<int, DIMN>& x,
                   const int mu ) const&;
      std::array<int, DIMN> rotate( const std::array<int, DIMN>& x, const int mu ) const&;

      // CONVERSIONS
      void global2unit( std::array<int, DIMN>& x_unit, std::array<int, DIMN>& x_block,
                        const std::array<int, DIMN>& x_global ) const&;
      void global2unit( int& x_idx, std::array<int, DIMN>& x_block,
                        const std::array<int, DIMN>& x_global ) const&;

      inline void unit2global( std::array<int, DIMN>& x_global,
                               const std::array<int, DIMN>& x_unit,
                               const std::array<int, DIMN>& x_block ) const&;
      void unit2global( std::array<int, DIMN>& x_global,
                        const int x_idx,
                        const std::array<int, DIMN>& x_block ) const&;
      std::array<int, DIMN> unit2global( const int x_idx,
                                         const std::array<int, DIMN>& x_block ) const;

      void unit2idx( int& x_idx, const std::array<int, DIMN>& x_unit ) const&;
      int unit2idx( const std::array<int, DIMN>& x_unit) const&;

      void idx2unit( std::array<int, DIMN>& x_unit, const int& x_idx ) const&;
      std::array<int, DIMN> idx2unit( const int& x_idx) const&;

      void quotient_remainder
      (int& quotient, int& remainder, const int x, const int y, const int safety=4) const&;

    private:
      int set_unit_size() const&;
      std::array<int, DIMN> set_n_unit() const&;
      std::array<int, DIMN> set_n_unit( const std::array<int, DIMN>& lattice_size ) const&;
    };


    inline void UnitInfo::shift(std::array<int, DIMN>& x_global, const int mu) const& {
      if(0<=mu && mu<DIMN) ++x_global[mu];
      else if(-DIMN<=mu && mu<0) --x_global[-mu-1];
      else assert(false);
    }

    int UnitInfo::rotate( const int nu, const int mu) const& {
      if(nu>=0 && mu>=0) return (nu+mu)%UnitInfo::DIMN;
      else if(nu>=0 && mu<0) {
        return (nu+mu+1+UnitInfo::DIMN)%UnitInfo::DIMN;
      }
      else if(nu<0 && mu>=0) {
        const int tmp = (-nu-1+mu+UnitInfo::DIMN)%UnitInfo::DIMN;
        return -tmp-1;
      }
      else { // nu<0 && mu<0
        const int tmp = (-nu-1+mu+1+UnitInfo::DIMN)%UnitInfo::DIMN;
        return -tmp-1;
      }
    }

    void UnitInfo::rotate
    (
     std::array<int, UnitInfo::DIMN>& x_prime,
     const std::array<int, UnitInfo::DIMN>& x,
     const int mu
     ) const&
    {
      if(mu>=0) {
        assert(0<=mu && mu<UnitInfo::DIMN);
        for(std::size_t rho=0; rho<x.size(); ++rho) x_prime[(rho+mu)%UnitInfo::DIMN] = x[rho];
      }
      else{
        assert(-UnitInfo::DIMN<=mu && mu<0);
        for(std::size_t rho=0; rho<x.size(); ++rho)
          x_prime[(rho+UnitInfo::DIMN+mu+1)%UnitInfo::DIMN] = x[rho];
      }
    }

    std::array<int, UnitInfo::DIMN> UnitInfo::rotate
    (
     const std::array<int, UnitInfo::DIMN>& x,
     const int mu
     ) const&
    {
      std::array<int, UnitInfo::DIMN> x_prime;
      rotate( x_prime, x, mu );
      return x_prime;
    }


    void UnitInfo::global2unit
    (
     std::array<int, UnitInfo::DIMN>& x_unit,
     std::array<int, UnitInfo::DIMN>& x_block,
     const std::array<int, UnitInfo::DIMN>& x_global
     ) const&
    {
      for(int mu=0; mu<DIMN; ++mu) {
        quotient_remainder(x_block[mu],x_unit[mu],x_global[mu],m_dimensions[mu]);
      }
    }

    void UnitInfo::global2unit
    (
     int& x_idx, std::array<int, UnitInfo::DIMN>& x_block,
     const std::array<int, UnitInfo::DIMN>& x_global
     ) const&
    {
      std::array<int, DIMN> x_unit;
      global2unit(x_unit,x_block,x_global);
      unit2idx(x_idx, x_unit);
    }

    inline void UnitInfo::unit2global
    (
     std::array<int, UnitInfo::DIMN>& x_global,
     const std::array<int, UnitInfo::DIMN>& x_unit,
     const std::array<int, UnitInfo::DIMN>& x_block
     ) const&
    {
      for(int mu=0; mu<DIMN; ++mu) x_global[mu] = x_block[mu]*m_dimensions[mu] + x_unit[mu];
    }

    void UnitInfo::unit2global
    (
     std::array<int, UnitInfo::DIMN>& x_global,
     const int x_idx,
     const std::array<int, UnitInfo::DIMN>& x_block
     ) const&
    {
      const std::array<int, DIMN> x_unit = idx2unit(x_idx);
      unit2global(x_global,x_unit,x_block);
    }

    std::array<int, UnitInfo::DIMN> UnitInfo::unit2global
    (
     const int x_idx,
     const std::array<int, UnitInfo::DIMN>& x_block
     ) const
    {
      std::array<int, DIMN> x_global;
      unit2global(x_global,x_idx,x_block);
      return x_global;
    }


    void UnitInfo::unit2idx( int& x_idx, const std::array<int, UnitInfo::DIMN>& x_unit ) const& {
      x_idx =0;
      for(int mu=0; mu<DIMN; ++mu) {
        assert(0<=x_unit[mu] && x_unit[mu]<m_dimensions[mu]);
        x_idx *= m_dimensions[mu];
        x_idx += x_unit[mu];
      }
    }

    int UnitInfo::unit2idx( const std::array<int, UnitInfo::DIMN>& x_unit) const& {
      int x_idx = 0;
      unit2idx(x_idx,x_unit);
      return x_idx;
    }

    void UnitInfo::idx2unit( std::array<int, UnitInfo::DIMN>& x_unit, const int& x_idx ) const& {
      int tmp1=x_idx;
      int tmp2=x_idx;
      for(int mu=DIMN-1; mu>=0; --mu) {
        quotient_remainder(tmp1,x_unit[mu],tmp2,m_dimensions[mu]);
        tmp2 = tmp1;
      }
    }

    std::array<int, UnitInfo::DIMN> UnitInfo::idx2unit( const int& x_idx) const& {
      std::array<int, DIMN> x_unit;
      idx2unit(x_unit,x_idx);
      return x_unit;
    }


    void UnitInfo::quotient_remainder
    (
     int& quotient,
     int& remainder,
     const int x,
     const int y, // >0
     const int safety // =4 >0
     ) const& {
      assert(y>0);
      assert(x+safety*y>=0);
      remainder = (x+safety*y)%y;
      quotient = (x-remainder)/y;
    }


    int UnitInfo::set_unit_size() const& {
      int res = 1;
      for(int mu=0; mu<DIMN; ++mu) res *= m_dimensions[mu];
      return res;
    }

    std::array<int, UnitInfo::DIMN> UnitInfo::set_n_unit() const& {
      std::array<int, DIMN> res;
      for(int mu=0; mu<DIMN; ++mu) {
        res[mu] = m_lattice_size[mu]/m_dimensions[mu];
        assert(m_lattice_size[mu] == res[mu]*m_dimensions[mu]);
      }
      return res;
    }

    std::array<int, UnitInfo::DIMN> UnitInfo::set_n_unit
    (
     const std::array<int, UnitInfo::DIMN>& lattice_size
     ) const&
    {
      std::array<int, DIMN> res;
      for(int mu=0; mu<DIMN; ++mu) {
        res[mu] = lattice_size[mu]/m_dimensions[mu];
        assert(lattice_size[mu] == res[mu]*m_dimensions[mu]);
      }
      return res;
    }

  }
}

#endif


// inline void unit2global_wrapped( std::array<int, DIMN>& x_global,
//                                  const std::array<int, DIMN>& x_unit,
//                                  const std::array<int, DIMN>& x_block,
//                                  const std::array<int, DIMN>& lattice_size) const {
//   std::array<int, DIMN> n_unit = this->set_n_unit(lattice_size);
//   for(int mu=0; mu<DIMN; ++mu) {
//     int x_block_mu = x_block[mu];
//     while(x_block_mu<0) x_block_mu += n_unit[mu];
//     x_block_mu %= n_unit[mu];
//     x_global[mu] = std::abs(x_block_mu)*m_dimensions[mu] + x_unit[mu];
//   }
// }

// inline std::array<int, DIMN> unit2global_wrapped( const std::array<int, DIMN>& x_unit,
//                                                   const std::array<int, DIMN>& x_block,
//                                                   const std::array<int, DIMN>& lattice_size
//                                                   ) const {
//   std::array<int, DIMN> x_global;
//   unit2global_wrapped(x_global,x_unit,x_block,lattice_size);
//   return x_global;
// }

// inline void unit2global_wrapped( std::array<int, DIMN>& x_global,
//                                  const int x_idx,
//                                  const std::array<int, DIMN>& x_block,
//                                  const std::array<int, DIMN>& lattice_size) const {
//   const std::array<int, DIMN> x_unit = idx2unit(x_idx);
//   unit2global_wrapped(x_global,x_unit,x_block,lattice_size);
// }

// inline std::array<int, DIMN> unit2global_wrapped( const int x_idx,
//                                                   const std::array<int, DIMN>& x_block,
//                                                   const std::array<int, DIMN>& lattice_size
//                                                   ) const {
//   std::array<int, DIMN> x_global;
//   unit2global_wrapped(x_global,x_idx,x_block,lattice_size);
//   return x_global;
// }
