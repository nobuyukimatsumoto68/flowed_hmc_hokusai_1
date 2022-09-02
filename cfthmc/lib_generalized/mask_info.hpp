#ifndef MASK_INFO_HPP
#define MASK_INFO_HPP

// mask_type: mask = 0, 1, 2, ..., n_masks-1;
#include <vector>
#include <array>
#include <cassert>
#include <functional>

#include "unit_info.hpp"

namespace qlat{
  namespace Generalized{

    class MaskInfo : public virtual UnitInfo {
    public:
      const int m_unit_size;
      const int m_n_masks; // number of mask types
    protected:
      std::vector<int> m_mask_pattern; // mask pattern in fund unit
      bool m_is_initialized; // true after properly initialized

    public:
      MaskInfo() = delete;
      MaskInfo
      (
       const std::array<int, DIMN>& dimensions_,
       const int unit_size_,
       const int n_masks_
       )
        : UnitInfo(dimensions_)
        , m_unit_size(unit_size_)
        , m_n_masks(n_masks_)
        , m_mask_pattern( m_unit_size )
        , m_is_initialized(false)
      {
        for(std::size_t i=0; i<m_mask_pattern.size(); ++i) m_mask_pattern[i] = n_masks_;
      }

      // SET
      void set_mask(const std::array<int, DIMN>& x_unit, const int mask_type) & {
        assert(0<= mask_type && mask_type<m_n_masks);
        m_mask_pattern[unit2idx(x_unit)] = mask_type;
      }

      void set_mask
      ( const std::function<int(const int, const int, const int, const int)>& masking_rule ) & {
        for(int x1=0; x1<m_dimensions[0]; ++x1){
          for(int x2=0; x2<m_dimensions[1]; ++x2){
            for(int x3=0; x3<m_dimensions[2]; ++x3){
              for(int x4=0; x4<m_dimensions[3]; ++x4){
                set_mask( std::array<int, DIMN>{{x1,x2,x3,x4}}, masking_rule(x1,x2,x3,x4) );
              }}}}
      }

      // PRINT
      const int& mask( const int x_idx ) const& {
        assert(is_initialized());
        return m_mask_pattern[x_idx];
      }
      const int& mask( const std::array<int, DIMN>& x_unit ) const& { return mask(unit2idx(x_unit)); }
      const bool& is_initialized() const& { return m_is_initialized; }
      //
      bool initialize() &;

    };

    bool MaskInfo::initialize() & {
      for(int mask=0; mask<m_n_masks; ++mask)
        if(m_mask_pattern[mask]==m_n_masks) return m_is_initialized;
      m_is_initialized = true;
      return m_is_initialized;
    }

  }
}

#endif
