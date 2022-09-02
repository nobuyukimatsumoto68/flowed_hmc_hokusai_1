#ifndef FLOW_TYPE_HPP
#define FLOW_TYPE_HPP

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

namespace qlat{
  namespace Generalized{

    class RowOfTable{
    public:
      static constexpr int kCol = 15;
      //
      static constexpr int kXIdx= 0;
      static constexpr int kXMask= 1;
      static constexpr int kYIdx= 2;
      static constexpr int kNu= 3;
      static constexpr int kYMask= 4;
      static constexpr int kYBlockBeg= 5; // 5--8 : block
      static constexpr int kIsPositive= 9;
      static constexpr int kL= 10;
      static constexpr int kS= 11;
      static constexpr int kMPrime= 12;
      static constexpr int kNext= 13;
      static constexpr int kWhichLoop= 14;
    private:
      std::array<int, kCol> row;

    public:
      RowOfTable()
        : row(std::array<int, kCol>{{}})
      {}
      RowOfTable
      (
       const int x_idx_,
       const int x_mask_,
       const int y_idx_,
       const int nu,
       const int y_mask_,
       const std::array<int, UnitInfo::DIMN>& y_block_,
       const int is_positive,
       const int l,
       const int s,
       const int which_loop = -1
       )
        : row (std::array<int, kCol>{{
            x_idx_, x_mask_, // 0-1
            y_idx_, nu, y_mask_, // 2-4
            y_block_[0], y_block_[1], y_block_[2], y_block_[3], // 5-8
            is_positive, l, s, // 9-11
            0, 0, // 12-13
            which_loop }} // 14
          )
      {}
      ~RowOfTable(){};

      bool operator< (const RowOfTable& another) const& ;

      int& operator[](const int i) & { return row[i]; }
      const int& operator[](const int i) const& { return row[i]; }

      const int& x_idx() const& { return row[kXIdx]; }
      const int& x_mask() const& { return row[kXMask]; }
      const int& y_idx() const& { return row[kYIdx]; }
      const int& nu() const& { return row[kNu]; }
      const int& y_mask() const& { return row[kYMask]; }
      const int& y_block( const int mu ) const& { return row[kYBlockBeg+mu]; }
      const int& is_positive() const& { return row[kIsPositive]; }
      const int& l() const& { return row[kL]; }
      const int& s() const& { return row[kS]; }
      int& m_prime() & { return row[kMPrime]; }
      const int& m_prime() const& { return row[kMPrime]; }
      int& next() & { return row[kNext]; }
      const int& next() const& { return row[kNext]; }
      const int& which_loop() const& { return row[kWhichLoop]; }

      void print(std::ostream& os) const& {
        os << "x: idx = " << std::setw(2) << x_idx() << ", "
           << "mask = " << std::setw(2)<< x_mask() << ", "
           << "y: idx = " << std::setw(2)<< y_idx() << ", "
           << "nu = " << std::setw(2)<< nu() << ", "
           << "mask = " << std::setw(2)<< y_mask() << ", "
           << "block = "
           << std::setw(2) << y_block(0) << ", "
           << std::setw(2) << y_block(1) << ", "
           << std::setw(2) << y_block(2) << ", "
           << std::setw(2) << y_block(3) << ", "
           << "is_positive = " << std::setw(2) << is_positive() << ", "
           << "l = " << std::setw(2) << l() << ", "
           << "s = " << std::setw(2) << s() << ", "
           << "m_prime = " << std::setw(2) << m_prime() << ", "
           << "next = " << std::setw(2) << next() << ", "
           << "which_loop = " << std::setw(2) << which_loop();
      }

      std::array<int, kCol>::const_iterator begin() const& { return row.begin(); }
      std::array<int, kCol>::const_iterator end() const& { return row.end(); }

    };


    class FlowType : public KernelInfo, public MaskInfo {
    public:
      const int unit_size; // total number of sites in a unit
      const int n_loops; // number of loops starting from a flowed link
      const int length_of_loop; // length of a loop (including a flowd link)
      const int n_masks; // number of mask types
      const std::string description;
      const int flow_size;
      const int multiplicative_const;
      const bool is_flowed_link_nonlinear; // up to two at present [ only affects "this->calculate_bound_inverse()" ]

    protected:
      // x in fund unit, y may not.
      std::vector<RowOfTable> m_table;
      std::vector<int> m_multiplicity_mask_yidx;
      std::vector<int> m_begin_mask_yidx;
      int m_multiplicity_max;
      std::vector<int> m_multiplicity_reduced_mask_yidx;
      int m_multiplicity_reduced_max;

    public:
      const double bound_inverse;

      FlowType
      (
       const std::array<int, DIMN>& dimensions_=std::array<int, DIMN>{{0,0,0,0}},
       const int n_loops_=0,
       const int length_of_loop_=0,
       const int n_masks_=0,
       const std::string description_="",
       const int flow_size_=0,
       const int multiplicative_const_=1, // default. meaningful value
       const bool is_flowed_link_nonlinear_=false,
       const int number_of_multiplied_loops_=0 // separate cases
       )
        : UnitInfo(dimensions_)
        , KernelInfo(dimensions_,n_loops_,length_of_loop_,number_of_multiplied_loops_)
        , MaskInfo(dimensions_,calculate_unit_size(dimensions_),n_masks_)
        , unit_size( calculate_unit_size(dimensions_) )
        , n_loops(n_loops_)
        , length_of_loop(length_of_loop_)
        , n_masks(n_masks_)
        , description( description_ )
        , flow_size( flow_size_ )
        , multiplicative_const( multiplicative_const_ )
        , is_flowed_link_nonlinear(is_flowed_link_nonlinear_)
        , m_table( calculate_unit_size(dimensions_)*length_of_loop_*n_loops_*(1+number_of_multiplied_loops_) )
        , m_multiplicity_mask_yidx( n_masks_*calculate_unit_size(dimensions_) )
        , m_begin_mask_yidx( n_masks_*calculate_unit_size(dimensions_) )
        , m_multiplicity_reduced_mask_yidx( n_masks_*calculate_unit_size(dimensions_) )
        , bound_inverse( this->calculate_bound_inverse() )
      {
        assert( typeid(dimensions_)==typeid(std::array<int,DIMN>) );
        assert( typeid(n_loops_)==typeid(int) );
        assert( typeid(length_of_loop_)==typeid(int) );
        assert( typeid(n_masks_)==typeid(int) );
        assert( typeid(description_)==typeid(std::string) );
        assert( typeid(flow_size_)==typeid(int) );
        assert( typeid(multiplicative_const_)==typeid(int) );
        assert( typeid(is_flowed_link_nonlinear_)==typeid(bool) );
        assert( typeid(number_of_multiplied_loops_)==typeid(int) );
      }

      std::vector<RowOfTable>::const_iterator begin() const& { return m_table.begin(); }
      std::vector<RowOfTable>::const_iterator end() const& { return m_table.end(); }

      // SET
      const RowOfTable& table(const int i) const& {
        assert(i<static_cast<int>(m_table.size()));
        return m_table[i];
      }
      bool initialize() &; // override

      // PRINT
      bool is_initialized() const& { // override
        return ( KernelInfo::is_initialized() && MaskInfo::is_initialized() );
      }
      const std::vector<RowOfTable>& table() const& {
        return m_table;
      }
      const int& table(const int col, const int row) const& { return m_table[col][row]; }
      const int& multiplicity(const int mask, const int y_idx) const& {
        return m_multiplicity_mask_yidx[mask*m_unit_size+y_idx];
      }
      const int& multiplicity(const int mask, const std::array<int, DIMN> y_unit) const& {
        return multiplicity(mask,unit2idx(y_unit));
      }
      const int& multiplicity_max() const& { return m_multiplicity_max; }
      const int& multiplicity_reduced(const int mask, const int y_idx) const& {
        return m_multiplicity_reduced_mask_yidx[mask*m_unit_size+y_idx];
      }
      const int& multiplicity_reduced(const int mask, const std::array<int, DIMN> y_unit) const& {
        return multiplicity_reduced(mask,unit2idx(y_unit));
      }
      const int& multiplicity_reduced_max() const& { return m_multiplicity_reduced_max; }
      const int& begin_mask_yidx(const int mask, const int y_idx) const& {
        return m_begin_mask_yidx[mask*m_unit_size+y_idx];
      }

      void print(std::ostream& os) const& {
        for(auto itr=this->begin(); itr!=this->end(); ++itr){
          itr->print(os);
          os << std::endl;
        }
      }

      double calculate_bound_inverse() const&;

    protected:
      int calculate_unit_size( const std::array<int, DIMN>& dimension ) const& {
        int res = 1;
        for(std::size_t i=0; i<dimension.size(); ++i) res *= dimension[i];
        return res;
      }
      void calculate_table() &;
      void calculate_multiplicity() &;
      void calculate_begin_mask_yidx() &;
      void supply_m_prime() &;
      void supply_next() &;
      void calculate_multiplicity_reduced() &;
    };


    bool RowOfTable::operator<
    (const RowOfTable& another) const& {

      int i=kXMask; // x_mask
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kYIdx; // y_idx
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kNu; // nu
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kXIdx; // x_idx
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      for(int rho=0; rho<UnitInfo::DIMN; ++rho){
        i=kYBlockBeg+rho; // y_blocks
        if( (*this)[i] < another[i] ) return true;
        else if( (*this)[i] > another[i] ) return false;
      }

      i=kIsPositive; // is_positive
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kYMask; // y_mask
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kMPrime; // m_prime
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kS; // s
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kL; // l
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kWhichLoop; // which_loop
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      i=kNext; // next
      if( (*this)[i] < another[i] ) return true;
      else if( (*this)[i] > another[i] ) return false;

      return false;
    }


    void FlowType::calculate_table() & {
      int counter = 0;

      for(int x_idx=0; x_idx<m_unit_size; ++x_idx){
        const std::array<int, DIMN> x_unit = this->idx2unit(x_idx);

        for(int l=0; l<m_n_loops; ++l){
          std::array<int, DIMN> y_global;
          for(int rho=0; rho<DIMN; rho++) y_global[rho] = x_unit[rho];

          for(int s=0; s<m_length_of_loop; ++s){
            // const int r = this->ls2r(l,s);
            const int nu = this->m_shapes[l][s];
            if( 0<=nu && nu<DIMN ) {
              int y_idx;
              std::array<int, DIMN> y_block;
              this->global2unit(y_idx, y_block, y_global);

              //m_table[x_idx*m_length_of_loop*m_n_loops + r]
              m_table[counter] = RowOfTable(x_idx,this->mask(x_idx),
                                            y_idx,nu,this->mask(y_idx),y_block,
                                            1,l,s);
              ++counter;
              this->shift( y_global, nu );
            }
            else if( -DIMN<=nu && nu<0 ){
              this->shift( y_global, nu );
              int y_idx;
              std::array<int, DIMN> y_block;
              this->global2unit(y_idx, y_block, y_global);

              //m_table[x_idx*m_length_of_loop*m_n_loops + r]
              m_table[counter] = RowOfTable(x_idx,this->mask(x_idx),
                                            y_idx,-nu-1,this->mask(y_idx),y_block,
                                            0,l,s);
              ++counter;
            }
            else assert(false);
          } // end for s

          for(int i=0; i<this->m_number_of_multiplied_loops; ++i){
            for(int rho=0; rho<DIMN; rho++) y_global[rho] = x_unit[rho] + m_multiplied_shapes_relative_coordinates[l][i][rho];

            for(int s=0; s<m_length_of_loop; ++s){
              // const int r = this->ls2r(l,s);
              const int nu = this->m_multiplied_shapes[l][i][s];
              if( 0<=nu && nu<DIMN ) {
                int y_idx;
                std::array<int, DIMN> y_block;
                this->global2unit(y_idx, y_block, y_global);

                m_table[counter] = RowOfTable(x_idx,this->mask(x_idx),
                                              y_idx,nu,this->mask(y_idx),y_block,
                                              1,l,s,i);
                ++counter;
                this->shift( y_global, nu );
              }
              else if( -DIMN<=nu && nu<0 ){
                this->shift( y_global, nu );
                int y_idx;
                std::array<int, DIMN> y_block;
                this->global2unit(y_idx, y_block, y_global);

                m_table[counter] = RowOfTable(x_idx,this->mask(x_idx),
                                              y_idx,-nu-1,this->mask(y_idx),y_block,
                                              0,l,s,i);
                ++counter;
              }
              else assert(false);

            } // end for s
          } // end for i
        } // end for l
      } // end for x_idx
      assert(counter==static_cast<int>(m_table.size()));
    }

    bool FlowType::initialize() & {
      const bool is_initialized_kernel = KernelInfo::initialize();
      const bool is_initialized_mask = MaskInfo::initialize();
      const bool flag = ( is_initialized_kernel && is_initialized_mask ); // 0210
      assert(flag);
      calculate_table();
      std::sort(m_table.begin(),m_table.end());
      calculate_multiplicity();
      calculate_begin_mask_yidx();
      supply_m_prime();
      supply_next();
      calculate_multiplicity_reduced();
      return flag;
    }

    void FlowType::calculate_multiplicity() & {
      for(std::size_t i=0; i<m_multiplicity_mask_yidx.size(); ++i){
        m_multiplicity_mask_yidx[i] = 0;
      }
      for(auto itr=m_table.begin(); itr!=m_table.end(); ++itr){
        m_multiplicity_mask_yidx[ itr->x_mask()*m_unit_size + itr->y_idx() ] += 1;
      }

      auto itr_max = max_element(m_multiplicity_mask_yidx.begin(),
                                 m_multiplicity_mask_yidx.end());
      m_multiplicity_max = *itr_max;
    }

    void FlowType::calculate_begin_mask_yidx() & {
      m_begin_mask_yidx[0] = 0;
      for(std::size_t i=0; i<m_multiplicity_mask_yidx.size()-1; ++i){
        m_begin_mask_yidx[i+1] = m_begin_mask_yidx[i] + m_multiplicity_mask_yidx[i];
      }
    }

    void FlowType::supply_m_prime() & {
      int m_prime=0;
      auto itr = m_table.begin();
      itr->m_prime() = m_prime;
      ++itr;
      for( ; itr!=m_table.end(); ++itr){
        if( std::vector<int>( (itr-1)->begin(), (itr-1)->begin()+RowOfTable::kYBlockBeg+4 ) ==
            std::vector<int>( itr->begin(), itr->begin()+RowOfTable::kYBlockBeg+4 ) ){
          itr->m_prime() = m_prime;
        }
        else{ // if two successive (x,y,nu) is NOT identical
          if( (itr-1)->y_idx()==itr->y_idx() ){ // y_idx is identical
            ++m_prime;
            itr->m_prime() = m_prime;
          }
          else{ // y is NOT identical
            m_prime=0;
            itr->m_prime() = m_prime;
          }
        }
      }
    }


    void FlowType::supply_next() & {
      for(auto itr = m_table.begin() ; itr!=m_table.end(); ++itr){
        int m_prime = itr->m_prime();
        auto itr_next = std::find_if(itr, m_table.end(),
                                     [&m_prime](const RowOfTable& row)
                                     { return row.m_prime()!=m_prime; });
        itr->next() = std::distance(itr,itr_next);
      }
    }


    void FlowType::calculate_multiplicity_reduced() & {
      for(int mask=0; mask<m_n_masks; ++mask){
        for(int index=0; index<m_unit_size; ++index){
          const int mul = multiplicity(mask, index);
          const int i = begin_mask_yidx(mask, index);
          m_multiplicity_reduced_mask_yidx[mask*m_unit_size+index] = m_table[i+mul-1][12]+1;
        }}
      auto itr_max = std::max_element(m_multiplicity_reduced_mask_yidx.begin(),
                                      m_multiplicity_reduced_mask_yidx.end());
      m_multiplicity_reduced_max = *itr_max;
    }

    double FlowType::calculate_bound_inverse() const& {
      double res = 0.0;
      const bool is_multiple = (this->number_of_multiplied_loops()==0);
      const bool is_nonlinear = this->is_flowed_link_nonlinear;

      if(is_multiple && is_nonlinear) res = 8.0 * this->number_of_multiplied_loops() * this->n_loops;
      else if(is_multiple && (!is_nonlinear)) res = 4.0 * this->number_of_multiplied_loops() * this->n_loops;
      else if((!is_multiple) && is_nonlinear) res = 4.0 * 2.0 * this->n_loops / 3.0; // assumed that the link has two loops
      else if((!is_multiple) && (!is_nonlinear)) res = 4.0 * this->n_loops / 3.0;
      else assert(false);

      return res;
    }
  }
}

#endif
