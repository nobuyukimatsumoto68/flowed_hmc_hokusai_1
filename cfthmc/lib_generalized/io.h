#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include <qlat/hmc.h>

namespace qlat{
  namespace Generalized{
    namespace IO{

      inline void wb( std::ofstream& ofs, const double& elem) {
        double d = elem;
        ofs.write(reinterpret_cast<char *>(&d), sizeof(double));
      }

      inline void wb( std::ofstream& ofs, const int& elem) {
        int i = elem;
        ofs.write(reinterpret_cast<char *>(&i), sizeof(int));
      }

      inline void wb
      (
       const std::string& f_out,
       const double& elem
       ){
        double d = elem;
        //
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          std::ofstream out(f_out, std::ios::binary | std::ios::app);
          assert(out.is_open());
          out.write(reinterpret_cast<char *>(&d), sizeof(double));
          out.close();
        }
      }

//       inline void wb
//       (
//        const std::string& f_out,
//        const int& elem
//        ){
//         int d = elem;
//         //
//         int id_node = 0;
// #ifdef USE_MULTI_NODE
//         MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
// #endif
//         if(id_node==0) {
//           std::ofstream out(f_out, std::ios::binary | std::ios::app);
//           assert(out.is_open());
//           out.write(reinterpret_cast<char *>(&i), sizeof(int));
//           out.close();
//         }
//       }

      template <class Stream>
      void setup_io( Stream& s ){
        s << std::scientific << std::setprecision(15);
      }

      void setup_file
      (
       const std::string& f_out,
       const std::ios_base::openmode& mode
       ){
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          std::ofstream out(f_out, mode);
          assert(out.is_open());
          out.close();
        }
      }

      template <class Out>
      void print_to( Out& out, const std::string& str )
      {
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          out << std::scientific << std::setprecision(15);
          out << str << std::flush;
        }
      }

      template <class Out>
      void print_to( Out& out, const double d )
      {
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          out << std::scientific << std::setprecision(15);
          out << d << std::flush;
        }
      }

      template <class Out>
      void print_to( Out& out, const int i )
      {
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          out << std::scientific << std::setprecision(15);
          out << i << std::flush;
        }
      }

      void print_to_file
      (
       const std::string& f_out,
       const std::string& str,
       const std::ios_base::openmode& mode
       ){
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          std::ofstream out(f_out, mode);
          assert(out.is_open());
          out << std::scientific << std::setprecision(15);
          out << str << std::flush;
          out.close();
        }
      }

      template <class M>
      void print_matrix( std::ostream& out, const M& m ){
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          const int size = std::sqrt(sizeof(m)/sizeof(m(0,0)));
          for(int i=0; i<size; ++i) for(int j=0; j<size; ++j) {
              double d = m(i,j).real();
              out << d << " + I* ";
              d = m(i,j).imag();
              out << d << ", " << std::flush;
            }
        }
      }


      template <class M>
      void print_matrix( std::ofstream& out, const M& m ){
        int id_node = 0;
#ifdef USE_MULTI_NODE
        MPI_Comm_rank(MPI_COMM_WORLD, &id_node);
#endif
        if(id_node==0) {
          const int size = std::sqrt(sizeof(m)/sizeof(m(0,0)));
          for(int i=0; i<size; ++i) for(int j=0; j<size; ++j) {
              double d = m(i,j).real();
              wb(out,d);
              d = m(i,j).imag();
              wb(out,d);
            }
        }
      }


      template <class Out>
      void print_matrix_field( Out& out, const GaugeMomentum& fm ){
        const Geometry geo = geo_reform(fm.geo(), 1);
        const Coordinate total_site = geo.total_site();
        const Geometry geo_ext = geo_reform(geo, 1, total_site, total_site);

        GaugeMomentum fm_ext;
        fm_ext.init(geo_ext);
        fm_ext = fm;
        refresh_expanded(fm_ext);

        for(int x0=0; x0<total_site[0]; ++x0){
          for(int x1=0; x1<total_site[1]; ++x1){
            for(int x2=0; x2<total_site[2]; ++x2){
              for(int x3=0; x3<total_site[3]; ++x3){
                const Coordinate xl = geo_ext.coordinate_l_from_g(Coordinate(x0,x1,x2,x3));
                assert(geo_ext.is_on_node(xl));
                Vector<ColorMatrix> fmx = fm_ext.get_elems(xl);
                for(int mu=0; mu<DIMN; ++mu) IO::print_matrix(out, fmx[mu]);
              }}}}
      }
    }
  }
}
