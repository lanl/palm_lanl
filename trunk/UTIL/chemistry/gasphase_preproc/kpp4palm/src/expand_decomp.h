#ifndef EXPDE
#define EXPDE 1

// ############################################################################
//
//     create mz_kpp_module                       
//
//     create code from .f90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################

#include <iostream>

#include <string>
#include <list>
#include <vector>

#include "fortran_file.h"
#include "program_line.h"

class expand_decomp {

  int                   NVAR;
  vector<int>           LU_IROW;
  vector<int>           LU_ICOL;
  vector<int>           LU_CROW;
  vector<int>           LU_DIAG;

  public:

  void create_sparse_info (vector<fortran_file> & include_list, string module_name);
  void create_routine (vector<fortran_file> & routine_list);

};

#endif
