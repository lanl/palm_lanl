#ifndef FFILE
#define FFILE 1

// ############################################################################
//
//     create mz_kpp_module                       
//
//     create code from .f90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################
//
//  Modifications:
//    RFo added global_subtolower

#include <iostream>

#include <string>
#include <list>
#include <vector>

#include "program_line.h"

class fortran_file {

  string                    name;

  public:
  vector<program_line>      pline;

  void   set_name(string s)       {name=s;return;};
  string get_name()                {return name;};
  void   add_line(program_line &s) {pline.push_back(s);return;};
  void   add_line(string s)        {program_line p; p.set_line(s); pline.push_back(p);return;};

  void read () ;
  void edit_fortran () ;
  void copy_to_subroutine_vector (vector<fortran_file> &subvec, fortran_file & header_var) ;
  void edit_Jac_SP () ;
  void edit_KppDecomp () ;
  void edit_KppSolve () ;
  void edit_Fun () ;
  void edit_FUNC () ;
  void edit_Update_RCONST (vector <Vvar> &var_list) ;
  void edit_inc (fortran_file & header_var) ;
  void create_species_list(vector <string> &species_list);
  void vector_variable_list(vector <Vvar> &var_list);
  void print () ;
  void write_file (ofstream & out);
  void copy_to_MZ_KPP (fortran_file & ka);
  void global_substitute(string &line, string old_s, string new_s);
  void global_subtolower(string &line);

// Routines for FORTRAN vector Mode
  void edit_inc_vec (vector<string> &gvl);
  void global_variables2vector (vector<string> &gvl);

};

#endif
