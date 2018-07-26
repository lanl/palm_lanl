#ifndef PLINE
#define PLINE 1

// ############################################################################
//
//     create_mz_kpp_module                       
//
//     create code from .90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################

#include <iostream>

#include <string>
#include <list>
#include <vector>

#include "utils.h"

using namespace std;

class program_line {

   string            line;
   string_token      tokens;

  public:

   string get_line()          {return line;};
   string get_token(int i)    {return tokens.get_token_by_index(i);};
   int    get_token_size()    {return tokens.size();};

   void   set_line(string s);
   void   substitute(string old_s, string new_s);
   void   global_substitute(string old_s, string new_s);
   int    get_token_number_from_string (string s);
   int    get_token_number_from_string_upper (string s);

   void   update_token                 (int i, string s);
   void   change_variable_to_vector    (string var);
   void   change_variable_to_vector_g  (Vvar  &var);
};

#endif
