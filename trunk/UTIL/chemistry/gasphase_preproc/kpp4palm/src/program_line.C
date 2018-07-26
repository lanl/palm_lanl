
// ############################################################################
//
//     create_mz_kpp_module                       
//
//     create vectorcode from .90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################

#include "program_line.h"

void  program_line::set_line(string s) {

   line = s;

   tokens.fill_token(line);

   return;
}
void  program_line::substitute(string old_s, string new_s) {

   int pos = line.find (old_s, 0);       // look for string

   if (pos != string::npos) {                // found
     line.replace(pos,old_s.size(),new_s);
   }   

   return;
}

void  program_line::global_substitute(string old_s, string new_s) {
   int         pos;

   int start = line.size()-1;

   while (1) {
     pos = line.rfind (old_s, start);       // look for string

     if (pos == string::npos) {
       break;
     }

     line.replace(pos,old_s.size(),new_s);

     start = pos-1;
   }

   return;
}

int  program_line::get_token_number_from_string (string s) {

   for (int i=0; i<tokens.size(); i++) {
     if(tokens.get_token_by_index(i) == s)  {
       return i;
     }
   }

   return -1;
}

int  program_line::get_token_number_from_string_upper (string s) {

   string s1 = my_to_upper(s);
   for (int i=0; i<tokens.size(); i++) {
     if(my_to_upper(tokens.get_token_by_index(i)) == s1)  {
       return i;
     }
   }

   return -1;
}

void  program_line::update_token (int i, string s) {

  string s_old=tokens.get_token_by_index(i);
  tokens.update (i,s);
  int pos=tokens.get_position(i);
  line.replace(pos,s_old.size(),s);

  return;
}

void  program_line::change_variable_to_vector (string var) {
  string   U_token;

  int     ind = tokens.size()-2;
  for (int i=ind; i >= 0; i--) {
    U_token.clear();
    U_token = my_to_upper(tokens.get_token_by_index(i));
    if(U_token == my_to_upper(var) && tokens.get_token_by_index(i+1) == "(" ) {
      string new_s="(1:VL,";
      tokens.update (i+1,new_s);
      int pos=tokens.get_position(i+1);
      line.replace(pos,1,new_s);
      tokens.reset (line);
    }
  }

  return;
}

void  program_line::change_variable_to_vector_g (Vvar &var) {
  string   U_token;
  string   new_s;

  tokens.reset (line);

  int     ind = tokens.size()-2;
  for (int i=ind; i >= 0; i--) {
    U_token.clear();
    U_token = tokens.get_token_by_index(i);
    if(var.nr_dim() == 0) {
      if(my_to_upper(U_token) == my_to_upper(var.name) ) {
        new_s=U_token + "(k)";
        tokens.update (i,new_s);
        int pos=tokens.get_position(i);
        line.replace(pos,U_token.size(),new_s);
        tokens.reset (line);
      }
    } else {
      if(my_to_upper(U_token) == my_to_upper(var.name) && tokens.get_token_by_index(i+1) == "(" ) {
        new_s="(k,";
        tokens.update (i+1,new_s);
        int pos=tokens.get_position(i+1);
        line.replace(pos,1,new_s);
        tokens.reset (line);
      }
    }
  }

  return;
}

