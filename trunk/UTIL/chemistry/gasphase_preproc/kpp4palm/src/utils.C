
// ############################################################################
//
//     create_mz_kpp_module 
//
//     create scalar code from .f90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################

//  mz_rs_20090111+
// stdlib is necessary to define abort():
#include <stdlib.h>
//  mz_rs_20090111-

#include "utils.h"

void string_token::fill_token (string s) {
   int      ip,ib;

   string buf=s+" ";
   size_val  = 0;

   tokens.clear();
   position.clear();

// extract fields from buffer

   int i=0;
   int pos=0;
   while( 1 ) {
     while (1) {
        ib = buf.find(seperator,0);
        if(ib == string::npos)   break;
        if(ib != 0)   break;
        buf.erase(0,1);
        pos++;
     }
     ip = buf.find(seperator,0);
     if(ip == string::npos)   break;
     tokens.push_back(buf.substr(0,ip));
     position.push_back(pos);
     size_val++;
     buf.erase(0,ip);
     pos += ib;

     i++;
   }

   is = tokens.begin();

   return;
};


void my_abort (string s) {
   cout << "*** ERROR: " << s  <<endl;
   abort();

   return;
}

string my_to_upper(string sinp) {
  string       s1,s2;

  const char   *c1;
  char         c2[2];
  int          i1,i2;

  for(int i=0; i<sinp.size(); i++) {
    s1.clear();
    s1 = sinp.substr(i,1);
    c1 = s1.c_str();

    i1 = *c1;
    i2 = toupper(i1);
    *c2 = i2;
    c2[1] = 0;

    s2 += c2;
  }

  return s2;
}
