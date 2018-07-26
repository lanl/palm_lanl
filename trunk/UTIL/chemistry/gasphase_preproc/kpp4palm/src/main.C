
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
// stdlib is necessary to define atoi:
#include <stdlib.h>
//  mz_rs_20090111-

#include "create_kpp_module.h"
control_switches kpp_switches;

main (int argc, char *argv [] ) {
   create_kpp_module         cm;
   string                    inp_name;
   char                      *lo_arg;

   cout << "####################################################" <<endl ;
   cout << "###             KP4 = KPP POST PROCESSOR         ###" <<endl ;
   cout << "###                   Version 1.0                ###" <<endl ;
   cout << "### (C) by Klaus Ketelsen and MPI-CH, April 2007 ###" <<endl ;
   cout << "####################################################" <<endl ;

   lo_arg = argv[1];

   if(argc > 1) {
     inp_name = lo_arg ;
   } else {
     inp_name = "kk_" ;
   }

   if(argc > 2) {
     lo_arg = argv[2];
     string vm = lo_arg;
     if(vm == "vector") {
       kpp_switches.set_vector_mode ();
     }
   }

   if(argc > 3) {
     lo_arg = argv[3];
     string vl = lo_arg;
     kpp_switches.set_vector_length (vl);
   }

   if(argc > 4) {
     lo_arg = argv[4];
     int de = atoi(lo_arg);
     kpp_switches.set_de_index_mode (de);
   }

   cm.do_work(inp_name);

   cout << "####################################################" <<endl ;
   cout << "###                END OF KP4                    ###" <<endl ;
   cout << "####################################################" <<endl ;

   return 0;
}
