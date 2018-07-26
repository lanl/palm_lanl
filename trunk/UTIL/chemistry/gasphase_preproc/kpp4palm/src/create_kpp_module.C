
// ############################################################################
//
//     create_kpp_module 
//
//     create scalar code from .f90 sources created by KPP
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################
//
//Current revisions:
//------------------
/ 
/ 
//Former revisions:
//-----------------
//$Id: create_kpp_module.C 2470 2017-09-14 13:56:42Z forkel $
// Removed preprocessor directive __chem again 
//
/ 2017-09-14 13:56:42Z forkel $
//
// 
//added phot
//change of some output to lowercase with uppercase Fortran
//
// Nov 2016: Intial version (Klaus Ketelsen)
//

#include <stdio.h>
//  mz_rs_20090111+
// stdlib is necessary to define getenv:
#include <stdlib.h>
//  mz_rs_20090111-

#include "create_kpp_module.h"
#include "utils.h"

void create_kpp_module::do_work (string s) {
   vector<fortran_file>::iterator  it;
   vector<string>::iterator        ic;
   vector<Vvar>::iterator          iv;

   expand_decomp                   exp_de;

   prefix = s;
   module_name = prefix;

   cout << "Create " << module_name << " from kpp Fortran sources" <<endl;
   cout << "Vector mode " << kpp_switches.is_vector() <<endl;

   create_fortran_files_and_read();
   cout << "## after create_fortran_files_and_read " <<endl;

// Generate first module lines 

     string first_line="MODULE " + module_name;
     mz_kpp.add_line(first_line);
     mz_kpp.add_line("!");
//   mz_kpp.add_line("#if defined( __chem )");
//   mz_kpp.add_line(" ");

//    string e5_line = first_line +"_e5";
//    e5_kpp.add_line(e5_line);
//    e5_line = "  USE             " + module_name;
//    e5_kpp.add_line(e5_line);
//    e5_kpp.add_line(" ");

// edit include files

   for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
     it->edit_inc(header_variables);

//   Create variable Species list and vector variable list

     if(it->get_name() == module_name + "_Parameters") {
       it->create_species_list(species_list);
     }
     if(it->get_name() == module_name + "_Global") {
       it->vector_variable_list(Vvar_list);
     }
   }

// Prepare expansion of decomposition subroutine

   if(kpp_switches.de_indexing () > 0 ) {
     exp_de.create_sparse_info (kpp_includes, module_name);
   }

// edit FORTRAN files

   for(it=kpp_files.begin();it!=kpp_files.end();it++) {
     it->edit_fortran ();
   }
   cout << "## after edit FORTRAN files            " <<endl;

// Generate a list of single subroutines from kpp-files
// kpp files are modules containing several subroutines

   copy_files_to_subroutines ();

// All header_variables to include list

   kpp_includes.push_back(header_variables);

// Create decomposition subroutine
   if(kpp_switches.de_indexing () > 0 ) {
     exp_de.create_routine (kpp_subroutines);
   }


   if(kpp_switches.is_vector()) {

//   Change header section
     for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
       it->edit_inc_vec(global_variable_list);
     }

//   Change global variables to vector
     
     for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
       it->global_variables2vector (global_variable_list);
     }

//   Edit individual subroutines 

     for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
       if(it->get_name() == "KppDecomp") {
         it->edit_KppDecomp();
       }
       if(it->get_name() == "KppSolve") {
         it->edit_KppSolve();
       }
       if(it->get_name() == "Jac_SP" ) {
         it->edit_Jac_SP();
       }
       if(it->get_name() == "Fun" ) {
         it->edit_Fun();
       }
     }
   cout << "## after Edit individual subroutines    " <<endl;

   }

// Update_RCONST has to be changed also in scalar mode

   for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
     if(it->get_name() == "update_rconst") {
       it->edit_Update_RCONST(Vvar_list);
     }
   }

// Add Solver template to subroutine list
   if(kpp_switches.is_vector()) {
     add_solver_to_subroutine_list ();
   }

// The module header will be taken from ../templates/module_header.
// Please edit if header has to be changed.

   generate_module_header();

// Create kpp_integrate subroutine (chem_gasphase_integrate) for skalar and vector mode

   create_kpp_integrate();
// Copy include files

   for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
     it->copy_to_MZ_KPP(mz_kpp);
   }

   mz_kpp.add_line(" ");
   mz_kpp.add_line("! Interface Block ");
   mz_kpp.add_line(" ");
   for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
     string          buf;

     string prefix = "  ";
     for(ic=interface_ignore.begin();ic!=interface_ignore.end();ic++) {
       if(it->get_name() == *ic) {
         prefix = "!interface not working  ";
         break;
       }
     }

     buf = prefix + "INTERFACE            " + it->get_name() ;
     mz_kpp.add_line(buf);
     buf = prefix + "  MODULE PROCEDURE   " + it->get_name();
     mz_kpp.add_line(buf);
     buf = prefix + "END INTERFACE        " + it->get_name();
     mz_kpp.add_line(buf);
     mz_kpp.add_line(" ");
   }
   for(iv=Vvar_list.begin();iv!=Vvar_list.end();iv++) {
     string          buf;

     string sub_name = "fill_" + iv->name;
     buf = "  INTERFACE            " + sub_name;
     mz_kpp.add_line(buf);
     buf = "    MODULE PROCEDURE   " + sub_name;
     mz_kpp.add_line(buf);
     buf = "  END INTERFACE        " + sub_name;
     mz_kpp.add_line(buf);
     buf = "  PUBLIC               " + sub_name;
     mz_kpp.add_line(buf);
     mz_kpp.add_line(" ");
   }

   mz_kpp.add_line(" ");

   for(iv=Vvar_list.begin();iv!=Vvar_list.end();iv++) {
     create_fill_routine(kpp_subroutines, *iv);
   }

// Copy FORTRAN subroutines to mz_kpp

   mz_kpp.add_line(" CONTAINS");
   
   for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
     mz_kpp.add_line(" ");
     it->copy_to_MZ_KPP(mz_kpp);
   }

// Finish module

//   mz_kpp.add_line("#endif");
     string last_line="END MODULE " + module_name;
     mz_kpp.add_line("");
     mz_kpp.add_line(last_line);

// Handle e5 module

//    for(it=e5_subroutines.begin();it!=e5_subroutines.end();it++) {
//      e5_kpp.add_line(" ");
//      it->copy_to_MZ_KPP(e5_kpp);
//    }

//    last_line = last_line + "_e5";
//    e5_kpp.add_line(" ");
//    e5_kpp.add_line(last_line);

// Write the complete module to file: mz_kpp.f

   write_module_file();

   return;
}

void create_kpp_module::create_fortran_files_and_read() {

   string                          name;
   ifstream                        in,in_c,in_b,in_i;
   fortran_file                    f_file;
   vector<fortran_file>::iterator  it;

// Open file with list of FORTRAN routines

   in.open("file_list");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("file_list");
   }
   
// Create kpp_fortran routines
   while ( 1 ) {
     in >> name;
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_1");
     f_file.set_name(name);
     kpp_files.push_back(f_file);
   }
   in.close();

// Read FORTRAN code

   for(it=kpp_files.begin();it!=kpp_files.end();it++) {
     it->read();
   }

// Open file with list of include files

   in_c.open("include_list");
   if( !in_c ) {
      cout << "cannot open " << endl; my_abort("include_list");
   }

// Create kpp_includes vector
   while ( 1 ) {
     in_c >> name;
     if( in_c.eof() ) break;
     if( in_c.bad() ) my_abort("ERROR_READ_3");
     f_file.set_name(name);
     kpp_includes.push_back(f_file);
   }
   in_c.close();

// Read include files

   for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
     it->read();
   }

// Read Ignore list

   in_i.open("interface_ignore_list");
   if( !in_i ) {
      cout << "cannot open " << endl; my_abort("include_list");
   }

// Create kpp_includes vector
   while ( 1 ) {
     in_i >> name;
     if( in_i.eof() ) break;
     if( in_i.bad() ) my_abort("ERROR_READ_4");
     interface_ignore.push_back(name);
   }
   in_c.close();

}

void create_kpp_module::copy_files_to_subroutines () {
   string                          name;
   ifstream                        in;
   fortran_file                    s_file;
   vector<fortran_file>::iterator  it;

// Open file with list of FORTRAN routines

   in.open("subroutine_list");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("subroutine_list");
   }

// Create vector kpp_subroutines

   while ( 1 ) {
     in >> name;
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_S1");
     s_file.set_name(name);
     kpp_subroutines.push_back(s_file);
   }
   in.close();

   header_variables.add_line(" ");
   header_variables.add_line("!  variable definations from  individual module headers ");
   header_variables.add_line(" ");

//  Loop over all FORTRAN Files

   for(it=kpp_files.begin();it!=kpp_files.end();it++) {
     it->copy_to_subroutine_vector(kpp_subroutines, header_variables);
   }
}

void create_kpp_module::add_solver_to_subroutine_list () {
   fortran_file                    s_file;

   string solver_name = getenv("KPP_SOLVER");
   cout << "KPP_SOLVER " <<solver_name <<endl;
   
   s_file.set_name(solver_name);
   s_file.read();
   kpp_subroutines.push_back(s_file);

   return;
}

void create_kpp_module::generate_module_header() {

   string                          buf;
   ifstream                        in;
   ifstream                        in_e5;
   program_line                    line;
   vector<fortran_file>::iterator  it;
   char                            distr[2];
   string                          diline;

// Read Modul Header from file $MZ_KPP_HOME/templates/module_header

   in.open("module_header");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("module_header");
   }

   while ( 1 ) {
     getline (in, buf);
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_4");
     line.set_line(buf);
     mz_kpp.add_line(line); 
   }
   mz_kpp.add_line("                                                                 "); 
   mz_kpp.add_line("! Variables used for vector mode                                 "); 
   mz_kpp.add_line("                                                                 "); 
   if(kpp_switches.is_vector()) {
       mz_kpp.add_line("  LOGICAL,  PARAMETER          :: l_vector = .TRUE.             ");
   } else {
       mz_kpp.add_line("  LOGICAL,  PARAMETER          :: l_vector = .FALSE.            ");
   }
//  mz_pj_20070531+
   sprintf(distr,"%i",kpp_switches.de_indexing());
   diline = distr ;
   mz_kpp.add_line("  INTEGER,  PARAMETER          :: i_lu_di = " + diline );
//  mz_pj_20070531-

   mz_kpp.add_line("  INTEGER,  PARAMETER          :: vl_dim = " 
                   + kpp_switches.get_vector_length() ); 
   mz_kpp.add_line("  INTEGER                     :: vl                              "); 
   mz_kpp.add_line("                                                                 "); 
   mz_kpp.add_line("  INTEGER                     :: vl_glo                          "); 
   mz_kpp.add_line("  INTEGER                     :: is,ie                           "); 
   mz_kpp.add_line("                                                                 "); 
   mz_kpp.add_line("  INTEGER,  DIMENSION(vl_dim)  :: kacc, krej                      "); 
   mz_kpp.add_line("  INTEGER,  DIMENSION(vl_dim)  :: ierrv                           "); 
   mz_kpp.add_line("  LOGICAL                     :: data_loaded = .false.             "); 
   
   in.close();

//    in_e5.open("module_header_e5");
//    if( !in_e5 ) {
//       cout << "cannot open " << endl; my_abort("module_header_e5");
//    }

//    while ( 1 ) {
//      getline (in_e5, buf);
//      if( in_e5.eof() ) break;
//      if( in_e5.bad() ) my_abort("ERROR_READ_4");
//      line.set_line(buf);
//      e5_kpp.add_line(line);
//    }
//    in_e5.close();

   return;
}

void create_kpp_module::write_module_file() {
   ofstream                    out;
   ofstream                    out_e5;

   string out_file  = "kk_kpp.f90";
   out.open(out_file.c_str(), ios::out);
   if( !out ) {
      cout << "cannot open " << endl; my_abort(out_file);
   }

   mz_kpp.write_file (out);

   out.close();
   
//    out_file  = "kk_mecca_kpp_e5.f90";
//    out_e5.open(out_file.c_str(), ios::out);
//    if( !out_e5 ) {
//       cout << "cannot open " << endl; my_abort(out_file);
//    }

//    e5_kpp.write_file (out_e5);

//    out_e5.close();

   return;
}

void create_kpp_module::create_kpp_integrate() {
   fortran_file          kppi;
   vector<Vvar>::iterator               iv;
   string                               line;

   kppi.set_name("chem_gasphase_integrate");

   kppi.add_line("SUBROUTINE chem_gasphase_integrate (time_step_len, conc, tempk, photo, ierrf, xnacc, xnrej, istatus, l_debug, pe )  ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  IMPLICIT NONE                                                     ");
   kppi.add_line("                                                                    ");

   kppi.add_line("  REAL(dp),  INTENT(IN)                    :: time_step_len           ");
   kppi.add_line("  REAL(dp),  DIMENSION(:,:),  INTENT(INOUT) :: conc                    ");
   kppi.add_line("  REAL(dp),  DIMENSION(:,:),  INTENT(INOUT) :: photo                   ");
   kppi.add_line("  REAL(dp),  DIMENSION(:),  INTENT(IN)      :: tempk                   ");
   kppi.add_line("  INTEGER,  INTENT(OUT),  OPTIONAL          :: ierrf(:)                ");
   kppi.add_line("  INTEGER,  INTENT(OUT),  OPTIONAL          :: xnacc(:)                ");
   kppi.add_line("  INTEGER,  INTENT(OUT),  OPTIONAL          :: xnrej(:)                ");
   kppi.add_line("  INTEGER,  INTENT(INOUT),  OPTIONAL        :: istatus(:)              ");
   kppi.add_line("  INTEGER,  INTENT(IN),  OPTIONAL           :: pe                      ");
   kppi.add_line("  LOGICAL,  INTENT(IN),  OPTIONAL           :: l_debug                 ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  INTEGER                                 :: k   ! loop variable     ");
   kppi.add_line("  REAL(dp)                                :: dt                      ");
   kppi.add_line("  INTEGER,  DIMENSION(20)                  :: istatus_u               ");
   kppi.add_line("  INTEGER                                 :: ierr_u                  ");
   kppi.add_line("                                                                    ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  if (present (istatus) )  istatus = 0                              ");
   kppi.add_line("                                                                    ");
// kppi.add_line("  vk_glo = size(tempk,1)                                            ");
// kppi.add_line("                                                                    ");

   kppi.add_line("  DO k=1,vl_glo,vl_dim                                              ");
   kppi.add_line("    is = k                                                          ");
   kppi.add_line("    ie = min(k+vl_dim-1,vl_glo)                                     ");
   kppi.add_line("    vl = ie-is+1                                                    ");

   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    c(1:vl,:) = conc(is:ie,:)                                     ");
   } else {
     kppi.add_line("    c(:) = conc(is,:)                                             ");
   }
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    temp(1:vl) = tempk(is:ie)                                     ");
   } else {
     kppi.add_line("    temp = tempk(is)                                              ");
   }
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    phot(1:vl,:) = photo(is:ie,:)                                     ");
   } else {
     kppi.add_line("    phot(:) = photo(is,:)                                             ");
   }
   kppi.add_line("                                                                    ");
   kppi.add_line("    CALL update_rconst                                              ");
   kppi.add_line("                                                                    ");
   kppi.add_line("    dt = time_step_len                                              ");
   kppi.add_line("                                                                    ");
   kppi.add_line("    ! integrate from t=0 to t=dt                                    ");
   kppi.add_line("    CALL integrate(0._dp, dt, icntrl, rcntrl, istatus_u = istatus_u, ierr_u=ierr_u)");
   kppi.add_line("                                                                    ");
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    conc(is:ie,:) = c(1:VL,:)                                     ");
   } else {
     kppi.add_line("    IF (PRESENT(l_debug) .AND. PRESENT(pe)) THEN                       ");
     kppi.add_line("       IF (l_debug) CALL error_output(conc(is,:),ierr_u, pe)           ");
     kppi.add_line("    ENDIF                                                              ");
     kppi.add_line("    conc(is,:) = c(:)                                                  ");
   }

   kppi.add_line("                                                                    ");
   kppi.add_line("    ! Return Diagnostic Information                                 ");
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    if(PRESENT(ierrf))    ierrf(is:ie) = ierrv(1:vl)              ");
     kppi.add_line("    if(PRESENT(xnacc))    xnacc(is:ie) = kacc(1:vl)               ");
     kppi.add_line("    if(PRESENT(xnrej))    xnrej(is:ie) = krej(1:vl)               ");
   } else {
     kppi.add_line("    if(PRESENT(ierrf))    ierrf(is) = ierr_u                      ");
     kppi.add_line("    if(PRESENT(xnacc))    xnacc(is) = istatus_u(4)                ");
     kppi.add_line("    if(PRESENT(xnrej))    xnrej(is) = istatus_u(5)                ");
   }
   kppi.add_line("                                                                    ");
   kppi.add_line("    if (PRESENT (istatus) )  then                                   ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("      istatus(4) =   istatus(4) + sum(kacc(1:vl))                  ");
     kppi.add_line("      istatus(5) =   istatus(5) + sum(krej(1:vl))                  ");
     kppi.add_line("      istatus(3) =   istatus(4) + istatus(5)                       ");
     kppi.add_line("      istatus(6) =   istatus(6) + istatus_u(6)                     ");
     kppi.add_line("      istatus(7) =   istatus(7) + istatus_u(7)                     ");
   } else {
     kppi.add_line("      istatus(1:8) = istatus(1:8) + istatus_u(1:8)                 ");
   }
   kppi.add_line("    ENDIF                                                          ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  END DO                                                            ");
   kppi.add_line(" ");

   kppi.add_line("                                                                    ");
   kppi.add_line("! Deallocate input arrays                                           ");
   kppi.add_line("                                                                    ");
   for(iv=Vvar_list.begin();iv!=Vvar_list.end();iv++) {
     kppi.add_line("  if (ALLOCATED("+ iv->name +"))   DEALLOCATE("+ iv->name +" )    ");
   }

   kppi.add_line("                                                                    ");
   kppi.add_line("  data_loaded = .false.                                             ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  RETURN                                                            ");
   kppi.add_line("END SUBROUTINE chem_gasphase_integrate                                        ");

//   e5_subroutines.push_back(kppi);
   kpp_subroutines.push_back(kppi);

   return;
}

void create_kpp_module::create_fill_routine(vector<fortran_file> &fi_list, Vvar & var) {
   fortran_file                         fi;
   vector<string>::iterator             is;
   string                               line;

   cout << "Generate fill subroutine for " << var.name << endl;

   fi.set_name(var.name);
   line = "  SUBROUTINE fill_" + var.name;
   fi.add_line(line + "(status, array) ");
     fi.add_line(" ");
     fi.add_line("    INTEGER,  INTENT(OUT)               :: status ");
   if(var.nr_dim() == 0) {
     fi.add_line("    REAL(dp),  INTENT(IN),  DIMENSION(:) :: array ");
     fi.add_line(" ");
     fi.add_line("    status = 0");
     fi.add_line("    IF (.not. ALLOCATED("+var.name+")) & ");
     fi.add_line("       ALLOCATE("+var.name+"(size(array))) ");
   } else if(var.nr_dim() == 1) {
     fi.add_line("    REAL (dp), INTENT(IN), DIMENSION(:,:) :: array ");
     fi.add_line(" ");
     fi.add_line("    status = 0 ");
     fi.add_line("    if (.not. ALLOCATED("+var.name+")) & ");
     fi.add_line("        ALLOCATE("+var.name+"(size(array,1),"+var.dim_var[0]+")) ");
   } else if(var.nr_dim() == 2) {
     fi.add_line(" ");
     fi.add_line("    REAL (dp), INTENT(IN), DIMENSION(:,:,:) :: array ");
     fi.add_line(" ");
     fi.add_line("    status = 0 ");
     fi.add_line("    if (.not. ALLOCATED("+var.name+")) & ");
     fi.add_line("        ALLOCATE("+var.name+"(size(array,1),"+var.dim_var[0]+var.dim_var[1]+")) ");
   } else {
     fi.add_line("    REAL (dp), INTENT(IN), DIMENSION(:,:,:,:) :: array ");
     fi.add_line(" ");
     fi.add_line("    status = 0 ");
     fi.add_line("    IF (.not. ALLOCATED("+var.name+")) & ");
     fi.add_line("        ALLOCATE("+var.name+"(size(array,1),"+var.dim_var[0]
                               +var.dim_var[1]+var.dim_var[3]+")) ");
   }

   fi.add_line(" ");
   fi.add_line("    IF (data_loaded .AND. (vl_glo /= size(array,1)) )  THEN ");
   fi.add_line("       status = 1 ");
   fi.add_line("       RETURN ");
   fi.add_line("    END IF ");
   fi.add_line(" ");
   fi.add_line("    vl_glo = size(array,1) ");
   fi.add_line("    "+var.name+ " = array ");
   fi.add_line("    data_loaded = .TRUE. ");
   fi.add_line(" ");
   fi.add_line("    RETURN");
   fi.add_line(" ");
   fi.add_line("  END " + line);

   fi_list.push_back(fi);

   return;
}
