
// ############################################################################
//
//     create_mz_kpp_module                       
//
//     create scalar code from .f90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################
//
//
//Current revisions:
//------------------
// 
// 
//Former revisions:
//-----------------
//$Id: fortran_file.C 2470 2017-09-14 13:56:42Z forkel $
//
// 
//changed KPP-generated code to lowercase with uppercase Fortran  expressions
//added photolysis variables
//
//
//
//
// Nov 2016: Intial version (Klaus Ketelsen)
//


#include <fstream>
#include "fortran_file.h"

#include "utils.h"
#include "ctype.h"
void fortran_file::read () {

  ifstream                        in;
  program_line                    pl;
  string                          line;

// Note: FORTRAN77 and include files are internally named .f90

  string file_name = name + ".f90";
  in.open(file_name.c_str() );
  if( !in ) {
    cout << "cannot open " << endl; my_abort(file_name);
  }

// Read kpp_fortran routines;
  while ( 1 ) {
     getline(in,line);
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_2");

// Remove trailing blanks
     while (1) {
       if(line.size() < 100) break;
       if(line.substr(line.size()-1,1) != " ")  break;
       line.erase(line.size()-1,1);
     }

// Now do some cosmetics to adapt the KPP generated output a bit o the looks of PALM, 
// i.e. add some blanks, convert all to lowercase except Fortran language elements, etc.
       if(line.find("'",0) == string::npos)  {     // No substitution in line with strings
        if(line.substr(0,2) =="! ") {
         global_substitute(line,"#","#");      // just a dummy so comments do not get lost
         } else {
         global_substitute(line,"("," ( ");
         global_substitute(line,")"," ) ");
         global_substitute(line,",",", ");
         global_substitute(line,",  ",", ");
         global_substitute(line,")  ",") ");
         global_substitute(line,"*","* ");
         global_substitute(line,"* *","**");
         global_substitute(line,"/JVS","/ JVS");
         global_substitute(line,"-","- ");
         global_substitute(line,"e- ","e-");
         global_substitute(line,"+","+ ");
         global_substitute(line,"+  ","+ ");
         global_substitute(line,"d- ","d-");
         global_substitute(line,"D- ","D-");
         global_substitute(line,"e+ ","e+");
         global_substitute(line,"E+ ","E+");
         global_substitute(line,"E- ","E-");
         global_substitute(line,")=",") =");
         global_subtolower(line);
         global_substitute(line,"call ","CALL ");
         global_substitute(line,"case","CASE");
         global_substitute(line,"character","CHARACTER");
         global_substitute(line,"contains","CONTAINS");
         global_substitute(line,"dimension","DIMENSION");
         global_substitute(line,"do ","DO ");
         global_substitute(line,"elseif","ELSEIF");
         global_substitute(line,"else","ELSE");
         global_substitute(line,"#ELSE","#else");
         global_substitute(line,"end do","ENDDO");
         global_substitute(line,"end if","ENDIF");
         global_substitute(line,"end ","END ");    // Modify "end" after all other strings containing "end..." are done!
         global_substitute(line,"#ENDIF","#endif");
         global_substitute(line,"function","FUNCTION");
         global_substitute(line,"if(","IF (");
         global_substitute(line," if "," IF ");
         global_substitute(line,"implicit","IMPLICIT");
         global_substitute(line,"include","INCLUDE");
         global_substitute(line,"intent","INTENT");
         global_substitute(line,"integer","INTEGER");
         global_substitute(line,"logical","LOGICAL");
         global_substitute(line,"module","MODULE");
         global_substitute(line,"none","NONE");
         global_substitute(line,"optional","OPTIONAL");
         global_substitute(line,"parameter","PARAMETER");
         global_substitute(line,"public","PUBLIC");
         global_substitute(line,"real","REAL");
         global_substitute(line,"return","RETURN");
         global_substitute(line,"use ","USE ");
         global_substitute(line,"save","SAVE");
         global_substitute(line,"subroutine","SUBROUTINE");
         global_substitute(line,"then","THEN");
         global_substitute(line,"while","WHILE");
       }
     }

     pl.set_line(line);
     pline.push_back(pl);
   }
   in.close();

  return;
}
void fortran_file::edit_fortran () {

  vector<program_line>::iterator     ip;
  string                             lo_line;
  bool                               deleted;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(ip->get_token(0) =="!")   continue;           // No editing of Comments
    deleted = false;
    lo_line = ip->get_line() ;
    if(ip->get_token(0) == "MODULE" && ip->get_token_size() >= 1) {
      lo_line.insert(0,"!DELETE ");
    }
    if(ip->get_token(0) == "USE" && ip->get_token_size() >= 1) {
      deleted = true;
      lo_line.insert(0,"!DELETE ");
    }
    if(ip->get_token(0) == "PUBLIC" && ip->get_token_size() >= 1) {
      deleted = true;
      lo_line.insert(0,"!DELETE ");
    }
//    if(ip->get_token(0) == "SAVE" && ip->get_token_size() >= 1) {
//      lo_line.insert(0,"!DELETE ");
//    }

//  Only IMPLICIT none, not IMPLICIT REAL (A-H,O-Z)
    if(ip->get_token(0) == "IMPLICIT" && ip->get_token_size() == 2) {
      lo_line.insert(0,"!DELETE ");
    }

//  Delete INCLUDE lines
    if(ip->get_token(0) == "INCLUDE" && ip->get_token_size() >= 1) {
      lo_line.insert(0,"!DELETE ");
    }

//  Update_RCONST has only to be called once per outer timeloop in KPP_FOR_PALM

//  if(ip->get_token(0) == "CALL" && ip->get_token(1) == "Update_RCONST" ) {
    if(ip->get_token(0) == "CALL" && ip->get_token(1) == "update_rconst" ) {
      lo_line.insert(0,"!DELETE ");
    cout << lo_line << endl;
    }

//  Update_SUN must not be called within in KPP_FOR_PALM

//  if(ip->get_token(0) == "CALL" && ip->get_token(1) == "Update_SUN" ) {
    if(ip->get_token(0) == "CALL" && ip->get_token(1) == "update_sun" ) {
      lo_line.insert(0,"!DELETE ");
    cout << lo_line << endl;
    }


    ip->set_line(lo_line);

//  Delete continuation lines

    if(deleted) {
      while( ip->get_token_number_from_string("&") > 1) {
        ip++;
        lo_line = ip->get_line() ;
        lo_line.insert(0,"!DELETE ");
        ip->set_line(lo_line);
      }
    }
  }

  return;
}

void fortran_file::copy_to_subroutine_vector (vector<fortran_file> &subvec, 
                                                         fortran_file & header_var) {
  vector<program_line>::iterator     ip;
  vector<fortran_file>::iterator     iv;
  bool                               active_subroutine;
  string                             active_name;
  int                                name_pos;

// loop over all lines in a fortran file


// First, copy Module variables into header_var 
// This variables wil later copied into the kpp-moduzle header

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(ip->get_token(0) == "CONTAINS")   break;          // Header done 

//  Special Treatment of defining variables in PARAMETER context
//  The intel compiler creates fatal error, if blanks are embedded between / /
//  The SX-6 Compiler gives a warning 

    ip->substitute("( / &","(/&");
    ip->substitute(",  & !",",& !");
    ip->substitute(" / )","/)");


    header_var.add_line(ip->get_line());
  }

//  look for SUBROUTINE statement
  active_subroutine = false;

  for (ip=pline.begin(); ip != pline.end(); ip++) {

//  look for SUBROUTINE statement
    if(ip->get_token(0) == "SUBROUTINE" && !active_subroutine) {

      for(iv=subvec.begin();iv!=subvec.end();iv++) {
        if(ip->get_token(1) == iv->get_name() ) { 
//        Subroutine is in list
          active_subroutine = true;
          active_name = ip->get_token(1);
          cout << "SUBROUTINE: " << active_name << " found in file " << name << endl;
          break;
        }
      }
    }

//  look for FUNCTION statement
    name_pos = ip->get_token_number_from_string("FUNCTION");
    if(name_pos != -1 && !active_subroutine) {
      name_pos++;
      for(iv=subvec.begin();iv!=subvec.end();iv++) {
        if(ip->get_token(name_pos) == iv->get_name() ) { 
//        Subroutine is in list
          active_subroutine = true;
          active_name = ip->get_token(name_pos);
          cout << "FUNCTION  : " << active_name << " found in file " << name << endl;
          break;
        }
      }
    }

    if(active_subroutine)  {
//    copy FORTRAN line from file to subroutine
      iv->add_line(ip->get_line());
    }
    if(ip->get_token(1) == "SUBROUTINE" && ip->get_token(2) == active_name) {
      cout << "SUBROUTINE: " << active_name << " done " << endl;
      active_subroutine = false;
      active_name = " ";
    }
    if(ip->get_token(1) == "FUNCTION" && ip->get_token(2) == active_name) {
      cout << "FUNCTION  : " << active_name << " done " << endl;
      active_subroutine = false;
      active_name = " ";
    }
  }

  return;
}

void fortran_file::edit_FUNC () {

  vector<program_line>::iterator     ip;
  string                             lo_line;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    ip->substitute("REAL*8   Y(NVAR), J(LU_NONZERO)","real (kind=8),dimension(:,:)  :: y,j");
    ip->substitute("REAL*8   Y(NVAR), P(NVAR)",      "real (kind=8),dimension(:,:)  :: y,p");
    lo_line = ip->get_line() ;
    if(ip->get_token(0) == "Told" && ip->get_token(1) == "=") {
      lo_line.insert(0,"!DELETE ");
    }
    if(ip->get_token(0) == "TIME" && ip->get_token(1) == "=") {
      lo_line.insert(0,"!DELETE ");
    }
    ip->set_line(lo_line);
  }

  return;
}

void fortran_file::edit_inc (fortran_file & header_var) {

  vector<program_line>::iterator     ip;
  string                             lo_line;
  string                             lo_line_2;
  string                             lo_line_3;
  bool                               deleted;
  int                                i,nr_var;
  string                             public_line;

// Delete module and end module lines

  header_var.add_line("! Automatic generated PUBLIC Statements for ip_ and ihs_ variables ");
  header_var.add_line(" ");

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    deleted = false;
    lo_line = ip->get_line() ;
    if(ip->get_token(0) == "MODULE" && ip->get_token_size() >= 1) {
      lo_line.insert(0,"!DELETE ");
    }
    if(ip->get_token(0) == "END" && ip->get_token_size() >= 1) {
      lo_line.insert(0,"!DELETE ");
    }
    if(ip->get_token(0) == "USE" && ip->get_token_size() >= 1) {
      deleted = true;
      lo_line.insert(0,"!DELETE ");
    }
    if(ip->get_token(0) == "PUBLIC" && ip->get_token_size() >= 1) {
      deleted = true;
      lo_line.insert(0,"!DELETE ");
    }
    if(ip->get_token(0) == "SAVE" && ip->get_token_size() >= 1) {
      lo_line.insert(0,"!DELETE ");
    }

//  Make ind_ variables public
    if(ip->get_token(3).substr(0,4) == "ind_") {
      global_substitute (lo_line,"PARAMETER","PARAMETER, PUBLIC");
    }
 
// Make ip_ variables public

//  mz_rs_20070907+
// ip_* are already public
// Make ip_ variables public
//     nr_var=0;
//     for(i=0;i<ip->get_token_size();i++) {
//       if(ip->get_token(i).substr(0,3) == "ip_" && ip->get_token(0).substr(0,1) != "!") {
//         nr_var++;
//         if(nr_var == 1)  {
//           public_line.clear();
//           public_line = "  PUBLIC  " + ip->get_token(i);

//         } else {
//           public_line += ", " + ip->get_token(i);
//         }
//       }
//     }
//     if(nr_var > 0) {
//       header_var.add_line(public_line);
//     }
//  mz_rs_20070907-

// Make ihs_ variables public

    nr_var=0;
    for(i=0;i<ip->get_token_size();i++) {
      if(ip->get_token(i).substr(0,4) == "ihs_" && ip->get_token(0).substr(0,1) != "!") {
        nr_var++;
        if(nr_var == 1)  {
          public_line.clear();
          public_line = "  PUBLIC  " + ip->get_token(i);

        } else {
          public_line += ", " + ip->get_token(i);
        }
      }
    }
    if(nr_var > 0) {
      header_var.add_line(public_line);
    }

    ip->set_line(lo_line);

//  Delete continuation lines

    if(deleted) {
      while( ip->get_token_number_from_string("&") > 1) {
        ip++;
        lo_line = ip->get_line() ;
        lo_line.insert(0,"!DELETE ");
        ip->set_line(lo_line);
      }
    }
  }

  for (ip=pline.begin(); ip != pline.end(); ip++) {
//  Special Treatment of defining variables in PARAMETER context
//  The intel compiler creates fatal error, if blanks are embedded between / /
//  The SX-6 Compiler gives a warning 

    if(ip->get_token(7) == "substep" && ip->get_token(9) == "nsubsteps") {
      ip->substitute("( / 0. / )","0.0");
    }
    ip->substitute("( /","(/");
    ip->substitute("/ )","/)");
  }

  return;
}

void fortran_file::create_species_list(vector <string> &species_list) {
  vector<program_line>::iterator     ip;
  string                             lo_line;
  string                             specname;
  string                             longname;

  bool to_do = false;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(ip->get_token(2) == "declaration" && ip->get_token(5) == "species") {
      to_do = true;
    }
    if(ip->get_token(5) == "0" && ip->get_token_size() > 5 ) {
       to_do = false;
    }
    if(ip->get_token(0) == "INTEGER," && to_do) {
      lo_line = ip->get_line() ;
      int pos = lo_line.find("ind_",0);
      int end = lo_line.find(" ",pos+1);
      if(pos > 1) {
        specname.clear();
        longname.clear();
        specname = lo_line.substr(pos+4,end-pos-4);
        longname = specname + "                 ";
	//        if(specname != "H2O")  {                     // NO idt_variable for H2O
          species_list.push_back(specname);
        //}
      }
    }
  }

  return;
}

void fortran_file::vector_variable_list(vector <Vvar> &var_list) {
  vector<program_line>::iterator     ip;
  string                             lo_line;
  Vvar                               vari;
  bool                               todo_1;

  bool to_do = false;
  for (ip=pline.begin(); ip != pline.end(); ip++) {
    todo_1 = false;
    if(ip->get_token(0) == "!KPPPP_DIRECTIVE" && ip->get_token(4) == "start") {
      to_do = true;
      ip++;
    }
    if(ip->get_token(0) == "!KPPPP_DIRECTIVE" && ip->get_token(4) == "end") {
       to_do = false;
    }
    if(ip->get_token(0).substr(0,1) != "!" && ip->get_token_number_from_string("temp") > 2) {
//  if(ip->get_token(0).substr(0,1) != "!" && ip->get_token_number_from_string("XXXX") > 2) {
      todo_1 = true;
    }
    if(to_do || todo_1) {
      if(ip->get_token(0).substr(0,1) == "!") continue;       // skip comment limes
      vari.clear();
      int pos1 = ip->get_token_number_from_string_upper ("DimensioN");
      int pos  = ip->get_token_number_from_string("::");
      vari.name = ip->get_token(pos+1);
      if(pos1 > 1) {
        for (int i=pos1+2; i< ip->get_token_size(); i++) {
          if(ip->get_token(i).substr(0,1) == "!" || ip->get_token(i) == ")" ) {
            break;
          } else {
            ip->substitute(","," "); 
            vari.dim_var.push_back(ip->get_token(i));
          }
        }
      } else {
        if(ip->get_token_size() > pos+1 && ip->get_token(pos+2).substr(0,1) != "!" 
                          && ip->get_token(pos+2).substr(0,1) != "=" )  {
          for (int i=pos+3; i< ip->get_token_size(); i++) {
            if(ip->get_token(i).substr(0,1) == "!" || ip->get_token(i) == ")" ) {
              break;
            } else {
              ip->substitute(","," "); 
              vari.dim_var.push_back(ip->get_token(i));
            }
          }
        }
      }
      var_list.push_back(vari);
      cout << "Vector variable " << ip->get_token(pos+1) <<" " << vari.nr_dim() <<endl;
      lo_line = ip->get_line() ;
      if(vari.nr_dim() == 0) {
        lo_line.clear();
        lo_line = "  REAL(dp),dimension(:),allocatable             :: " + vari.name;
      }
      if(vari.nr_dim() == 1) {
        lo_line.clear();
        lo_line = "  REAL(dp),dimension(:,:),allocatable           :: " + vari.name;
      }
      if(vari.nr_dim() == 2) {
        lo_line.clear();
        lo_line = "  REAL(dp),dimension(:,:,:),allocatable         :: " + vari.name;
      }
      if(vari.nr_dim() == 3) {
        lo_line.clear();
        lo_line = "  REAL(dp),dimension(:,:,:,:),allocatable       :: " + vari.name;
      }
      ip->set_line(lo_line);
    }
  }

  return;
}

void fortran_file::print () {

  vector<program_line>::iterator     ip;

  cout << " " <<endl;
  cout << "FORTRAN file " << name << endl;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    cout << ip->get_line() <<endl;
  }

  cout << " " <<endl;

  return;
}

void fortran_file::write_file (ofstream & out) {

  vector<program_line>::iterator     ip;
  string                             lo_line;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    ip->global_substitute(" ( ","(");
    ip->global_substitute(" ) ",")");
    ip->global_substitute(" )",")");
    ip->global_substitute("A (","A(");
    ip->global_substitute("B (","B(");
    ip->global_substitute("V (","V(");
    ip->global_substitute("Vdot (","Vdot(");
    ip->global_substitute("JVS (","JVS(");
    ip->global_substitute("RCT (","RCT(");
    ip->global_substitute("* ","*");
    ip->global_substitute("* ","*");
    ip->global_substitute("d- ","d-");
    ip->global_substitute("d+ ","d+");
    ip->global_substitute(", ",",");
    ip->global_substitute(")=",") =");
    ip->global_substitute("dp,","dp, ");
    ip->global_substitute(", - ",", -");

//  line break if more than 130 character

    lo_line = ip->get_line();
    if( lo_line.size() < 130 ) {
      out << ip->get_line() <<endl;
    } else {
      int cp  = lo_line.rfind("!",129);
      int pos = lo_line.rfind(" ",129);
      out << lo_line.substr(0,pos) << " &"<<endl;
      lo_line.erase (0,pos);
      if(ip->get_token (0).substr(0,1)  == "!" || cp != string::npos ) {
        out << "!   " << lo_line   <<endl;                 // comment also in next line
      } else {
        out << "                    " << lo_line   <<endl;
      }
    }
  }

  out << " " <<endl;

  return;
}

void fortran_file::copy_to_MZ_KPP (fortran_file & ka) {

  vector<program_line>::iterator     ip;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
//  Do not copy lines marked for delete
    if(ip->get_token(0) != "!DELETE") {
      ka.add_line( ip->get_line() );
    }
  }

  return;
}
void  fortran_file::global_substitute(string &line, string old_s, string new_s) {
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

void  fortran_file::global_subtolower(string &line) {

   int start = line.size()-1;
   char c;

    int i = 0;
    while (line[i])
    {
      c = line[i];
      line[i] = tolower(c);
      i++;
    }
   return;
}
