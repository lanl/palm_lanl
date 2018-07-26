
// ############################################################################
//
//     create_mz_kpp_module                       
//
//     create scalar code from .f90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################

#include <fstream>
#include <stdio.h>

 
#include "expand_decomp.h"
#include "utils.h"

void expand_decomp::create_sparse_info (vector<fortran_file> & include_list, string module_name) {
  vector<fortran_file>::iterator  it;
  vector<program_line>::iterator  ip;
  bool                            done;
  int                             number;
  
// Generate a dummy vector entry [0] to be compatable with FORTRAN indexing
  number = 0;
  LU_IROW.push_back(number); 
  LU_ICOL.push_back(number); 
  LU_CROW.push_back(number); 
  LU_DIAG.push_back(number); 
  
  cout <<endl;
  cout << "Sparse matrix: " ;

// Get NVAR

  for(it=include_list.begin();it!=include_list.end();it++) {
    if(it->get_name() == module_name + "_Parameters") {
      for (ip=it->pline.begin(); ip != it->pline.end(); ip++) {
        if(ip->get_token(0).substr(0,7) == "INTEGER" && ip->get_token(3) == "NVAR") {
          sscanf (ip->get_token(5).c_str(),"%i",&NVAR);
          cout << "      NVAR = "<< NVAR <<endl;
        }
      }
    }
  }
 
// Get LU_IROW

  for(it=include_list.begin();it!=include_list.end();it++) {
    if(it->get_name() == module_name + "_JacobianSP") {
      for (ip=it->pline.begin(); ip != it->pline.end(); ip++) {
        if(ip->get_token(0).substr(0,7) == "INTEGER" && 
                                   ip->get_token(7).substr(0,7) == "LU_IROW") {
          done = false;
          while(1) {
            ip++;
            for(int i=0; i<ip->get_token_size (); i++) {
              if(ip->get_token(i).substr(0,1) == "/") {
                done = true;
                break;
              } else if(ip->get_token(i).substr(0,1) != "&")  {
		  if ( sscanf (ip->get_token(i).c_str(),"%i",&number) == 1)
		      LU_IROW.push_back(number); 
              }
            }
            if(done)   break;
          }
        }
      }
    }
  }
  cout << "extracted LU_IROW:   " << LU_IROW.size()-1 << " values " <<endl;

// Get LU_ICOL

  for(it=include_list.begin();it!=include_list.end();it++) {
    if(it->get_name() == module_name + "_JacobianSP") {
      for (ip=it->pline.begin(); ip != it->pline.end(); ip++) {
        if(ip->get_token(0).substr(0,7) == "INTEGER" &&
                                   ip->get_token(7).substr(0,7) == "LU_ICOL") {
          done = false;
          while(1) {
            ip++;
            for(int i=0; i<ip->get_token_size (); i++) {
              if(ip->get_token(i).substr(0,1) == "/") {
                done = true;
                break;
              } else if(ip->get_token(i).substr(0,1) != "&")  {
		  if (sscanf (ip->get_token(i).c_str(),"%i",&number) == 1)
		      LU_ICOL.push_back(number);
              }
            }
            if(done)   break;
          }
        }
      }
    }
  }
  cout << "extracted LU_ICOL:   " << LU_ICOL.size()-1 << " values " <<endl;

// Get LU_CROW
// In case of large system, for LU_CROW and LU_DIAG might be the same logic
// necessary like LU_IROW and LU_ICOL

  for(it=include_list.begin();it!=include_list.end();it++) {
    if(it->get_name() == module_name + "_JacobianSP") {
      for (ip=it->pline.begin(); ip != it->pline.end(); ip++) {
        if(ip->get_token(0).substr(0,7) == "INTEGER" &&
                                   ip->get_token(7).substr(0,7) == "LU_CROW") {
          done = false;
          while(1) {
            ip++;
            for(int i=0; i<ip->get_token_size (); i++) {
              if(ip->get_token(i).substr(0,1) == "/") {
                done = true;
                break;
              } else if(ip->get_token(i).substr(0,1) != "&")  {
		  if (sscanf (ip->get_token(i).c_str(),"%i",&number) == 1)
		      LU_CROW.push_back(number);
              }
            }
            if(done)   break;
          }
        }
      }
    }
  }
  cout << "extracted LU_CROW:   " << LU_CROW.size()-1 << " values " <<endl;

// Get LU_DIAG

  for(it=include_list.begin();it!=include_list.end();it++) {
    if(it->get_name() == module_name + "_JacobianSP") {
      for (ip=it->pline.begin(); ip != it->pline.end(); ip++) {
        if(ip->get_token(0).substr(0,7) == "INTEGER" &&
                                   ip->get_token(7).substr(0,7) == "LU_DIAG") {
          done = false;
          while(1) {
            ip++;
            for(int i=0; i<ip->get_token_size (); i++) {
              if(ip->get_token(i).substr(0,1) == "/") {
                done = true;
                break;
              } else if(ip->get_token(i).substr(0,1) != "&")  {
		  if (sscanf (ip->get_token(i).c_str(),"%i",&number) == 1)
		      LU_DIAG.push_back(number);
              }
            }
            if(done)   break;
          }
        }
      }
    }
  }
  cout << "extracted LU_DIAG:   " << LU_DIAG.size()-1 << " values " <<endl;

  cout <<endl;

  return;
}

void expand_decomp::create_routine (vector<fortran_file> & routine_list) {

  vector<fortran_file>::iterator  it;
  fortran_file                    de;
  int                             k,kk,j,jj,i,ii;
  string                          line;
  char                            cline[80];

// Delete original KppDecomp from subroutine list

  for(it=routine_list.begin();it!=routine_list.end();it++) {
    if(it->get_name() == "KppDecomp") {
      routine_list.erase(it);
    }
  }

  de.set_name("KppDecomp");

  de.add_line ("  SUBROUTINE KppDecomp( JVS, IER)                                     ");
  de.add_line ("  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ");
  de.add_line ("  !        Sparse LU factorization                                    ");
  de.add_line ("  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ");
  de.add_line ("  ! Loop expansion generated by kp4                                   ");
  de.add_line ("                                                                      ");
  de.add_line ("    INTEGER  :: IER                                                   ");
//  de.add_line ("    REAL ( kind=dp ) :: JVS (:), W ( NVAR ), a                        ");
  de.add_line ("    REAL ( kind=dp ) :: JVS ( LU_NONZERO ), W ( NVAR ), a             ");
  de.add_line ("    INTEGER  :: k, kk, j, jj                                          ");
  de.add_line ("                                                                      ");
  de.add_line ("    a = 0.                                                            ");
  de.add_line ("    IER = 0                                                           ");
  de.add_line ("                                                                      ");

  if( kpp_switches.de_indexing () == 1 ) {
    cout << " *** de-indexing of deomposition algorithm: 1 "  <<endl;
    for (k=1; k<= NVAR; k++) {
      if(LU_CROW[k] <= LU_DIAG[k]-1)   {
        sprintf(cline,"!   k = %i",k);  
        line = cline;    de.add_line(line);
        for (kk=LU_CROW[k]; kk<=LU_CROW[k+1]-1; kk++) {
          sprintf(cline,"    W ( %i )= JVS ( %i )",LU_ICOL[kk],kk);  
          line = cline;    de.add_line(line);
        }
        for (kk=LU_CROW[k]; kk<=LU_DIAG[k]-1; kk++) {
          j = LU_ICOL[kk];
          sprintf(cline,"    a = - W ( %i ) / JVS ( %i )",j,LU_DIAG[j]);  
          line = cline;    de.add_line(line);
          sprintf(cline,"    W ( %i ) = -a",j);  
          line = cline;    de.add_line(line);
          for (jj=LU_DIAG[j]+1; jj<=LU_CROW[j+1]-1; jj++) {
            sprintf(cline,"      W ( %i ) = W ( %i ) + a * JVS ( %i )",LU_ICOL[jj],LU_ICOL[jj],jj);  
            line = cline;    de.add_line(line);
          }
        }
        for (kk=LU_CROW[k]; kk<=LU_CROW[k+1]-1; kk++) {
          sprintf(cline,"    JVS ( %i )= W ( %i )",kk,LU_ICOL[kk]);  
          line = cline;    de.add_line(line);
        }
      } else {
        sprintf(cline,"!   k = %i     Nothing to do    LU_CROW[k] > LU_DIAG[k]-1   (%i %i)",
                                                                    k,LU_CROW[k],LU_DIAG[k]-1);  
        line = cline;    de.add_line(line);
      }
    }
  }

  if( kpp_switches.de_indexing () == 2 ) {

//  No W array required
//  No data copying
//  better data reuse

    int           flag  [NVAR+1] [NVAR+1];
    int           index [NVAR+1] [NVAR+1];
    int           ij,ji,jk,ki;
    int           count;
    bool          term;

    cout << " *** de-indexing of decomposition algorithm: 2 (fast) " <<endl;

    for (i=1; i<= NVAR; i++) {
      for (j=1; j<= NVAR; j++) {
        flag[i][j] = 0;
      }
    }
    for (k=1; k<= LU_ICOL.size()-1; k++) {
      i = LU_IROW[k];
      j = LU_ICOL[k];
      flag [j][i] = 1;
      index[j][i] = k;
    }

//    for (i=1; i<= NVAR; i++) {
//      cout << "  " <<i << "    ";
//      for (j=1; j<= NVAR; j++) {
//        cout  << flag[i][j] << " ";
//      }
//      cout << endl;
//    }

    for (i=1; i<= NVAR; i++) {
      sprintf(cline,"!   i = %i",i);  
      line = cline;    de.add_line(line);
      for (j=1; j<= i-1; j++) {
        if(flag[j][i] == 1)  {
          sprintf(cline,"    a = 0.0; a = a ");  
          line = cline;
          term = false;
          count=0;
          for (k=1; k<= j-1; k++) {
            if(flag[j][k] == 1 && flag[k][i] == 1) {
              term = true;
              count++;
              if(count == 5) {
                line += "&";
                de.add_line(line);
                count = 0;
                line="         ";
              }
              jk = index [j][k];
              ki = index [k][i];
              sprintf(cline," - JVS ( %i ) * JVS ( %i )",jk,ki);  
              line += cline;
            }
          }
          ji = index [j][i];
          jj = index [j][j];
          if(term) {
            de.add_line(line);
            sprintf(cline,"    JVS ( %i ) =  ( JVS ( %i ) +a ) / JVS ( %i )  ",ji,ji,jj);  
          } else {
            sprintf(cline,"    JVS ( %i ) =  ( JVS ( %i ) ) / JVS ( %i )  ",ji,ji,jj);  
          }
          line = cline; de.add_line(line);
        }
      }
      for (j=i; j<= NVAR; j++) {
        if(flag[j][i] == 1)  {
          ij = index [j][i];
          sprintf(cline,"    JVS ( %i )= JVS ( %i )",ij,ij);  
          line = cline;
          term = false;
          count=0;
          for (k=1; k<= i-1; k++) {
            if(flag[j][k] == 1 && flag[k][i] == 1) {
              term = true;
              count++;
              if(count == 5) {
                line += "&";
                de.add_line(line);
                count = 0;
                line="         ";
              }
              jk = index [j][k];
              ki = index [k][i];
              sprintf(cline," - JVS ( %i ) * JVS ( %i )",jk,ki);  
              line += cline;
            }
          }
          if(term) {
            de.add_line(line);
          }
        }
      }
    }
  }

  if( kpp_switches.de_indexing () == 3 ) {
      cout << " *** de-indexing of decomposition algorithm: 3 (old method)"
           <<endl;
      for (k=1; k<= NVAR; k++) {
	  for (kk=LU_CROW[k]; kk<=LU_CROW[k+1]-1; kk++) {
	      sprintf(cline," W ( %i ) = JVS ( %i )", LU_ICOL[kk], kk);
	      line = cline;    de.add_line(line);
	  }

	  for (kk=LU_CROW[k]; kk<=LU_DIAG[k]-1; kk++) {
	      j = LU_ICOL[kk] ;
	      sprintf(cline, " a = - W ( %i ) / JVS ( %i )", j, LU_DIAG[j]);
	      line = cline;    de.add_line(line);

	      sprintf(cline, " W ( %i ) = -a ", j);
	      line = cline;    de.add_line(line);

	      for (jj=LU_DIAG[j]+1; jj<=LU_CROW[j+1]-1; jj++) {
		  sprintf(cline," W ( %i ) = W ( %i ) + a * JVS ( %i )",
			  LU_ICOL[jj], LU_ICOL[jj], jj);
		  line = cline;    de.add_line(line);
	      }
	  }

	  if ( LU_DIAG[k]-1 >= LU_CROW[k] ) {
	      // removes unnecessary statements
	      for (kk=LU_CROW[k]; kk<=LU_CROW[k+1]-1; kk++) {
		  sprintf(cline, " JVS ( %i ) = W ( %i )",kk, LU_ICOL[kk]);
		  line = cline;    de.add_line(line);

	      }
	  }
      }
  }

  de.add_line ("    return                                                            ");
  de.add_line ("                                                                      ");
  de.add_line ("  End SUBROUTINE KppDecomp                                            ");
  routine_list.push_back(de);

  return;
}
