#ifndef UTIL
#define UTIL 1

// ############################################################################
//
//     create_mz_kpp_module 
//
//     create scalar code from .f90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################

#include <string>
#include <list>
#include <vector>

#include <iostream>
#include <string>

using namespace std;

extern void my_abort(string s);
extern string my_to_upper(string sinp);

class Vvar {
 public:
  string                name;
  vector<string>        dim_var;

  int nr_dim() {return dim_var.size(); };
  void clear() {name.clear(); dim_var.clear(); return; };
};

class control_switches {
 private:
  bool   vector_mode;            // false     Generate scalar code
                                 // true      Generate vector code

  string vector_length;

  int de_index_mode;              // 0: no de-indexing
                                  // 1: normal de-indexing
                                  // 2: 'fast' de-indexing  (LU)
                                  // 3: 'old' de-indexing

 public:

  void   set_vector_mode()              {vector_mode = true; return;};
  void   set_vector_length (string s)   {vector_length=s; return;};
  void   set_de_index_mode(int i)       {de_index_mode = i; return;};


  bool   is_vector()           {return vector_mode;};
  string get_vector_length ()  {return vector_length;};
  int de_indexing ()           {return de_index_mode;};

  control_switches () {vector_mode = false; vector_length="1"; de_index_mode=0;} 
};

extern control_switches kpp_switches;

class string_token {
 private:
  int                         size_val;
  string                      seperator;
  vector <string>             tokens;
  vector <string>::iterator   is;
  vector <int>                position;

 public:
    void set_separator (string s) {seperator.clear();seperator=s; return; };

    void    fill_token (string s);
    bool    get_next_token(string &s) { s=*is++; return (is == tokens.end()); };
    string  get_token_by_index (int index) { if(index <= size_val-1) {
                                               return tokens[index]; 
                                             } else {
                                               return " "; }; };
    int     size()                         { return size_val; };
    void    update (int i, string s)       { tokens[i].clear(); tokens[i] = s;return;};
    int     get_position(int i)            { return position[i]; };
    void    reset (string s)               { tokens.clear(); fill_token (s); return; };


    string_token ()         {size_val=0;seperator=" ";};
    string_token (string s) {size_val=0;seperator=s;};

};

#endif
