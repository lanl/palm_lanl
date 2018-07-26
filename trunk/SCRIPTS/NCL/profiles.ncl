load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;***************************************************
; Checking the kind of the script
;***************************************************

function which_script()
local script 
begin 
   script = "profiles"
   return(script)
end

;***************************************************
; Retrieving the NCL version used
;***************************************************
   
ncl_version_ch = systemfunc("ncl -V") 
ncl_version    = stringtofloat(ncl_version_ch)

;***************************************************
; Function for checking file existence in dependence
; on NCL version
;***************************************************

function file_exists(version:string,file_name:string)
begin
   if( version .GE. "6.2.1" ) then
      existing = fileexists(file_name)
   else
      existing = isfilepresent(file_name)
   end if
   return(existing)
end
  
;***************************************************
; load .ncl.config or .ncl.config.default
;***************************************************

existing_file = file_exists(ncl_version_ch,file_config)
if (existing_file) then
   loadscript(file_config)
else
   palm_bin_path = getenv("PALM_BIN")
   print(" ")
   print("Neither the personal configuration file '.ncl.config' exists in")
   print("~/palm/current_version")
   print("nor the default configuration file '.ncl.config.default' "+\
         "exists in")
   print(palm_bin_path + "/NCL")
   print(" ")
   exit
end if

    
begin

   ;***************************************************
   ; Retrieving the double quote character
   ;***************************************************
   
   dq=str_get_dq()

   ;***************************************************
   ; set up default parameter values and strings 
   ;***************************************************

   if (no_files .LT. 1 .OR. no_files .GT. 6) then
      print(" ")
      print("Assign 'no_files' between 1 and 6")  
      print(" ")
      exit
   end if
   
   file_in   = new(no_files,string)
   file_in_1 = new(no_files,logical)
   start_f_t = new(no_files,integer)
   end_f_t   = new(no_files,integer)
   
   if (file_1 .EQ. "File in") then
      print(" ")
      print("Declare 1st input file 'file_1=' in '.ncl.config' or prompt")
      print(" ")
      exit 
   else
      file_in(0) = file_1   
   end if      
   file_in_1(0) = False
   if (isStrSubset(file_in(0), ".nc"))then
      start_f = -2
      end_f = -2 
      file_in_1(0) = True     
   end if  
   if (start_f .EQ. -1)then
      print(" ")
      print("'start_f' must be one of the cyclic numbers (at least 0) "+\
            "of your input file(s)")
      print(" ")  
      exit
   end if 
   if (end_f .EQ. -1)then
      print(" ")
      print("'end_f' must be one of the cyclic numbers (at least 0) of "+\
            "your input file(s)")
      print(" ")  
      exit
   end if            
   start_f_t(0) = start_f 
   end_f_t(0) = end_f   
 
   if (no_files .GT. 1) then
      if (file_2 .EQ. "File in") then
         print(" ")
         print("Declare 2nd input file 'file_2=' in '.ncl.config' or prompt")
         print(" ")
         exit 
      else
         file_in(1) = file_2   
      end if      
      file_in_1(1) = False
      if (isStrSubset(file_in(1), ".nc"))then
         start_f_2 = -2
         end_f_2 = -2 
         file_in_1(1) = True     
      end if
      if (start_f_2 .EQ. -1)then
         print(" ")
         print("'start_f_2' must be one of the cyclic numbers (at least 0) "+\
              "of your input file(s)")
         print(" ")  
         exit
      end if 
      if (end_f_2 .EQ. -1)then
         print(" ")
         print("'end_f_2' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if      
      start_f_t(1) = start_f_2 
      end_f_t(1) = end_f_2   
   end if    
    
   if (no_files .GT. 2) then
      if (file_3 .EQ. "File in") then
         print(" ")
         print("Declare 3rd input file 'file_3=' in '.ncl.config' or prompt")
         print(" ")
         exit 
      else
         file_in(2) = file_3   
      end if      
      file_in_1(2) = False
      if (isStrSubset(file_in(2), ".nc"))then
         start_f_3 = -2
         end_f_3 = -2 
         file_in_1(2) = True     
      end if
      if (start_f_3 .EQ. -1)then
         print(" ")
         print("'start_f_3' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if 
      if (end_f_3 .EQ. -1)then
         print(" ")
         print("'end_f_3' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if      
      start_f_t(2) = start_f_3 
      end_f_t(2) = end_f_3  
   end if 
   
   if (no_files .GT. 3) then
      if (file_4 .EQ. "File in") then
         print(" ")
         print("Declare 4th input file 'file_4=' in '.ncl.config' or prompt")
         print(" ")
         exit 
      else
         file_in(3) = file_4   
      end if      
      file_in_1(3) = False
      if (isStrSubset(file_in(3), ".nc"))then
         start_f_4 = -2
         end_f_4 = -2 
         file_in_1(3) = True     
      end if
      if (start_f_4 .EQ. -1)then
         print(" ")
         print("'start_f_4' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if 
      if (end_f_4 .EQ. -1)then
         print(" ")
         print("'end_f_4' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if      
      start_f_t(3) = start_f_4 
      end_f_t(3) = end_f_4 
   end if 
   
   if (no_files .GT. 4) then
      if (file_5 .EQ. "File in") then
         print(" ")
         print("Declare 5th input file 'file_5=' in '.ncl.config' or prompt")
         print(" ")
         exit 
      else
         file_in(4) = file_5   
      end if      
      file_in_1(4) = False
      if (isStrSubset(file_in(4), ".nc"))then
         start_f_5 = -2
         end_f_5 = -2 
         file_in_1(4) = True     
      end if
      if (start_f_5 .EQ. -1)then
         print(" ")
         print("'start_f_5' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if 
      if (end_f_5 .EQ. -1)then
         print(" ")
         print("'end_f_5' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if      
      start_f_t(4) = start_f_5 
      end_f_t(4) = end_f_5    
   end if 
   
   if (no_files .GT. 5) then
      if (file_6 .EQ. "File in") then
         print(" ")
         print("Declare 6th input file 'file_6=' in '.ncl.config' or prompt")
         print(" ")
         exit 
      else
         file_in(5) = file_6   
      end if      
      file_in_1(5) = False
      if (isStrSubset(file_in(5), ".nc"))then
         start_f_6 = -2
         end_f_6 = -2 
         file_in_1(5) = True     
      end if
      if (start_f_6 .EQ. -1)then
         print(" ")
         print("'start_f_6' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if 
      if (end_f_6 .EQ. -1)then
         print(" ")
         print("'end_f_6' must be one of the cyclic numbers (at least 0) "+\
               "of your input file(s)")
         print(" ")  
         exit
      end if      
      start_f_t(5) = start_f_6 
      end_f_t(5) = end_f_6    
   end if 

   
   if (format_out .NE. "x11" .AND. format_out .NE. "pdf" .AND.  \
      format_out .NE. "eps" .AND. format_out .NE. "ps" .AND.    \
      format_out .NE. "epsi" .AND. format_out .NE. "ncgm" .AND. \
      format_out .NE. "png")then
      print(" ")
      print("'format_out = "+format_out+"' is invalid and set to'x11'")
      print(" ")
      format_out="x11"
   end if  

   if (ncl_version .LE. 5.1 .AND. format_out .EQ. "png") then
      print(" ")
      print("Output of png files not available")
      print("png output is avaiable with NCL version 5.2.0 and higher ")
      print("NCL version used: " + ncl_version_ch)
      print(" ")
      exit
   end if

   if (over .NE. 0 .AND. over .NE. 1) then
      print(" ")
      print("'over'= "+over+" is invalid and set to 0")
      print(" ")
      over = 0
   end if
   if (no_files .GT. 1) then
      over = 0
      print(" ")
      print("If you have more than one input file - you cannot overlay "+\
            "variables: over is set to 0")
      print(" ")
   end if
   if (no_files .GT. 1 .AND. var .EQ. ",cross_profiles,") then
      var = "all"
      print(" ")
      print("If you have more than one input file - you cannot combine "+\
            "profiles with cross_profiles ")
      print(" ")
   end if

   ;************************************************************************
   ; check if cross_profiles is defined in file_1
   ;************************************************************************

    if (isStrSubset(file_in(0), ".nc")) then
      firstfile=addfile(file_1,"r")
    else
      firstfile=addfile(file_in(0) + "." + start_f + ".nc","r")
    end if

   cross_check=isatt(firstfile,"cross_profiles")

   if (var .EQ. "all") then
      VAR_LIST = firstfile@VAR_LIST
      no_plots = str_fields_count(VAR_LIST, ";")
   else
      no_plots = str_fields_count(var, ",")
   end if

   if ((.not. cross_check) .AND. (var .EQ. ",cross_profiles,")) then
       print(" ")
       print("cross_profiles is not defined in the NetCDF header."+\
             "Using var='all' instead.")
       print(" ") 
       var="all"
   end if

   wks_gsn_log = True

   ;************************************************************************
   ; check if cross_profiles can be used
   ;************************************************************************
   
   VAR_LIST = firstfile@VAR_LIST
   if ( var .EQ. ",cross_profiles," .AND. str_fields_count(VAR_LIST,"_") \
        .GT. str_fields_count(VAR_LIST,";") ) then
       var = "all"
       print(" ")
       print("It is not possible to use var='cross_profiles' for"+\
             " statistic regions. Using var='all' instead.")
       print(" ")
   end if

   if ( var .EQ. ",cross_profiles," .AND. cross_check .AND. (over .EQ. 0) \
       .AND. (no_files .LT. 2)) then 
      c_var_log = 1
   else
      c_var_log = 0
   end if

   if ((no_rows .EQ. -1) .AND. isatt(firstfile,"no_rows")) then
      no_rows = firstfile@no_rows
   else if ((no_rows .EQ. -1) .AND. .not. isatt(firstfile,"no_rows")) then
       print(" ")
       print("no_rows is not defined in the NetCDF header."+\
             "Using no_rows=3 instead.")
       print(" ")
      no_rows = 3
   end if
   end if

   if ((no_columns .EQ. -1) .AND. isatt(firstfile,"no_columns")) then
      no_columns = firstfile@no_columns
   else if ((no_columns .EQ. -1) .AND. .not. isatt(firstfile,"no_columns")) then
       print(" ")
       print("no_columns is not defined in the NetCDF header."+\
             "Using no_columns=2 instead.")
       print(" ")
      no_columns = 2
   end if
   end if

   if (c_var_log .EQ. 1) then 
      c_var_str = firstfile@cross_profiles

      c_var_all = str_split(c_var_str,";")

      c_var_all = ","+c_var_all+","

      ;************************************************************************
      ; check if an additional combined plot is demanded by the user
      ;************************************************************************

      if ((c_var .EQ. "c_variables") .AND. (number_comb .EQ. -1) .AND. \
             (combine .EQ. 0)) then
         max_com_i = max(dimsizes(c_var_all)) - 1
      else
         max_com_i = max(dimsizes(c_var_all))
         c_var_all_temp = c_var_all
         delete(c_var_all)
         c_var_all = new(max_com_i+1,string)
         c_var_all(0:max_com_i-1) = c_var_all_temp
         delete(c_var_all_temp)
         c_var_all(max_com_i) = c_var
      end if
      number_comb_all = str_fields_count(c_var_all,",")

      ;************************************************************************
      ; currently, more than 3 plots cannot be drawn together. the profiles 
      ; are split to a maximum of 3 plots per panel.
      ;***********************************************************************

      do com_i = 0, max_com_i
         if(number_comb_all(com_i) .GT. 3) then
            print(" ")
            print("Currently, it is not possible to plot more than three "+\
                  "profiles together ("+c_var_all(com_i)+"). Hence, they "+\
                  "are split to a maximum of three profiles per "+\
                  "coordinate plane.")
            print(" ") 
	    c_var_all_temp = c_var_all
            number_comb_all_temp = number_comb_all
	    delete(c_var_all)
            delete(number_comb_all)
            max_com_i = max_com_i + 1 
	    c_var_all = new(max_com_i+1,string)
            number_comb_all = new(max_com_i+1, integer)
            if (com_i-1 .GE. 0) then
	      c_var_all(0:com_i-1) = c_var_all_temp(0:com_i-1)
	      number_comb_all(0:com_i-1) = number_comb_all_temp(0:com_i-1)
            end if
	    c_var_long = str_split(c_var_all_temp(com_i),",")
            c_var_all(com_i) = "," + c_var_long(0) + "," + c_var_long(1) + "," + c_var_long(2) + ","
            number_comb_all(com_i) = 3
            c_var_all(com_i + 1) = ","
            do com_j = 3, max(dimsizes(c_var_long)) - 1 
               c_var_all(com_i + 1) = c_var_all(com_i + 1) + c_var_long(com_j) + ","
            end do
            number_comb_all(com_i+1) = max(dimsizes(c_var_long)) - 3 
            if (com_i+2 .le. max_com_i) then
	       c_var_all(com_i+2:max_com_i) = c_var_all_temp(com_i+1:max_com_i-1)
               number_comb_all(com_i+2:max_com_i) = number_comb_all_temp(com_i+1:max_com_i-1)
            end if
	    delete(c_var_all_temp)
            delete(number_comb_all_temp)
            delete(c_var_long)
	 end if
      end do

      c_var_plot = new(max_com_i+1,graphic)
      c_var_end_time_step = end_time_step
      c_var_start_time_step = start_time_step
      number_comb_1_log = 0
      number_comb_23_log = 0
   else
      max_com_i = 0
   end if

   do com_i = 0,max_com_i

      ;************************************************************************
      ; changes basic properties of the plots according to cross_profiles
      ;************************************************************************

      if (c_var_log .EQ. 1) then
         if (number_comb_all(com_i) .EQ. 1) then
            number_combine= -1
	    combine = 0
	    c_var = "c_variables"
	    var = c_var_all(com_i)
	    delete(end_time_step)
	    delete(start_time_step)
	    end_time_step = c_var_end_time_step
            start_time_step= c_var_start_time_step
	    if (number_comb_1_log .EQ. 1) then
		delete(var_char)
	    end if
	    number_comb_1_log=1
	 else
	    var = "all"
	    number_comb = number_comb_all(com_i)
	    combine = 1
	    c_var = c_var_all(com_i)
	    delete(end_time_step)
	    delete(start_time_step)
	    end_time_step = c_var_end_time_step
            start_time_step= c_var_start_time_step
            if ((number_comb_1_log .EQ. 1) .OR. (number_comb_23_log .EQ. 1)) then
		delete(label)
	    end if
	    if (number_comb_23_log .EQ. 1) then
		delete(plot_o)
;		delete(label)
		delete(color_o)
		delete(mini)
		delete(maxi)
		delete(plot)
	    end if
	    number_comb_23_log=1
	 end if
      end if

      if (combine .NE. 0 .AND. combine .NE. 1)then
	  print(" ")
	  print("Your 'combine'= "+combine+" is invalid and set to 0")
	  print(" ")
	  combine = 0
      end if
      if (no_files .GT. 1) then
	  combine = 0
	  print(" ")
	  print("If you have more than one input file you cannot combine "+\
		"variables: combine is set to 0")
	  print(" ")
      end if   
      if (combine .EQ. 1 .AND. number_comb .EQ. -1) then
	  print(" ") 
	  print("Set 'number_comb' to 2 or 3 or combine to 0")
	  print(" ")
	  exit
      end if 
      if (combine .EQ. 1)then
	  if (number_comb .EQ. -1) then
	    print(" ") 
	    print("Set 'number_comb' to 2 or 3 or combine to 0")
	    print(" ")
	    exit
	  end if
	  if (number_comb .NE. 2 .AND. number_comb .NE. 3) then
	    print(" ") 
	    print("Set 'number_comb' to 2 or 3 or combine to 0")
	    print(" ")
	    exit
	  end if
      end if  
      if (combine .EQ. 1 .AND. c_var .EQ. "c_variables") then
	  print(" ") 
	  print("Select variables for overlaying ('c_var') or set combine to 0")
	  print(" ")
	  exit
      end if

      if (log_z .NE. 0 .AND. log_z .NE. 1)then
	  print(" ")
	  print("'log_z'= "+log_z+" is invalid and set to 0")
	  print(" ")
	  log_z = 0
      end if    
      
      if (norm_z .EQ. 0) then 
	  print(" ")
	  print("Normalising with 0 is not allowed, 'norm_z' is set to 1.0")
	  print(" ")
	  norm_z = 1.0
      end if

      ;***************************************************
      ; open input file
      ;***************************************************

      do nof=0,no_files-1

      if (com_i .EQ. 0)
      files=new(end_f_t(nof)-start_f_t(nof)+1,string)
      if (file_in_1(nof))then
	  if (isfilepresent(file_in(nof)))then
	    files(0)=file_in(nof)
	  else
	    print(" ")
	    print("1st input file: '"+file_in(nof)+"' does not exist")
	    print(" ")
	    exit
	  end if
      else
	  if (start_f_t(nof) .EQ. 0)then
	    if (isfilepresent(file_in(nof)+".nc"))then
		files(0)=file_in(nof)+".nc"
		do i=1,end_f_t(nof)
		  if (isfilepresent(file_in(nof)+"."+i+".nc"))then    
		      files(i)=file_in(nof)+"."+i+".nc"
		  else
		      print(" ")
		      print("Input file: '"+file_in(nof)+"."+i+".nc' does not "+\
			    "exist")
		      print(" ")
		      exit  
		  end if	    
		end do         
	    else
		print(" ")
		print("Input file: '"+file_in(nof)+".nc' does not exist")
		print(" ")
		exit
	    end if
	  else
	    do i=start_f_t(nof),end_f_t(nof)
		if (isfilepresent(file_in(nof)+"."+i+".nc"))then    
		  files(i-start_f_t(nof))=file_in(nof)+"."+i+".nc"
		else
		  print(" ")
		  print("Input file: '"+file_in(nof)+"."+i+".nc' does not exist")
		  print(" ")
		  exit  
		end if
	    end do
	  end if
      end if

      f=addfiles(files,"r")
      f_att=addfile(files(0),"r")
      ListSetType(f,"cat")
      
      vNam  = getfilevarnames(f_att)
      if( ncl_version .GE. 6.1 ) then
         vNam = vNam(::-1)
      end if
      vType = getfilevartypes(f_att,vNam) 

      if ((all(vType .eq. "double"))) then ;distinction if data is double or float
	  check_vType = True
      else
	  check_vType = False
      end if
      
      if (nof .EQ. 0)then
	  vNam0=vNam
      end if
      if (nof .NE. 0)then
	  if (dim0 .NE. dim)then
	    print(" ")
	    print("There are 'no_files'="+no_files+" input files but they do "+\
		  "not contain the same variables")
	    print(" ")
	    exit
	  else
	    do i=0,dim0-1
		if (vNam0(i) .NE. vNam(i))then
		  print(" ")
		  print("There are 'no_files'="+no_files+" input files but "+\
			"they do not contain the same variables")
		  print(" ")
		  exit
		end if
	    end do
	  end if
      end if
      nof=nof+1
      if (com_i .EQ. 0) then
	  print(" ")
	  print("Variables in input file "+nof+":")
	  print("- "+ vNam)
	  print(" ")
      end if
      nof=nof-1
      dim = dimsizes(vNam)
      dim0=dim
      
      if (dim .EQ. 0) then
	  print(" ")
	  print("There is no data on file")
	  print(" ")
      end if
      
      prof3d = 0
      do varn = dim-1,0,1
	  if ( isStrSubset( vNam(varn), "zu_3d") .OR. \
	      isStrSubset( vNam(varn), "zw_3d")) then
	    prof3d = 1
	    break
	  end if
      end do
      
      end if

      if (var .NE. "all" ) then
	  ;rearrange the order of the variables in vNam so that the variables 
	  ;specified by "var" are in the top of vNam

	  vNam_static = new((/dim/),string)
	  vNam_temp   = new((/dim/),string)

	  var_char    = stringtocharacter(var)
	  no_char     = dimsizes(var_char)
	  comma       = 0
	  
	  do j=0,no_char-1
	    if(var_char(j) .eq. ",")
		comma = comma + 1
	    end if
	  end do

	  if(comma .le. 1)
	      print(" ")
	      print("The variables 'var="+var+"'" )
	      print("do not exist on your input file;")
	      print("be sure to have one comma before and after each variable")
	      print(" ")
	      exit
	  end if

	  if(comma .gt. dim)
	      print(" ")
	      print("The variables 'var="+var+"'" )
	      print("do not exist on your input file;")
	      print("be sure to have one comma before and after each variable")
	      print(" ")
	      exit
	  end if

	  indices = new((/comma/),integer)
	  comma   = 0
	    
	  do j=0,no_char-1
	    if(var_char(j) .eq. ",")
	      indices(comma) = j
	      comma          = comma + 1
	    end if
	  end do

	  do j=0,comma-2
	    vNam_temp(j) = charactertostring(\
				      var_char(indices(j)+1:indices(j+1)-1))
	  end do

	  do j=0,comma-2
	    count_check = 0
	    do varn=0,dim-1 
	      if(vNam_temp(j) .ne. vNam(varn)) 
		count_check=count_check+1
	      end if
	    end do    

	    if(count_check .eq. dim)
		print(" ")
		print("The variables 'var="+var+"'" )
		print("do not exist on your input file;")
		print("be sure to have one comma before and after each variable")
		print(" ")
		exit
	    end if

	    vNam_static(j) = vNam_temp(comma - 2 - j)
	  end do

	  vNam_temp   = vNam_static
	  vNam_static = ""

	  if (prof3d .EQ. 0) then
	    i = 0
	    do j=0,dim-1,2
		vNam_static(j)   = vNam_temp(i)
		vNam_static(j+1) = "z" + vNam_static(j)
		i = i+1
		if(i .eq. comma-1)
		    break
		end if
	    end do

	    delete(vNam_temp)
	    count_variable = 2*(comma-1) 
	  else
	    i = 0
	      do j=0,dim-1,1
		vNam_static(j) = vNam_temp(i)
		i = i+1
		if(i .eq. comma-1)
		    break
		end if
	      end do

	      delete(vNam_temp)
	      count_variable = comma-1 
	  end if

	  counter  = 0
	  counter2 = 0

	  do j=0,dim-1
	    counter = 0
	    do i = 0, count_variable - 1        
		if(vNam_static(i) .ne. vNam(j))
		  counter = counter + 1
		end if
	    end do
	  
	    if (counter .eq. count_variable)
		vNam_static(count_variable + counter2) = vNam(j)
		counter2 = counter2 + 1
	    end if
	  end do
	
	  vNam = vNam_static 
	  delete(vNam_static)
	end if

      ;---------below steps only for first file -> nof=0
      if (nof .EQ. 0) then

	  plot = new(dim,graphic)
	  plot_ = new(dim,graphic)
	  if (no_files .GT. 1) then
	    multi_plot = new((/no_files,dim/),graphic)
	    if (check_vType) then
	      max_nof = new((/no_files,dim/),double)
	      min_nof = new((/no_files,dim/),double)
	    else
	      max_nof = new((/no_files,dim/),float)
	      min_nof = new((/no_files,dim/),float)
	    end if
	    name    = new((/no_files,dim/),string)
	    unit_   = new((/no_files,dim/),string)
	  end if

      if (prof3d .EQ. 0) then
    
      do varn=0,dim-1
	  if ( isStrSubset( vNam(varn), "time") .OR. \
	      isStrSubset( vNam(varn), "NORM")) then
	    continue
	  end if
	  if (vNam(varn) .EQ. "u" .OR. isStrSubset(vNam(varn), "u_"))then
	    z_u = f_att->$vNam(varn+1)$
	    break
	  else
	    if (vNam(varn) .EQ. "v" .OR. isStrSubset(vNam(varn), "v_"))then  
		z_u = f_att->$vNam(varn+1)$
		break
	    else
		if(vNam(varn) .EQ. "pt" .OR. isStrSubset(vNam(varn), "pt_"))then
		  z_u = f_att->$vNam(varn+1)$
		  break
		else
		  if(vNam(varn) .EQ. "vpt" .OR. isStrSubset(vNam(varn), "vpt_"))then
		      z_u = f_att->$vNam(varn+1)$
		      break
		  else
		      if(vNam(varn) .EQ. "lpt" .OR. isStrSubset(vNam(varn), "lpt_"))then
			z_u = f_att->$vNam(varn+1)$ 
			break
		      else   
			if(vNam(varn) .EQ. "q" .OR. isStrSubset(vNam(varn), "q_") )then
			    z_u = f_att->$vNam(varn+1)$ 
			    break
			else
			    if(vNam(varn) .EQ. "qv" .OR. isStrSubset(vNam(varn), "qv_"))then
			      z_u = f_att->$vNam(varn+1)$ 
			      break
			    else
			      if(vNam(varn) .EQ. "ql" .OR. isStrSubset(vNam(varn), "ql_"))then
				  z_u = f_att->$vNam(varn+1)$ 
				  break
			      else
				  if(vNam(varn) .EQ. "rho" .OR. isStrSubset(vNam(varn), "rho_"))then
				    z_u = f_att->$vNam(varn+1)$ 
				    break
				  else
				    if(vNam(varn) .EQ. "s" .OR. isStrSubset(vNam(varn), "s_") )then
					z_u = f_att->$vNam(varn+1)$ 
					break
				    else
					if(vNam(varn) .EQ. "sa" .OR. isStrSubset(vNam(varn), "sa_"))then
					  z_u = f_att->$vNam(varn+1)$ 
					  break
					else
					  if(vNam(varn) .EQ. "e" .OR. isStrSubset(vNam(varn), "e_"))then
					    
					      z_u = f_att->$vNam(varn+1)$ 
					      break
					  else
					      if(vNam(varn) .EQ. "es" .OR. \
					      isStrSubset(vNam(varn), "es_") .OR. vNam(varn) .EQ. "e*" .OR. isStrSubset(vNam(varn), "e*_"))then
						z_u = f_att->$vNam(varn+1)$ 
						break
					      else
						if(vNam(varn) .EQ. "km" .OR. isStrSubset(vNam(varn), "km_"))then
						    z_u = f_att->$vNam(varn+1)$ 
						    break
						else
						    if(vNam(varn) .EQ. "kh" .OR. isStrSubset(vNam(varn), "kh_"))then
						      z_u = f_att->$vNam(varn+1)$ 
						      break
						    else
						      if(vNam(varn) .EQ. "l" .OR. isStrSubset(vNam(varn), "l_"))then
							  z_u = f_att->$vNam(varn+1)$ 
							  break
						      else
							  if(vNam(varn) .EQ. "us2" .OR. isStrSubset(vNam(varn), "us2_")\
							    .OR. vNam(varn) .EQ. "u*2" .OR. isStrSubset(vNam(varn), "u*2_"))then
							    z_u = f_att->$vNam(varn+1)$ 
							    break
							  else
							    if(vNam(varn) .EQ. "vs2" .OR. isStrSubset(vNam(varn), "vs2_")\
								.OR. vNam(varn) .EQ. "v*2" .OR. isStrSubset(vNam(varn), "v*2_"))then
								z_u = f_att->$vNam(varn+1)$ 
								break
							    else
								if(vNam(varn) .EQ. "pts2" .OR. isStrSubset(vNam(varn), "pts2_")\
								  .OR. vNam(varn) .EQ. "pt*2" .OR. isStrSubset(vNam(varn), "pt*2_"))then
								  z_u = f_att->$vNam(varn+1)$ 
								  break
								else
								  if(vNam(varn) .EQ. "wsususodz" .OR. isStrSubset(vNam(varn), "wsususodz_")\
								      .OR. vNam(varn) .EQ. "w*u*u*:dz" .OR. isStrSubset(vNam(varn), "w*u*u*:dz_"))then
								      z_u = f_att->$vNam(varn+1)$ 
								      break
								  else
								      if(vNam(varn) .EQ. "wspsodz" .OR. isStrSubset(vNam(varn), "wspsodz_")\
									.OR. vNam(varn) .EQ. "w*p*:dz" .OR. isStrSubset(vNam(varn), "w*p*:dz_"))then
									z_u = f_att->$vNam(varn+1)$ 
									break
								      else
									if(vNam(varn) .EQ. "wpeodz" .OR. isStrSubset(vNam(varn), "wpeodz_")\
									    .OR. vNam(varn) .EQ. "w"+dq+"e:dz" .OR. isStrSubset(vNam(varn), "w"+dq+"e:dz_"))then
									    z_u = f_att->$vNam(varn+1)$
									    break   
									else
									    if(vNam(varn) .EQ. "qs2" .OR. isStrSubset(vNam(varn), "qs2_")\
									      .OR. vNam(varn) .EQ. "q*2" .OR. isStrSubset(vNam(varn), "q*2_"))then
									      z_u = f_att->$vNam(varn+1)$
									      break  
									  end if               
									end if
								      end if
								  end if
								end if
							    end if
							  end if
						      end if
						    end if
						end if
					      end if
					  end if
					end if
				    end if
				  end if
			      end if
			    end if
			end if
		      end if
		  end if
		end if 
	    end if
	  end if
      end do

      do varn=0,dim-1
	  if ( isStrSubset( vNam(varn), "time") .OR. isStrSubset( vNam(varn), "NORM")) then
	    continue
	  end if 
	  if (vNam(varn) .EQ. "w" .OR. isStrSubset(vNam(varn), "w_"))then
	    z_w = f_att->$vNam(varn+1)$
	    break
	  else
	    if (vNam(varn) .EQ. "wpup" .OR. isStrSubset(vNam(varn), "wpup_")\
		.OR. vNam(varn) .EQ. "w"+dq+"u"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"u"+dq+"_"))then
		z_w = f_att->$vNam(varn+1)$
		break
	    else
		if(vNam(varn) .EQ. "wsus" .OR. isStrSubset(vNam(varn), "wsus_")\
		  .OR. vNam(varn) .EQ. "w*u*" .OR. isStrSubset(vNam(varn), "w*u*_"))then
		  z_w = f_att->$vNam(varn+1)$
		  break
		else
		  if(vNam(varn) .EQ. "wu" .OR. isStrSubset(vNam(varn), "wu_"))then
		      z_w = f_att->$vNam(varn+1)$ 
		      break
		  else
		      if(vNam(varn) .EQ. "wpvp" .OR. isStrSubset(vNam(varn), "wpvp_")\
			.OR. vNam(varn) .EQ. "w"+dq+"v"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"v"+dq+"_"))then
			z_w = f_att->$vNam(varn+1)$ 
			break
		      else   
			if(vNam(varn) .EQ. "wsvs" .OR. isStrSubset(vNam(varn), "wsvs_")\
			    .OR. vNam(varn) .EQ. "w*v*" .OR. isStrSubset(vNam(varn), "w*v*_"))then
			    z_w = f_att->$vNam(varn+1)$ 
			    break
			else
			    if(vNam(varn) .EQ. "wv" .OR. isStrSubset(vNam(varn), "wv_"))then
			      z_w = f_att->$vNam(varn+1)$ 
			      break
			    else
			      if(vNam(varn) .EQ. "wptpp" .OR. isStrSubset(vNam(varn), "wptpp_")\
				  .OR. vNam(varn) .EQ. "w"+dq+"pt"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"pt"+dq+"_"))then
				  z_w = f_att->$vNam(varn+1)$ 
				  break
			      else
				  if(vNam(varn) .EQ. "wspts" .OR. isStrSubset(vNam(varn), "wspts_")\
				    .OR. vNam(varn) .EQ. "w*pt*" .OR. isStrSubset(vNam(varn), "w*pt*_"))then
				    z_w = f_att->$vNam(varn+1)$ 
				    break
				  else
				    if(vNam(varn) .EQ. "wpt" .OR. isStrSubset(vNam(varn), "wpt_"))then
					z_w = f_att->$vNam(varn+1)$
					break
				    else
					if(vNam(varn) .EQ. "wsptsBC" .OR. isStrSubset(vNam(varn), "wsptsBC_")\
					  .OR. vNam(varn) .EQ. "w*pt*BC" .OR. isStrSubset(vNam(varn), "w*pt*BC_"))then
					  z_w = f_att->$vNam(varn+1)$ 
					  break
					else
					  if(vNam(varn) .EQ. "wptBC" .OR. isStrSubset(vNam(varn), "wptBC_"))then
					      z_w = f_att->$vNam(varn+1)$ 
					      break
					  else
					      if(vNam(varn) .EQ. "wpvptp" .OR. isStrSubset(vNam(varn), "wpvptp_")\
						  .OR. vNam(varn) .EQ. "w"+dq+"vpt"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"vpt"+dq+"_"))then
						z_w = f_att->$vNam(varn+1)$ 
						break
					      else
						if(vNam(varn) .EQ. "wsvpts" .OR. isStrSubset(vNam(varn), "wsvpts_")\
						  .OR. vNam(varn) .EQ. "w*vpt*" .OR. isStrSubset(vNam(varn), "w*vpt*_"))then
						    z_w = f_att->$vNam(varn+1)$
						    break
						else
						    if(vNam(varn) .EQ. "wvpt" .OR. isStrSubset(vNam(varn), "wvpt_"))then
						      z_w = f_att->$vNam(varn+1)$ 
						      break
						    else
						      if(vNam(varn) .EQ. "wpqp" .OR. isStrSubset(vNam(varn), "wpqp_")\
							  .OR. vNam(varn) .EQ. "w"+dq+"q"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"q"+dq+"_"))then
							  z_w = f_att->$vNam(varn+1)$ 
							  break
						      else
							  if(vNam(varn) .EQ. "wsqs" .OR. isStrSubset(vNam(varn), "wsqs_")\
							    .OR. vNam(varn) .EQ. "w*q*" .OR. isStrSubset(vNam(varn), "w*q*_"))then
							    z_w = f_att->$vNam(varn+1)$ 
							    break
							  else
							    if(vNam(varn) .EQ. "wq" .OR. isStrSubset(vNam(varn), "wq_"))then
								z_w = f_att->$vNam(varn+1)$ 
								break
							    else
								if(vNam(varn) .EQ. "wpqvp" .OR. isStrSubset(vNam(varn), "wpqvp_")\
								  .OR. vNam(varn) .EQ. "w"+dq+"qv"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"qv"+dq+"_"))then
								  z_w = f_att->$vNam(varn+1)$ 
								  break
								else
								  if(vNam(varn) .EQ. "wsqvs" .OR. isStrSubset(vNam(varn), "wsqvs_")\
								      .OR. vNam(varn) .EQ. "w*qv*" .OR. isStrSubset(vNam(varn), "w*qv*_"))then
								      z_w = f_att->$vNam(varn+1)$ 
								      break
								  else
								      if(vNam(varn) .EQ. "wqv" .OR. isStrSubset(vNam(varn), "wqv_"))then
									z_w = f_att->$vNam(varn+1)$ 
									break
								      else
									if(vNam(varn) .EQ. "wpsp" .OR. isStrSubset(vNam(varn), "wpsp_")\
									    .OR. vNam(varn) .EQ. "w"+dq+"s"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"s"+dq+"_"))then
									    z_w = f_att->$vNam(varn+1)$
									    break 
									else
									    if(vNam(varn) .EQ. "wsss" .OR. isStrSubset(vNam(varn), "wsss_")\
										.OR. vNam(varn) .EQ. "w*s*" .OR. isStrSubset(vNam(varn), "w*s*_"))then
									      z_w = f_att->$vNam(varn+1)$
									      break 
									    else
									      if(vNam(varn) .EQ. "ws" .OR. isStrSubset(vNam(varn), "ws_"))then
										  z_w = f_att->$vNam(varn+1)$
										  break 
									      else
										  if(vNam(varn) .EQ. "wpsap" .OR. isStrSubset(vNam(varn), "wpsap_")\
										      .OR. vNam(varn) .EQ. "w"+dq+"sa"+dq .OR. isStrSubset(vNam(varn), "w"+dq+"sa"+dq+"_"))then
										    z_w = f_att->$vNam(varn+1)$
										    break 
										  else
										    if(vNam(varn) .EQ. "wssas" .OR. isStrSubset(vNam(varn), "wssas_")\
											.OR. vNam(varn) .EQ. "w*sa*" .OR. isStrSubset(vNam(varn), "w*sa*_"))then
											z_w = f_att->$vNam(varn+1)$
											break 
										    else
											if(vNam(varn) .EQ. "wsa" .OR. isStrSubset(vNam(varn), "wsa_"))then
											  z_w = f_att->$vNam(varn+1)$
											  break 
											else
											  if(vNam(varn) .EQ. "wses" .OR. isStrSubset(vNam(varn), "wses_")\
											      .OR. vNam(varn) .EQ. "w*e*" .OR. isStrSubset(vNam(varn), "w*e*_"))then
											      z_w = f_att->$vNam(varn+1)$
											      break 
											  else
											      if(vNam(varn) .EQ. "ws2" .OR. isStrSubset(vNam(varn), "ws2_")\
												.OR. vNam(varn) .EQ. "w*2" .OR. isStrSubset(vNam(varn), "w*2_"))then
												z_w = f_att->$vNam(varn+1)$
												break 
											      else
												if(vNam(varn) .EQ. "ws3" .OR. isStrSubset(vNam(varn), "ws3_")\
												    .OR. vNam(varn) .EQ. "w*3" .OR. isStrSubset(vNam(varn), "w*3_"))then
												    z_w = f_att->$vNam(varn+1)$
												    break 
												else
												    if(vNam(varn) .EQ. "Sw" .OR. isStrSubset(vNam(varn), "Sw_"))then
												      z_w = f_att->$vNam(varn+1)$
												      break 
												    else
												      if(vNam(varn) .EQ. "ws2pts".OR. isStrSubset(vNam(varn), "ws2pts_")\
													  .OR. vNam(varn) .EQ. "w*2pt*" .OR. isStrSubset(vNam(varn), "w*2pt*_"))then
													  z_w = f_att->$vNam(varn+1)$
													  break 
												      else
													  if(vNam(varn) .EQ. "wspts2" .OR. isStrSubset(vNam(varn), "wspts2_")\
													    .OR. vNam(varn) .EQ. "w*pt*2" .OR. isStrSubset(vNam(varn), "w*pt*2_"))then
													    z_w = f_att->$vNam(varn+1)$
													    break                                            
													  end if
												      end if 
												    end if
												end if
											      end if
											  end if
											end if
										    end if
										  end if
									      end if
									    end if    
									end if
								      end if
								  end if
								end if
							    end if
							  end if
						      end if
						    end if
						end if
					      end if
					  end if
					end if
				    end if
				  end if
			      end if
			    end if
			end if
		      end if
		  end if
		end if 
	    end if
	  end if
      end do
      
      if (.not. isvar("z_u") .AND. .not. isvar ("z_w")) then
	  co=0
	  do varn=0,dim-1      
	    if ( isStrSubset( vNam(varn), "time") .OR. \
		  isStrSubset( vNam(varn), "NORM")) then
		check = False
	    else
		check = True
		if (var .NE. "all") then
		  check = isStrSubset( var,","+vNam(varn)+"," )
		end if          
	    end if
	    if (check)then
		co=co+1
		z = f_att->$vNam(varn+1)$
		if (getvardims(z) .EQ. "zu")then
		  z_u = z
		else
		  if (getvardims(z) .EQ. "zw")then
		      z_w = z
		  end if
		end if
		dimz = dimsizes(z)
		break
	    end if   
	  end do
	  if (co .EQ. 0) then
	    print(" ")
	    print("The variables 'var="+var+"'" )
	    print("do not exist on your input file;")
	    print("be sure to have one comma before and after each variable")
	    print(" ")
	    exit           
	  end if
      end if

      if (isvar("z_u") ) then 
	  dimz  = dimsizes(z_u)
      else
	  if (isvar("z_w"))then 
	    dimz  = dimsizes(z_w)
	  end if
      end if

      else

	  do varn = dim-1,0,1
	    if (vNam(varn) .EQ. "zu_3d")then
		z_u = f_att->zu_3d  
		dimz  = dimsizes(z_u)         
	    else 
		if (vNam(varn) .EQ. "zw_3d")then
		  z_w = f_att->zw_3d
		  dimz  = dimsizes(z_w)
		end if
	    end if
	    if (vNam(varn) .EQ. "x")then
		x = f_att->x
		dimx=dimsizes(x)
	    else
		if (vNam(varn) .EQ. "xu")then
		  x = f_att->xu
		  dimx=dimsizes(x)
		end if
	    end if
	    if (vNam(varn) .EQ. "y")then
		y = f_att->y
		dimy=dimsizes(y)
	    else
		if (vNam(varn) .EQ. "yv")then
		  y = f_att->yv
		  dimy=dimsizes(y)
		end if
	    end if
	  end do

      end if
;
;--   Since no other z-variables have been found, a z-variable has been derived 
;--   from the user profiles 
      if (.not. isvar("z_u") .AND. .not. isvar ("z_w")) then
         do varn=0,dim-2  
            if( isStrSubset( vNam(varn+1), vNam(varn)) ) then
               z_u = f_att->$vNam(varn+1)$
               break
            end if
         end do
      end if

      if(.not. isvar("z_u") .AND. .not. isvar("z_w"))then
	  print(" ")
	  print("Program aborts - there are no z-variables available")
	  print("Be sure if 'plot_3d' is set correctly")
	  print(" ")
	  exit
      end if
      

      t_all = f[:]->time
      nt    = dimsizes(t_all)

      
      if (nt .EQ. 1)then
	  delta_t=t_all(nt-1)/nt
      else
	  delta_t=(t_all(nt-1)-t_all(0))/(nt-1)
      end if

      ; ****************************************************	
      ; start of time step and different types of mistakes that could be done
      ; ****************************************************
      
      if (start_time_step .EQ. -1.) then   
	  delete(start_time_step)	
	  start_time_step=t_all(0)/3600    
      else 
	  if (start_time_step .GT. t_all(nt-1)/3600)then
	    print(" ")
	    print("'start_time_step' = "+ start_time_step +"h is greater "+\
		  "than last time step = " \
		  + t_all(nt-1)+"s = "+t_all(nt-1)/3600+"h")
	    print(" ")
	    print("Select another 'start_time_step'")
	    print(" ")
	    exit
	  end if
	  if (start_time_step .LT. t_all(0)/3600)then
	    print(" ")
	    print("'start_time_step' = "+ start_time_step +"h is lower "+\
		  "than first time step = " + t_all(0)+"s = "+t_all(0)/3600+"h")
	    print(" ")
	    exit
	  end if
      end if

      do i=0,nt-1    
	  if (start_time_step .GE. (t_all(i)-delta_t/2)/3600 .AND.\
	      start_time_step .LT. (t_all(i)+delta_t/2)/3600)then
	    st=i
	    break
	  else
	    st=0
	  end if
      end do
      
      ; ****************************************************
      ; end of time step and different types of mistakes that could be done
      ; ****************************************************

      if (end_time_step .EQ. -1.) then     	
	  delete(end_time_step)
	  end_time_step = t_all(nt-1)/3600  
      else
	  if (end_time_step .GT. t_all(nt-1)/3600)then
	    print(" ")
	    print("'end_time_step' = "+ end_time_step +"h is greater "+\
		  "than last time step = " +\
		    t_all(nt-1)+"s = "+t_all(nt-1)/3600+"h")
	    print(" ")
	    print("Select another 'end_time_step'")  
	    print(" ")
	    exit
	  end if
	  if (end_time_step .LT. start_time_step/3600)then
	    print(" ")
	    print("'end_time_step' = "+ end_time_step +"h is lower "+\
		  "than 'start_time_step' = "+start_time_step+"h")
	    print(" ")
	    print("Select another 'start_time_step' or 'end_time_step'")
	    print(" ")
	    exit
	  end if
      end if

      do i=0,nt-1     
	  if (end_time_step .GE. (t_all(i)-delta_t/2)/3600 .AND.\
	      end_time_step .LT. (t_all(i)+delta_t/2)/3600)then
	    et=i
	    break
	  else
	    et=0
	  end if
      end do 
      
      delete(start_time_step)
      start_time_step=round(st,3)
      delete(end_time_step)
      end_time_step=round(et,3)

      no_time = end_time_step-start_time_step+1
      
      ; ****************************************************
      ; time_stride and different types of mistakes that could be done
      ; ****************************************************
    
      if (time_stride .LT. 1) then
	  print(" ")
	  print("'time_stride' has to be positive and is set to 1")
	  print(" ")
	  time_stride = 1
      end if

      if (time_stride .GT. no_time) then
	  print(" ")
	  print("'time_stride' is greater than number of available "+\
	      "time steps,")
	  print("=> 'time_stride' is set to 1")
	  time_stride = 1
      end if

      ti_in = ispan(start_time_step,end_time_step,time_stride) ;ti_in contents 
								;the time indices 
								;to plot
      np    = dimsizes(ti_in)  

      if (com_i .EQ. max_com_i) then
	  print(" ")
	  print("Output of time steps from "+t_all(start_time_step)/3600+\
	    " h = "+t_all(start_time_step)+" s => index = "+start_time_step)
	  print("                     till "+t_all(ti_in(np-1))/3600+" h = "\
	    +t_all(ti_in(np-1))+" s => index = "+end_time_step)
	  print("                     with temporal stride = "+time_stride)
	  print(" ")
      end if

      if (com_i .EQ. 0) then

      ; ***************************************************
      ; set up minimum and maximum height
      ; *************************************************** 

      if (log_z .EQ. 1)then
	  if (min_z .EQ. -1)then
	    if (isvar("z_u"))then
		min_z=z_u(1)
	    else
		min_z=z_w(1)
	    end if       
	  else
	    if (isvar("z_w"))then
		if (min_z .GE. max(z_w) ) then
		  print(" ")
		  print("Minimum of height ('min_z'="+min_z+") is greater "+\
			"than available heights (="+max(z_w)+")")
		  print(" ")
		  exit
		end if	    
	    else
		if (min_z .GE. max(z_u) ) then
		  print(" ")
		  print("Minimum of height ('min_z'="+min_z+") is greater "+\
			"than available heights (="+max(z_u)+")")
		  print(" ")
		  exit
		end if 
	    end if
	    if (isvar("z_u"))then   
		if (min_z .LT. z_u(1) ) then
		  print(" ")
		  print("Begin height 'min_z' at least at level k=1 (="+\
			z_u(1)+"m) due to the logarithmic scale of the y-axis")
		  print(" ")
		  exit
		end if
	    else
		if (min_z .LT. z_w(1) ) then
		  print(" ")
		  print("Begin height 'min_z' at least at level k=1 (="+\
			z_w(1)+"m) due to the logarithmic scale of the y-axis")
		  print(" ")
		  exit
		end if   
	    end if
	  end if
      else
	  if (isvar("z_u"))then
	    if (min_z .EQ. -1)then
		min_z=z_u(0)
	    end if 
	  else
	    if (min_z .EQ. -1)then
		min_z=z_w(0)
	    end if    
	  end if
      
	  if (isvar("z_w"))then
	    if (min_z .GE. max(z_w) ) then
		print(" ")
		print("Minimum of height ('min_z'="+min_z+") is greater than "+\
		      "available heights (="+max(z_w)+")")
		print(" ")
		exit
	    end if
	  else
	    if (min_z .EQ. -1)then
		min_z=z_u(0)
	    end if
	    if (min_z .GE. max(z_u) ) then
		print(" ")
		print("Minimum of height ('min_z'="+min_z+") is greater than "+\
		      "available heights (="+max(z_u)+")")
		print(" ")
		exit
	    end if 
	  end if
      end if
      
      if (isvar("z_w"))then
	  if (max_z .EQ. -1)then
	    max_z=max(z_w)
	  end if 
      else
	  if (max_z .EQ. -1)then
	    max_z=max(z_u)
	  end if    
      end if
      
      if (isvar("z_w"))then
	  if (max_z .GT. max(z_w) ) then
	    print(" ")
	    print("Maximum of height ('max_z'="+max_z+") is greater than "+\
		  "available heights (="+max(z_w)+")")
	    print(" ")
	    exit
	  end if
      end if

      min_z=min_z/norm_z
      max_z=max_z/norm_z

      end if

      ; ****************************************************
      ; set up legend and colors 
      ; ****************************************************
      
      legend_label=new(np,string)
      do p=0, np-1
	  legend_label(p)=sprintf("%6.2f", t_all(ti_in(p))/3600)
      end do

      ; ***************************************************
      ; set up recourses
      ; ***************************************************

      res                         = True
      res@gsnDraw                 = False
      res@gsnFrame                = False
      res@txFont                  = "helvetica"
      res@tiMainFont              = "helvetica"
      res@tiXAxisFont             = "helvetica"
      res@tiYAxisFont             = "helvetica"
      res@tmXBLabelFont           = "helvetica"
      res@tmYLLabelFont           = "helvetica"
      res@lgLabelFont             = "helvetica"
      res@tmLabelAutoStride       = True


      if (legend .EQ. 1)then
	  res@pmLegendDisplayMode  = "Always"
      end if

      res@pmLegendSide            = "Left";"Right"
      res@xyExplicitLegendLabels  = legend_label
      res@pmLegendWidthF          = 0.12
      res@pmLegendHeightF         = 0.05*np
      res@pmLegendParallelPosF    = 1.0

      if (max_z .GE. 1000.) then
         res@pmLegendOrthogonalPosF  = -1.227
      else if ((max_z .LT. 1000.) .AND. (max_z .GE. 100.)) then
         res@pmLegendOrthogonalPosF  = -1.202
      else if ((max_z .LT. 100.) .AND. (max_z .GE. 10.)) then
         res@pmLegendOrthogonalPosF  = -1.177
      else if ((max_z .LT. 10.) .AND. (max_z .GE. 1.)) then
         res@pmLegendOrthogonalPosF  = -1.190
      else if ((max_z .LT. 1.)) then
         res@pmLegendOrthogonalPosF  = -1.215
      end if
      end if
      end if
      end if
      end if
      res@lgLabelFontHeightF      = font_size_legend
      res@lgTitleString           = "Time (h)"
      res@lgTitleFontHeightF      = font_size
      res@txFontHeightF           = font_size
      res@tiXAxisFontHeightF      = font_size
      res@tiYAxisFontHeightF      = font_size
      res@tmXBLabelFontHeightF    = font_size
      res@tmYLLabelFontHeightF    = font_size
      res@tiXAxisString           = " " 
      res@lgJustification         = "TopLeft"

      if ( black .eq. 0 ) then  
	  res@xyLineColors = -(ispan(-237,-2,235/np))
      end if
      if (norm_z .EQ. 1)then
	  res@tiYAxisString = "Height (m)"
      else   
	  res@tiYAxisString = "Height / "+norm_z+" (m)"
      end if
	
      if (log_z .EQ. 1) then
	  res@trYLog = True 
      end if

      if (dash .EQ. 0 ) then
	  res@xyMonoDashPattern = True 
      else
	  res@xyMonoDashPattern = False
	  if (no_files .GT. 1)
	    res@xyMonoDashPattern = True  
	    print(" ") 
	    print("If you use more than one file, patterns for different "+\
		  "timesteps cannot be used")
	    print(" ") 
	  end if       
      end if 

      res@tmXBMinorPerMajor = 4
      res@tmYLMinorPerMajor = 4

      resP                            = True
      resP@gsnMaximize                = True
      if ( no_files .EQ. 1 )  then
         resP@gsnPanelXWhiteSpacePercent = 6.0
         resP@gsnPanelYWhiteSpacePercent = 6.0
      else
         resP@gsnPanelXWhiteSpacePercent = 6.0
         resP@gsnPanelYWhiteSpacePercent = 0.0
      end if
      resP@txFont                     = "helvetica"
      resP@txString                   = f_att@title
      resP@txFuncCode                 = "~"
      resP@txFontHeightF              = 0.0105

      ; ***************************************************
      ; set up graphics for plot
      ; ***************************************************

      if (combine .EQ. 1) then
	  plot_o = new(number_comb,graphic)  
          label=new(number_comb,string)
	  color_o=new(number_comb,integer)

	  if (check_vType) then
	    mini=new(number_comb,double)
	    maxi=new(number_comb,double)
	  else
	    mini=new(number_comb,float)
	    maxi=new(number_comb,float)
	  end if
      end if
    
      if ( format_out .EQ. "pdf" .OR. format_out .EQ. "ps" ) then
	  format_out@wkPaperSize = "A4"
      end if
      if ( format_out .EQ. "png" ) then
	  format_out@wkWidth  = 1000
	  format_out@wkHeight = 1000
      end if

      if (wks_gsn_log) then
	  wks=gsn_open_wks(format_out,file_out)
	  gsn_define_colormap(wks,"rainbow+white") 
	  wks_gsn_log = False
      end if

      ; ***************************************************
      ; read data and create plots
      ; ***************************************************
	  
      do ti =0,np-1
	  if( t_all(ti_in(ti)) .lt. 10^36) then
	    start_time_step = ti
	    break
	  end if
      end do  
      
      if (log_z .EQ. 1) then      
	  if (check_vType) then
	    data   = new((/dim,np,dimz-1/),double) 
	    data_0 = new((/np,dimz-1/),double)
	  else
	    data   = new((/dim,np,dimz-1/),float) 
	    data_0 = new((/np,dimz-1/),float)
	  end if
	  data@_FillValue=9.96921e+36
	  data_0 = 0.1
	  t      = new((/np,dimz-1/),float) 
	  t      = 0.0
	  unit   = new(dim,string)
	  if (isvar("z_u"))then
	    if (typeof(z_u) .EQ. "double")then
		z_v    = new((/dim,dimz/),double)
		z_     = new((/dim,dimz-1/),double)
	    else
		if (typeof(z_u) .EQ. "float")then 
		  z_v    = new((/dim,dimz/),float)
		  z_     = new((/dim,dimz-1/),float)
		end if
	    end if
	  else
	    if (isvar("z_w"))then
		if (typeof(z_w) .EQ. "double")then
		  z_v    = new((/dim,dimz/),double)
		  z_     = new((/dim,dimz-1/),double)
		else
		  if (typeof(z_w) .EQ. "float")then 
		      z_v    = new((/dim,dimz/),float)
		      z_     = new((/dim,dimz-1/),float)
		  end if
		end if
	    end if
	  end if
      else
	  if (check_vType) then
	    data   = new((/dim,np,dimz/),double) 
	    data_0 = new((/np,dimz/),double)
	  else
	    data   = new((/dim,np,dimz/),float) 
	    data_0 = new((/np,dimz/),float)
	  end if
	  data@_FillValue=9.96921e+36      
	  data_0 = 0.0
	  t      = new((/np,dimz/),float) 
	  t      = 0.0
	  unit   = new(dim,string)
	  if (isvar("z_u"))then
	    if (typeof(z_u) .EQ. "double")then
		z_v    = new((/dim,dimz/),double)
		z_     = new((/dim,dimz/),double)
	    else
		if (typeof(z_u) .EQ. "float")then 
		  z_v    = new((/dim,dimz/),float)
		  z_     = new((/dim,dimz/),float)
		end if
	    end if
	  else
	    if (isvar("z_w"))then
		if (typeof(z_w) .EQ. "double")then
		  z_v    = new((/dim,dimz/),double)
		  z_     = new((/dim,dimz/),double)
		else
		  if (typeof(z_w) .EQ. "float")then 
		      z_v    = new((/dim,dimz/),float)
		      z_     = new((/dim,dimz/),float)
		  end if
		end if
	    end if
	  end if
      end if

      end if
      ;-------above steps only for first file

      ; ***************************************************
      ; indicate plot number 
      ; ***************************************************
      
      if (combine .EQ. 1) then
	  n = 1
      else
	  n = 0
      end if

      if (over .EQ. 1) then 
	  plot_u         = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_v         = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_w         = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_pt        = gsn_csm_xy(wks,t,data_0(:,:),res)      
	  plot_vpt       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_lpt       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_q         = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_qv        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_ql        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_rho       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_s         = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_sa        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_e         = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_es        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_km        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_kh        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_l         = gsn_csm_xy(wks,t,data_0(:,:),res)     
	  plot_wpup      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsus      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wu        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpvp      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsvs      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wv        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpptp     = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wspts     = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpt       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsptsBC   = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wptBC     = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpvptp    = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsvpts    = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wvpt      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpqp      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsqs      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wq        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpqvp     = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsqvs     = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wqv       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpsp      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsss      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_ws        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpsap     = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wssas     = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsa       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wses      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_us2       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_vs2       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_ws2       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_pts2      = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_ws3       = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_Sw        = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_ws2pts    = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wspts2    = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wsususodz = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wspsodz   = gsn_csm_xy(wks,t,data_0(:,:),res)
	  plot_wpeodz    = gsn_csm_xy(wks,t,data_0(:,:),res)
      
	  if (check_vType) then
	    miniu         =  100000.d
	    maxiu         = -100000.d
	    miniv         =  100000.d
	    maxiv         = -100000.d
	    miniw         =  100000.d
	    maxiw         = -100000.d
	    minipt        =  100000.d
	    maxipt        = -100000.d
	    minivpt       =  100000.d
	    maxivpt       = -100000.d
	    minilpt       =  100000.d
	    maxilpt       = -100000.d
	    miniq         =  100000.d
	    maxiq         = -100000.d
	    miniqv        =  100000.d
	    maxiqv        = -100000.d
	    miniql        =  100000.d
	    maxiql        = -100000.d
	    minie         =  100000.d
	    maxie         = -100000.d
	    minies        =  100000.d
	    maxies        = -100000.d
	    minikm        =  100000.d
	    maxikm        = -100000.d
	    minikh        =  100000.d
	    maxikh        = -100000.d
	    miniwpup      =  100000.d
	    maxiwpup      = -100000.d
	    miniwsus      =  100000.d
	    maxiwsus      = -100000.d
	    miniwu        =  100000.d
	    maxiwu        = -100000.d
	    miniwpvp      =  100000.d
	    maxiwpvp      = -100000.d
	    miniwsvs      =  100000.d
	    maxiwsvs      = -100000.d
	    miniwv        =  100000.d
	    maxiwv        = -100000.d
	    miniwpptp     =  100000.d
	    maxiwpptp     = -100000.d
	    miniwspts     =  100000.d
	    maxiwspts     = -100000.d
	    miniwpt       =  100000.d
	    maxiwpt       = -100000.d
	    miniwsptsBC   =  100000.d
	    maxiwsptsBC   = -100000.d
	    miniwptBC     =  100000.d
	    maxiwptBC     = -100000.d
	    miniwpvptp    =  100000.d
	    maxiwpvptp    = -100000.d
	    miniwsvpts    =  100000.d
	    maxiwsvpts    = -100000.d
	    miniwvpt      =  100000.d
	    maxiwvpt      = -100000.d
	    miniwpqp      =  100000.d
	    maxiwpqp      = -100000.d
	    miniwsqs      =  100000.d
	    maxiwsqs      = -100000.d
	    miniwq        =  100000.d
	    maxiwq        = -100000.d
	    miniwpqvp     =  100000.d
	    maxiwpqvp     = -100000.d
	    miniwsqvs     =  100000.d
	    maxiwsqvs     = -100000.d
	    miniwqv       =  100000.d
	    maxiwqv       = -100000.d
	    miniwpsp      =  100000.d
	    maxiwpsp      = -100000.d
	    miniwsss      =  100000.d
	    maxiwsss      = -100000.d
	    miniws        =  100000.d
	    maxiws        = -100000.d
	    miniwpsap     =  100000.d
	    maxiwpsap     = -100000.d
	    miniwssas     =  100000.d
	    maxiwssas     = -100000.d
	    miniwsa       =  100000.d
	    maxiwsa       = -100000.d
	    minius2       =  100000.d
	    maxius2       = -100000.d
	    minivs2       =  100000.d
	    maxivs2       = -100000.d
	    miniws2       =  100000.d
	    maxiws2       = -100000.d
	    miniwsususodz =  100000.d
	    maxiwsususodz = -100000.d
	    miniwspsodz   =  100000.d
	    maxiwspsodz   = -100000.d
	    miniwpeodz    =  100000.d
	    maxiwpeodz    = -100000.d
	  else
	    miniu         =  100000.
	    maxiu         = -100000.
	    miniv         =  100000.
	    maxiv         = -100000.
	    miniw         =  100000.
	    maxiw         = -100000.
	    minipt        =  100000.
	    maxipt        = -100000.
	    minivpt       =  100000.
	    maxivpt       = -100000.
	    minilpt       =  100000.
	    maxilpt       = -100000.
	    miniq         =  100000.
	    maxiq         = -100000.
	    miniqv        =  100000.
	    maxiqv        = -100000.
	    miniql        =  100000.
	    maxiql        = -100000.
	    minie         =  100000.
	    maxie         = -100000.
	    minies        =  100000.
	    maxies        = -100000.
	    minikm        =  100000.
	    maxikm        = -100000.
	    minikh        =  100000.
	    maxikh        = -100000.
	    miniwpup      =  100000.
	    maxiwpup      = -100000.
	    miniwsus      =  100000.
	    maxiwsus      = -100000.
	    miniwu        =  100000.
	    maxiwu        = -100000.
	    miniwpvp      =  100000.
	    maxiwpvp      = -100000.
	    miniwsvs      =  100000.
	    maxiwsvs      = -100000.
	    miniwv        =  100000.
	    maxiwv        = -100000.
	    miniwpptp     =  100000.
	    maxiwpptp     = -100000.
	    miniwspts     =  100000.
	    maxiwspts     = -100000.
	    miniwpt       =  100000.
	    maxiwpt       = -100000.
	    miniwsptsBC   =  100000.
	    maxiwsptsBC   = -100000.
	    miniwptBC     =  100000.
	    maxiwptBC     = -100000.
	    miniwpvptp    =  100000.
	    maxiwpvptp    = -100000.
	    miniwsvpts    =  100000.
	    maxiwsvpts    = -100000.
	    miniwvpt      =  100000.
	    maxiwvpt      = -100000.
	    miniwpqp      =  100000.
	    maxiwpqp      = -100000.
	    miniwsqs      =  100000.
	    maxiwsqs      = -100000.
	    miniwq        =  100000.
	    maxiwq        = -100000.
	    miniwpqvp     =  100000.
	    maxiwpqvp     = -100000.
	    miniwsqvs     =  100000.
	    maxiwsqvs     = -100000.
	    miniwqv       =  100000.
	    maxiwqv       = -100000.
	    miniwpsp      =  100000.
	    maxiwpsp      = -100000.
	    miniwsss      =  100000.
	    maxiwsss      = -100000.
	    miniws        =  100000.
	    maxiws        = -100000.
	    miniwpsap     =  100000.
	    maxiwpsap     = -100000.
	    miniwssas     =  100000.
	    maxiwssas     = -100000.
	    miniwsa       =  100000.
	    maxiwsa       = -100000.
	    minius2       =  100000.
	    maxius2       = -100000.
	    minivs2       =  100000.
	    maxivs2       = -100000.
	    miniws2       =  100000.
	    maxiws2       = -100000.
	    miniwsususodz =  100000.
	    maxiwsususodz = -100000.
	    miniwspsodz   =  100000.
	    maxiwspsodz   = -100000.
	    miniwpeodz    =  100000.
	    maxiwpeodz    = -100000.
	  end if

      end if

      if (prof3d .EQ. 1)then

      if (end_x .EQ. -1) then 
	  end_x=dimx-2
      end if
      if (end_y .EQ. -1)then
	  end_y=dimy-2
      end if
      if (start_x .LT. 0)then
	  print(" ")
	  print("'start_x' is lower than 0 and set to 0")
	  print(" ")
	  start_x=0
      end if
      if (start_x .GT. dimx-1)then
	  print(" ")
	  print("'start_x' is greater than available x-range and set to "+\
		"maximum of x-range (excluding ghostpoint)")
	  print(" ")
	  start_x=dimx-2
      end if
      if (end_x .EQ. dimx-1)then
	  print(" ")
	  print("'end_x' = "+end_x+" and includes the ghostpoint")
	  print(" ")
      end if
      if (end_x .GT. dimx-1)then
	  print(" ")
	  print("'end_x' = "+end_x+" is greater than available x-range and set "+\
		"to maximum of x-range (excluding ghostpoint)")
	  print(" ")
	  end_x=dimx-2
      end if
      if (end_x .LT. start_x)then
	  print(" ")
	  print("'end_x' = "+end_x+" is lower than 'start_x' = "+start_x+\
		" and set to maximum of x-range (excluding ghostpoint)")
	  print(" ")
	  end_x=dimx-2
      end if
      if (start_y .LT. 0)then
	  print(" ")
	  print("'start_y' is lower than 0 and set to 0")
	  print(" ")
	  start_y=0
      end if
      if (start_y .GT. dimy-1)then
	  print(" ")
	  print("'start_y' is greater than available y-range and set to "+\
		"maximum of y-range (excluding ghostpoint)")
	  print(" ")
	  start_x=dimy-2
      end if
      if (end_y .EQ. dimy-1)then
	  print(" ")
	  print("'end_y' = "+end_y+" and includes the ghostpoint")
	  print(" ")
      end if
      if (end_y .GT. dimy-1)then
	  print(" ")
	  print("'end_y' = "+end_y+" is greater than available y-range and "+\
		"set to maximum of y-range (excluding ghostpoint)")
	  print(" ")
	  end_x=dimy-2
      end if
      if (end_y .LT. start_y)then
	  print(" ")
	  print("'end_y' = "+end_y+" is lower than 'start_y' = "+start_y+\
		" and set to maximum of y-range (excluding ghostpoint)")
	  print(" ")
	  end_y=dimy-2
      end if
      
      end if
      
      n_o=0
      count_var=0

      res@xyDashPattern = 1*nof

      over_remind = False
      if ( over .eq. 1)then
	over_remind = True
      end if 

      do varn = 0,dim-1 

	  check = True
	  
	  if (prof3d .EQ. 0) then
	    if ( isStrSubset( vNam(varn), "time") .OR. \
		  isStrSubset( vNam(varn), "NORM")) then
		check = False
	    end if
	  else
	    if ( isStrSubset( vNam(varn), "time") .OR.  \
		  isStrSubset( vNam(varn), "zusi") .OR.  \
		  isStrSubset( vNam(varn), "zwwi") .OR.  \
		  isStrSubset( vNam(varn), "x") .OR.     \
		  isStrSubset( vNam(varn), "xu") .OR.    \
		  isStrSubset( vNam(varn), "y") .OR.     \
		  isStrSubset( vNam(varn), "yv") .OR.    \ 
		  isStrSubset( vNam(varn), "zu_3d") .OR. \
		  isStrSubset( vNam(varn), "zw_3d")) then
		check = False
	    end if
	  end if

	  if (var .NE. "all") then
	    check = isStrSubset( var,","+vNam(varn)+"," )
	  end if

	  if (combine .EQ. 1) then         
	    com=isStrSubset(c_var,","+vNam(varn)+"," )      
	    if (com) then      
		if (prof3d .EQ. 0) then	      	     
		  temp     = f[:]->$vNam(varn)$	      
		  temp_att = f_att->$vNam(varn)$ 
		  if (log_z .EQ. 1) then
		      do j=0,np-1
			data(varn,j,:) = temp(ti_in(j),1:dimz-1)
		      end do
		  else
		      do j=0,np-1
			data(varn,j,:) = temp(ti_in(j),0:dimz-1)
		      end do
		  end if
		else
		  if (log_z .EQ. 1) then
		      do i=1,dimz-1
			do j=0,np-1
			    temp = f[:]->$vNam(varn)$
			    temp_att = f_att->$vNam(varn)$ 
			    data_temp = temp(ti_in(j),i,\
					    start_y:end_y,start_x:end_x)
			    data(varn,j,i-1) = dim_avg_Wrap(\
						      dim_avg_Wrap(data_temp))
			end do
		      end do
		  else
		      do i=0,dimz-1
			do j=0,np-1
			    temp = f[:]->$vNam(varn)$
			    temp_att = f_att->$vNam(varn)$ 
			    data_temp = temp(ti_in(j),i,\
						start_y:end_y,start_x:end_x)
			    data(varn,j,i) = dim_avg_Wrap(dim_avg_Wrap(data_temp))
			end do
		      end do
		  end if
		  print(" ")
		  print("Variable for combine '"+vNam(varn)+"' is read")
		  print(" ")
		end if                 
		unit(varn) = temp_att@units
		if (n_o .GT. number_comb-1) then
		  print(" ")
		  print("Set 'number_comb' to the number of overlaying "+\
			"variables ('c_var' = "+c_var+")")
		  print(" ")
		  exit
		end if 
		if (abs(min(data(varn,:,:))) .GT. 10)then
		  min_value = abs(0.01*min(data(varn,:,:)))
		  max_value = abs(0.01*max(data(varn,:,:)))
		else
		  if (abs(min(data(varn,:,:))) .LT. 0.01 .AND. \
		      abs(max(data(varn,:,:))) .GT. 0.01)then
		      min_value = abs(0.1*max(data(varn,:,:)))
		      max_value = abs(0.1*max(data(varn,:,:)))
		  else
		      if (abs(max(data(varn,:,:))) .LT. 0.01 .AND. \
			  abs(min(data(varn,:,:))) .GT. 0.01)then
			min_value = abs(0.1*min(data(varn,:,:)))
			max_value = abs(0.1*min(data(varn,:,:)))
		      else
			min_value = abs(0.1*min(data(varn,:,:)))
			max_value = abs(0.1*max(data(varn,:,:)))
		      end if
		  end if
		end if 
		if (min(data(varn,:,:)) .EQ. 0 .AND. \
		    max(data(varn,:,:)) .EQ. 0)then
		  min_value = 0.1
		  max_value = 0.1
		end if
		mini(n_o)=min(data(varn,:,:))
		maxi(n_o)=max(data(varn,:,:)) 
		n_o=n_o+1
	    end if
	  end if

	  if(check) then
	    
	    count_var=count_var+1

	    if (prof3d .EQ. 0) then	  
		temp = f[:]->$vNam(varn)$
		temp_att = f_att->$vNam(varn)$ 
	    else
		if (log_z .EQ. 1) then
		  do i=1,dimz-1
		      do j=0,np-1
			temp= f[:]->$vNam(varn)$
			temp_att = f_att->$vNam(varn)$ 
			data_temp = temp(ti_in(j),i,start_y:end_y,start_x:end_x)
			data(varn,j,i-1) = dim_avg_Wrap(dim_avg_Wrap(data_temp))
			delete(data_temp)
		      end do
		  end do
		else
		  do i=0,dimz-1
		      do j=0,np-1
			temp= f[:]->$vNam(varn)$
			temp_att = f_att->$vNam(varn)$ 
			data_temp = temp(ti_in(j),i,start_y:end_y,start_x:end_x)
			data_temp!0 = "t"
			data_temp!1 = "z"
			data(varn,j,i) = dim_avg_Wrap(dim_avg_Wrap(data_temp))
			delete(data_temp)
		      end do
		  end do
		end if
		print(" ")
		print("Variable '"+vNam(varn)+"' is read")
		print(" ")
		unit(varn) = temp_att@units
		a=getvardims(temp_att)
		b=dimsizes(a)
	    end if

	    if (prof3d .EQ. 0) then
		if (log_z .EQ. 1) then
		  z = f_att->$vNam(varn+1)$(1:dimz-1) 
		  unit(varn) = temp_att@units
		  do j=0,np-1
		      data(varn,j,:) = temp(ti_in(j),1:dimz-1) 
		  end do
		else
		  z = f_att->$vNam(varn+1)$ 
		  unit(varn) = temp_att@units
		  do j=0,np-1
		      data(varn,j,:) = temp(ti_in(j),:) 
		  end do
		end if
	    else 
		do i=0,b-1           
		  if (isStrSubset( a(i),"zu_3d" ))then
		      z_v(varn,:) = z_u
		      if (log_z .EQ. 1) then
			z = z_v(varn,1:dimz-1)
		      else
			z = z_v(varn,:)
		      end if 
		  else
		      if (isStrSubset( a(i),"zw_3d" ))then
			z_v(varn,:) = z_w
			if (log_z .EQ. 1) then
			    z = z_v(varn,1:dimz-1)
			else
			    z = z_v(varn,:)
			end if
		      end if                   
		  end if
		end do            
	    end if

	    if (isStrSubset(data@long_name," SR " ) ) then
	      over = 0
	    end if 
	    
	    ;if (nof .EQ. 0) then
		;z_(n,:)=z/norm_z
		;z    = z_(n,:)
	    ;else
		z=z/norm_z
	    ;end if

	    delta_z = z(2) - z(1)

	    max_z_int=doubletoint(max_z/delta_z)
	    min_z_int=doubletoint(min_z/delta_z)

	    if(max_z_int .eq. min_z_int)
		print(" ")
		print("Please increase 'max_z' or decrease 'min_z' so that "+\
		      "there are")
		print("at least two layers for the z-axis to plot")
		print(" ")
		exit
	    end if

	    if(min_z_int .gt. max_z_int)
		print(" ")
		print("'min_z' is greater than 'max_z',")
		print("please change this")
		print(" ")
		exit
	    end if

	    if (max_z_int .ge. dimz-1)
		max_z_int = dimz-1
		if (log_z .EQ. 1) then
		max_z_int = max_z_int - 1
		end if
	    end if

	    if (min_z_int .lt. 0)
		min_z_int = 0
		if (log_z .EQ. 1) then
		  min_z_int = min_z_int + 1
		end if
	    end if         

	    ;data can contain missing values in case of output of t=0h (inital profiles)
	    ;where no output is possible
	    n_not_ismissing = num(.not.ismissing(data(varn,:,:)))

	    if (n_not_ismissing .GT. 0 ) then    
    
	      if (abs(min(data(varn,:,min_z_int:max_z_int))) .GT. 10)then
		  min_value = abs(0.001*min(data(varn,:,min_z_int:max_z_int)))
		  max_value = abs(0.001*max(data(varn,:,min_z_int:max_z_int)))
	      else
		  if (abs(min(data(varn,:,min_z_int:max_z_int))) .LT. 0.01 .AND. \
		      abs(max(data(varn,:,min_z_int:max_z_int))) .GT. 0.01)then
		    min_value = abs(0.1*max(data(varn,:,min_z_int:max_z_int)))
		    max_value = abs(0.1*max(data(varn,:,min_z_int:max_z_int)))
		  else
		    if (abs(max(data(varn,:,min_z_int:max_z_int))) .LT. 0.01 .AND.\
			abs(min(data(varn,:,min_z_int:max_z_int))) .GT. 0.01)then
			min_value = abs(0.1*min(data(varn,:,min_z_int:max_z_int)))
			max_value = abs(0.1*min(data(varn,:,min_z_int:max_z_int)))
		    else
			min_value = abs(0.1*min(data(varn,:,min_z_int:max_z_int)))
			max_value = abs(0.1*max(data(varn,:,min_z_int:max_z_int)))
		    end if
		  end if
	      end if 
	    
	      if (min(data(varn,:,min_z_int:max_z_int)) .EQ. 0 .AND. \
		  max(data(varn,:,min_z_int:max_z_int)) .EQ. 0)then
		  min_value = 0.1
		  max_value = 0.1
	      end if

	    else     
	      print(" ")
	      print(vNam(varn) + " contains only missing values")
	      print(" ")
	    end if   

	    if (over .EQ. 0) then  
		res@gsnLeftString      = " "
		res@tiXAxisString      = "("+unit(varn)+")"
		res@gsnRightString     = " "
		res@trYMinF            = min_z
		res@trYMaxF            = max_z 
		if (xs .EQ. -1) then 
		  res@trXMinF            = min(data(varn,:,min_z_int:max_z_int))-\
									  min_value
		else 
		  res@trXMinF            = xs     
		end if
		if (xe .EQ. -1) then 
		  res@trXMaxF            = max(data(varn,:,min_z_int:max_z_int))+\
									  max_value
		else 
		  res@trXMaxF            = xe  
		end if         

		if (c_var_log .EQ. 1) then
		  if ((mod(com_i, no_rows * no_columns) .EQ. 0) .AND. (no_files .LT. 2))then
		    res@pmLegendDisplayMode  = "Always"
		  else
		    res@pmLegendDisplayMode  = "Never"
		  end if
		else
		  if ((mod(no_plots-1-n+combine, no_rows * no_columns) .EQ. 0) .AND. (no_files .LT. 2))then
		    res@pmLegendDisplayMode  = "Always"
		  else
		    res@pmLegendDisplayMode  = "Never"
		  end if
		end if

		plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)  

	    end if
      
	
	    if (vNam(varn) .EQ. "u" ) then
		  miniu=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiu=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_u = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "                 
		      if (xs .EQ. -1) then  
			res@trXMinF            = miniu
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiu 
		      else 
			res@trXMaxF            = xe  
		      end if                
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)  
		  end if
		end if
		if (vNam(varn) .EQ. "v") then 
		  miniv=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiv=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyMonoDashPattern = True
		      res@xyDashPattern  = 1
		      plot_v = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "                  
		      if (xs .EQ. -1) then
			res@trXMinF            = miniv
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF            = maxiv
		      else 
			res@trXMaxF            = xe  
		      end if                
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if  
		end if
		if (vNam(varn) .EQ. "w") then
		  miniw=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiw=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern = 2
		      plot_w = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF            = miniw
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF            = maxiw
		      else 
			res@trXMaxF            = xe  
		      end if           
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if   
		end if

		if (vNam(varn) .EQ. "pt") then
		  minipt=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxipt=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_pt = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = minipt
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then       
			res@trXMaxF            = maxipt 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if  
		end if
		if (vNam(varn) .EQ. "vpt") then
		  minivpt=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxivpt=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_vpt = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = minivpt 
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then       
			res@trXMaxF            = maxivpt 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if 
		end if
		if (vNam(varn) .EQ. "lpt") then
		  minilpt=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxilpt=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2
		      plot_lpt = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = minilpt
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxilpt
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if

		if (vNam(varn) .EQ. "q") then
		  miniq=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiq=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_q = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniq
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiq
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if
		if (vNam(varn) .EQ. "qv") then
		  miniqv=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiqv=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then 
		      res@xyDashPattern  = 1
		      plot_qv = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniqv
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then         
			res@trXMaxF            = maxiqv 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "ql") then
		  miniql=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiql=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2
		      plot_ql = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniql
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then       
			res@trXMaxF            = maxiql 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "e") then
		  minie=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxie=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_e = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = minie
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then       
			res@trXMaxF            = maxie 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if
		if (vNam(varn) .EQ. "es" .OR. vNam(varn) .EQ. "e*") then
		  minies=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxies=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_es = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF            = minies
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxies
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if

		if (vNam(varn) .EQ. "km") then
		  minikm=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxikm=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_km = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF            = minikm
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then      
			res@trXMaxF            = maxikm
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if
		if (vNam(varn) .EQ. "kh") then
		  minikh=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxikh=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_kh = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = minikh
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then     
			res@trXMaxF            = maxikh
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if

		if (vNam(varn) .EQ. "wpup" .OR. vNam(varn) .EQ. "w"+dq+"u"+dq) then
		  miniwpup=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpup=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wpup = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpup
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF            = maxiwpup 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if
		if (vNam(varn) .EQ. "wsus" .OR. vNam(varn) .EQ. "w*u*") then
		  miniwsus=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsus=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wsus = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsus
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then      
			res@trXMaxF            = maxiwsus 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "wu") then
		  miniwu=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwu=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2 
		      plot_wu = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF            = miniwu
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then 
			res@trXMaxF            = maxiwu 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "wpvp" .OR. vNam(varn) .EQ. "w"+dq+"v"+dq) then
		  miniwpvp=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpvp=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wpvp = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpvp
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then     
			res@trXMaxF            = maxiwpvp 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if
		if (vNam(varn) .EQ. "wsvs" .OR. vNam(varn) .EQ. "w*v*") then
		  miniwsvs=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsvs=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wsvs = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF            = miniwsvs
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then      
			res@trXMaxF            = maxiwsvs 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "wv") then
		  miniwv=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwv=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2 
		      plot_wv = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwv
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwv 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "wpptp" .OR. vNam(varn) \
		  .EQ. "w"+dq+"pt"+dq) then
		  miniwpptp=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpptp=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wpptp = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpptp
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpptp 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if
		if (vNam(varn) .EQ. "wspts" .OR. vNam(varn) .EQ. "w*pt*") then
		  miniwspts=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwspts=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wspts = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwspts
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwspts 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "wpt") then       
		  miniwpt=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpt=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2 
		      plot_wpt = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpt
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpt 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "wsptsBC".OR. vNam(varn) .EQ. "w*pt*BC" ) then
		  miniwsptsBC=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsptsBC=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wsptsBC = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsptsBC
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwsptsBC 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if             
		if (vNam(varn) .EQ. "wptBC") then
		  miniwptBC=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwptBC=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wptBC = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwptBC
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwptBC 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if

		if (vNam(varn) .EQ. "wpvptp" .OR. vNam(varn) \
		    .EQ. "w"+dq+"vpt"+dq) then
		  miniwpvptp=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpvptp=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wpvptp = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpvptp
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpvptp 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if
		if (vNam(varn) .EQ. "wsvpts" .OR. vNam(varn) .EQ. "w*vpt*") then
		  miniwsvpts=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsvpts=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wsvpts = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsvpts
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwsvpts 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "wvpt") then
		  miniwvpt=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwvpt=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2 
		      plot_wvpt = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwvpt
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwvpt 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "wpqp" .OR. vNam(varn) .EQ. "w"+dq+"q"+dq) then
		  miniwpqp=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpqp=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wpqp = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpqp
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpqp 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if
		if (vNam(varn) .EQ. "wsqs".OR. vNam(varn) .EQ. "w*s*" ) then
		  miniwsqs=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsqs=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wsqs = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsqs
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwsqs 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "wq") then
		  miniwq=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwq=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2
		      plot_wq = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwq
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwq 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "wpqvp" .OR. \
		  vNam(varn) .EQ. "w"+dq+"qv"+dq) then
		  miniwpqvp=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpqvp=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0 
		      plot_wpqvp = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpqvp
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpqvp 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if
		if (vNam(varn) .EQ. "wsqvs" .OR. vNam(varn) .EQ. "w*qv*") then
		  miniwsqvs=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsqvs=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wsqvs = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsqvs
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwsqvs 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "wqv") then
		  miniwqv=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwqv=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2
		      plot_wqv = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwqv
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwqv 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "wpsp" .OR. vNam(varn) .EQ. "w"+dq+"s"+dq) then
		  miniwpsp=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpsp=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wpsp = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpsp
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpsp 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if 
		if (vNam(varn) .EQ. "wsss" .OR. vNam(varn) .EQ. "w*s*" ) then
		  miniwsss=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsss=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1 
		      plot_wsss = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsss
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwsss 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "ws") then
		  miniws=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiws=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2
		      plot_ws = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniws
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiws 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if

		if (vNam(varn) .EQ. "wpsap" .OR. \
		    vNam(varn) .EQ. "w"+dq+"sa"+dq) then
		  miniwpsap=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpsap=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wpsap = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpsap
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpsap 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if
		if (vNam(varn) .EQ. "wssas" .OR. vNam(varn) .EQ. "w*sa*") then
		  miniwssas=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwssas=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wssas = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwssas
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwssas 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "wsa") then
		  miniwsa=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsa=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2
		      plot_wsa = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsa
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwsa 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "us2" .OR. vNam(varn) .EQ. "u*2") then
		  minius2=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxius2=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_us2 = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = minius2
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxius2 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if 
		end if
		if (vNam(varn) .EQ. "vs2" .OR. vNam(varn) .EQ. "v*2") then
		  minivs2=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxivs2=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_vs2 = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = minivs2
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxivs2 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if
		end if
		if (vNam(varn) .EQ. "ws2" .OR. vNam(varn) .EQ. "w*2") then
		  miniws2=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiws2=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2 
		      plot_ws2 = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniws2
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiws2 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if

		if (vNam(varn) .EQ. "wsususodz" .OR. \
		    vNam(varn) .EQ. "w*u*u*:dz") then
		  miniwsususodz=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwsususodz=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 0
		      plot_wsususodz = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwsususodz
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwsususodz 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if  
		end if
		if (vNam(varn) .EQ. "wspsodz" .OR. vNam(varn) .EQ. "w*p*:dz") then
		  miniwspsodz=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwspsodz=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 1
		      plot_wspsodz = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwspsodz
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwspsodz 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  end if 
		end if
		if (vNam(varn) .EQ. "wpeodz" .OR. \
		    vNam(varn) .EQ. "w"+dq+"p:dz") then
		  miniwpeodz=min(data(varn,:,min_z_int:max_z_int))-min_value
		  maxiwpeodz=max(data(varn,:,min_z_int:max_z_int))+max_value
		  if (over .EQ. 1) then
		      res@xyDashPattern  = 2 
		      plot_wpeodz = gsn_csm_xy(wks,data(varn,:,:),z,res)
		  else
		      res@gsnLeftString      = " "
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then 
			res@trXMinF            = miniwpeodz
		      else 
			res@trXMinF            = xs
		      end if
		      if (xe .EQ. -1) then        
			res@trXMaxF            = maxiwpeodz 
		      else 
			res@trXMaxF            = xe  
		      end if
		      plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  end if
		end if 

		; ***************************************************
		; legend for each plot
		; ***************************************************
		
		if (over .EQ. 0) then

		   label = " " + vNam(varn)

		   lgres1                    = True
		   lgMonoDashIndex           = False
		   lgres1@lgDashIndexes      = (/0,1,2/)
		   lgres1@lgLabelFont        = "helvetica"    
		   lgres1@lgLabelFontHeightF = font_size_legend * 10.0 * 0.7
		   lgres1@vpWidthF           = 0.14           
		   lgres1@vpHeightF          = 0.1 / 3.0 + 0.01         
		
		   lbid = gsn_create_legend(wks,1,label,lgres1)       

		   amres = True
		   amres@amParallelPosF   = -0.5		    
		   amres@amOrthogonalPosF = -0.5   
		   amres@amJust           = "TopLeft"

		   annoid1 = gsn_add_annotation(plot(n),lbid,amres) 
      
		end if

	    if (no_files .GT. 1) then
;		print("nof="+nof+" und n="+n)
		multi_plot(nof,n)=plot(n) 
		max_nof(nof,n)=max(data(varn,:,min_z_int:max_z_int))
		min_nof(nof,n)=min(data(varn,:,min_z_int:max_z_int))
		name(nof,n)   =vNam(varn)
		unit_(nof,n)  =unit(varn)
	    end if
	    if (over .EQ. 0) then
		n=n+1  
	    end if  
	    if (prof3d .EQ. 0)then    
		varn=varn+1
	    end if
	    delete(temp) 
	    delete(temp_att)
	  end if          
      end do
      if (no_files .GT. 1) then
	  delete(vNam)
	  delete(files)
      end if
      
      end do

      ;#########ENDE DO LOOP FOR no_files GT 1#############

      if (isStrSubset(data@long_name," SR " ) .and. over_remind) then
	  print(" ")
	  print("If you have outputs of statistic regions you cannot overlay "+\
		"variables;")
	  print("'over' is set to 0" )
	  print(" ")
	  over = 0
      end if 

      if (count_var .EQ. 0) then
	  print(" ")
	  print("The variables 'var="+var+"'" )
	  print("do not exist on your input file;")
	  print("be sure to have one comma before and after each variable")
	  print(" ")
	  exit        
      end if
      
      if (no_files .GT. 1) then 
	  multi_legend=new(6,string)
	  string_len=new(6,integer)
	  multi_dash=new(no_files,string)
	  multi_legend(0)="  "+name_legend_1
	  string_len(0)=strlen(multi_legend(0))
	  multi_legend(1)="  "+name_legend_2 
	  string_len(1)=strlen(multi_legend(1))
	  multi_legend(2)="  "+name_legend_3
	  string_len(2)=strlen(multi_legend(2))
	  multi_legend(3)="  "+name_legend_4 
	  string_len(3)=strlen(multi_legend(3))
	  multi_legend(4)="  "+name_legend_5 
	  string_len(4)=strlen(multi_legend(4))
	  multi_legend(5)="  "+name_legend_6 
	  string_len(5)=strlen(multi_legend(5))
	  do ml=1,no_files
	    multi_dash(ml-1)=ml-1
	  end do
	  delete(plot)
	  plot = new(dim,graphic)
	  do pl=0,n-1
	    plot0 = new(1,graphic)
	    res@trXMinF = min(min_nof(:,pl))
	    res@trXMaxF = max(max_nof(:,pl))
            res@xyDashPattern=0
	    res@gsnLeftString  = " "; name(0,pl)
	    res@gsnRightString = " ";unit_(0,pl) 
	    res@tiXAxisString  = "(" + unit_(0,pl) + ")"
	    data_0(:,:) = min(min_nof(:,pl))

	    if ( (mod(no_plots - pl, no_rows * no_columns) .EQ. 1) .OR. (no_rows * no_columns .EQ. 1) ) then
		res@pmLegendDisplayMode  = "Always"
	    else
		res@pmLegendDisplayMode  = "Never"
	    end if

	    plot0 = gsn_csm_xy(wks,data_0(:,:),z_(pl,:),res)

	    ; ***************************************************
	    ; legend for combined plot (File 1, File 2, ...)
	    ; ***************************************************

	    if ( (mod(no_plots - pl, no_rows * no_columns) .EQ. 1) .OR. (no_rows * no_columns .EQ. 1) )then
		lgres                    = True
		lgMonoDashIndex          = False
		lgres@lgLabelFont        = "helvetica"    
		lgres@lgLabelFontHeightF = font_size_legend * 10.0       
		lgres@vpWidthF           = max(string_len)*0.02;0.015           
		lgres@vpHeightF          = 0.035*no_files    
		lgres@lgDashIndexes      = multi_dash(no_files-1:0)
                lgres@lgLineDashSegLenF  = 0.075
		lbid = gsn_create_legend(\
				    wks,no_files,multi_legend(no_files-1:0),lgres)

		amres = True
		amres@amParallelPosF   = -0.5	    
		amres@amOrthogonalPosF = -0.55
		amres@amJust           = "BottomLeft"
		annoid1 = gsn_add_annotation(plot0,lbid,amres) 
	    end if

	    do plo=0,no_files-1
		overlay(plot0,multi_plot(plo,pl))
		plot(pl)=plot0
	    end do
	    delete(plot0)
	  end do
      end if

      if (count_var .EQ. 0) then
	  print(" ")
	  print("Select a variable 'var=' or use the default value")
	  print(" ")
	  print("Your selection '"+var+"' does not exist on the input file")
	  print(" ")
	  exit
      end if

      if (over .EQ. 1 ) then

	  overlay(plot_u,plot_v)
	  overlay(plot_u,plot_w)
	  u=0
	  overlay(plot_pt,plot_vpt)
	  overlay(plot_pt,plot_lpt)
	  pt=0
	  overlay(plot_q,plot_qv)
	  overlay(plot_q,plot_ql)
	  q=0
	  overlay(plot_e,plot_es)
	  e=0
	  overlay(plot_km,plot_kh)
	  km=0
	  overlay(plot_wpup,plot_wsus)
	  overlay(plot_wpup,plot_wu)
	  wpup=0
	  overlay(plot_wpvp,plot_wsvs)
	  overlay(plot_wpvp,plot_wv)
	  wpvp=0
	  overlay(plot_wpptp,plot_wspts)
	  overlay(plot_wpptp,plot_wpt)
	  wpptp=0
	  overlay(plot_wsptsBC,plot_wptBC)
	  wsptsBC=0
	  overlay(plot_wpvptp,plot_wsvpts)
	  overlay(plot_wpvptp,plot_wvpt)
	  wpvptp=0
	  overlay(plot_wpqp,plot_wsqs)
	  overlay(plot_wpqp,plot_wq)
	  wpqp=0
	  overlay(plot_wpqvp,plot_wsqvs)
	  overlay(plot_wpqvp,plot_wqv)
	  wpqvp=0
	  overlay(plot_wpsp,plot_wsss)
	  overlay(plot_wpsp,plot_ws)
	  wpsp=0
	  overlay(plot_wpsap,plot_wssas)
	  overlay(plot_wpsap,plot_wsa)
	  wpsap=0
	  overlay(plot_us2,plot_vs2)
	  overlay(plot_us2,plot_ws2)
	  us2=0
	  overlay(plot_wsususodz,plot_wspsodz)
	  overlay(plot_wsususodz,plot_wpeodz)
	  wsususodz=0

      end if

      if (over .EQ. 1) then 
      
	  do varn = 0,dim-1    
	    
	    check = True
	  
	    if (prof3d .EQ. 0) then
		if ( isStrSubset( vNam(varn), "time") .OR. \
		    isStrSubset( vNam(varn), "NORM")) then
		  check = False
		end if
	    else
		if ( isStrSubset( vNam(varn), "time") .OR.  \
		    isStrSubset( vNam(varn), "zusi") .OR.  \
		    isStrSubset( vNam(varn), "zwwi") .OR.  \
		    isStrSubset( vNam(varn), "x") .OR.     \
		    isStrSubset( vNam(varn), "xu") .OR.    \
		    isStrSubset( vNam(varn), "y") .OR.     \
		    isStrSubset( vNam(varn), "yv") .OR.    \
		    isStrSubset( vNam(varn), "zu_3d") .OR. \
		    isStrSubset( vNam(varn), "zw_3d")) then
		  check = False
		end if
	    end if

	    if (var .NE. "all") then      
		check = isStrSubset( var,","+vNam(varn)+"," )
	    end if      

	    if (check)then

		if (prof3d .EQ. 0) then
		  if (log_z .EQ. 1) then
		      z = f_att->$vNam(varn+1)$(1:dimz-1) 
		  else
		      z = f_att->$vNam(varn+1)$               
		  end if
		else 
		  do i=0,b-1            
		      if (isStrSubset( a(i),"zu_3d" ))then
			z_v(varn,:) = z_u
			if (log_z .EQ. 1) then
			    z = z_v(varn,1:dimz-1)
			else
			    z = z_v(varn,:)
			end if 
		      else
			if (isStrSubset( a(i),"zw_3d" ))then
			    z_v(varn,:) = z_w
			    if (log_z .EQ. 1) then
			      z = z_v(varn,1:dimz-1)
			    else
			      z = z_v(varn,:)
			    end if
			end if                   
		      end if
		  end do            
		end if

		z=z/norm_z
		
		if (abs(min(data(varn,:,min_z_int:max_z_int))) .GT. 10)then
		  min_value = abs(0.01*min(data(varn,:,min_z_int:max_z_int)))
		  max_value = abs(0.01*max(data(varn,:,min_z_int:max_z_int)))
		else
		  if (abs(min(data(varn,:,min_z_int:max_z_int))) .LT. 0.01 .AND. \
		      abs(max(data(varn,:,min_z_int:max_z_int))) .GT. 0.01)then
		      min_value = abs(0.1*max(data(varn,:,min_z_int:max_z_int)))
		      max_value = abs(0.1*max(data(varn,:,min_z_int:max_z_int)))
		  else
		      if (abs(max(data(varn,:,:))) .LT. 0.01 .AND. \
			  abs(min(data(varn,:,min_z_int:max_z_int))) .GT. 0.01)then
			min_value = abs(0.1*min(data(varn,:,min_z_int:max_z_int)))
			max_value = abs(0.1*min(data(varn,:,min_z_int:max_z_int)))
		      else
			min_value = abs(0.1*min(data(varn,:,min_z_int:max_z_int)))
			max_value = abs(0.1*max(data(varn,:,min_z_int:max_z_int)))
		      end if
		  end if
		end if 
		if (min(data(varn,:,min_z_int:max_z_int)) .EQ. 0 .AND.\
		    max(data(varn,:,min_z_int:max_z_int)) .EQ. 0)then
		  min_value = 0.1
		  max_value = 0.1
		end if

		res@gsnLeftString      = vNam(varn)
		res@tiXAxisString      = "("+unit(varn)+")"
		res@gsnRightString     = " "
		res@trYMinF            = min_z
		res@trYMaxF            = max_z  
		res@xyDashPattern = 0

		if (xs .EQ. -1) then 
		  res@trXMinF = min(data(varn,:,min_z_int:max_z_int))-min_value
		else 
		  res@trXMinF = xs
		end if
		if (xe .EQ. -1) then 
		  res@trXMaxF = max(data(varn,:,min_z_int:max_z_int))+max_value
		else 
		  res@trXMaxF = xe  
		end if 
		plot(n) = gsn_csm_xy(wks,data(varn,:,:),z,res)  
		
		if (vNam(varn) .EQ. "u" .OR. vNam(varn) .EQ. "v" .OR. \
		    vNam(varn) .EQ. "w") then
		  if (u .EQ. 0) then
		      res@gsnLeftString      = "u, v and w"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniu,miniv,miniw/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiu,maxiv,maxiw/))
		      else 
			res@trXMaxF = xe  
		      end if  
		      if (vNam(varn) .EQ. "v")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "w")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.07           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(wks,3,(/"u","v","w"/),lgres)       

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_u)
		      u=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue
		  end if        
		end if  
	
		if (vNam(varn) .EQ. "pt" .OR. vNam(varn) .EQ. "vpt" .OR. \
		    vNam(varn) .EQ. "lpt") then
		  if (pt .EQ. 0) then
		      res@gsnLeftString      = "pt, vpt and lpt"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/minipt,minivpt,minilpt/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxipt,maxivpt,maxilpt/))
		      else 
			res@trXMaxF = xe  
		      end if
		      if (vNam(varn) .EQ. "vpt")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "lpt")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0          
		      lgres@vpWidthF           = 0.07           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(wks,3,(/"pt","vpt","lpt"/),lgres) 
		      amres = True
		      amres@amParallelPosF   = 0.88	    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_pt)
		      pt=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue       
		  end if
		end if            
		if (vNam(varn) .EQ. "q" .OR. vNam(varn) .EQ. "qv" .OR. \
		    vNam(varn) .EQ. "ql") then
		  if (q .EQ. 0) then
		      res@gsnLeftString      = "q, qv and ql"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniq,miniqv,miniql/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiq,maxiqv,maxiql/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "qv")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "ql")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0          
		      lgres@vpWidthF           = 0.07          
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(wks,3,(/"q","qv","ql"/),lgres)

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_q)
		      q=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if    
	      
		if (vNam(varn) .EQ. "e" .OR. vNam(varn) .EQ. "es" .OR. \
		    vNam(varn) .EQ. "e*" ) then
		  if (e .EQ. 0) then
		      res@gsnLeftString      = "e and e*"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/minie,minies/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxie,maxies/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "es")then
			res@xyDashPattern = 1
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0           
		      lgres@vpWidthF           = 0.07           
		      lgres@vpHeightF          = 0.08          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(wks,2,(/"e","e*"/),lgres)       

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.365            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_e)
		      e=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if            
		if (vNam(varn) .EQ. "km" .OR. vNam(varn) .EQ. "kh") then
		  if (km .EQ. 0) then
		      res@gsnLeftString      = "km and kh"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/minikm,minikh/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxikm,maxikh/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "kh")then
			res@xyDashPattern = 1
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0           
		      lgres@vpWidthF           = 0.07           
		      lgres@vpHeightF          = 0.08          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(wks,2,(/"km","kh"/),lgres)       

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.365            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_km) 
		      km=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if           
	      
		if (vNam(varn) .EQ. "wpup" .OR. vNam(varn) .EQ. "wsus" .OR.      \
		    vNam(varn) .EQ. "wu" .OR. vNam(varn) .EQ. "w"+dq+"u"+dq .OR. \
		    vNam(varn) .EQ. "w*u*") then
		  if (wpup .EQ. 0) then
		      res@gsnLeftString      = "w"+dq+"u"+dq+", w*u* and wu"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpup,miniwsus,miniwu/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpup,maxiwsus,maxiwu/))
		      else 
			res@trXMaxF = xe  
		      end if  
		      if (vNam(varn) .EQ. "wsus")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wu")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.08           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
				    wks,3,(/"w"+dq+"u"+dq,"w*u*","wu"/),lgres)

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpup)
		      wpup=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue    
		  end if
		end if
		if (vNam(varn) .EQ. "wpvp" .OR. vNam(varn) .EQ. "wsvs" .OR.     \
		  vNam(varn) .EQ. "wv" .OR. vNam(varn) .EQ. "w"+dq+"v"+dq .OR. \
		  vNam(varn) .EQ. "w*v*") then
		  if (wpvp .EQ. 0) then
		      res@gsnLeftString      = "w"+dq+"v"+dq+", w*v* and wv"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpvp,miniwsvs,miniwv/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpvp,maxiwsvs,maxiwv/))
		      else 
			res@trXMaxF = xe  
		      end if
		      if (vNam(varn) .EQ. "wsvs")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wv")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.08          
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
				wks,3,(/"w"+dq+"v"+dq,"w*v*","wv"/),lgres)       

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpvp)
		      wpvp=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if
		if (vNam(varn) .EQ. "wpptp" .OR. vNam(varn) .EQ. "wspts" .OR.     \
		    vNam(varn) .EQ. "wpt" .OR. vNam(varn) .EQ. "w"+dq+"pt"+dq .OR.\
		    vNam(varn) .EQ. "w*pt*") then
		  if (wpptp .EQ. 0) then                 
		      res@gsnLeftString      = "w"+dq+"pt"+dq+", w*pt* and wpt"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpptp,miniwspts,miniwpt/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpptp,maxiwspts,maxiwpt/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "wspts" .OR. vNam(varn) .EQ. "w*pt*")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wpt")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.09           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)              
		      lbid = gsn_create_legend(\
				  wks,3,(/"w"+dq+"pt"+dq,"w*pt*","wpt"/),lgres)

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpptp)
		      wpptp=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue    
		  end if
		end if
		if (vNam(varn) .EQ. "wsptsBC" .OR. vNam(varn) .EQ. "wptBC" .OR.\
		    vNam(varn) .EQ. "w*pt*BC") then
		  if (wsptsBC .EQ. 0) then
		      res@gsnLeftString      = "w*pt*BC and wptBC"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwsptsBC,miniwptBC/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwsptsBC,maxiwptBC/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "wptBC")then
			res@xyDashPattern = 1
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.1           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
					    wks,3,(/"w*pt*BC","wptBC"/),lgres)

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wsptsBC)
		      wsptsBC=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if  
		end if             
		if (vNam(varn) .EQ. "wpvptp" .OR. vNam(varn) .EQ. "wsvpts" .OR. \
		    vNam(varn) .EQ. "wvpt" .OR. vNam(varn) .EQ. \
		    "w"+dq+"vpt"+dq .OR. vNam(varn) .EQ. "w*vpt*") then
		  if (wpvptp .EQ. 0) then
		      res@gsnLeftString      = "w"+dq+"vpt"+dq+", w*vpt* and wvpt"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpvptp,miniwsvpts,miniwvpt/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpvptp,maxiwsvpts,maxiwvpt/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "wsvpts")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wvpt")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.1           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
				wks,3,(/"w"+dq+"vpt"+dq,"w*vpt*","wvpt"/),lgres)
		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpvptp)
		      wpvptp=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if
		if (vNam(varn) .EQ. "wpqp" .OR. vNam(varn) .EQ. "wsqs" .OR.      \
		    vNam(varn) .EQ. "wq" .OR. vNam(varn) .EQ. "w"+dq+"q"+dq .OR. \
		    vNam(varn) .EQ. "w*q*") then
		  if (wpqp .EQ. 0) then 
		      res@gsnLeftString      = "w"+dq+"q"+dq+", w*q* and wq"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpqp,miniwsqs,miniwq/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpqp,maxiwsqs,maxiwq/))
		      else 
			res@trXMaxF = xe  
		      end if  
		      if (vNam(varn) .EQ. "wsqs")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wq")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0           
		      lgres@vpWidthF           = 0.08           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
				wks,3,(/"w"+dq+"q"+dq,"w*q*","wq"/),lgres)       

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpqp)
		      wpqp=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if 
		      continue   
		  end if
		end if
		if (vNam(varn) .EQ. "wpqvp" .OR. vNam(varn) .EQ. "wsqvs" .OR.     \
		    vNam(varn) .EQ. "wqv" .OR. vNam(varn) .EQ. "w"+dq+"qv"+dq .OR.\
		    vNam(varn) .EQ. "w*qv*") then
		  if (wpqvp .EQ. 0) then
		      res@gsnLeftString      ="w"+dq+"qv"+dq+" , w*qv* and wqv"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpqp,miniwsqvs,miniwqv/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpqp,maxiwsqvs,maxiwqv/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "wsqvs")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wqv")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.09           
		      lgres@vpHeightF          = 0.12         
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
				    wks,3,(/"w"+dq+"qv"+dq,"w*qv*","wqv"/),lgres)

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpqvp) 
		      wpqvp=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if
		if (vNam(varn) .EQ. "wpsp" .OR. vNam(varn) .EQ. "wsss" .OR.     \
		    vNam(varn) .EQ. "ws" .OR. vNam(varn) .EQ. "w"+dq+"s"+dq .OR.\
		    vNam(varn) .EQ. "w*s*") then
		  if (wpsp .EQ. 0) then
		      res@gsnLeftString      = "w"+dq+"s"+dq+", w*s* and ws"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpsp,miniwsss,miniws/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpsp,maxiwsss,maxiws/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "wsss")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "ws")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.08           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
			      wks,3,(/"w"+dq+"s"+dq,"w*s*","ws"/),lgres)       

		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpsp)
		      wpsp=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if 
		if (vNam(varn) .EQ. "wpsap" .OR.vNam(varn) .EQ. "wssas" .OR.      \
		    vNam(varn) .EQ. "wsa" .OR. vNam(varn) .EQ. "w"+dq+"sa"+dq .OR.\
		    vNam(varn) .EQ. "w*sa*") then
		  if (wpsap .EQ. 0) then
		      res@gsnLeftString      = "w"+dq+"sa"+dq+", w*sa* and wsa"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwpsap,miniwssas,miniwsa/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwpsap,maxiwssas,maxiwsa/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "wssas")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wsa")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.09           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
				wks,3,(/"w"+dq+"sa"+dq,"w*sa*","wsa"/),lgres)
		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wpsap)
		      wpsap=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if
	      
		if (vNam(varn) .EQ. "us2" .OR. vNam(varn) .EQ. "vs2" .OR. \
		    vNam(varn) .EQ. "ws2" .OR. vNam(varn) .EQ. "u*2" .OR. \
		    vNam(varn) .EQ. "v*2" .OR. vNam(varn) .EQ. "w*2" ) then
		  if (us2 .EQ. 0) then
		      res@gsnLeftString      = "u*2, v*2 and w*2"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/minius2,minivs2,miniws2/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxius2,maxivs2,maxiws2/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "vs2")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "ws2")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0            
		      lgres@vpWidthF           = 0.07           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(wks,3,(/"u*2","v*2","w*2"/),lgres)
		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_us2)
		      us2=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if
	      
		if (vNam(varn) .EQ. "wsususodz" .OR. \
		    vNam(varn) .EQ. "wspsodz" .OR.   \
		    vNam(varn) .EQ. "wpeodz" .OR.    \
		    vNam(varn) .EQ. "w*u*u*:dz" .OR. \
		    vNam(varn) .EQ. "w*p*:dz" .OR.   \
		    vNam(varn) .EQ. "w"+dq+"e:dz") then
		  if (wsususodz .EQ. 0) then 
		      res@gsnLeftString      = "w*u*u*:dz, w*p*:dz and w"+dq+"e:dz"
		      res@tiXAxisString      = "("+unit(varn)+")"
		      res@gsnRightString     = " "
		      if (xs .EQ. -1) then
			res@trXMinF = min((/miniwsususodz,\
						      miniwspsodz,miniwpeodz/))
		      else 
			res@trXMinF = xs
		      end if
		      if (xe .EQ. -1) then
			res@trXMaxF = max((/maxiwsususodz,maxiwspsodz,\
								  maxiwpeodz/))
		      else 
			res@trXMaxF = xe  
		      end if 
		      if (vNam(varn) .EQ. "wspsodz")then
			res@xyDashPattern = 1
		      end if
		      if (vNam(varn) .EQ. "wpeodz")then
			res@xyDashPattern = 2
		      end if
		      plot(n) =  gsn_csm_xy(wks,data(varn,:,:),z,res) 

		      ; ***************************************************
		      ; legend for overlaid plot
		      ; ***************************************************
	  
		      lgres                    = True
		      lgMonoDashIndex          = False
		      lgres@lgLabelFont        = "helvetica"    
		      lgres@lgLabelFontHeightF = font_size_legend*10.0           
		      lgres@vpWidthF           = 0.12           
		      lgres@vpHeightF          = 0.12          
		      lgres@lgDashIndexes      = (/0,1,2/)
		      lbid = gsn_create_legend(\
			      wks,3,(/"w*u*u*:dz","w*p*:dz","w"+dq+"e:dz"/),lgres)
		      amres = True
		      amres@amParallelPosF   = 0.88		    
		      amres@amOrthogonalPosF = 0.33            
		      annoid1 = gsn_add_annotation(plot(n),lbid,amres)
		      overlay(plot(n),plot_wsususodz)
		      wsususodz=1
		  else
		      if (prof3d .EQ. 0)then
			varn=varn+1
		      end if
		      continue   
		  end if
		end if     
		n=n+1 
		if (prof3d .EQ. 0)then
		  varn=varn+1
		end if   
	    end if
	  end do 
      end if 
      
      if (com_i .NE. 0) then
	  delete(com_var_avail)
      end if
      com_var_avail=new(count_var,string)    

      if (combine .EQ. 1) then
	  co=0     
	  n_o=0
	  do varn = 0,dim-1
    
	    check = True
	  
	    if (prof3d .EQ. 0) then
		if ( isStrSubset( vNam(varn), "time") .OR. \
		    isStrSubset( vNam(varn), "NORM")) then
		  check = False
		end if
	    else
		if ( isStrSubset( vNam(varn), "time") .OR.  \
		    isStrSubset( vNam(varn), "zusi") .OR.  \
		    isStrSubset( vNam(varn), "zwwi") .OR.  \
		    isStrSubset( vNam(varn), "x") .OR.     \
		    isStrSubset( vNam(varn), "xu") .OR.    \
		    isStrSubset( vNam(varn), "y") .OR.     \
		    isStrSubset( vNam(varn), "yv") .OR.    \
		    isStrSubset( vNam(varn), "zu_3d") .OR. \
		    isStrSubset( vNam(varn), "zw_3d")) then
		  check = False
		end if
	    end if

	    if (var .NE. "all") then         
		check = isStrSubset( var,","+vNam(varn)+"," )
	    end if      

	    if (check)then
		
		if (prof3d .EQ. 0) then
		  if (log_z .EQ. 1) then
		      z = f_att->$vNam(varn+1)$(1:dimz-1) 
		  else
		      z = f_att->$vNam(varn+1)$               
		  end if
		else 
		  do i=0,b-1            
		      if (isStrSubset( a(i),"zu_3d" ))then
			z_v(varn,:) = z_u
			if (log_z .EQ. 1) then
			    z = z_v(varn,1:dimz-1)
			else
			    z = z_v(varn,:)
			end if 
		      else
			if (isStrSubset( a(i),"zw_3d" ))then
			    z_v(varn,:) = z_w
			    if (log_z .EQ. 1) then
			      z = z_v(varn,1:dimz-1)
			    else
			      z = z_v(varn,:)
			    end if
			end if                   
		      end if
		  end do            
		end if

		z=z/norm_z
		
		com_var_avail(n_o)=vNam(varn)
		
		com=isStrSubset( c_var,","+vNam(varn)+"," )
		if (com)then 
		  co = co+1           
		  if (n_o .EQ. 1) then
		      res@xyDashPattern  = 1                                   
		  else            
		      if (n_o .EQ. 2) then
			res@xyDashPattern  = 2
		      else
			res@xyDashPattern  = 0
			res@gsnLeftString  = " " ; "Combined Plot of "+c_var_name
			res@tiXAxisString  = "("+unit(varn)+")"
			res@gsnRightString = " "
			if (xs .EQ. -1) then 
			    res@trXMinF = min(mini)
			else 
			    res@trXMinF = xs
			end if
			if (xe .EQ. -1) then 
			    res@trXMaxF = max(maxi)
			else 
			    res@trXMaxF = xe  
			end if 
		      end if
		  end if
		  label(n_o)=" " + vNam(varn)
		  color_o(n_o)=237
		  plot_o(n_o)=gsn_csm_xy(wks,data(varn,:,:),z,res) 
		  n_o=n_o+1
		end if 
		if (prof3d .EQ. 0)then
		  varn=varn+1
		end if           
	    end if
	  end do
      
	  if(number_comb .EQ. 2)then
	    if (co .EQ. 2)then 
		overlay(plot_o(0),plot_o(1))
	    else
		print(" ")
		print("combining is not possible,")
		print("'c_var'(= "+c_var+") must include two variables of "+\
		      "the general plots = ")
		print("- "+com_var_avail)
		print("be sure to have one comma before and after the variable")
		print(" ")
		exit  
	    end if 
	  end if
	  if(number_comb .EQ. 3)then
	    if (co .EQ. 3)then 
		overlay(plot_o(0),plot_o(1))
		overlay(plot_o(0),plot_o(2))
	    else
		print(" ")
		print("combining is not possible,")
		print("'c_var'(= "+c_var+") must include three variables of "+\
		      "the general plots = ")
		print("- "+com_var_avail)
		print("be sure to have one comma before and after the variable")
		print(" ")
		exit
	    end if 
	  end if

	  ; ***************************************************
	  ; legend for combined plot
	  ; ***************************************************
	  
	  lgres                    = True
	  lgMonoDashIndex          = False
          lgres@lgLineDashSegLenF  = 0.075
	  lgres@lgDashIndexes      = (/0,1,2/)
	  lgres@lgLabelFont        = "helvetica"    
	  lgres@vpWidthF           = 0.14           
	  lgres@vpHeightF          = 0.1 * number_comb / 3.0 + (number_comb + 1) * 0.01        
          lgres@lgLabelFontHeightF = 10.0 * font_size_legend * ((number_comb \
				     - 2.0) * 0.35 + 0.65)

	  lbid = gsn_create_legend(wks,number_comb,label,lgres)       

	  amres = True
	  amres@amParallelPosF   = -0.5		    
	  amres@amOrthogonalPosF = -0.5   
	  amres@amJust           = "TopLeft"

	  annoid1 = gsn_add_annotation(plot_o(0),lbid,amres) 
	
	  plot(0) = plot_o(0)

      end if

      if (c_var_log .EQ. 1) then
	  c_var_plot(com_i) = plot(0)
      end if

   end do

   ; ***************************************************
   ; merge plots onto one page 
   ; ***************************************************

   if (c_var_log .EQ. 1) then
      delete(plot_)
      plot_=new(max_com_i+1,graphic)
      plot_=c_var_plot
      n = max_com_i + 1 
   else
      do m=0,n-1
         plot_(m)=plot(n-1-m)
      end do
   end if

   no_frames = 0

   if ((format_out .EQ. "eps" .OR. format_out .EQ. "epsi") .AND. \
       n .gt. no_rows*no_columns) then
      gsn_panel(wks,plot_,(/n,1/),resP)
      print(" ")
      print("Outputs to .eps or .epsi have only one frame")
      print(" ")
   else    
      do i = 0,n-1, no_rows*no_columns
         if( (i+no_rows*no_columns) .gt. (n-1)) then
            gsn_panel(wks,plot_(i:n-1),(/no_rows,no_columns/),resP)
            no_frames = no_frames + 1  
         else
            gsn_panel(wks,plot_(i:i+no_rows*no_columns-1),\
                                                (/no_rows,no_columns/),resP)
            no_frames = no_frames + 1  
         end if
      end do
   end if


   if (format_out .EQ. "png" ) then
     png_output = new((/no_frames/), string)
     j = 0

     if (no_frames .eq. 1) then
        if (ncl_version .GE. 6 ) then
           png_output(0) = file_out+".png"
        else
           png_output(0) = file_out+".00000"+1+".png"
        end if
        ;using imagemagick's convert for reducing the white 
        ;space around the plot
        cmd = "convert -geometry 1000x1000 -density 300 -trim " +  \
              png_output(0) + " " + png_output(0)
       system(cmd)
     else

       do i=0, no_frames-1
         j = i + 1
         if (j .LE. 9) then
           png_output(i) = file_out+".00000"+j+".png"
         end if
         if (j .GT. 9 .AND. j .LE. 99) then
           png_output(i) = file_out+".0000"+j+".png"
         end if
         if (j .GT. 99 .AND. j .LE. 999) then
           png_output(i) = file_out+".000"+j+".png"
         end if
         if (j .GT. 999) then
           png_output(i) = file_out+".00"+j+".png"
         end if

         ;using imagemagick's convert for reducing the white 
         ;space around the plot
         cmd = "convert -geometry 1000x1000 -density 300 -trim " +  \
                png_output(i) + " " + png_output(i)
         system(cmd)
       end do
     end if

     print(" ")
     print("Output to: "+ png_output)
     print(" ")
   else
     print(" ")
     print("Output to: " + file_out +"."+ format_out)
     print(" ")
   end if

end
