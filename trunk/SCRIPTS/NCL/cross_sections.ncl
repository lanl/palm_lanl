load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;***************************************************
; Checking the kind of the script
;***************************************************

function which_script()
local script 
begin 
   script = "cross_section"
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

   ; ***************************************************
   ; Retrieving the double quote character
   ; ***************************************************
   
   dq=str_get_dq()

   ; ***************************************************
   ; set default parameter values and strings if not assigned 
   ; in prompt or parameter list
   ; ***************************************************

   ; Selection of type of cross section is done in palmplot
   if (.not. isvar("xyc"))then        
      xyc = 0	
   end if
   if (.not. isvar("xzc"))then  
      xzc = 0	
   end if
   if (.not. isvar("yzc"))then     
      yzc = 0
   end if
   
   if (file_1 .EQ. "File in") then
      print(" ")
      print("Declare input file 'file_1=' in '.ncl.config' or prompt")
      print(" ")
      exit 
   else
      file_in = file_1   
   end if

   if (format_out .NE. "x11" .AND. format_out .NE. "pdf" .AND.   \
       format_out .NE. "eps" .AND. format_out .NE. "ps" .AND.    \
       format_out .NE. "epsi" .AND. format_out .NE. "ncgm" .AND. \
       format_out .NE. "png") then
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

   if (sort .NE. "layer" .AND. sort .NE. "time") then  
      print(" ")
      print("'sort'= "+sort+" is invalid and set to 'layer'")
      print(" ")
      sort = "layer"  
   end if
   
   if (mode .NE. "Fill" .AND. mode .NE. "Line" .AND. mode .NE. "Both") then
      print(" ")
      print("'mode'= "+mode+" is invalid and set to 'Fill'")
      print(" ")
      mode = "Fill"
   end if

   if (fill_mode .NE. "AreaFill" .AND. fill_mode .NE. "RasterFill" .AND.\
       fill_mode .NE. "CellFill") then
      print(" ")
      print("'fill_mode'= "+fill_mode+" is invalid and set to 'AreaFill'")
      print(" ")
      fill_mode = "AreaFill"
   end if
   
   if (shape .NE. 0 .AND. shape .NE. 1) then
      print(" ")
      print("'shape'= "+shape+" is invalid and set to 1")
      print(" ")
      shape = 1
   end if 

   if (xyc .NE. 0 .AND. xyc .NE. 1)then
      print(" ")
      print("'xyc'= "+xyc+" is invalid and set to 0")
      print(" ")
      xyc = 0  
   end if
   
   if (xzc .NE. 0 .AND. xzc .NE. 1)then
      print(" ")
      print("'xzc'= "+xzc+" is invalid and set to 0")
      print(" ")
      xyc = 0  
   end if
   
   if (yzc .NE. 0 .AND. yzc .NE. 1)then
      print(" ")
      print("'yzc'= "+yzc+" is invalid and set to 0")
      print(" ")
      xyc = 0  
   end if
   
   if (xyc .EQ. 0 .AND. xzc .EQ. 0 .AND. yzc .EQ. 0) then
      print(" ")
      print("Select one crossection (xyc=1 or xzc=1 or yzc=1)")
      print(" ")
      exit
   end if
   if (xyc .EQ. 1 ) then 
      if (xzc .EQ. 1 .OR. yzc .EQ. 1) then
         print(" ")
         print("Select just one crossection (xyc=1 or xzc=1 or yzc=1)")
         print(" ") 
         exit
      end if
   end if
   if (xzc .EQ. 1) then
      if (xyc .EQ. 1 .OR. yzc .EQ. 1) then
         print(" ")
         print("Select just one crossection (xyc=1 or xzc=1 or yzc=1)")
         print(" ")
         exit
      end if
   end if
   if (yzc .EQ. 1) then
      if (xyc .EQ. 1 .OR. xzc .EQ. 1) then
         print(" ")
         print("Select just one crossection (xyc=1 or xzc=1 or yzc=1)")
         print(" ")
         exit
      end if
   end if

   if (vector .NE. 0 .AND. vector .NE. 1) then
      print(" ")
      print("Set 'vector' to 0 or 1")
      print(" ")
      exit
   end if  
   
   if (norm_x .EQ. 0) then 
      print(" ")
      print("Normalising with 0 is not allowed, 'norm_x' is set to 1.0")
      print(" ")
      norm_x = 1.0
   end if
   if (norm_y .EQ. 0) then 
      print(" ")
      print("Normalising with 0 is not allowed, 'norm_y' is set to 1.0")
      print(" ")
      norm_y = 1.0
   end if 
   if (norm_z .EQ. 0) then 
      print(" ")
      print("Normalising with 0 is not allowed, 'norm_z' is set to 1.0")
      print(" ")
      norm_z = 1.0
   end if 
 
   ; ***************************************************
   ; open input file
   ; ***************************************************
   
   file_in_1 = False
   if (isStrSubset(file_in, ".nc"))then
      start_f = -2
      end_f = -2 
      file_in_1 = True     
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

   files=new(end_f-start_f+1,string)

   if (file_in_1)then
      if (isfilepresent(file_in))then
         files(0)=file_in
      else
         print(" ")
         print("1st input file: '"+file_in+"' does not exist")
         print(" ")
         exit
      end if
   else   
      if (start_f .EQ. 0)then
         if (isfilepresent(file_in+".nc"))then
            files(0)=file_in+".nc"
            do i=1,end_f
	       if (isfilepresent(file_in+"."+i+".nc"))then    
                  files(i)=file_in+"."+i+".nc"
	       else
	          print(" ")
                  print("Input file: '"+file_in+"."+i+".nc' does not exist")
                  print(" ")
                  exit  
	       end if	    
            end do         
         else
            print(" ")
            print("Input file: '"+file_in+".nc' does not exist")
            print(" ")
            exit
         end if
      else
         do i=start_f,end_f
            if (isfilepresent(file_in+"."+i+".nc"))then    
               files(i-start_f)=file_in+"."+i+".nc"
	    else
	       print(" ")
               print("Input file: '"+file_in+"."+i+".nc' does not exist")
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

   print(" ")
   print("Variables in input file:")
   print("- "+ vNam)
   print(" ")
   dim   = dimsizes(vNam)

   ; ***************************************************
   ; check for kind of input file
   ; ***************************************************

   check_3d = True
   do varn=0,dim-1
      checkxy = isStrSubset( vNam(varn),"xy")
      checkxz = isStrSubset( vNam(varn),"xz")
      checkyz = isStrSubset( vNam(varn),"yz")
      if (checkxy .OR. checkxz .OR. checkyz) then
         check_3d = False
         break
      end if
   end do
   if (.not. check_3d) then
      if (xyc .EQ. 1)then
         if (.not. checkxy) then
            print(" ")
            print("Your input file doesn't have values for xy cross-sections")
            if (checkxz)then
               print("Select another input file or xz cross-sections")
            else
               print("Select another input file or yz cross-sections")
            end if
            print(" ")
            exit
         else
            print(" ")
            print("Your input file contains xy data")
            print(" ") 
	 end if
      end if
      if (xzc .EQ. 1)then
         if (.not. checkxz) then
            print(" ")
            print("Your input file doesn't have values for xz cross-sections")
            if (checkxy)then
               print("Select another input file or xy cross-sections")
            else
               print("Select another input file or yz cross-sections")
            end if
            print(" ")
            exit
         else
            print(" ")
            print("Your input file contains xz data")
            print(" ")
	 end if
      end if
      if (yzc .EQ. 1)then
         if (.not. checkyz) then
            print(" ")
            print("Your input file doesn't have values for yz cross-sections")
            if (checkxy)then
               print("Select another input file or xy cross-sections")
            else
               print("Select another input file or xz cross-sections")
            end if
            print(" ")
            exit
         else
            print(" ")
            print("Your input file contains yz data")
            print(" ")
	 end if
      end if
   else
      print(" ")
      print("Your input file contains 3d or other data")
      print(" ")
   end if
   
   if (dim .EQ. 0) then
      print(" ")
      print("There is no data on file")
      print(" ")
      exit
   end if

   ; ***************************************************
   ; set up recourses
   ; ***************************************************

   cs_res                          = True
   vecres                          = True
   cs_res@gsnYAxisIrregular2Linear = True  
   vecres@gsnYAxisIrregular2Linear = True
  
   cs_res@gsnDraw                 = False
   cs_res@gsnFrame                = False
   cs_res@gsnMaximize                = True
   
   cs_res@tmXBLabelFontHeightF   = font_size
   cs_res@tmYLLabelFontHeightF   = font_size
   cs_res@tiXAxisFontHeightF     = font_size
   cs_res@tiYAxisFontHeightF     = font_size
   cs_res@tmXBMode                ="Automatic"
   cs_res@tmYLMode                ="Automatic"
  
   cs_res@cnLevelSelectionMode    = "AutomaticLevels"
   cs_res@lbLabelFontHeightF = font_size_legend
   cs_res@lbLabelStride = legend_label_stride
   

   cs_resP = True
   cs_resP@gsnMaximize                = True
   cs_resP@gsnPanelXWhiteSpacePercent = 4.0
   cs_resP@gsnPanelYWhiteSpacePercent = 4.0
   cs_resP@txFont                     = "helvetica"
   cs_resP@txString                   = f_att@title
   cs_resP@txFuncCode                 = "~"         	
   cs_resP@txFontHeightF              = 0.0105        
  
   if ( mode .eq. "Fill" ) then
      cs_res@cnFillOn             = True
      cs_res@gsnSpreadColors      = True
      cs_res@cnFillMode           = fill_mode	
      cs_res@cnLinesOn            = False	
      cs_res@cnLineLabelsOn       = False
   end if

   if ( mode .eq. "Both" ) then
      cs_res@cnFillOn             = True
      cs_res@gsnSpreadColors      = True
      cs_res@cnFillMode           = fill_mode 
      cs_res@cnLinesOn            = True
      cs_res@cnLineLabelsOn       = True
   end if

   ; *********************************************
   ; set up of start_time_step and end_time_step
   ; *********************************************

   t_all         = f[:]->time
   nt            = dimsizes(t_all)

   if (nt .EQ. 1)then
      delta_t=t_all(nt-1)/nt
   else
      delta_t_array = new(nt-1,double)

      do i=0,nt-2
         delta_t_array(i) = t_all(i+1)-t_all(i)
      end do

      delta_t = min(delta_t_array)
      delete(delta_t_array)
   end if


   ; ****************************************************	
   ; start of time step and different types of mistakes that could be done
   ; ****************************************************

   if (start_time_step .EQ. -1.) then   	
      start_time_step=t_all(0)/3600      
   else 
      if (start_time_step .GT. t_all(nt-1)/3600)
         print(" ")
         print("'start_time_step' = "+ start_time_step +"h is greater than "+\
               "last time step = " + t_all(nt-1)+"s = "+t_all(nt-1)/3600+"h")
         print(" ")
         print("Select another 'start_time_step'")
         print(" ")
         exit
      end if
      if (start_time_step .LT. t_all(0)/3600)
         print(" ")
         print("'start_time_step' = "+ start_time_step +"h is lower than "+\
               "first time step = " + t_all(0)+"s = "+t_all(0)/3600+"h")
         print(" ")
         print("Select another 'start_time_step'")
         print(" ")
         exit
      end if
   end if
   do i=0,nt-1    
      if (start_time_step .GE. (t_all(i)-delta_t/2)/3600 .AND. \
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
      end_time_step = t_all(nt-1)/3600  
   else  
      if (end_time_step .LT. t_all(0)/3600)
         print(" ")
         print("'end_time_step' = "+end_time_step+ "h is lower than "+\
               "first time step = " + t_all(0)+"s = "+t_all(0)/3600+"h")
         print(" ")
         print("Select another 'end_time_step'")
         print(" ")
         exit
      end if
      if (end_time_step .GT. t_all(nt-1)/3600)
         print(" ")
         print("'end_time_step' = "+ end_time_step +"h is greater than "+\
               "last time step = " + t_all(nt-1)+"s = "+t_all(nt-1)/3600+"h")
         print(" ")
         print("Select another 'end_time_step'")  
         print(" ")
         exit
      end if
      if (end_time_step .LT. start_time_step/3600)
         print(" ")
         print("'end_time_step' = "+ end_time_step +"h is lower than "+\
               "'start_time_step' = "+start_time_step+"h")
         print(" ")
         print("Select another 'start_time_step' or 'end_time_step'")
         print(" ")
         exit
      end if
   end if
   do i=0,nt-1     
      if (end_time_step .GE. (t_all(i)-delta_t/2)/3600 .AND. \
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

   print(" ")
   print("Output of time steps from "+t_all(start_time_step)/3600+" h = "+\
                t_all(start_time_step)+" s => index = "+start_time_step)
   print("                     till "+t_all(end_time_step)/3600+" h = "  +\
                t_all(end_time_step)+" s => index = "+end_time_step)
   print(" ")
  
   no_time=(end_time_step-start_time_step)+1

   ; ****************************************************
   ; set up variables for vector plot if required
   ; ****************************************************   

   if (vector .EQ. 1) then
      if (vec1 .EQ. "component1") then
         print(" ")
         print("Indicate Vector 1 ('vec1') for Vector-Plot or "+\
               "set 'vector' to 0")
         print(" ")
         exit
      end if
      if (vec2 .EQ. "component2") then
         print(" ")
         print("Indicate Vector 2 ('vec2') for Vector-Plot or "+\
               "set 'vector' to 0")
         print(" ")
         exit
      end if            
   end if 

   check_vec1 = False
   check_vec2 = False
   check_vecp = False

   ; ****************************************************
   ; get data for plots
   ; ****************************************************

   if (xyc .EQ. 1) then
      do varn=0,dim-1
         if (vNam(varn) .eq. "xu" .OR. vNam(varn) .eq. "x")then    
            x_d     = f_att->$vNam(varn)$
            xdim    = dimsizes(x_d)-1
            delta_x = x_d(1)-x_d(0)
            break
         end if
      end do 
      do varn=0,dim-1
         if (vNam(varn) .eq. "yv" .or. vNam(varn) .eq. "y")then   
            y_d     = f_att->$vNam(varn)$
            ydim    = dimsizes(y_d)-1
            delta_y = y_d(1)-y_d(0)
            break
         end if
      end do
      do varn=0,dim-1
         if (vNam(varn) .eq. "zu_3d" .OR. vNam(varn) .eq. "zw_3d")then
            z_d     = f_att->$vNam(varn)$
            zdim    = dimsizes(z_d)-1
            delta_z = 0
            break
         else
            if (vNam(varn) .eq. "zu_xy" .OR. vNam(varn) .eq. "zw_xy") then
               z_d     = f_att->$vNam(varn)$
               zdim    = dimsizes(z_d)-1 
               delta_z = -1.d
               break
            else
               zdim    = 0
               delta_z = -1.d
            end if
         end if
         if (vNam(varn) .eq. "zu1_xy") then
            zu1    = f_att->$vNam(varn)$
         end if
      end do
   end if
   if (xzc .EQ. 1) then
      do varn=0,dim-1
         if (vNam(varn) .eq. "xu"  .OR. vNam(varn) .eq. "x") then 
            x_d     = f_att->$vNam(varn)$
            xdim    = dimsizes(x_d)-1
            delta_x = x_d(1)-x_d(0)
            break
         end if
      end do
      do varn=0,dim-1   
         if (vNam(varn) .eq. "zw" .or. vNam(varn) .eq. "zu" .OR. \
             vNam(varn) .eq. "zu_3d" .or. vNam(varn) .eq. "zw_3d") then
            z_d     = f_att->$vNam(varn)$
            zdim    = dimsizes(z_d)-1 
            delta_z = z_d(1)-z_d(0)
            break
         end if
      end do
      do varn=0,dim-1 
         if (vNam(varn) .eq. "yv" .or. vNam(varn) .eq. "y") then
            y_d     = f_att->$vNam(varn)$
            ydim    = dimsizes(y_d)-1
            delta_y = y_d(1)-y_d(0)
            break
         else
            if (vNam(varn) .eq. "y_xz" .or. vNam(varn) .eq. "yv_xz") then
               y_d     = f_att->$vNam(varn)$
               ydim    = dimsizes(y_d)-1 
               delta_y = -1.d
               break
            else
               ydim    = 0
               delta_y = -1.d
            end if
         end if 
      end do
   end if
   if (yzc .EQ. 1) then
      do varn=0,dim-1  
         if (vNam(varn) .eq. "yv" .or. vNam(varn) .eq. "y") then
            y_d     = f_att->$vNam(varn)$
            ydim    = dimsizes(y_d)-1
            delta_y = y_d(1)-y_d(0)
            break
         end if
      end do
      do varn=0,dim-1 
         if (vNam(varn) .eq. "zw" .or. vNam(varn) .eq. "zu" .OR. \
             vNam(varn) .eq. "zu_3d" .or. vNam(varn) .eq. "zw_3d") then
            z_d     = f_att->$vNam(varn)$
            zdim    = dimsizes(z_d)-1 
            delta_z = z_d(1)-z_d(0)
            break
         end if
      end do
      do varn=0,dim-1 
         if (vNam(varn) .eq. "xu" .or. vNam(varn) .eq. "x")
            x_d     = f_att->$vNam(varn)$
            xdim    = dimsizes(x_d)-1 
            delta_x = x_d(1)-x_d(0)
            break
         else
            if (vNam(varn) .eq. "xu_yz" .or. vNam(varn) .eq. "x_yz" ) then
               x_d     = f_att->$vNam(varn)$
               xdim    = dimsizes(x_d)-1 
               delta_x = -1.d
               break
            else
               xdim    = 0
               delta_x = -1.d
            end if
         end if 
      end do
   end if
 
   ; ****************************************************
   ; set up ranges of x-, y- and z-coordinates
   ; ****************************************************         
                    
   if (xs .EQ. -1.d) then     		
      xs = x_d(0)
      if (delta_x .EQ. -1) then
         print(" ") 
         print("You cannot choose a start value for x, "+\
               "there are preset layers for x")
         print(" ")
         xstart=0
      end if
   else 
      if (delta_x .EQ. -1) then
         print(" ") 
         print("You cannot choose a start value for x, "+\
               "there are preset layers for x")
         print(" ")
         xstart=0
      else
         if (xs .LT. 0-delta_x/2) then
            print(" ")
            print("range start for x-coordinate = "+\
                   xs+"m is lower than first value = "+(0-delta_x/2)+"m")
            print(" ")
            exit
         end if
         if (xyc .EQ. 1 .OR. xzc .EQ. 1) then
            if (xs .GE. xdim*delta_x) then
               print(" ")
               print("range start for x-coordinate = "+\
                     xs+"m is equal or greater than last value = "+\
                     xdim*delta_x+"m")
               print(" ")
               exit
            end if 
         else
            if (xs .GT. xdim*delta_x) then
               print(" ")
               print("range start for x-coordinate = "+\
                     xs+"m is greater than last value = "+\
                     xdim*delta_x+"m")
               print(" ")
               exit
            end if 
         end if
      end if
   end if 

   if (ys .EQ. -1.d) then    
      ys = y_d(0)
      if (delta_y .EQ. -1) then
         print(" ") 
         print("You cannot choose a start value for y, "+\
               "there are preset layers for y")
         print(" ")
         ystart=0
      end if
   else
      if (delta_y .EQ. -1) then
         print(" ") 
         print("You cannot choose a start value for y, "+\
                "there are preset layers for y")
         print(" ")
         ystart=0
      else
         if (ys .LT. 0-delta_y/2) then
            print(" ")
            print("range start for y-coordinate = "+\
                  ys+"m is lower than first value = "+0-delta_y/2+"m")
            print(" ")
            exit
         end if
         if (xyc .EQ. 1 .OR. yzc .EQ. 1) then
            if (ys .GE. ydim*delta_y) then
               print(" ")
               print("range start for y-coordinate = "+\
                     ys+"m is equal or greater than last value = "+\
                     ydim*delta_y+"m")
               print(" ")
               exit
            end if
         else
            if (ys .GT. ydim*delta_y) then
               print(" ")
               print("range start for y-coordinate = "+\
                     ys+"m is greater than last value = "+ydim*delta_y+"m")
               print(" ")
               exit
            end if
         end if
      end if
   end if
  
   if (zs .EQ. -1) then     			
      zs = 0
      if (delta_z .EQ. -1) then
         print(" ") 
         print("You cannot choose a start value for z, "+\
               "there are preset layers for z")
         print(" ")
      end if
   else
      if (delta_z .EQ. -1) then
         print(" ") 
         print("You cannot choose a start value for z, "+\
               "there are preset layers for z")
         print(" ")
         zs = 0
      else
         if (zs .LT. 0) then
            print(" ")
            print("range start for z-coordinate = "+\
                  zs+" is lower than first gridpoint = 0")
            print(" ")
            exit
         end if
         if (xzc .EQ. 1 .OR. yzc .EQ. 1) then
            if (zs .GE. zdim) then
               print(" ")
               print("range start for z-coordinate = "+\
                     zs+" is equal or greater than last gridpoint = "+zdim)
               print(" ")
               exit
            end if
         else
            if (zs .GT. zdim) then
               print(" ")
               print("range start for z-coordinate = "+\
                     zs+" is greater than last gridpoint = "+zdim)
               print(" ")
               exit
            end if
         end if
      end if
   end if  

   if (xe .EQ. -1) then    
      xe = x_d(xdim)
      if (delta_x .EQ. -1) then
         print(" ") 
         print("You cannot choose an end value for x, "+\
               "there are preset layers for x")
         print(" ")
         xend=xdim
      end if
   else
      if (delta_x .EQ. -1) then
         print(" ") 
         print("You cannot choose an end value for x, "+\
               "there are preset layers for x")
         print(" ")
         xend=xdim
      else
         if (xe .GT. xdim*delta_x) then
            print(" ")
            print("range end for x-coordinate = "+\
                  xe+"m is greater than last value = "+xdim*delta_x+"m")
            print(" ") 
            exit 
         end if 
         if (xyc .EQ. 1 .OR. xzc .EQ. 1) then
            if (xe .LE. 0-delta_x/2)
               print(" ")
               print("range end for x-coordinate = "+\
                     xe+"m is equal or lower than first value = "+\
                     (0-delta_x/2)+"m")
               print(" ")
               exit
            end if
            if (xe .LE. xs) then
               print(" ")
               print("range end for x-coordinate = "+\
                     xe+"m is equal or lower than start range = "+xs+"m")
               print(" ")
               exit
            end if
            if (xe .EQ. xs+1) then
               print(" ")
               print("range end for x-coordinate = "+\
                     xe+"m must be at least two more gridpoints "+\
                     "greater than start range = "+xs+"m")
               print(" ")
               exit
            end if
         else
            if (xe .LT. 0-delta_x/2)
               print(" ")
               print("range end for x-coordinate = "+\
                     xe+"m is lower than first value = "+(0-delta_x/2)+"m")
               print(" ")
               exit
            end if
            if (xe .LT. xs) then
               print(" ")
               print("range end for x-coordinate = "+\
                     xe+"m is lower than start range = "+xs+"m")
               print(" ")
               exit
            end if
         end if 
      end if               
   end if
   
   if (ye .EQ. -1) then  
      ye = y_d(ydim)
      if (delta_y .EQ. -1) then
         print(" ") 
         print("You cannot choose an end value for y, "+\
               "there are preset layers for y")
         print(" ")
         yend=ydim
      end if
   else
      if (delta_y .EQ. -1) then
         print(" ") 
         print("You cannot choose an end value for y, "+\
               "there are preset layers for y")
         print(" ")
         yend=ydim
      else
         if (ye .GT. ydim*delta_y) then
            print(" ")
            print("range end for y-coordinate = "+\
                  ye+"m is greater than last value = "+ydim*delta_y+"m")
            print(" ") 
            exit 
         end if 
         if (xyc .EQ. 1 .OR. yzc .EQ. 1) then
            if (ye .LE. 0-delta_y/2)
               print(" ")
               print("range end for y-coordinate = "+\
                     ye+"m is equal or lower than first value = "+\
                     (0-delta_y/2)+"m")
               print(" ")
               exit
            end if
            if (ye .LE. ys) then
               print(" ")
               print("range end for y-coordinate = "+\
                    ye+"m is equal or lower than start range = "+ys+"m")
               print(" ")
               exit
            end if
            if (ye .EQ. ys+1) then
               print(" ")
               print("range end for y-coordinate = "+\
                     ye+"m must be at least two more gridpoints greater "+\
                     "than start range = "+ys+"m")
               print(" ")
               exit
            end if
         else
            if (ye .LT. 0-delta_y/2)
               print(" ")
               print("range end for y-coordinate = "+\
                     ye+"m is lower than first value = "+(0-delta_y/2)+"m")
               print(" ")
               exit
            end if
            if (ye .LT. ys) then
               print(" ")
               print("range end for y-coordinate = "+\
                     ye+"m is lower than start range = "+ys+"m")
               print(" ")
               exit
            end if
         end if
      end if
   end if 
 
   if (ze .EQ. -1) then  
      ze = zdim
      if (delta_z .EQ. -1) then
         print(" ") 
         print("You cannot choose an end value for z, "+\
               "there are preset layers for z")
         print(" ")
         ze = zdim
      end if
   else
      if (delta_z .EQ. -1) then
         print(" ") 
         print("You cannot choose an end value for z, "+\
               "there are preset layers for z")
         print(" ")
         ze = zdim
      else
         if (ze .GT. zdim) then
            print(" ")
            print("range end for z-coordinate = "+\
                  ze+" is greater than last gridpoint = "+zdim)
            print(" ") 
            exit 
         end if 
         if (xzc .EQ. 1 .OR. yzc .EQ. 1) then
            if (ze .LE. 0)
               print(" ")
               print("range end for z-coordinate = "+\
                     ze+" is equal or lower than first gridpoint = 0")
               print(" ")
               exit
            end if
            if (ze .LE. zs) then
               print(" ")
               print("range end for z-coordinate = "+\
                     ze+" is equal or lower than start range = "+zs)
               print(" ")
               exit
            end if
            if (ze .EQ. zs+1) then
               print(" ")
               print("range end for z-coordinate = "+\
                     ze+" must be at least two more gridpoints "+\
                     "greater than start range = "+zs)
               print(" ")
               exit
            end if
         else
            if (ze .LT. 0)
               print(" ")
               print("range end for z-coordinate = "+\
                     ze+" is lower than first gridpoint = 0")
               print(" ")
               exit
            end if
            if (ze .LT. zs) then
               print(" ")
               print("range end for z-coordinate = "+\
                     ze+" is lower than start range = "+zs)
               print(" ")
               exit
            end if
         end if
      end if
   end if 

   if (delta_x .NE. -1) then
      do i=0,xdim    
         if (xs .GT. x_d(i)-delta_x/2 .AND. xs .LE. x_d(i)+delta_x/2)then
            xstart=i
            break
         end if
      end do 
      do i=0,xdim    
         if (xe .GT. x_d(i)-delta_x/2 .AND. xe .LE. x_d(i)+delta_x/2)then
            xend=i
            break
         end if
      end do
   end if
   if (delta_y .NE. -1) then
      do i=0,ydim    
         if (ys .GT. y_d(i)-delta_y/2 .AND. ys .LE. y_d(i)+delta_y/2)then
            ystart=i
            break
         end if
      end do 
      do i=0,ydim    
         if (ye .GT. y_d(i)-delta_y/2 .AND. ye .LE. y_d(i)+delta_y/2)then
            yend=i
            break
         end if
      end do 
   end if

   if( shape .eq. 1 ) then 
      if (xyc .EQ. 1)then
         cs_res@vpWidthF    = (xe-xs)/(ye-ys)        
         cs_res@vpHeightF   = 1
	 vecres@vpWidthF    = (xe-xs)/(ye-ys)        
         vecres@vpHeightF   = 1
	 if (xe-xs .GT. ye-ys)then
	    if (no_rows .LE. no_columns)then
	       cs_res@gsnPaperOrientation     = "landscape"
	       vecres@gsnPaperOrientation     = "landscape"
	       cs_res@lbOrientation = "Horizontal" 
	    else 
	       cs_res@gsnPaperOrientation     = "portrait"
	       vecres@gsnPaperOrientation     = "portrait"
	       cs_res@lbOrientation = "Vertical"
	    end if  
	 else
	    if (no_rows .GE. no_columns)then
	       cs_res@gsnPaperOrientation     = "portrait"
	       vecres@gsnPaperOrientation     = "portrait"
	       cs_res@lbOrientation = "Vertical" 
	    else
	       cs_res@gsnPaperOrientation     = "landscape"
	       vecres@gsnPaperOrientation     = "landscape"
	       cs_res@lbOrientation = "Horizontal"
	    end if
	 end if
      end if
      if (xzc .EQ. 1)then
         cs_res@vpWidthF    = (xe-xs)/(delta_x*(ze-zs))        
         cs_res@vpHeightF   = 1
	 vecres@vpWidthF    = (xe-xs)/(delta_x*(ze-zs))        
         vecres@vpHeightF   = 1
	 if (xe-xs .GT. (delta_x*(ze-zs)))then
	    if (no_rows .LE. no_columns)then
	       cs_res@gsnPaperOrientation     = "landscape"
	       vecres@gsnPaperOrientation     = "landscape"
	       cs_res@lbOrientation = "Horizontal" 
	    else 
	       cs_res@gsnPaperOrientation     = "portrait"
	       vecres@gsnPaperOrientation     = "portrait"
	       cs_res@lbOrientation = "Vertical"
	    end if  
	 else
	    if (no_rows .GE. no_columns)then
	       cs_res@gsnPaperOrientation     = "portrait"
	       vecres@gsnPaperOrientation     = "portrait"
	       cs_res@lbOrientation = "Vertical" 
	    else
	       cs_res@gsnPaperOrientation     = "landscape"
	       vecres@gsnPaperOrientation     = "landscape"
	       cs_res@lbOrientation = "Horizontal"
	    end if
	 end if
      end if
      if (yzc .EQ. 1)then
         cs_res@vpWidthF    = (ye-ys)/(delta_y*(ze-zs))        
         cs_res@vpHeightF   = 1
	 vecres@vpWidthF    = (ye-ys)/(delta_y*(ze-zs))        
         vecres@vpHeightF   = 1
	 if (ye-ys .GT. (delta_y*(ze-zs)))then
            if (no_rows .LE. no_columns)then
	       cs_res@gsnPaperOrientation     = "landscape"
	       vecres@gsnPaperOrientation     = "landscape"
	       cs_res@lbOrientation = "Horizontal" 
	    else 
	       cs_res@gsnPaperOrientation     = "portrait"
	       vecres@gsnPaperOrientation     = "portrait"
	       cs_res@lbOrientation = "Vertical"
	    end if  
	 else
	    if (no_rows .GE. no_columns)then
	       cs_res@gsnPaperOrientation     = "portrait"
	       vecres@gsnPaperOrientation     = "portrait"
	       cs_res@lbOrientation = "Vertical" 
	    else
	       cs_res@gsnPaperOrientation     = "landscape"
	       vecres@gsnPaperOrientation     = "landscape"
	       cs_res@lbOrientation = "Horizontal"
	    end if
	 end if
      end if
   end if

   delete(xs)
   delete(xe)   
   delete(ys)
   delete(ye)

   xs=xstart
   xe=xend
   ys=ystart
   ye=yend

   if (xyc .EQ. 1) then
      d= (ye-ys+1)/(major_ticks_y-1)
      e=(xe-xs+1)/(major_ticks_x-1)
      array_yl       =new(major_ticks_y,integer)
      array_yl_labels=new(major_ticks_y,double)
      array_minor_yl =new((major_ticks_y-1)*4,double)
      array_xb       =new(major_ticks_x,integer)
      array_xb_labels=new(major_ticks_x,double)
      array_minor_xb =new((major_ticks_x-1)*4,double)
      array_yl(0)=ys
      array_xb(0)=xs
      array_yl_labels(0)=y_d(ys)
      array_xb_labels(0)=x_d(xs)
      do ar=1,major_ticks_y-1
         array_yl(ar)=d*(ar-1)+d-1
	 array_yl_labels(ar) = y_d(array_yl(ar))
	 if (ar .GT. 0)
	    do min_ar=0,3
                array_minor_yl(4*(ar-1)+min_ar)= y_d(array_yl(ar-1))+\
                         (y_d(array_yl(ar))-y_d(array_yl(ar-1)))/5*(min_ar+1)
	    end do 
	 end if	 
      end do 
      do br=1,major_ticks_x-1
         array_xb(br)=e*(br-1)+e-1 
	 array_xb_labels(br) = x_d(array_xb(br)) 
	 if (br .GT. 0)
	    do min_br=0,3
	       array_minor_xb(4*(br-1)+min_br)= x_d(array_xb(br-1))+\
                         (x_d(array_xb(br))-x_d(array_xb(br-1)))/5*(min_br+1)
	    end do 
	 end if	
      end do 
   end if
 
   if (xzc .EQ. 1) then
      d=(ze-zs+1)/(major_ticks_y-1)
      e=(xe-xs+1)/(major_ticks_x-1)
      array_yl       =new(major_ticks_y,integer)
      array_yl_labels=new(major_ticks_y,double)
      array_minor_yl =new((major_ticks_y-1)*4,double)
      array_xb       =new(major_ticks_x,integer)
      array_xb_labels=new(major_ticks_x,double)
      array_minor_xb =new((major_ticks_x-1)*4,double)
      array_yl(0)=zs
      array_xb(0)=xs
      array_yl_labels(0)=z_d(zs)
      array_xb_labels(0)=x_d(xs)
      do ar=1,major_ticks_y-1
         array_yl(ar)=d*(ar-1)+d-1
	 array_yl_labels(ar) = z_d(array_yl(ar))
	 if (ar .GT. 0)
	    do min_ar=0,3
	       array_minor_yl(4*(ar-1)+min_ar)= z_d(array_yl(ar-1))+\
                          (z_d(array_yl(ar))-z_d(array_yl(ar-1)))/5*(min_ar+1)
	    end do 
	 end if	 
      end do 
      do br=1,major_ticks_x-1
         array_xb(br)=e*(br-1)+e-1 
	 array_xb_labels(br) = x_d(array_xb(br)) 
	 if (br .GT. 0)
	    do min_br=0,3
	       array_minor_xb(4*(br-1)+min_br)= x_d(array_xb(br-1))+\
                         (x_d(array_xb(br))-x_d(array_xb(br-1)))/5*(min_br+1)
	    end do 
	 end if	
      end do 
   end if

   if (yzc .EQ. 1) then
      d=(ze-zs+1)/(major_ticks_y-1)
      e=(ye-ys+1)/(major_ticks_x-1)
      array_yl       =new(major_ticks_y,integer)
      array_yl_labels=new(major_ticks_y,double)
      array_minor_yl =new((major_ticks_y-1)*4,double)
      array_xb       =new(major_ticks_x,integer)
      array_xb_labels=new(major_ticks_x,double)
      array_minor_xb =new((major_ticks_x-1)*4,double)
      array_yl(0)=zs
      array_xb(0)=ys
      array_yl_labels(0)=z_d(zs)
      array_xb_labels(0)=y_d(ys)
      do ar=1,major_ticks_y-1
         array_yl(ar)=d*(ar-1)+d-1
	 array_yl_labels(ar) = z_d(array_yl(ar))
	 if (ar .GT. 0)
	    do min_ar=0,3
               array_minor_yl(4*(ar-1)+min_ar)= z_d(array_yl(ar-1))+\
                         (z_d(array_yl(ar))-z_d(array_yl(ar-1)))/5*(min_ar+1)
	    end do 
	 end if	 
      end do 
      do br=1,major_ticks_x-1
         array_xb(br)=e*(br-1)+e-1 
	 array_xb_labels(br) = y_d(array_xb(br)) 
	 if (br .GT. 0)
	    do min_br=0,3
	       array_minor_xb(4*(br-1)+min_br)= y_d(array_xb(br-1))+\
                         (y_d(array_xb(br))-y_d(array_xb(br-1)))/5*(min_br+1)
	    end do 
	 end if	
      end do
   end if
 
   if (axes_explicit .EQ. 1)then
      cs_res@tmYLMode = "Explicit"
      cs_res@tmXBMode = "Explicit"
      if (xyc .EQ. 1)then 
         cs_res@tmYLValues = y_d(array_yl) 
         cs_res@tmXBValues = x_d(array_xb)      
         cs_res@tmYLLabels = array_yl_labels/norm_y
	 cs_res@tmXBLabels = array_xb_labels/norm_x
	 if (norm_x .NE. 1.)then
	    cs_res@tiXAxisString = "x ("+unit_x+")"
	 else
	    cs_res@tiXAxisString = "x (m)"
	 end if 
	 if (norm_y .NE. 1.)then
            cs_res@tiYAxisString = "y ("+unit_y+")"
	 else
	    cs_res@tiYAxisString = "y (m)"
	 end if   
      end if
      if (xzc .EQ. 1)then
         cs_res@tmYLValues = z_d(array_yl) 
         cs_res@tmXBValues = x_d(array_xb)  
         cs_res@tmYLLabels = array_yl_labels/norm_z
	 cs_res@tmXBLabels = array_xb_labels/norm_x
	 if (norm_x .NE. 1.)then
	    cs_res@tiXAxisString = "x ("+unit_x+")"
	 else
	    cs_res@tiXAxisString = "x (m)"
	 end if 
	 if (norm_z .NE. 1.)then
            cs_res@tiYAxisString = "z ("+unit_z+")"
	 else
	    cs_res@tiYAxisString = "z (m)"
	 end if
      end if
      if (yzc .EQ. 1)then
         cs_res@tmYLValues = z_d(array_yl) 
         cs_res@tmXBValues = y_d(array_xb)  
         cs_res@tmYLLabels = array_yl_labels/norm_z
	 cs_res@tmXBLabels = array_xb_labels/norm_y
	 if (norm_y .NE. 1.)then
	    cs_res@tiXAxisString = "y ("+unit_y+")"
	 else
	    cs_res@tiXAxisString = "y (m)"
	 end if 
	 if (norm_z .NE. 1.)then
            cs_res@tiYAxisString = "z ("+unit_z+")"
	 else
	    cs_res@tiYAxisString = "z (m)"
	 end if
      end if
      cs_res@tmYLMinorValues = array_minor_yl
      cs_res@tmXBMinorValues = array_minor_xb
   else 
      if (axes_explicit .EQ. 0)then  
         cs_res@tmYLMinorPerMajor = 4
         cs_res@tmXBMinorPerMajor = 4
      else
         print(" ")
	 print("'axes_explicit' is invalid and set to 0")
	 print(" ")
	 axes_explicit = 0
	 cs_res@tmYLMinorPerMajor = 4
         cs_res@tmXBMinorPerMajor = 4
      end if
      if (xyc .EQ. 1)then
         cs_res@tiXAxisString = "x (m)"
         cs_res@tiYAxisString = "y (m)"
      end if
      if (xzc .EQ. 1)then
         cs_res@tiXAxisString = "x (m)"
         cs_res@tiYAxisString = "z (m)"
      end if
      if (yzc .EQ. 1)then
         cs_res@tiXAxisString = "y (m)"
         cs_res@tiYAxisString = "z (m)"
      end if
   end if

   if (xyc .EQ. 1 .OR. xzc .EQ.1)then
      if (xe .EQ. xs+1) then
         print(" ")
         print("range end for x-coordinate="+xe*delta_x+\
               "m must be at least two") 
         print("more gridpoints(="+2*delta_x+"m) greater than start range="+\
               xs*delta_x+"m)")
         print(" ")
         exit
      end if 
   end if
   if (xyc .EQ. 1 .OR. yzc .EQ.1)then
      if (ye .EQ. ys+1) then
         print(" ")
         print("range end for y-coordinate="+ye*delta_y+\
               "m must be at least two") 
         print("more gridpoints(="+2*delta_y+"m greater than start range="+\
               ys*delta_y+"m)")
         print(" ")
         exit
      end if 
   end if
   if (xzc .EQ. 1 .OR. yzc .EQ.1)then
      if (ze .EQ. zs+1) then
         print(" ")
         print("range end for z-coordinate="+ze+" must be at least two") 
         print("more gridpoints greater than start range="+zs+")")
         print(" ")
         exit
      end if 
   end if
  
   if (xyc .EQ. 1) then
      no_layer = (ze-zs)+1
      if (check_vType) then
         data = new((/dim,nt,(ze-zs)+1,(ye-ys)+1,(xe-xs)+1/),double)
      else
         data = new((/dim,nt,(ze-zs)+1,(ye-ys)+1,(xe-xs)+1/),float)
      end if
   end if
   if (xzc .EQ. 1) then
      no_layer = (ye-ys)+1
      if (check_vType) then
         data = new((/dim,nt,(ze-zs)+1,(ye-ys)+1,(xe-xs)+1/),double)
      else
         data = new((/dim,nt,(ze-zs)+1,(ye-ys)+1,(xe-xs)+1/),float)
      end if
   end if
   if (yzc .EQ. 1) then
      no_layer = (xe-xs)+1
      if (check_vType) then
         data = new((/dim,nt,(ze-zs)+1,(ye-ys)+1,(xe-xs)+1/),double)
      else
         data = new((/dim,nt,(ze-zs)+1,(ye-ys)+1,(xe-xs)+1/),float)
      end if
   end if

   if (check_vType) then
      MinVal = new(dim,double)
      MaxVal = new(dim,double)
   else
      MinVal = new(dim,float)
      MaxVal = new(dim,float)
   end if
   unit   = new(dim,string) 

   ; ****************************************************
   ; define inner and outer loops depending on "sort"
   ; ****************************************************       
    
   if ( xyc .eq. 1 ) then
      lays = zs
      laye = ze
   end if
   if ( xzc .eq. 1 ) then
      lays = ys
      laye = ye
   end if
   if ( yzc .eq. 1) then
      lays = xs
      laye = xe
   end if

   if ( sort .eq. "time" ) then
      lis = start_time_step
      lie = end_time_step
      los = lays
      loe = laye
   else
      lis = lays
      lie = laye
      los = start_time_step 
      loe = end_time_step
   end if 

   check  = True
   v1     = 0
   v2     = 0
   no_var = 0
   n      = 0
   no_zu1 = 0
   no_topo= 0

   ;****************************************************
   ; Preparation of vector plots
   ; ****************************************************       

   do varn=dim-1,0,1 

     if (vector .EQ. 1) then    
        check_vec1 = isStrSubset( vec1,","+vNam(varn)+",")
        if (check_vec1) then
            temp = f[:]->$vNam(varn)$
	    data_att = f_att->$vNam(varn)$
            vect1=temp(:,zs:ze,ys:ye,xs:xe) 
            delete(temp)
         end if

         check_vec2 = isStrSubset( vec2,","+vNam(varn)+",")
         if (check_vec2) then
            temp = f[:]->$vNam(varn)$
	    data_att = f_att->$vNam(varn)$
            vect2=temp(:,zs:ze,ys:ye,xs:xe) 
            delete(temp)
         end if
    
         if (","+vNam(varn)+"," .EQ. vec1) then
            v1=v1+1
         end if
         if (","+vNam(varn)+"," .EQ. vec2) then
             v2=v2+1
          end if
     end if

   end do



   do varn=dim-1,0,1     
           
      if ( vNam(varn) .eq. "zw" .or. vNam(varn) .eq. "yv" .or.          \
           vNam(varn) .eq. "xu" .or. vNam(varn) .eq. "zu" .or.          \
           vNam(varn) .eq. "time" .or. vNam(varn) .eq. "zu_3d" .or.     \
           vNam(varn) .eq. "zw_3d" .or. vNam(varn) .eq. "x".or.         \
           vNam(varn) .eq. "y" .or. vNam(varn) .eq. "zu_xy" .or.        \
           vNam(varn) .eq. "zw_xy" .or. vNam(varn) .eq. "zu1_xy" .or.   \
           vNam(varn) .eq. "ind_z_xy" .or. vNam(varn) .eq. "y_xz" .or.  \
           vNam(varn) .eq. "yv_xz" .or. vNam(varn) .eq. "ind_y_xz" .or. \
           vNam(varn) .eq. "x_yz" .or. vNam(varn) .eq. "xu_yz" .or.     \
           vNam(varn) .eq. "ind_x_yz") then
         check = False
      else
         check = True
      end if  

      if (var .NE. "all") then
         check = isStrSubset( var,","+vNam(varn)+"," ) 
      end if

      if(check) then
      
         no_var=no_var+1

         if (xyc .EQ. 1) then
	    temp = f[:]->$vNam(varn)$
            if (vNam(varn) .eq. "zwwi" .or. vNam(varn) .eq. "zusi")
	       dummy=0
            else
               data_att = f_att->$vNam(varn)$
            end if
            if(vNam(varn) .eq. "lwps_xy" .or. vNam(varn) .eq. "pras_xy"     \
              .or. vNam(varn) .eq. "prrs_xy" .or. vNam(varn) .eq. "qsws_xy" \
              .or. vNam(varn) .eq. "shfs_xy" .or. vNam(varn) .eq. "ts_xy"   \
              .or. vNam(varn) .eq. "us_xy" .or. vNam(varn) .eq. "z0s_xy"    \
              .or. vNam(varn) .eq. "lwp*_xy" .or. vNam(varn) .eq. "pra*_xy" \
              .or. vNam(varn) .eq. "prr*_xy" .or. vNam(varn) .eq. "qsws*_xy"\
              .or. vNam(varn) .eq. "shf*_xy" .or. vNam(varn) .eq. "t*_xy"   \
              .or. vNam(varn) .eq. "u*_xy" .or. vNam(varn) .eq. "z0*_xy") then
              ;these variables depend on zu1_xy and that's why they have 
              ;only one z-layer
              data(varn,:,0,:,:)=temp(:,0,ys:ye,xs:xe) 
              no_zu1=no_zu1+1
            else
              if (vNam(varn) .eq. "zwwi" .or. vNam(varn) .eq. "zusi")
                  ;these variables depend on x and y
                  data(varn,0,0,:,:)=doubletofloat(temp(ys:ye,xs:xe))
                  no_topo=no_topo+1
               else
                  data(varn,:,:,:,:)=temp(:,zs:ze,ys:ye,xs:xe) 
               end if
            end if
            delete(temp)	    	
         end if
         if ( xzc .eq. 1 ) then
	    temp = f[:]->$vNam(varn)$
	    data_att = f_att->$vNam(varn)$
            data(varn,:,:,:,:)=temp(:,zs:ze,ys:ye,xs:xe) 
            delete(temp)
         end if
         if ( yzc .eq. 1) then
	    temp = f[:]->$vNam(varn)$
	    data_att = f_att->$vNam(varn)$
            data(varn,:,:,:,:)=temp(:,zs:ze,ys:ye,xs:xe) 
            delete(temp)
         end if

         data!0 = "var"
   	 data!1 = "t"
   	 data!2 = "z"
   	 data!3 = "y"
   	 data!4 = "x"  

         MinVal(varn) = min(data(varn,start_time_step:end_time_step,\
                                               0:(ze-zs),0:(ye-ys),0:(xe-xs)))
         MaxVal(varn) = max(data(varn,start_time_step:end_time_step,\
                                               0:(ze-zs),0:(ye-ys),0:(xe-xs)))
	 
         if (vNam(varn) .eq. "zwwi" .or. vNam(varn) .eq. "zusi")
	    unit(varn) = "meters"
          else
            unit(varn) = data_att@units
            delete(data_att)
          end if

      end if

   end do

   if (no_var .EQ. 0) then
      print(" ")
      print("The variables var='"+var+"' do not exist on your input file;")
      print("be sure to have one comma before and after each variable")
      print(" ")
      exit
   end if

   var_input=new(no_var,string)
   numb=0
   no_var=0
   
   do varn=dim-1,0,1   
           
      if ( vNam(varn) .eq. "zw" .or. vNam(varn) .eq. "yv" .or.          \
           vNam(varn) .eq. "xu" .or. vNam(varn) .eq. "zu" .or.          \
           vNam(varn) .eq. "time" .or. vNam(varn) .eq. "zu_3d" .or.     \
           vNam(varn) .eq. "zw_3d" .or. vNam(varn) .eq. "x".or.         \
           vNam(varn) .eq. "y" .or. vNam(varn) .eq. "zu_xy" .or.        \
           vNam(varn) .eq. "zw_xy" .or. vNam(varn) .eq. "zu1_xy" .or.   \
           vNam(varn) .eq. "ind_z_xy" .or. vNam(varn) .eq. "y_xz" .or.  \
           vNam(varn) .eq. "yv_xz" .or. vNam(varn) .eq. "ind_y_xz" .or. \
           vNam(varn) .eq. "x_yz" .or. vNam(varn) .eq. "xu_yz" .or.     \
           vNam(varn) .eq. "ind_x_yz") then
         check = False
      else
         check = True
      end if  

      if (var .NE. "all") then
         check = isStrSubset( var,","+vNam(varn)+"," ) 
      end if

      if(check) then     
         var_input(numb)=vNam(varn) 
	 numb=numb+1     
      end if
      
      if (check) then
         no_var = no_var+1
      end if
      
   end do

   if (no_var .EQ. 0) then
      print(" ")
      print("The variables var='"+var+"' do not exist on your input file;")
      print("be sure to have one comma before and after each variable")
      print(" ")
      exit
   end if

   if (vector .EQ. 1) then
      if (v1 .EQ. 0)then
         print(" ")
         print("Component 1 for the vector-plot ('vec1') must be one "+\
               "of the variables on the input file:")
	 print("- "+var_input)
	 print("be sure to have one comma before and after the variable")
         print(" ")
         exit
      end if

      if (v2 .EQ. 0)then
         print(" ")
         print("Component 2 for the vector-plot ('vec2') must be one "+\
               "of the variables on the input file:")
	 print("- "+var_input)
	 print("be sure to have one comma before and after the variable")
         print(" ")
         exit
      end if
   end if

   ; ***************************************************
   ; open workstation(s) 
   ; ***************************************************

   if ( format_out .EQ. "pdf" .OR. format_out .EQ. "ps" ) then
      format_out@wkPaperSize = "A4"
   end if
   if ( format_out .EQ. "png" ) then
      format_out@wkWidth  = 1000
      format_out@wkHeight = 1000
   end if

   wks_ps  = gsn_open_wks(format_out,file_out)
   gsn_define_colormap(wks_ps,"rainbow+white")

   if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
      plot=new((/no_time*no_layer/),graphic)
   else
      plot=new((/no_time*no_layer*(no_var-no_topo) + no_topo - \
                       no_time*(no_layer-1)*no_zu1/),graphic)
   end if
   dim_plot=dimsizes(plot)

   page = 0
   n=0
   print(" ")
   print("Plots sorted by '"+sort+"'")
   print(" ")

   ; ***************************************************
   ; create plots 
   ; ***************************************************


   if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
      do lo = los, loe           				
         do li = lis, lie
	    if (xyc .EQ. 1)then
	       if (sort .EQ. "time")then
                  if(z_d(lo) .eq. -1.d) then
                    level = "z-average"
                  else
                    level = "z=" + z_d(lo) + "m"
                  end if
	       else
                  if(z_d(li) .eq. -1.d) then
                    level = "z-average"
                  else
	            level = "z=" + z_d(li) + "m"
                  end if
	       end if
	    end if
	    if (xzc .EQ. 1)then
	       if (sort .EQ. "time")then
                 if(y_d(lo) .eq. -1.d) then
                   level = "y-average"
                 else
                   level = "y=" + y_d(lo) + "m" 
                 end if
	       else
                  if(y_d(li) .eq. -1.d) then
                    level = "y-average"
                  else
	            level = "y=" + y_d(li) + "m"
                  end if
	       end if
	    end if
	    if (yzc .EQ. 1)then
	       if (sort .EQ. "time")then
                 if(x_d(lo) .eq. -1.d) then
                    level = "x-average"
                 else
                    level = "x=" + x_d(lo) + "m" 
                 end if
	       else
                 if(x_d(li) .eq. -1.d) then
                    level = "x-average"
                 else
	            level = "x=" + x_d(li) + "m" 
                 end if
	       end if
	    end if	    	  
            vecres                  = True            ; vector only resources
            vecres@gsnDraw          = False           ; don't draw
  	    vecres@gsnFrame         = False           ; don't advance frame
            vecres@vcGlyphStyle     = "CurlyVector"   ; curly vectors
  	    vecres@vcRefMagnitudeF  = ref_mag         ; define vector ref mag
  	    vecres@vcRefLengthF     = 0.05            ; define length of 
                                                      ; vec ref
  	    vecres@gsnRightString   = "Vector plot of "+vec1+" and "+vec2
	    vecres@tmXBLabelFontHeightF   = font_size
            vecres@tmYLLabelFontHeightF   = font_size
            vecres@tiXAxisFontHeightF     = font_size
            vecres@tiYAxisFontHeightF     = font_size
	    if (sort .EQ. "time")then
               vecres@gsnLeftString = "t=" + \
                             decimalPlaces(t_all(li)/3600,2,True) +"h   "+level
	    else
	       vecres@gsnLeftString = "t=" + \
                             decimalPlaces(t_all(lo)/3600,2,True) +"h   "+level
	    end if
	    if (xyc .EQ. 1) then  
               vecres@tiXAxisString = "x (m)"
	       vecres@tiYAxisString = "y (m)"   
	       if (sort .EQ. "time")then
                  plot(n) = gsn_csm_vector(wks_ps,vect1(li,lo-los,:,:),\
                                               vect2(li,lo-los,:,:),vecres)
	       else
	          plot(n) = gsn_csm_vector(wks_ps,vect1(lo,li-lis,:,:),\
                                                vect2(lo,li-lis,:,:),vecres)
	       end if
	    end if
	    if (xzc .EQ. 1) then
               vecres@tiXAxisString = "x (m)"
	       vecres@tiYAxisString = "z (m)" 
	       if (sort .EQ. "time")then 
	          plot(n) = gsn_csm_vector(wks_ps,vect1(li,:,lo-los,:),\
                                               vect2(li,:,lo-los,:),vecres)
	       else
	          plot(n) = gsn_csm_vector(wks_ps,vect1(lo,:,li-lis,:),\
                                                vect2(lo,:,li-lis,:),vecres)
	       end if
	    end if
	    if (yzc .EQ. 1) then
               vecres@tiXAxisString = "y (m)"
	       vecres@tiYAxisString = "z (m)" 
	       if (sort .EQ. "time")then 
	          plot(n) = gsn_csm_vector(wks_ps,vect1(li,:,:,lo-los),\
                                                vect2(li,:,:,lo-los),vecres)
	       else
	          plot(n) = gsn_csm_vector(wks_ps,vect1(lo,:,:,li-lis),\
                                                vect2(lo,:,:,li-lis),vecres)
	       end if
	    end if
            n=n+1
         end do
      end do
   end if
 
   do varn=dim-1,0,1 

      if (vector .EQ. 1 ) then    
         check_vecp = isStrSubset( plotvec,","+vNam(varn)+",")
      end if 
     
      if ( vNam(varn) .eq. "zw" .or. vNam(varn) .eq. "yv" .or.          \
           vNam(varn) .eq. "xu" .or. vNam(varn) .eq. "zu" .or.          \
           vNam(varn) .eq. "time" .or. vNam(varn) .eq. "zu_3d" .or.     \
           vNam(varn) .eq. "zw_3d" .or. vNam(varn) .eq. "x".or.         \
           vNam(varn) .eq. "y" .or. vNam(varn) .eq. "zu_xy" .or.        \
           vNam(varn) .eq. "zw_xy" .or. vNam(varn) .eq. "zu1_xy" .or.   \
           vNam(varn) .eq. "ind_z_xy" .or. vNam(varn) .eq. "y_xz" .or.  \
           vNam(varn) .eq. "yv_xz" .or. vNam(varn) .eq. "ind_y_xz" .or. \
           vNam(varn) .eq. "x_yz" .or. vNam(varn) .eq. "xu_yz" .or.     \
           vNam(varn) .eq. "ind_x_yz") then
         check = False 
      else
         check = True
      end if   
  
      if (var .NE. "all") then
         check = isStrSubset( var,","+vNam(varn)+"," ) 
      end if 
    
      if(check) then

         space=(decimalPlaces(MaxVal(varn),3,True)-decimalPlaces(MinVal(varn),3,True))/24
  
         ;cs_res@cnMinLevelValF = decimalPlaces(MinVal(varn),3,True)
         ;cs_res@cnMaxLevelValF = decimalPlaces(MaxVal(varn),3,True)

         ;cs_res@cnLevelSpacingF = space
         cs_res@cnMaxLevelCount = 26
      
         ; ****************************************************
         ; loops over time and layer
         ; **************************************************** 

         
         do lo = los, loe           				
            do li = lis, lie

               ; ****************************************************
               ; xy cross section
               ; ****************************************************

               if ( xyc .eq. 1 ) then 
          
                  cs_res@gsnLeftString = "Plot of "+vNam(varn)
		  cs_res@gsnRightString = unit(varn)
                  
                  if ( sort .eq. "time" ) then
                     if ( z_d(lo) .eq. -1)then
		        if (delta_z .EQ. -1) then
                           level = "-average"
			else
			   level = "=" + z_d(lo) + "m"
			end if
                     else 
                        level = "=" + z_d(lo) + "m"
                     end if

                     if(vNam(varn) .eq. "zwwi"  .or. \
                        vNam(varn) .eq. "zusi") then
                         loe = 0
                         los = 0
                         lie  = 0
                         lis = 0
                         level = ""
                     end if

                     if(vNam(varn) .eq. "lwps_xy"  .or. \
                        vNam(varn) .eq. "pras_xy"  .or. \
                        vNam(varn) .eq. "prrs_xy"  .or. \
                        vNam(varn) .eq. "qsws_xy"  .or. \
                        vNam(varn) .eq. "shfs_xy"  .or. \
                        vNam(varn) .eq. "ts_xy"    .or. \
                        vNam(varn) .eq. "us_xy"    .or. \
                        vNam(varn) .eq. "z0s_xy"   .or. \
                        vNam(varn) .eq. "lwp*_xy"  .or. \
                        vNam(varn) .eq. "pra*_xy"  .or. \
                        vNam(varn) .eq. "prr*_xy"  .or. \
                        vNam(varn) .eq. "qsws*_xy" .or. \
                        vNam(varn) .eq. "shf*_xy"  .or. \
                        vNam(varn) .eq. "t*_xy"    .or. \
                        vNam(varn) .eq. "u*_xy"    .or. \
                        vNam(varn) .eq. "z0*_xy") then
                         loe = 0
                         level = "=" + zu1(0) + "m"
                      else
                         lis = start_time_step
                         lie = end_time_step
                         los = lays
                         loe = laye
                     end if                 

                     if(vNam(varn) .eq. "zwwi"  .or. \
                        vNam(varn) .eq. "zusi") then
                         cs_res@gsnCenterString = ""
                     else
                         cs_res@gsnCenterString = "t=" + \
                          decimalPlaces(t_all(li)/3600,2,True) +"h  z"+level
                     end if

                     if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
                        ;nothing to be done
                     else
                         plot(n) = gsn_csm_contour(wks_ps,\
                                           data(varn,li,lo-los,:,:),cs_res) 
                     end if


                     if (vector .EQ. 1 .AND. check_vecp) then
                        vecres                 = True           ; vector only
                                                                ; resources
  			vecres@gsnDraw         = False          ; don't draw
  	                vecres@gsnFrame        = False          ; don't 
                                                                ; advance frame
  	                vecres@vcGlyphStyle    = "CurlyVector"  ; curly 
                                                                ; vectors
  		        vecres@vcRefMagnitudeF = ref_mag        ; define 
                                                                ; vector ref 
                                                                ; mag
  		        vecres@vcRefLengthF    = 0.05           ; define 
                                                                ; length of 
                                                                ; vec ref  
                        vecres@gsnRightString  = " "
                        vecres@gsnLeftString   = " " 
  			vecres@tiXAxisString   = " "   
                        plot_vec=gsn_csm_vector(wks_ps,\
                              vect1(li,lo-los,:,:),vect2(li,lo-los,:,:),vecres)
                        overlay(plot(n), plot_vec)
                     end if                         
                  end if
                 
                  if ( sort .eq. "layer" ) then
		     if (z_d(li) .eq. -1) then
		        if (delta_z .EQ. -1) then
                           level = "-average"
			else
			   level = "=" + z_d(li) + "m"
			end if
                     else 
                        level = "=" + z_d(li) + "m"
                     end if 

                     if(vNam(varn) .eq. "zwwi"  .or. \
                        vNam(varn) .eq. "zusi") then
                         lie = 0
                         lis = 0
                         loe = 0
                         los = 0
                         level = ""
                      end if

                     if(vNam(varn) .eq. "lwps_xy"  .or. \
                        vNam(varn) .eq. "pras_xy"  .or. \
                        vNam(varn) .eq. "prrs_xy"  .or. \
                        vNam(varn) .eq. "qsws_xy"  .or. \
                        vNam(varn) .eq. "shfs_xy"  .or. \
                        vNam(varn) .eq. "ts_xy"    .or. \
                        vNam(varn) .eq. "us_xy"    .or. \
                        vNam(varn) .eq. "z0s_xy"   .or. \
                        vNam(varn) .eq. "lwp*_xy"  .or. \
                        vNam(varn) .eq. "pra*_xy"  .or. \
                        vNam(varn) .eq. "prr*_xy"  .or. \
                        vNam(varn) .eq. "qsws*_xy" .or. \
                        vNam(varn) .eq. "shf*_xy"  .or. \
                        vNam(varn) .eq. "t*_xy"    .or. \
                        vNam(varn) .eq. "u*_xy"    .or. \
                        vNam(varn) .eq. "z0*_xy") then
                         lie = 0
                         level = "=" + zu1(0) + "m"
                     else
                         lis = lays
                         lie = laye
                         los = start_time_step 
                         loe = end_time_step
                     end if

                     if(vNam(varn) .eq. "zwwi"  .or. \
                        vNam(varn) .eq. "zusi") then
                         cs_res@gsnCenterString = ""
                     else
                         cs_res@gsnCenterString = "t=" + \
                          decimalPlaces(t_all(lo)/3600,2,True) +"h  z"+level
                     end if

                     ;if (topo .EQ. 0) then
                        ;nothing to be done
                     ;else
                     ;   print("'topo' is set to " + topo)

                        ;plot(n) = gsn_csm_contour(wks_ps,\
                                            ;data(varn,lo,li-lis,:,:),cs_res)
                        ;if(vNam(varn) .EQ. "zwwi" .OR. \
                         ;  vNam(varn) .EQ. "zusi") then
                          ; to_res@cnFillPalette = "GMT_gray"
                           ;plot_topo = gsn_csm_contour(wks_ps,\
                            ;  data(varn,lo,li-lis,:,:),to_res)
                         ;end if
                        ;overlay(plot(n), plot_topo) 
                     ;end if

                     if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
                        ;nothing to be done
                     else
                        plot(n) = gsn_csm_contour(wks_ps,\
                                            data(varn,lo,li-lis,:,:),cs_res)
                     end if

                     if (vector .EQ. 1 .AND. check_vecp) then
                        vecres                 = True           ; vector only
                                                                ; resources
  			vecres@gsnDraw         = False          ; don't draw
  			vecres@gsnFrame        = False          ; don't 
                                                                ; advance frame
  			vecres@vcGlyphStyle    = "CurlyVector"  ; curly vectors
  			vecres@vcRefMagnitudeF = ref_mag        ; define 
                                                                ; vector ref 
                                                                ; mag
  			vecres@vcRefLengthF    = 0.05           ; define 
                                                                ; length of 
                                                                ; vec ref
  			vecres@gsnRightString  = " "            ; turn off 
                                                                ; right string
  			vecres@gsnLeftString   = " "            ; turn off 
                                                                ; left string
  			vecres@tiXAxisString   = " "   
                        plot_vec=gsn_csm_vector(wks_ps,vect1(lo,li-lis,:,:),\
                                                   vect2(lo,li-lis,:,:),vecres)
                        overlay(plot(n), plot_vec)
                     end if
                  end if
               end if

               ; ****************************************************
               ; xz cross section
               ; ****************************************************

               if ( xzc .eq. 1 ) then
                  
		  cs_res@gsnLeftString = "Plot of "+vNam(varn)
               
                  if ( sort .eq. "time" ) then
                     if ( y_d(lo) .eq. -1 ) then
                        level = "-average"
                     else 
                        level = "=" + y_d(lo) + "m"
                     end if
                     cs_res@gsnCenterString = "t=" + \
                           decimalPlaces(t_all(li)/3600,2,True) + "h  y"+ level

                     if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
                        ;nothing to be done
                     else
                         plot(n) = gsn_csm_contour(wks_ps,\
                                              data(varn,li,:,lo-los,:),cs_res)
                     end if

                     if (vector .EQ. 1 .AND. check_vecp) then
                        vecres                 = True           ; vector only 
                                                                ; resources
  			vecres@gsnDraw         = False          ; don't draw
  			vecres@gsnFrame        = False          ; don't 
                                                                ; advance frame
  			vecres@vcGlyphStyle    = "CurlyVector"  ; curly vectors
  			vecres@vcRefMagnitudeF = ref_mag        ; define 
                                                                ; vector ref 
                                                                ; mag
  			vecres@vcRefLengthF    = 0.05           ; define 
                                                                ; length of 
                                                                ; vec ref
  			vecres@gsnRightString  = " "            ; turn off 
                                                                ; right string
  			vecres@gsnLeftString   = " "            ; turn off 
                                                                ; left string
  			vecres@tiXAxisString   = " "   
                        plot_vec=gsn_csm_vector(wks_ps,vect1(li,:,lo-los,:),\
                                                   vect2(li,:,lo-los,:),vecres)
                        overlay(plot(n), plot_vec)
                     end if
                  end if

                  if ( sort .eq. "layer" ) then
                     if ( y_d(li) .eq. -1 ) then
                        level = "-average"
                     else 
                        level = "=" + y_d(li) + "m"
                     end if
                     cs_res@gsnCenterString = "t=" + \
                           decimalPlaces(t_all(lo)/3600,2,True) + "h  y"+ level

                     if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
                        ;nothing to be done
                     else
                         plot(n) = gsn_csm_contour(wks_ps,\
                                              data(varn,lo,:,li-lis,:),cs_res)  
                     end if

                     if (vector .EQ. 1 .AND. check_vecp) then
                        vecres                 = True           ; vector only 
                                                                ; resources
  			vecres@gsnDraw         = False          ; don't draw
  			vecres@gsnFrame        = False          ; don't 
                                                                ; advance frame
  			vecres@vcGlyphStyle    = "CurlyVector"  ; curly vectors
  			vecres@vcRefMagnitudeF = ref_mag        ; define 
                                                                ; vector ref 
                                                                ; mag
  			vecres@vcRefLengthF    = 0.05           ; define 
                                                                ; length of 
                                                                ; vec ref
  			vecres@gsnRightString  = " "            ; turn off 
                                                                ; right string
  			vecres@gsnLeftString   = " "            ; turn off 
                                                                ; left string
  			vecres@tiXAxisString   = " "   
                        plot_vec=gsn_csm_vector(wks_ps,vect1(lo,:,li-lis,:),\
                                                   vect2(lo,:,li-lis,:),vecres)
                        overlay(plot(n), plot_vec)
                     end if
                  end if                  
               end if

               ; ****************************************************
               ; yz cross section
               ; ****************************************************

               if ( yzc .eq. 1 ) then
                                  
		  cs_res@gsnLeftString = "Plot of "+vNam(varn)

                  if ( sort .eq. "time" ) then
                     if ( x_d(lo) .eq. -1 ) then
                        level = "-average"
                     else 
                        level = "=" + x_d(lo) + "m"
                     end if
                     cs_res@gsnCenterString = "t=" + \
                          decimalPlaces(t_all(li)/3600,2,True) + "h  x"+ level

                     if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
                        ;nothing to be done
                     else
                        plot(n) = gsn_csm_contour(wks_ps,\
                                          data(varn,li,:,:,lo-los),cs_res)
                     end if

    	             if (vector .EQ. 1 .AND. check_vecp) then
                        vecres                 = True           ; vector only 
                                                                ; resources
  			vecres@gsnDraw         = False          ; don't draw
  			vecres@gsnFrame        = False          ; don't 
                                                                ; advance frame
  			vecres@vcGlyphStyle    = "CurlyVector"  ; curly vectors
  			vecres@vcRefMagnitudeF = ref_mag        ; define 
                                                                ; vector ref 
                                                                ; mag
  			vecres@vcRefLengthF    = 0.05           ; define 
                                                                ; length of 
                                                                ; vec ref
  			vecres@gsnRightString  = " "            ; turn off 
                                                                ; right string
  			vecres@gsnLeftString   = " "            ; turn off 
                                                                ; left string
  			vecres@tiXAxisString   = " "   
                        plot_vec=gsn_csm_vector(wks_ps,vect1(li,:,:,lo-los),\
                                                   vect2(li,:,:,lo-los),vecres)
                        overlay(plot(n), plot_vec)
                     end if
                  end if

                  if ( sort .eq. "layer" ) then
                     if ( x_d(li) .eq. -1 ) then
                        level = "-average"
                     else 
                        level = "=" + x_d(li) + "m"
                     end if
                     cs_res@gsnCenterString = "t=" + \
                           decimalPlaces(t_all(lo)/3600,2,True) + "h  x"+ level

                     if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then
                        ;nothing to be done
                     else
                        plot(n) = gsn_csm_contour(wks_ps,\
                                          data(varn,lo,:,:,li-lis),cs_res)
                     end if

                     if (vector .EQ. 1 .AND. check_vecp)then
                        vecres                 = True           ; vector only 
                                                                ;resources
  			vecres@gsnDraw         = False          ; don't draw
  			vecres@gsnFrame        = False          ; don't 
                                                                ; advance frame
  			vecres@vcGlyphStyle    = "CurlyVector"  ; curly vectors
  			vecres@vcRefMagnitudeF = ref_mag        ; define 
                                                                ; vector ref 
                                                                ; mag
  			vecres@vcRefLengthF    = 0.05           ; define 
                                                                ; length of 
                                                                ; vec ref
  			vecres@gsnRightString  = " "            ; turn off 
                                                                ; right string
  			vecres@gsnLeftString   = " "            ; turn off 
                                                                ; left string
  			vecres@tiXAxisString   = " "   
                        plot_vec=gsn_csm_vector(wks_ps,vect1(lo,:,:,li-lis),\
                                                   vect2(lo,:,:,li-lis),vecres)
                        overlay(plot(n), plot_vec)
                     end if 
                  end if
               end if	      
               n=n+1  
            end do 	
         end do  
      end if
   end do

   ; ***************************************************
   ; merge plots onto one page 
   ; ***************************************************   

   no_frames = 0

   if (vector .EQ. 1 .AND. plotvec .EQ. "plotvec") then 
      if ((format_out .EQ. "eps" .OR. format_out .EQ. "epsi") .AND. \
           no_time*no_layer*(no_var+1) .gt. no_rows*no_columns) then
         gsn_panel(wks_ps,plot(0:(no_time*no_layer*(no_var+1))-1),\
                                       (/no_var+1,no_layer*no_time/),cs_resP)
         print(" ")
         print("Outputs to .eps or .epsi have only one frame")
         print(" ")
      else
         do np = 0,dim_plot-1,no_rows*no_columns   
            if ( np + no_rows*no_columns .gt. dim_plot-1) then 
               gsn_panel(wks_ps, plot(np:dim_plot-1),(/no_rows,no_columns/),\
                                                                      cs_resP)
               no_frames = no_frames + 1
            else
               gsn_panel(wks_ps, plot(np:np+no_rows*no_columns-1),\
                                              (/no_rows,no_columns/),cs_resP)
               no_frames = no_frames + 1
            end if
         end do
      end if
   else       
      if ((format_out .EQ. "eps" .OR. format_out .EQ. "epsi") \
                               .AND. dim_plot .gt. no_rows*no_columns) then 
         gsn_panel(wks_ps,plot(0:dim_plot-1),(/dim_plot,1/),cs_resP)
         print(" ")
         print("Outputs to .eps or .epsi have only one frame")
         print(" ")
      else
         do np = 0,dim_plot-1,no_rows*no_columns  
            if ( np + no_rows*no_columns .gt. dim_plot-1) then
               gsn_panel(wks_ps, plot(np:dim_plot-1),(/no_rows,no_columns/),\
                                                                     cs_resP)
               no_frames = no_frames + 1
            else
               gsn_panel(wks_ps, plot(np:np+no_rows*no_columns-1),\
                                              (/no_rows,no_columns/),cs_resP)
               no_frames = no_frames + 1
            end if
         end do
      end if
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
