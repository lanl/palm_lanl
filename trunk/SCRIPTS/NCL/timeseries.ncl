load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;***************************************************
; Checking the kind of the script
;***************************************************

function which_script()
local script 
begin 
   script = "timeseries"
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
   
   if (norm_t .EQ. 0) then 
      print(" ")
      print("Normalising with 0 is not allowed, 'norm_t' is set to 1.0")
      print(" ")
      norm_t = 1.0
   end if 
 

   ;***************************************************
   ; open input file
   ;***************************************************
   
   file_in_1 = False
   if (isStrSubset(file_in, ".nc"))then
      start_f = -2
      end_f = -2 
      file_in_1 = True     
   end if

   if (start_f .EQ. -1)then
      print(" ")
      print("'start_f' must be one of the cyclic numbers (at least 0) of "+\
            "your input file(s)")
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

   vNam  =getfilevarnames(f_att)
   if( ncl_version .GE. 6.1 ) then
       vNam = vNam(::-1)
   end if
   vType = getfilevartypes(f_att,vNam) 
   vNam_static = vNam

   if ((all(vType .eq. "double"))) then ;distinction if data is double or float
      check_vType = True
   else
      check_vType = False
   end if

   print(" ")
   print("Variables in input file:")
   print("- "+ vNam)
   print(" ")
   dim = dimsizes(vNam)
   if (dim .EQ. 0) then
      print(" ")
      print("There is no data on file")
      print(" ")
   end if

   t_all   = f[:]->time
   nt      = dimsizes(t_all)

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

   if (nt .LE. 1) then
      print(" ")
      print("Input file contains only one time step -> " +\
            "plot of timeseries is not possible ")
      print(" ")
      exit
   end if

   ;****************************************************	
   ; start of time step and different types of mistakes that could be done
   ;****************************************************

   if (start_time_step .EQ. -1.) then   	
      start_time_step=t_all(0)/3600      
   else 
      if (start_time_step .GE. t_all(nt-1)/3600)
         print(" ")
         print("'start_time_step' = "+ start_time_step +\
               "h is equal or greater than last time step = " + \
               t_all(nt-1)+"s = "+t_all(nt-1)/3600+"h")
         print(" ")
         print("Select another 'start_time_step'")
         print(" ")
         exit
      end if
      if (start_time_step .LT. t_all(0)/3600)
         print(" ")
         print("'start_time_step' = "+ start_time_step +\
               "h is lower than first time step = " + t_all(0)+"s = "+\
               t_all(0)/3600+"h")
         print(" ")
         print("Select another 'start_time_step'")
         print(" ")
         exit
      end if
   end if
   do i=0,nt-2     
      if (start_time_step .GE. (t_all(i)-delta_t/2)/3600 .AND. \
          start_time_step .LT. (t_all(i)+delta_t/2)/3600)then
         st=i
         break
       else
         st=0
      end if
   end do
   if (start_time_step .GE. t_all(nt-1)-delta_t/2 .AND. \
       start_time_step .LT. t_all(nt-1)) then
      st=nt-2   
    else
      st=0
   end if
      
   ;****************************************************
   ; end of time step and different types of mistakes that could be done
   ;****************************************************

   if (end_time_step .EQ. -1.) then     	
      end_time_step = t_all(nt-1)/3600  
   else
      if (end_time_step .LE. t_all(0)/3600)
         print(" ")
         print("'end_time_step' = "+end_time_step+ \
              "h is lower or equal than first time step = " + \
               t_all(0)+"s = "+t_all(0)/3600+"h")
         print(" ")
         print("Select another 'end_time_step'")
         print(" ")
         exit
      end if
      if (end_time_step .GT. t_all(nt-1)/3600)
         print(" ")
         print("'end_time_step' = "+ end_time_step +\
               "h is greater than last time step = " + t_all(nt-1)+"s = "+\
               t_all(nt-1)/3600+"h")
         print(" ")
         print("Select another 'end_time_step'")  
         print(" ")
         exit
      end if
      if (end_time_step .LE. start_time_step/3600)
         print(" ")
         print("'end_time_step' = "+ end_time_step +\
               "h is equal or lower than 'start_time_step' = "+\
               start_time_step+"h")
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
   print("                     till "+t_all(end_time_step)/3600+" h = "+\
         t_all(end_time_step)+" s => index = "+end_time_step)
   print(" ")

   t = t_all(start_time_step:end_time_step)/norm_t
  
   ; ***************************************************
   ; set up recourses
   ; ***************************************************

   res                         = True
   res@gsnDraw                 = False
   res@gsnFrame                = False
   res@tmXBMode                = True
   res@tmYLMode                = True
   res@txFont                  = "helvetica"
   res@tiMainFont              = "helvetica"
   res@tiXAxisFont             = "helvetica"
   res@tiYAxisFont             = "helvetica"
   res@tmXBLabelFont           = "helvetica"
   res@tmYLLabelFont           = "helvetica" 
   res@xyLineColors            = (/237/) 
   res@trXMaxF                 = t(start_time_step)
   res@trXMinF                 = t(end_time_step)
   
   res@lgLabelFontHeightF     = 0.02

   resP                            = True
   resP@gsnMaximize                = True
   resP@gsnPanelXWhiteSpacePercent = 4.0
   resP@gsnPanelYWhiteSpacePercent = 4.0
   resP@txFont                     = "helvetica"
   resP@txString                   = f_att@title+" time series "
   resP@txFuncCode                 = "~" 
   if (format_out .NE. "x11") then
      if (no_rows .GE. 5) then
         resP@txFontHeightF        = 0.0105
      else
         resP@txFontHeightF        = 0.015
      end if
   else
      resP@txFontHeightF           = 0.015
   end if
   res@tmXBMinorPerMajor = 4
   res@tmYLMinorPerMajor = 4

   res@vpWidthF=4

   txres = True

   ; ***************************************************
   ; read data and create plots
   ; ***************************************************
   
   if ( format_out .EQ. "pdf" .OR. format_out .EQ. "ps" ) then
      format_out@wkPaperSize = "A4"
   end if
   if ( format_out .EQ. "png" ) then
      format_out@wkWidth  = 1000
      format_out@wkHeight = 1000
   end if

   wks_ps = gsn_open_wks(format_out,file_out)		
   gsn_define_colormap(wks_ps,"rainbow+white") 
   plot_ps=new(dim,graphic) 
                           	 
   n=0 
   if (check_vType) then
      minE  =  100000.d
      maxE  = -100000.d
      minus =  100000.d
      maxus = -100000.d
      minu  =  100000.d
      maxu  = -100000.d
      minz  =  100000.d
      maxz  = -100000.d
      minw  =  100000.d
      maxw  = -100000.d
      minp  =  100000.d
      maxp  = -100000.d
      mins  =  100000.d
      maxs  = -100000.d
   else
      minE  =  100000.
      maxE  = -100000.
      minus =  100000.
      maxus = -100000.
      minu  =  100000.
      maxu  = -100000.
      minz  =  100000.
      maxz  = -100000.
      minw  =  100000.
      maxw  = -100000.
      minp  =  100000.
      maxp  = -100000.
      mins  =  100000.
      maxs  = -100000.
   end if

   if (check_vType) then
      data   = new((/dim,(end_time_step-start_time_step)+1/),double)
      data_0 = new((end_time_step-start_time_step)+1,double) 
      mini   = new(dim,double)
      maxi   = new(dim,double)
   else
      data   = new((/dim,(end_time_step-start_time_step)+1/),float) 
      data_0 = new((end_time_step-start_time_step)+1,float) 
      mini   = new(dim,float)
      maxi   = new(dim,float)
   end if 
   unit   = new(dim,string)   
   data_0 = 0.0
  
   
   if (over .EQ. 1) then 
      plot_E       = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_Es      = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_us      = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_ws      = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_umax    = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_vmax    = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_wmax    = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_z_i_wpt = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_z_i_pt  = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_wpptp0  = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_wpptp   = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_wpt     = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_pt_0_   = gsn_csm_xy(wks_ps,t,data_0(:),res)
      plot_pt_zp_  = gsn_csm_xy(wks_ps,t,data_0(:),res)
   end if


   if (var .NE. "all") then
      vNam_temp = new((/dim/),string)
      comma     = 0
      var_char  = stringtocharacter(var)
      no_char   = dimsizes(var_char)
       
      do j = 0, no_char-1
         if(var_char(j) .eq. ",")
            comma = comma + 1
         end if
      end do

      if(comma .eq. 0 .or. comma .eq. 1)
          print(" ")
          print("The variables 'var="+var+"'" )
          print("do not exist on your input file;")
          print("be sure to have one comma before and after each variable")
          print(" ")
          exit
      end if

      indices = new((/comma/), integer)
      comma   = 0
        
      do j = 0, no_char-1
         if(var_char(j) .eq. ",")
           indices(comma) = j
           comma = comma + 1
         end if
      end do

      do j = 0, comma -2
         vNam_temp(j) = charactertostring(\
                                   var_char(indices(j)+1:indices(j+1)-1))
      end do
     
      delete(vNam)
      vNam = new((/comma-1/),string)
      do varn = 0, comma - 2
         vNam(varn) = vNam_temp(comma - 2 - varn)
      end do
     
      do j = 0, comma-2
        count_check = 0
        do varn = 0, dim-1 
           if(vNam(j) .ne. vNam_static(varn)) 
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
      end do

      dim = comma-1
    end if


   do varn = dim-1, 0, 1 

      if( isStrSubset (vNam(varn), "time") )
         check = False
      else
         check = True
      end if
     
      if(check) then

         data_all = f[:]->$vNam(varn)$ 
	 data_att = f_att->$vNam(varn)$
         unit(varn) = data_att@units
  
      
         data(varn,:)=data_all(start_time_step:end_time_step)
         
         if (over .EQ. 1) then

            mini(varn) = min(data(varn,:))
            maxi(varn) = max(data(varn,:))
           
            if (vNam(varn) .EQ. "E" .OR. vNam(varn) .EQ. "Es" .OR. \
                vNam(varn) .EQ. "E*") then
               if (mini(varn) .EQ. maxi(varn)) then
                  if (min(data(varn,:)) .EQ. 0)then
                     mini(varn)= mini(varn)-0.1
                     maxi(varn)= maxi(varn)+0.1
                  end if
                  if (min(data(varn,:)) .LT. 0)then
                     mini(varn)= mini(varn)-1.+(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.-(maxi(varn))/2
                  end if
                  if (min(data(varn,:)) .GT. 0)then
                     mini(varn)= mini(varn)-1.-(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.+(maxi(varn))/2
                  end if
               end if
               minE=min((/minE,mini(varn)/))
               maxE=max((/maxE,maxi(varn)/))
            end if

            if (vNam(varn) .EQ. "us" .OR. vNam(varn) .EQ. "ws"\
                .OR. vNam(varn) .EQ. "u*" .OR. vNam(varn) .EQ. "w*") then
               if (mini(varn) .EQ. maxi(varn)) then
                  if (min(data(varn,:)) .EQ. 0)then
                     mini(varn)= mini(varn)-0.1
                     maxi(varn)= maxi(varn)+0.1
                  end if
                  if (min(data(varn,:)) .LT. 0)then
                     mini(varn)= mini(varn)-1.+(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.-(maxi(varn))/2
                  end if
                  if (min(data(varn,:)) .GT. 0)then
                     mini(varn)= mini(varn)-1.-(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.+(maxi(varn))/2
                  end if
               end if
               minus=min((/minus,mini(varn)/))
               maxus=max((/maxus,maxi(varn)/))
            end if

            if (vNam(varn) .EQ. "umax" .OR. vNam(varn) .EQ. "vmax" .OR. \
                vNam(varn) .EQ. "wmax") then
               if (mini(varn) .EQ. maxi(varn)) then
                  if (mini(varn) .EQ. 0)then
                     mini(varn)= mini(varn)-0.1
                     maxi(varn)= maxi(varn)+0.1
                  end if
                  if (mini(varn) .LT. 0)then
                     mini(varn)= mini(varn)-1.+(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.-(maxi(varn))/2
                  end if
                  if (mini(varn) .GT. 0)then
                     mini(varn)= mini(varn)-1.-(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.+(maxi(varn))/2
                  end if
               end if
               minu=min((/minu,mini(varn)/))
               maxu=max((/maxu,maxi(varn)/))
            end if

            if (vNam(varn) .EQ. "z_i_wpt" .OR. vNam(varn) .EQ. "z_i_pt") then
               if (mini(varn) .EQ. maxi(varn)) then
                  if (min(data(varn,:)) .EQ. 0)then
                     mini(varn)= mini(varn)-0.1
                     maxi(varn)= maxi(varn)+0.1
                  end if
                  if (min(data(varn,:)) .LT. 0)then
                     mini(varn)= mini(varn)-1.+(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.-(maxi(varn))/2
                  end if
                  if (min(data(varn,:)) .GT. 0)then
                     mini(varn)= mini(varn)-1.-(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.+(maxi(varn))/2
                  end if
               end if
               minz=min((/minz,mini(varn)/))
               maxz=max((/maxz,maxi(varn)/))
            end if

            if (vNam(varn) .EQ. "wpptp0" .OR. vNam(varn) .EQ. "wpptp" .OR. \
                vNam(varn) .EQ. "wpt" .OR. vNam(varn) .EQ. \
                "w"+dq+"pt"+dq+"0" .OR. vNam(varn) .EQ. "w"+dq+"pt"+dq ) then
               if (mini(varn) .EQ. maxi(varn)) then
                  if (min(data(varn,:)) .EQ. 0)then
                     mini(varn)= mini(varn)-0.1
                     maxi(varn)= maxi(varn)+0.1
                  end if
                  if (min(data(varn,:)) .LT. 0)then
                     mini(varn)= mini(varn)-1.+(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.-(maxi(varn))/2
                  end if
                  if (min(data(varn,:)) .GT. 0)then
                     mini(varn)= mini(varn)-1.-(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.+(maxi(varn))/2
                  end if
               end if
               minw=min((/minw,mini(varn)/))
               maxw=max((/maxw,maxi(varn)/))
            end if

            if (vNam(varn) .EQ. "pt_0_" .OR. vNam(varn) .EQ. "pt_zp_"\
               .OR. vNam(varn) .EQ. "pt(0)" .OR. vNam(varn) .EQ. "pt(zp)") then
               if (mini(varn) .EQ. maxi(varn)) then
                  if (min(data(varn,:)) .EQ. 0)then
                     mini(varn)= mini(varn)-0.1
                     maxi(varn)= maxi(varn)+0.1
                  end if
                  if (min(data(varn,:)) .LT. 0)then
                     mini(varn)= mini(varn)-1.+(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.-(maxi(varn))/2
                  end if
                  if (min(data(varn,:)) .GT. 0)then
                     mini(varn)= mini(varn)-1.-(mini(varn))/2
                     maxi(varn)= maxi(varn)+1.+(maxi(varn))/2
                  end if
               end if
               minp=min((/minp,mini(varn)/))
               maxp=max((/maxp,maxi(varn)/))
            end if
            
         end if
      end if
   end do
 
  if (isStrSubset(data@long_name," SR " ))then
      print(" ")
      print("If you have outputs of statistic regions you cannot overlay "+\
            "variables;")
      print("'over' is set to 0" )
      print(" ")
      over = 0
   end if

   do varn = dim-1, 0, 1 
      
      if( isStrSubset (vNam(varn), "time") )
         check = False
      else
         check = True
      end if     
    
      if (var .NE. "all") then
         check = isStrSubset( var,","+vNam(varn)+"," ) 
      end if 
  
      if(check) then
       
        if (over .EQ. 1) then
 
            res@gsnLeftString       = "overlayed plot"
            res@gsnRightString      = unit(varn)
	    if (norm_t .NE. 1.)then
               res@tiXAxisString        = "t ("+unit_t+")"
	    else
	       res@tiXAxisString        = "t (s)"
	    end if
	    res@tiYAxisString = " "
            res@tiXAxisFontHeightF   = font_size
            res@txFontHeightF        = font_size
            res@tiYAxisFontHeightF   = font_size            

            if (vNam(varn) .EQ. "E")
               E=0
               res@xyLineColors     = (/237/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 237
               res@trYMaxF           = maxE
               res@trYMinF           = minE
               plot_E = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if
            if (vNam(varn) .EQ. "Es" .OR. vNam(varn) .EQ. "E*")
               Es=0
               res@xyLineColors            = (/144/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 144
               plot_Es = gsn_csm_xy(wks_ps,t,data(varn,:),res) 
            end if 
             
            if (vNam(varn) .EQ. "us" .OR. vNam(varn) .EQ. "u*")
               us=0
               res@xyLineColors            = (/237/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 237
               res@trYMaxF           = maxus
               res@trYMinF           = minus
               plot_us = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if             
            if (vNam(varn) .EQ. "ws" .OR. vNam(varn) .EQ. "w*")
               ws=0
               res@xyLineColors            = (/144/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 144
               plot_ws = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if
             
            if (vNam(varn) .EQ. "umax")
               u=0
               res@xyLineColors            = (/237/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 237
               res@trYMaxF           = maxu
               res@trYMinF           = minu
               plot_umax = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if
            if (vNam(varn) .EQ. "vmax")
               v=0
               res@xyLineColors            = (/144/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 144
               plot_vmax = gsn_csm_xy(wks_ps,t,data(varn,:),res) 
            end if
            if (vNam(varn) .EQ. "wmax")
               w=0
               res@xyLineColors            = (/80/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 80
               plot_wmax = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if 
   
            if (vNam(varn) .EQ. "z_i_wpt")
               zw=0
               res@xyLineColors            = (/237/) 
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 237
               res@trYMaxF           = maxz
               res@trYMinF           = minz 
               plot_z_i_wpt = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if
            if (vNam(varn) .EQ. "z_i_pt")
               z=0
               res@xyLineColors            = (/144/)  
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 144 
               plot_z_i_pt = gsn_csm_xy(wks_ps,t,data(varn,:),res) 
            end if

            if (vNam(varn) .EQ. "wpptp0" .OR. \
                vNam(varn) .EQ. "w"+dq+"pt"+dq+"0" )
               w0=0
               res@xyLineColors            = (/237/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 237
               res@trYMaxF           = maxw
               res@trYMinF           = minw
               plot_wpptp0 = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if
            if (vNam(varn) .EQ. "wpptp" .OR. vNam(varn) .EQ. "w"+dq+"pt"+dq)
               wp=0
               res@xyLineColors            = (/144/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 144
               plot_wpptp = gsn_csm_xy(wks_ps,t,data(varn,:),res)  
            end if
            if (vNam(varn) .EQ. "wpt")
               wt=0
               res@xyLineColors            = (/80/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 80
               plot_wpt = gsn_csm_xy(wks_ps,t,data(varn,:),res)
            end if

            if (vNam(varn) .EQ. "pt_0_" .OR. vNam(varn) .EQ. "pt(0)")
               p=0
               res@xyLineColors            = (/237/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 237
               res@trYMaxF           = maxp
               res@trYMinF           = minp
               plot_pt_0_ = gsn_csm_xy(wks_ps,t,data(varn,:),res)  
            end if
            if (vNam(varn) .EQ. "pt_zp_" .OR. vNam(varn) .EQ. "pt(zp)")
               pz=0
               res@xyLineColors            = (/144/)
               res@xyLineLabelFontHeightF = 0.05 
               res@xyLineLabelFontColor   = 144
               plot_pt_zp_ = gsn_csm_xy(wks_ps,t,data(varn,:),res) 
            end if

         end if
      end if
   end do


   do varn = dim-1, 0, 1 
   
      if( isStrSubset (vNam(varn), "time") )
         check = False
      else 
         check = True
      end if    
  
      if (var.NE. "all") then
         check = isStrSubset( var,","+vNam(varn)+"," ) 
      end if 
    
      if(check) then
        
         if (over .EQ. 1) then       
   	
            if (vNam(varn) .EQ. "E" .AND. Es .NE. 1) then 
               E=1   
               overlay(plot_E,plot_Es)
               n=n+1
               plot_ps(n) = plot_E    

               ; ***************************************************
               ; legend for overlaid plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1            
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"E","E*"/),lgres)       

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)            
            end if
            if ((vNam(varn) .EQ. "Es" .OR. vNam(varn) .EQ. "E*") .AND. \
                E .NE. 1) then
               Es=1
               overlay(plot_E,plot_Es)
               n=n+1
               plot_ps(n) = plot_E  

               ; ***************************************************
               ; legend for overlaid plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1            
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"E","E*"/),lgres)       

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)         
            end if

            if ((vNam(varn) .EQ. "us" .OR. vNam(varn) .EQ. "u*") .AND. \
                 ws .NE. 1) then
               us=1
               overlay(plot_us,plot_ws)
               n=n+1
               plot_ps(n) = plot_us

               ; ***************************************************
               ; legend for overlaid plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1           
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"u*","w*"/),lgres)       

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
            end if
            if ((vNam(varn) .EQ. "ws" .OR. vNam(varn) .EQ. "w*") .AND. \
                 us .NE. 1) then
               ws=1
               overlay(plot_us,plot_ws)
               n=n+1
               plot_ps(n) = plot_us

               ; ***************************************************
               ; legend for overlaid  plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1            
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"u*","w*"/),lgres)       

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
            end if
         
            if (vNam(varn) .EQ. "umax" .AND. v .NE. 1)
               if (w .NE. 1) then 
                  u=1          
                  overlay(plot_umax,plot_vmax)
                  overlay(plot_umax,plot_wmax)
                
                  n=n+1
                  plot_ps(n) = plot_umax

                  ; ***************************************************
                  ; legend for overlaid plot
                  ; ***************************************************
      
                  lgres                    = True
                  lgMonoDashIndex          = False
                  lgres@lgLabelFont        = "helvetica"    
                  lgres@lgLabelFontHeightF = .1            
                  lgres@vpWidthF           = 0.4           
                  lgres@vpHeightF          = 0.4          
                  lgres@lgDashIndexes      = (/0,0,0/)
                  lgres@lgLineColors       = (/237,144,80/)
                  lbid = gsn_create_legend(\
                                  wks_ps,3,(/"umax","vmax","wmax"/),lgres)

                  amres = True
                  amres@amParallelPosF   = 0.6		    
                  amres@amOrthogonalPosF = -0.2            
                  annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
               end if
            end if
            if (vNam(varn) .EQ. "vmax" .AND. u .NE. 1)
               if (w .NE. 1) then
                  v=1  
                  overlay(plot_umax,plot_vmax)
                  overlay(plot_umax,plot_wmax)
                  n=n+1
                  plot_ps(n) = plot_umax

                  ; ***************************************************
                  ; legend for overlaid plot
                  ; ***************************************************
      
                  lgres                    = True
                  lgMonoDashIndex          = False
                  lgres@lgLabelFont        = "helvetica"    
                  lgres@lgLabelFontHeightF = .1          
                  lgres@vpWidthF           = 0.4           
                  lgres@vpHeightF          = 0.4          
                  lgres@lgDashIndexes      = (/0,0,0/)
                  lgres@lgLineColors       = (/237,144,80/)
                  lbid = gsn_create_legend(\
                                      wks_ps,3,(/"umax","vmax","wmax"/),lgres)
                  amres = True
                  amres@amParallelPosF   = 0.6		    
                  amres@amOrthogonalPosF = -0.2            
                  annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
               end if
            end if
            if (vNam(varn) .EQ. "wmax" .AND. v .NE. 1)
               if(u .NE. 1) then  
                  w=1        
                  overlay(plot_umax,plot_vmax)
                  overlay(plot_umax,plot_wmax)
                  n=n+1
                  plot_ps(n) = plot_umax

                  ; ***************************************************
                  ; legend for overlaid plot
                  ; ***************************************************
      
                  lgres                    = True
                  lgMonoDashIndex          = False
                  lgres@lgLabelFont        = "helvetica"    
                  lgres@lgLabelFontHeightF = .1            
                  lgres@vpWidthF           = 0.4           
                  lgres@vpHeightF          = 0.4          
                  lgres@lgDashIndexes      = (/0,0,0/)
                  lgres@lgLineColors       = (/237,144,80/)
                  lbid = gsn_create_legend(\
                                       wks_ps,3,(/"umax","vmax","wmax"/),lgres)
                  amres = True
                  amres@amParallelPosF   = 0.6		    
                  amres@amOrthogonalPosF = -0.2            
                  annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
               end if
            end if
            
            if (vNam(varn) .EQ. "z_i_wpt" .AND. z .NE. 1) then
               zw=1       
               overlay(plot_z_i_wpt,plot_z_i_pt)
               n=n+1
               plot_ps(n) = plot_z_i_wpt 
       
               ; ***************************************************
               ; legend for overlaid plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1            
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"z_i_wpt","z_i_pt"/),lgres)

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)          
            end if 
            if (vNam(varn) .EQ. "z_i_pt" .AND. zw .NE. 1) then
               z=1      
               overlay(plot_z_i_wpt,plot_z_i_pt)
               n=n+1
               plot_ps(n) = plot_z_i_wpt

               ; ***************************************************
               ; legend for overlaid plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1            
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"z_i_wpt","z_i_pt"/),lgres) 

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)           
            end if    
         
            if ((vNam(varn) .EQ. "wpptp0" .OR. vNam(varn) .EQ. "w"+dq+"pt"+dq+"0") .AND. wp .NE. 1)
               if (wt .NE. 1) then
                  w0=1
                  overlay(plot_wpptp0,plot_wpptp)
                  overlay(plot_wpptp0,plot_wpt) 
                  n=n+1
                  plot_ps(n) = plot_wpptp0 

                  ; ***************************************************
                  ; legend for overlaid plot
                  ; ***************************************************
      
                  lgres                    = True
                  lgMonoDashIndex          = False
                  lgres@lgLabelFont        = "helvetica"    
                  lgres@lgLabelFontHeightF = .1            
                  lgres@vpWidthF           = 0.4           
                  lgres@vpHeightF          = 0.4          
                  lgres@lgDashIndexes      = (/0,0,0/)
                  lgres@lgLineColors       = (/237,144,80/)
                  lbid = gsn_create_legend(\
                                    wks_ps,3,(/"wpptp0","wpptp","wpt"/),lgres)
                  amres = True
                  amres@amParallelPosF   = 0.6		    
                  amres@amOrthogonalPosF = -0.2            
                  annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
               end if            
            end if
            if ((vNam(varn) .EQ. "wpptp" .OR. \
                 vNam(varn) .EQ. "w"+dq+"pt"+dq) .AND. w0 .NE. 1)
               if (wt .NE. 1) then
                 wp=1
                 overlay(plot_wpptp0,plot_wpptp)
                 overlay(plot_wpptp0,plot_wpt) 
                 n=n+1
                 plot_ps(n) = plot_wpptp0
                 
                 ; ***************************************************
                  ; legend for overlaid plot
                  ; ***************************************************
      
                  lgres                    = True
                  lgMonoDashIndex          = False
                  lgres@lgLabelFont        = "helvetica"    
                  lgres@lgLabelFontHeightF = .1            
                  lgres@vpWidthF           = 0.4           
                  lgres@vpHeightF          = 0.4          
                  lgres@lgDashIndexes      = (/0,0,0/)
                  lgres@lgLineColors       = (/237,144,80/)
                  lbid = gsn_create_legend(\
                                     wks_ps,3,(/"wpptp0","wpptp","wpt"/),lgres)

                  amres = True
                  amres@amParallelPosF   = 0.6		    
                  amres@amOrthogonalPosF = -0.2            
                  annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)  
               end if           
            end if
            if (vNam(varn) .EQ. "wpt" .AND. wp .NE. 1)
               if (w0 .NE. 1) then
                  wt=1
                  overlay(plot_wpptp0,plot_wpptp)
                  overlay(plot_wpptp0,plot_wpt) 
                  n=n+1
                  plot_ps(n) = plot_wpptp0

                  ; ***************************************************
                  ; legend for overlaid plot
                  ; ***************************************************
      
                  lgres                    = True
                  lgMonoDashIndex          = False
                  lgres@lgLabelFont        = "helvetica"    
                  lgres@lgLabelFontHeightF = .1            
                  lgres@vpWidthF           = 0.4           
                  lgres@vpHeightF          = 0.4          
                  lgres@lgDashIndexes      = (/0,0,0/)
                  lgres@lgLineColors       = (/237,144,80/)
                  lbid = gsn_create_legend(\
                                  wks_ps,3,(/"wpptp0","wpptp","wpt"/),lgres)
                  amres = True
                  amres@amParallelPosF   = 0.6		    
                  amres@amOrthogonalPosF = -0.2            
                  annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
               end if
            end if

            if ((vNam(varn) .EQ. "pt_0_" .OR. vNam(varn) .EQ. "pt(0)") .AND. \
                pz .NE. 1) then
               p=1      
               overlay(plot_pt_0_,plot_pt_zp_)
               n=n+1 
               plot_ps(n) = plot_pt_0_  
      
               ; ***************************************************
               ; legend for overlaid plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1            
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"pt_0_","pt_zp_"/),lgres)

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
            end if 
            if ((vNam(varn) .EQ. "pt_zp_" .OR. vNam(varn) .EQ. "pt(zp)") .AND.\
                 p .NE. 1) then
               pz=1        
               overlay(plot_pt_0_,plot_pt_zp_)
               n=n+1 
               plot_ps(n) = plot_pt_0_   

               ; ***************************************************
               ; legend for overlaid plot
               ; ***************************************************
      
               lgres                    = True
               lgMonoDashIndex          = False
               lgres@lgLabelFont        = "helvetica"    
               lgres@lgLabelFontHeightF = .1            
               lgres@vpWidthF           = 0.4           
               lgres@vpHeightF          = 0.4          
               lgres@lgDashIndexes      = (/0,0,0/)
               lgres@lgLineColors       = (/237,144,80/)
               lbid = gsn_create_legend(wks_ps,2,(/"pt_0_","pt_zp_"/),lgres)

               amres = True
               amres@amParallelPosF   = 0.6		    
               amres@amOrthogonalPosF = -0.2            
               annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
            end if

            if(vNam(varn) .NE. "pt_zp_" .AND. vNam(varn) .NE. "pt(zp)" .AND. \
               vNam(varn) .NE. "pt_0_" .AND. vNam(varn) .NE. "pt(0)" .AND.   \
               vNam(varn) .NE. "wpt" .AND. vNam(varn) .NE. "wpptp" .AND.     \
               vNam(varn) .NE. "w"+dq+"pt"+dq .AND.                          \
               vNam(varn) .NE. "wpptp0" .AND.                                \
               vNam(varn) .NE. "w"+dq+"pt"+dq+"0" .AND.                      \
               vNam(varn) .NE. "z_i_pt" .AND.                                \
               vNam(varn) .NE. "z_i_wpt" .AND. vNam(varn) .NE. "wmax" .AND.  \
               vNam(varn) .NE. "vmax" .AND. vNam(varn) .NE. "umax" .AND.     \
               vNam(varn) .NE. "ws" .AND. vNam(varn) .NE. "w*" .AND.         \
               vNam(varn) .NE. "us" .AND. vNam(varn) .NE. "u*" .AND.         \
               vNam(varn) .NE. "Es" .AND. vNam(varn) .NE. "E*" .AND.         \
               vNam(varn) .NE. "E") then

               n=n+1
               res@xyLineColors   = (/237/)
               res@xyLabelMode    = False
               res@gsnLeftString  = vNam(varn)
               res@gsnRightString = unit(varn) 
               res@trYMaxF        = max(data(varn,:))
               res@trYMinF        = min(data(varn,:))
               if (min(data(varn,:)) .EQ. max(data(varn,:))) then
                  if (min(data(varn,:)) .EQ. 0)then
                     res@trYMaxF = max(data(varn,:))-0.1
                     res@trYMinF = min(data(varn,:))+0.1
                  end if
                  if (min(data(varn,:)) .LT. 0)then
                     res@trYMaxF = max(data(varn,:))+(min(data(varn,:)))/2
                     res@trYMinF = min(data(varn,:))-(max(data(varn,:)))/2
                  end if
                  if (min(data(varn,:)) .GT. 0)then
                     res@trYMaxF = max(data(varn,:))-(min(data(varn,:)))/2
                     res@trYMinF = min(data(varn,:))+(max(data(varn,:)))/2
                  end if
               end if
               plot_ps(n) = gsn_csm_xy(wks_ps,t,data(varn,:),res)  
            end if

            ; ***************************************************
            ; transparent legend as placeholder for single plots
            ; ***************************************************
            lgresmono                    = True
            lgMonoDashIndex              = False
            lgresmono@lgMonoLineColor    = True
            lgresmono@lgLabelFont        = "helvetica"    
            lgresmono@lgLabelFontHeightF = .1            
            lgresmono@vpWidthF           = 0.4           
            lgresmono@vpHeightF          = 0.4          
            lgresmono@lgLineColor        = 0
            lgresmono@lgPerimOn          = False
            lbid = gsn_create_legend(wks_ps,1,(/""/),lgresmono)       

            amres = True
            amres@amParallelPosF   = 0.6		    
            amres@amOrthogonalPosF = -0.2            
            annoid1 = gsn_add_annotation(plot_ps(n),lbid,amres)
       
         else
            
            print("plot of " + vNam(varn)) 
            
            n=n+1
            res@xyLineColors            = (/237/)
            res@gsnLeftString       = vNam(varn)
            res@gsnRightString      = unit(varn)
        
	    if (norm_t .NE. 1.)then
               res@tiXAxisString        = "t ("+unit_t+")"
	    else
	       res@tiXAxisString        = "t (s)"
	    end if
	    res@tiYAxisString = " "
            res@tiXAxisFontHeightF   = font_size
            res@txFontHeightF        = font_size
            res@tiYAxisFontHeightF   = font_size
	    res@tmXBLabelFontHeightF = font_size
            res@tmYLLabelFontHeightF = font_size
            res@trYMaxF        = max(data(varn,:))
            res@trYMinF        = min(data(varn,:)) 
            if (min(data(varn,:)) .EQ. max(data(varn,:))) then
               if (min(data(varn,:)) .EQ. 0)then
                     res@trYMaxF = max(data(varn,:))-0.1
                     res@trYMinF = min(data(varn,:))+0.1
                  end if
                  if (min(data(varn,:)) .LT. 0)then
                     res@trYMaxF = max(data(varn,:))+(min(data(varn,:)))/2
                     res@trYMinF = min(data(varn,:))-(max(data(varn,:)))/2
                  end if
                  if (min(data(varn,:)) .GT. 0)then
                     res@trYMaxF = max(data(varn,:))-(min(data(varn,:)))/2
                     res@trYMinF = min(data(varn,:))+(max(data(varn,:)))/2
                  end if
            end if
            plot_ps(n) = gsn_csm_xy(wks_ps,t,data(varn,:),res)
           
         end if     
      end if
   end do
  
   ; ***************************************************
   ; merge plots onto one page 
   ; ***************************************************

   no_frames = 0
 
   if ((format_out .EQ. "eps" .OR. format_out .EQ. "epsi") .AND. \
        n .gt. no_rows*no_columns) then
      gsn_panel(wks_ps,plot_ps(1:n),(/n,1/),resP)
      print(" ")
      print("Outputs to .eps or .epsi have only one frame")
      print(" ")
   else
      do np = 1,n,no_rows*no_columns    
         if ( np + no_rows*no_columns .gt. n) then   
            gsn_panel(wks_ps, plot_ps(np:n),(/no_rows,no_columns/),resP)
            no_frames = no_frames + 1            
         else
            gsn_panel(wks_ps, plot_ps(np:np+no_rows*no_columns-1),\
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
