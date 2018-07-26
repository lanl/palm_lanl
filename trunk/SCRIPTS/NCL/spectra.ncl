load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;***************************************************
; Checking the kind of the script
;***************************************************

function which_script()
local script 
begin 
   script = "spectra"
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
   
   if (log_x .NE. 0 .AND. log_x .NE. 1)then
      print(" ")
      print("'log_x'= "+log_x+" is invalid and set to 1")
      print(" ")
      log_x = 1
   end if  
   
   if (log_y .NE. 0 .AND. log_y .NE. 1)then
      print(" ")
      print("'log_y'= "+log_y+" is invalid and set to 1")
      print(" ")
      log_y = 1
   end if    
  
   if (norm_y .EQ. 0.) then 
      print(" ")
      print("Normalising with 0 is not allowed, 'norm_y' is set to 1.0")
      print(" ")
      norm_y = 1.0
   end if
   if (norm_x .EQ. 0.) then 
      print(" ")
      print("Normalising with 0 is not allowed, 'norm_x' is set to 1.0")
      print(" ")
      norm_x= 1.0
   end if
   
   if (sort .NE. "height" .AND. sort .NE. "time") then  
      print(" ")
      print("'sort'= "+sort+" is invalid and set to 'height'")
      print(" ")
      sort = "height"  
   end if
   
   if (black .NE. 0 .AND. black .NE. 1)then
      print(" ")
      print("'black'= "+black+" is invalid and set to 0")
      print(" ")
      black = 0
   end if
 
   if (dash .NE. 0 .AND. dash .NE. 1)then
      print(" ")
      print("'dash'= "+dash+" is invalid and set to 0")
      print(" ")
      dash = 0
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
   
   vNam = getfilevarnames(f_att) 
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
   dim = dimsizes(vNam)
   vDim = getfiledimsizes(f_att) 
 
   t_all = f[:]->time
   nt    = dimsizes(t_all)

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
   
   k_x=f_att->k_x
   dimx=dimsizes(k_x)
   k_y=f_att->k_y
   dimy=dimsizes(k_y)
   
   
   dim_level=dimsizes(height_level)

   do i=0,dim-1
      if (vNam(i) .EQ. "zu_sp")then
         zu=f_att->zu_sp	 
	 if (height_level(0) .EQ. -1)then
	    dimz=dimsizes(zu)
	 else
	    if (dim_level .GT. dimsizes(zu))then
               print(" ")
               print("'height_level' has more elements than available "+\
                     "height levels in input file (= "+dimsizes(zu)+")")
               print(" ")
               exit
	    else
	       zuh=new(dim_level,double)
	       do le=0,dim_level-1
                  if (height_level(le) .GE. dimsizes(zu)) then
                     no_levels=dimsizes(zu)-1
                     print(" ")
                     print("Element "+le+" of 'height_level' is greater "  +\
                          "than the maximum available index in input file "+\
                          "which is "+no_levels+". Note that the first "   +\
                          "element has the index 0.")
                     print(" ")
                     exit
                  end if
	          zuh(le)=zu(height_level(le))
	       end do
	       dimz=dim_level
            end if   
	 end if	 
      else 
         if (vNam(i) .EQ. "zw_sp")then
	    zw=f_att->zw_sp
	    if (height_level(0) .EQ. -1)then              
               dimz=dimsizes(zw)
	    else
	       if (dim_level .GT. dimsizes(zw))then
                  print(" ")
                  print("'height_level' has more elements than available "+\
                        "height levels in input file (= "+dimsizes(zw)+")")
                  print(" ")
                  exit
	       else
	          zwh=new(dim_level,double)
	          do le=0,dim_level-1
                     if (height_level(le) .GE. dimsizes(zw)) then
                        no_levels=dimsizes(zw)-1
                        print(" ")
                        print("Element "+le+" of 'height_level' is greater "+\
                              "than the maximum available index in input "  +\
                              "file which is "+no_levels+". Note that the " +\
                              "first element has the index 0.")
                        print(" ")
                        exit
                     end if
	             zwh(le)=zw(height_level(le))
	          end do
	          dimz=dim_level
               end if   
	    end if
         end if
      end if
   end do

   ;****************************************************	
   ; start of time step and different types of mistakes that could be done
   ;****************************************************
   
   if (start_time_step .EQ. -1.d) then   	
      start_time_step=t_all(0)/3600
   else
      if (start_time_step .GT. t_all(nt-1)/3600)then
         print(" ")
         print("'start_time_step' = "+ start_time_step +"h is greater than "+\
               "last time step = "+ t_all(nt-1)+"s = "+t_all(nt-1)/3600+"h")
         print(" ")
         print("Select another 'start_time_step'")
         print(" ")
         exit
      end if
      if (start_time_step .LT. t_all(0)/3600)then
         print(" ")
         print("'start_time_step' = "+ start_time_step +"h is lower than "+\
               "first time step = "+ t_all(0)+"s = "+t_all(0)/3600+"h")        
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
   
   ;****************************************************
   ; end of time step and different types of mistakes that could be done
   ;****************************************************

   if (end_time_step .EQ. -1.d) then     	
      end_time_step = t_all(nt-1)/3600
   else
      if (end_time_step .GT. t_all(nt-1)/3600)then
         print(" ")
         print("'end_time_step' = "+ end_time_step +"h is greater than "+\
               "last time step = " + t_all(nt-1)+"s = "+t_all(nt-1)/3600+"h")
         print(" ")
         print("Select another 'end_time_step'")  
         print(" ")
         exit
      end if
      if (end_time_step .LT. start_time_step/3600)then
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
   print("                     till "+t_all(end_time_step)/3600+" h = "+\
         t_all(end_time_step)+" s => index = "+end_time_step)
   print(" ")

   dimt = end_time_step-start_time_step+1
 
   ;***************************************************
   ; set up recourses
   ;***************************************************
     
   res = True
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
   
   res@lgLabelFontHeightF     = font_size_legend 
   res@lgTitleFontHeightF     = font_size
   res@txFontHeightF      = font_size
   res@tiXAxisFontHeightF = font_size
   res@tiYAxisFontHeightF = font_size
   res@tmXBLabelFontHeightF = font_size
   res@tmYLLabelFontHeightF = font_size
   
   res@tmXBMinorPerMajor = 4
   res@tmYLMinorPerMajor = 4
   
   if (log_x .EQ. 1) then
      res@trXLog = True 
   else
      res@trXLog = False 
   end if
   if (log_y .EQ. 1)then
      res@trYLog = True 
   else  
      res@trYLog = False
   end if

   legend_label=new(dimt,string)
   legend_label_zu=new(dimz,double)
   legend_label_zw=new(dimz,double)
   legend_label_z=new(dimz,double)
   do p=start_time_step,end_time_step
      legend_label(p-start_time_step)=sprintf("%6.2f", t_all(p)/3600)
   end do  
   if (sort .EQ. "time")
      plot  = new(dim*dimz,graphic)
      np=dimt
      res@lgTitleString = "Time (h)"
   else
      plot  = new(dim*dimt,graphic)
      np=dimz
      
      do p=0,dimz-1
         if (height_level(0) .EQ. -1)then
            legend_label_zu(p)=round(zu(p),3)
            legend_label_zw(p)=round(zw(p),3)
	 else
	    legend_label_zu(p)=round(zuh(p),3)
	    legend_label_zw(p)=round(zwh(p),3)
	 end if
      end do
   end if

   if (black .eq. 0 ) then 
      if (np .EQ. 1)then
         color = 237
      else   
         step=round(235/(np-1),3) 
         color = ispan(2,237,step)
      end if
   else
      color = 2
   end if
   if ( dash .eq. 0 ) then
      res@xyMonoDashPattern = True 
   end if

   if ( format_out .EQ. "pdf" .OR. format_out .EQ. "ps" ) then
      format_out@wkPaperSize = "A4"
   end if
   if ( format_out .EQ. "png" ) then
      format_out@wkWidth  = 1000
      format_out@wkHeight = 1000
   end if

   wks=gsn_open_wks(format_out,file_out)
   gsn_define_colormap(wks,"rainbow+white") 

   n=0
   do varn =dim-1,0,1 
      
      check = True

      if ( isStrSubset( vNam(varn), "time") .OR.  \
           isStrSubset( vNam(varn), "zu_sp") .OR. \
           isStrSubset( vNam(varn), "zw_sp") .OR. \
           isStrSubset( vNam(varn), "k_x") .OR.   \
           isStrSubset( vNam(varn), "k_y")) then
            check = False
      end if

      if (var .NE. "all") then
         check = isStrSubset( var,","+vNam(varn)+"," )
      end if

      if(check) then 

         temp = f[:]->$vNam(varn)$
	 data = temp(start_time_step:end_time_step,0:dimz-1,:) 

         temp_att = f_att->$vNam(varn)$
         a=getvardims(temp_att)
         b=dimsizes(a)
         delete(temp_att)

         if (height_level(0) .NE. -1)then
	    do te=0,dimz-1
	       data(:,te,:) = temp(start_time_step:end_time_step,\
                                                    height_level(te),:)	       
	    end do
	 end if

         data=data/(norm_y*norm_x) 
            
         do i=0,b-1            
            if (isStrSubset( a(i),"zu_sp" ))then
               legend_label_z=legend_label_zu
	       if (height_level(0) .NE. -1)then
	          z=zuh
	       else
	          z=zu
	       end if
            else
               if (isStrSubset( a(i),"zw_sp" ))then
                  legend_label_z=legend_label_zw 
		  if (height_level(0) .NE. -1)then
	             z=zwh
	          else
	             z=zw
	          end if   
               end if
            end if
         end do  
	 
         if (check_vType) then
            min_y=new(dimz,double)
	    max_y=new(dimz,double)
         else
            min_y=new(dimz,float)
	    max_y=new(dimz,float)
         end if 
         min_x=new(dimz,double)
	 max_x=new(dimz,double)
	 
	 plot_h  = new(dimz,graphic)
	 
         if (isStrSubset(vNam(varn),"x"))then
	    x_axis = new((/dimz,dimx/),double)
	    do q=0,dimz-1
               x_axis(q,:) = f_att->k_x
	       x_axis = x_axis/norm_x
	    end do
	    if (norm_x .NE. 1.)then
               res@tiXAxisString = "k~B~x~N~ ("+unit_x+")"
	    else
	       if (norm_height .EQ. 1)then
	          res@tiXAxisString = "k~B~x~N~ x z (1)"
	       else
	          res@tiXAxisString = "k~B~x~N~ (1/m)"
	       end if
	    end if
	    dim_r=dimx
         else
	    x_axis=new((/dimz,dimy/),double)
	    do q=0,dimz-1
               x_axis(q,:) = f_att->k_y 
	       x_axis = x_axis/norm_x
	    end do
	    if (norm_x .NE. 1.)then
               res@tiXAxisString = "k~B~y~N~ ("+unit_x+")"
	    else
	       if (norm_height .EQ. 1)then
	          res@tiXAxisString = "k~B~y~N~ x z (1)"
	       else
	          res@tiXAxisString = "k~B~y~N~ (1/m)"
	       end if 
	    end if
	    dim_r=dimy
         end if
        
         if (sort .EQ. "time")
	    res@xyLineColors = color
	    res@pmLegendDisplayMode     = "Always"
            res@pmLegendSide            = "Top"
            res@pmLegendParallelPosF    = 1.2
            res@pmLegendOrthogonalPosF  = -1.0
            res@pmLegendWidthF          = 0.12
            res@pmLegendHeightF         = 0.04*\
                                          (end_time_step-start_time_step+1)
            do p=dimz-1,0,1
	       if (log_y .EQ. 1)then  
                  do q=0,dimt-1
                     do r=0,dim_r-1 
                        if (data(q,p,r) .EQ. 0)then
                           st=p+start_time_step
                           print(" ")
                           print("'"+vNam(varn)+"("+st+","+q+","+r+")' is "+\
                                 "equal 0. Logarithmic scale for y-axis "+\
                                 "and height "+z(p)+" cannot be used")
                           print(" ")
                           res@trYLog = False
                        end if
                     end do
                  end do
	       end if
               res@gsnLeftString      = vNam(varn)
               res@gsnRightString     = "Height = "+z(p)+"m"               
               res@tiYAxisString      = "("+unit_y+")"           
               res@xyExplicitLegendLabels  = legend_label	       
	       if (norm_height .EQ. 1)then
	          data(:,p,:)=data(:,p,:)*doubletofloat(z(p))
		  x_axis(p,:) = x_axis(p,:)*z(p)
	       end if
	       res@trXMinF = min(x_axis(p,:))
               res@trXMaxF = max(x_axis(p,:))
               plot(n)  = gsn_csm_xy(wks,x_axis(p,:),data(:,p,:),res)
               n=n+1
            end do
         else
            if (sort .EQ. "height")
               do p=0,dimt-1 	       
                  do q=0,dimz-1
                     do r=0,dim_r-1 
                        if (data(p,q,r) .EQ. 0)then
                           st=p+start_time_step
                           print(" ")
                           print("'"+vNam(varn)+"("+st+","+q+","+r+")' "+\
                                 "is equal 0. Logarithmic scale for y-axis "+\
                                 "and time "+legend_label(p)+" h cannot be used")
                           print(" ")
                           res@trYLog = False
                        end if
                     end do
		     if (norm_height .EQ. 1 .AND. p .EQ. 0)then
		        data(p,q,:) = data(p,q,:)*\
                                           doubletofloat(legend_label_z(q))
			x_axis(q,:) = x_axis(q,:)*\
                                           doubletofloat(legend_label_z(q))
		     end if
	             max_y(q)=max(data(p,q,:))
		     min_y(q)=min(data(p,q,:))
		     min_x(q)=min(x_axis(q,:))
		     max_x(q)=max(x_axis(q,:))
                  end do 
		  do q=0,dimz-1
		     res@xyLineColor = color(q)
		     if (dash .EQ. 1)then
		        res@xyDashPattern = q
	             end if
		     if (q .EQ. 0)then
		        res@tiYAxisString      = "("+unit_y+")"
		        res@gsnLeftString      = vNam(varn)
                        res@gsnRightString     = "Time = "+legend_label(p)+"h"
		        res@trXMinF = min(min_x)
                        res@trXMaxF = max(max_x)
                        res@trYMinF = min(min_y)
                        res@trYMaxF = max(max_y*10)
			
			plot_h(q)  = gsn_csm_xy(wks,x_axis(q,:),\
                                                         data(p,q,:),res)
			
			lgres = True
			if (dash .EQ. 0)then
                           lgres@lgMonoDashIndex = True
			else
                           lgres@lgDashIndexes   = ispan(0,dimz-1,1)
			end if
			if (black .EQ. 1)then
			   lgres@lgMonoLineColors = True
			else
			   lgres@lgLineColors = color
			end if
			lgres@lgTitleString      = "Height (m)"  
                        lgres@lgLabelFont        = "helvetica"
			lgres@lgLabelFontHeightF = font_size_legend*6 
                        lgres@lgTitleFontHeightF = font_size      
                        lgres@vpWidthF           = 0.12           
                        lgres@vpHeightF          = font_size_legend*(dimz+3)
  
                        lbid = gsn_create_legend(wks,dimz,legend_label_z,lgres)

                        amres = True
                        amres@amParallelPosF   = 0.75		    
                        amres@amOrthogonalPosF = 0.15            
                        annoid1 = gsn_add_annotation(plot_h(q),lbid,amres)
		     else
		        plot_h(q)  = gsn_csm_xy(wks,x_axis(q,:),\
                                                        data(p,q,:),res) 
			overlay(plot_h(0),plot_h(q))
		     end if
		  end do              
		  plot(n)=plot_h(0)
		  n=n+1
               end do
            end if
	 end if
         delete(data)
	 delete(temp)
         delete(x_axis)
	 delete(min_x)
	 delete(max_x)
	 delete(min_y)
	 delete(max_y)
	 delete(plot_h)
      end if
   end do

   if (n .EQ. 0) then
      print(" ")
      print("The variables 'var="+var+"' do not exist on your input file;")
      print("be sure to have one comma berfore and after each variable")
      print(" ")
      exit
   end if
  
   ; ***************************************************
   ; merge plots onto one page 
   ; ***************************************************

   resP                            = True
   resP@gsnMaximize                = True 
   resP@gsnPanelXWhiteSpacePercent = 4.0
   resP@gsnPanelYWhiteSpacePercent = 4.0
   resP@txFont                     = "helvetica"
   resP@txString                   = f_att@title
   resP@txFuncCode                 = "~"
   resP@txFontHeightF              = 0.0105

   no_frames = 0

   if ((format_out .EQ. "eps" .OR. format_out .EQ. "epsi") .AND. \
                                          n .gt. no_rows*no_columns) then
      gsn_panel(wks,plot,(/n,1/),resP)
      print(" ")
      print("Outputs to .eps or .epsi have only one frame")
      print(" ")
   else    
      do i = 0,n-1, no_rows*no_columns
         if( (i+no_rows*no_columns) .gt. (n-1)) then
            gsn_panel(wks,plot(i:n-1),(/no_rows,no_columns/),resP)
            no_frames = no_frames + 1
         else
            gsn_panel(wks,plot(i:i+no_rows*no_columns-1),\
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
