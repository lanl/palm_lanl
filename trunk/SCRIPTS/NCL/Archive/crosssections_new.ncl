load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; cross-sections
; last change: $Id: crosssections_new.ncl 127 2007-10-23 11:05:25Z letzel $

begin
;
; set default value(s) for shell script variables assigned on command line
  if ( .not. isvar("cm") ) then     ; colormap
    cm = "ncview_default"
  end if
  if ( .not. isvar("di") ) then     ; input directory (with final /)
    di = ""
  end if
  if ( .not. isvar("d") ) then      ; output directory (with final /)
    d = di
  end if
  if ( .not. isvar("fi") ) then     ; base name of input file (without suffix)
    fi = "example_xy"
  end if
  if ( .not. isvar("fill_mode") ) then ; "AreaFill", "RasterFill" or "CellFill"
    fill_mode = "AreaFill"
  end if
  if ( .not. isvar("fo") ) then     ; base name of output files (without suffix)
    fo = ""
  end if
  if ( .not. isvar("mode") ) then   ; output mode ("Fill" or "Line")
    mode = "Fill"
  end if
  if ( .not. isvar("t") ) then      ; output time step
    t = 0
  end if
  if ( .not. isvar("var") ) then    ; variable to be output
    var = "u_xy"
  end if
  if ( .not. isvar("xs") ) then     ; output x-coordinate range start (in m)
    xs = -1e+38
  end if
  if ( .not. isvar("xe") ) then     ; output x-coordinate range end (in m)
    xe = 1e+38
  end if
  if ( .not. isvar("ys") ) then     ; output y-coordinate range start (in m)
    ys = -1e+38
  end if
  if ( .not. isvar("ye") ) then     ; output y-coordinate range end (in m)
    ye = 1e+38
  end if
  if ( .not. isvar("zs") ) then     ; output z-coordinate range start (in m)
    zs = -1e+38
  end if
  if ( .not. isvar("ze") ) then     ; output z-coordinate range end (in m)
    ze = 1e+38
  end if
;
; open input file
  f  = addfile( di + fi  + ".nc", "r" )
;
; open workstation(s) and set colormap
  wks_x11 = gsn_open_wks("x11","cross-section") ; X11 workstation
  gsn_define_colormap(wks_x11,cm)
  if ( isvar("fo") ) then
    wks_pdf = gsn_open_wks("pdf",d+fo)  ; optional workstations
    gsn_define_colormap(wks_pdf,cm)
    wks_eps = gsn_open_wks("eps",d+fo)  ; for output on file
    gsn_define_colormap(wks_eps,cm)
    wks_ps  = gsn_open_wks("ps",d+fo)
    gsn_define_colormap(wks_ps,cm)
  end if
;
; read input data using 'coordinate subscripting'
; NCL uses the closest corresponding values (in case of two equally distant 
; values NCL chooses the smaller one)
  raw_data = f->$var$(t:t,{zs:ze},{ys:ye},{xs:xe})
  raw_data!0 = "t"
  raw_data!1 = "z"
  raw_data!2 = "y"
  raw_data!3 = "x"
  time = raw_data&t
;
; reduce variable dimensions from 4D to 2D according to output ranges
  if ( zs .eq. ze ) then
    data = raw_data(0,0,:,:)
    x_axis = "x"
    y_axis = "y"
    plane  = "z"

    dx =  f->x(1)
    dy =  f->y(1)

    if ( raw_data&z .eq. -1 ) then
      level = "-average"
    else
      level = "=" + raw_data&z + "m"
    end if
  else
    if ( ys .eq. ye ) then
      data = raw_data(0,:,0,:)
      x_axis = "x"
      y_axis = "z"
      plane  = "y"

      dx =  f->x(1)
      if (isfilevar(f, "zw_3d")) then
        dz =  f->zw_3d(1)
      else
        dz =  f->zw(1)
      end if

      if ( raw_data&y .eq. -1 ) then
        level = "-average"
      else
        level = "=" + raw_data&y + "m"
      end if
    else
      if ( xs .eq. xe ) then
        data = raw_data(0,:,:,0)
        x_axis = "y"
        y_axis = "z"
        plane  = "x"

        dy =  f->y(1)
        if (isfilevar(f, "zw_3d")) then
          dz =  f->zw_3d(1)
        else
          dz =  f->zw(1)
        end if

        if ( raw_data&x .eq. -1 ) then
          level = "-average"
        else
          level = "=" + raw_data&x + "m"
        end if
      end if
    end if
  end if
  delete( raw_data )
;
; set up resources
  cs_res                         = True
  cs_res@gsnDraw                 = False  
  cs_res@gsnFrame                = False 
  cs_res@gsnMaximize             = True 
  cs_res@gsnPaperOrientation     = "portrait"
  cs_res@gsnPaperWidth           = 8.27
  cs_res@gsnPaperHeight          = 11.69
  cs_res@gsnPaperMargin          = 0.79
  cs_res@tiMainFuncCode          = "~"
  cs_res@tiMainFontHeightF       = 0.015
  cs_res@tiMainString            = f@title
  cs_res@tiXAxisString           = x_axis + " [m]"
  cs_res@tiYAxisString           = y_axis + " [m]"
;  cs_res@gsnLeftString           = ""   ; gsn_csm_* scripts use default data
  cs_res@gsnCenterString         = "t=" + time + "s  " + plane + level
;  cs_res@gsnRightString          = ""   ; gsn_csm_* scripts use default data
  cs_res@tmXBMode                ="Automatic"
  cs_res@tmYLMode                ="Automatic"
  if ( mode .eq. "Fill" ) then
    cs_res@cnFillOn                = True
    cs_res@gsnSpreadColors         = True
    cs_res@cnFillMode              = fill_mode
    cs_res@lbOrientation           = "Vertical"  ; vertical label bar
    cs_res@cnLinesOn               = False
    cs_res@cnLineLabelsOn          = False
  end if

  vector     = True
  co_overlay = True

  if ( .not. isvar("wv1") ) then
    vector = False
  end if
  if ( .not. isvar("wv2") ) then
    vector = False
  end if
  if ( .not. isvar("co") ) then
    co_overlay = False
  end if

  if ( vector .eq. True ) then

    if ( wv1 .eq. "u_xy" .or.  wv1 .eq. "u_xz" .or. wv1 .eq. "u_yz" ) then

      raw_data = f->$wv1$(t:t,{zs:ze},{ys:ye},{xs-(0.5*dx):xe-(0.5*dx)}) 
      raw_data!0 = "t"
      raw_data!1 = "z"
      raw_data!2 = "y"
      raw_data!3 = "x"
      time = raw_data&t

    else

      if ( wv1 .eq. "v_xy" .or.  wv1 .eq. "v_xz" .or. wv1 .eq. "v_yz" ) then

        raw_data = f->$wv1$(t:t,{zs:ze},{ys-(0.5*dy):ye-(0.5*dy)},{xs:xe}) 
        raw_data!0 = "t"
        raw_data!1 = "z"
        raw_data!2 = "y"
        raw_data!3 = "x"
        time = raw_data&t

      else

        if ( wv1 .eq. "w_xy" .or.  wv1 .eq. "w_xz" .or. wv1 .eq. "w_yz" ) then

          raw_data = f->$wv1$(t:t,{zs+(0.5*dz):ze+(0.5*dz)},{ys:ye},{xs:xe}) 
          raw_data!0 = "t"
          raw_data!1 = "z"
          raw_data!2 = "y"
          raw_data!3 = "x"
          time = raw_data&t

        end if
      end if
    end if

    if ( zs .eq. ze ) then
      data_u = raw_data(0,0,:,:)
    else
      if ( ys .eq. ye ) then
        data_u = raw_data(0,:,0,:)
      else
        if ( xs .eq. xe ) then
          data_u = raw_data(0,:,:,0)
        end if
      end if
    end if
    delete( raw_data )

    if ( wv2 .eq. "u_xy" .or.  wv2 .eq. "u_xz" .or. wv2 .eq. "u_yz" ) then

      raw_data = f->$wv2$(t:t,{zs:ze},{ys:ye},{xs-(0.5*dx):xe-(0.5*dx)}) 
      raw_data!0 = "t"
      raw_data!1 = "z"
      raw_data!2 = "y"
      raw_data!3 = "x"
      time = raw_data&t

    else

      if ( wv2 .eq. "v_xy" .or.  wv2 .eq. "v_xz" .or. wv2 .eq. "v_yz" ) then

        raw_data = f->$wv2$(t:t,{zs:ze},{ys-(0.5*dy):ye-(0.5*dy)},{xs:xe}) 
        raw_data!0 = "t"
        raw_data!1 = "z"
        raw_data!2 = "y"
        raw_data!3 = "x"
        time = raw_data&t

      else

        if ( wv2 .eq. "w_xy" .or.  wv2 .eq. "w_xz" .or. wv2 .eq. "w_yz" ) then

          raw_data = f->$wv2$(t:t,{zs+(0.5*dz):ze+(0.5*dz)},{ys:ye},{xs:xe}) 
          raw_data!0 = "t"
          raw_data!1 = "z"
          raw_data!2 = "y"
          raw_data!3 = "x"
          time = raw_data&t

        end if
      end if
    end if


    if ( zs .eq. ze ) then
      data_v = raw_data(0,0,:,:)
    else
      if ( ys .eq. ye ) then
        data_v = raw_data(0,:,0,:)
      else
        if ( xs .eq. xe ) then
          data_v = raw_data(0,:,:,0)
        end if
      end if
    end if
    delete( raw_data )

    wv_res                     = True
    wv_res@gsnDraw             = False  
    wv_res@gsnFrame            = False
    wv_res@gsnRightString      = ""
    wv_res@gsnLeftString       = ""
    wv_res@gsnCenterString     = ""
    wv_res@vpClipOn            = True

  end if

  if ( co_overlay .eq. True ) then

    raw_data = f->$co$(t:t,{zs:ze},{ys:ye},{xs:xe}) 
    raw_data!0 = "t"
    raw_data!1 = "z"
    raw_data!2 = "y"
    raw_data!3 = "x"
    time = raw_data&t

    if ( zs .eq. ze ) then
      data2 = raw_data(0,0,:,:)
    else
      if ( ys .eq. ye ) then
        data2 = raw_data(0,:,0,:)
      else
        if ( xs .eq. xe ) then
          data2 = raw_data(0,:,:,0)
        end if
      end if
    end if
    delete( raw_data )

    co_res                              = True
    co_res@cnLineThicknessF             = 1.0
    co_res@gsnDraw                      = False
    co_res@gsnFrame                     = False
    co_res@gsnContourZeroLineThicknessF = 2.0
    co_res@gsnContourNegLineDashPattern = 4
    co_res@cnLineColor                  = "Black"
    co_res@gsnRightString   = ""
    co_res@gsnLeftString    = ""
    co_res@gsnCenterString  = ""

;    co_res1@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;    co_res1@cnMinLevelValF       = -5.                ; set min contour level
;    co_res1@cnMaxLevelValF       =  5.                ; set max contour level
;    co_res1@cnLevelSpacingF      =  0.5               ; set contour spacing
     co_res@cnLineLabelsOn        = True
     co_res@cnLineLabelPlacementMode   = "constant"
;    co_res@cnLineLabelBackgroundColor = "white"

  end if

  if ( co_overlay .eq. True .and. vector .eq. True ) then 

    if ( isvar("fo") ) then

      plot_pdf2 = gsn_csm_contour(wks_pdf,data2,co_res)
      plot_pdf1 = gsn_csm_contour(wks_pdf,data,cs_res)
      plot_pdf3 = gsn_csm_vector(wks_pdf,data_u,data_v,wv_res)
      overlay(plot_pdf1,plot_pdf2)
      overlay(plot_pdf1,plot_pdf3)
      draw(plot_pdf1)
      frame(wks_pdf)

      plot_eps2 = gsn_csm_contour(wks_eps,data2,co_res)
      plot_eps1 = gsn_csm_contour(wks_eps,data,cs_res)
      plot_eps3 = gsn_csm_vector(wks_eps,data_u,data_v,wv_res)
      overlay(plot_eps1,plot_eps2)
      overlay(plot_eps1,plot_eps3)
      draw(plot_eps1)
      frame(wks_eps)

      plot_ps2 = gsn_csm_contour(wks_ps,data2,co_res)
      plot_ps1 = gsn_csm_contour(wks_ps,data,cs_res)
      plot_ps3 = gsn_csm_vector(wks_ps,data_u,data_v,wv_res)
      overlay(plot_ps1,plot_ps2)
      overlay(plot_ps1,plot_ps3)
      draw(plot_ps1)
      frame(wks_ps)

    end if

    plot1 = gsn_csm_contour(wks_x11,data,cs_res)
    plot2 = gsn_csm_contour(wks_x11,data2,co_res)
    plot3 = gsn_csm_vector(wks_x11,data_u,data_v,wv_res)

    overlay(plot1,plot2)
    overlay(plot1,plot3)

    draw(plot1)
    frame(wks_x11)

  else

    if ( co_overlay .eq. True .and. vector .ne. True ) then

      if ( isvar("fo") ) then

        plot_pdf2 = gsn_csm_contour(wks_pdf,data2,co_res)
        plot_pdf1 = gsn_csm_contour(wks_pdf,data,cs_res)
        overlay(plot_pdf1,plot_pdf2)
        draw(plot_pdf1)       
        frame(wks_pdf)

        plot_eps2 = gsn_csm_contour(wks_eps,data2,co_res)
        plot_eps1 = gsn_csm_contour(wks_eps,data,cs_res)
        overlay(plot_eps1,plot_eps2)
        draw(plot_eps1)
        frame(wks_eps)

        plot_ps2 = gsn_csm_contour(wks_ps,data2,co_res)
        plot_ps1 = gsn_csm_contour(wks_ps,data,cs_res)
        overlay(plot_ps1,plot_ps2)
        draw(plot_ps1)
        frame(wks_ps)

      end if 

      plot1 = gsn_csm_contour(wks_x11,data,cs_res)
      plot2 = gsn_csm_contour(wks_x11,data2,co_res)

      overlay(plot1,plot2)

      draw(plot1)
      frame(wks_x11)

    else

    if ( co_overlay .ne. True .and. vector .eq. True ) then 

      if ( isvar("fo") ) then

        plot_pdf1 = gsn_csm_contour(wks_pdf,data,cs_res)
        plot_pdf3 = gsn_csm_vector(wks_pdf,data_u,data_v,wv_res)
        overlay(plot_pdf1,plot_pdf3)
        draw(plot_pdf1)
        frame(wks_pdf)

        plot_eps1 = gsn_csm_contour(wks_eps,data,cs_res)
        plot_eps3 = gsn_csm_vector(wks_eps,data_u,data_v,wv_res)
        overlay(plot_eps1,plot_eps3)
        draw(plot_eps1)
        frame(wks_eps)

        plot_ps1 = gsn_csm_contour(wks_ps,data,cs_res)
        plot_ps3 = gsn_csm_vector(wks_ps,data_u,data_v,wv_res)
        overlay(plot_ps1,plot_ps3)
        draw(plot_ps1)
        frame(wks_ps)

      end if 

      plot1 = gsn_csm_contour(wks_x11,data,cs_res)
      plot3 = gsn_csm_vector(wks_x11,data_u,data_v,wv_res)

      overlay(plot1,plot3)

      draw(plot1)
      frame(wks_x11)

    else

    if ( co_overlay .ne. True .and. vector .ne. True) then 

      if ( isvar("fo") ) then

        plot_pdf1 = gsn_csm_contour(wks_pdf,data,cs_res)
        draw(plot_pdf1)
        frame(wks_pdf)

        plot_eps1 = gsn_csm_contour(wks_eps,data,cs_res)
        draw(plot_eps1)
        frame(wks_eps)

        plot_ps1 = gsn_csm_contour(wks_ps,data,cs_res)
        draw(plot_ps1)
        frame(wks_ps)

      end if 

      plot1 = gsn_csm_contour(wks_x11,data,cs_res)

      draw(plot1)
      frame(wks_x11)

    end if
    end if
    end if

  end if

end
