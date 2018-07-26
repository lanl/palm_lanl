load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; cross-sections
; last change: $Id: crosssections.ncl 127 2007-10-23 11:05:25Z letzel $

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
;
; data output
  plot_x11 = gsn_csm_contour(wks_x11,data,cs_res)
  if ( isvar("fo") ) then
    plot_pdf = gsn_csm_contour(wks_pdf,data,cs_res)
    plot_eps = gsn_csm_contour(wks_eps,data,cs_res)
    plot_ps  = gsn_csm_contour(wks_ps, data,cs_res)
  end if
end
