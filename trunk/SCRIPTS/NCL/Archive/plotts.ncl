load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; time series
; last change: $Id: plotts.ncl 127 2007-10-23 11:05:25Z letzel $

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
  if ( .not. isvar("fo") ) then     ; base name of output files (without suffix)
    fo = ""
  end if
  if ( .not. isvar("ts") ) then      ; output time step
    ts = 0
  end if
  if ( .not. isvar("te") ) then     ; output x-coordinate range start (in m)
    te = 0
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

  if ( .not. isvar("var1") .and. .not. isvar("var2") .and. .not. isvar("var3")  .and. .not. isvar("var4") .and. .not. isvar("var5") .and. .not. isvar("var6")) then
    vNam = getfilevarnames(f)
    n2 = dimsizes(vNam)-1

    print(vNam)

  else

    vNam = new((/6/), string)

    vNam(0) = var1
    vNam(1) = var2
    vNam(2) = var3
    vNam(3) = var4
    vNam(4) = var5
    vNam(5) = var6
    n2 = 6
  end if

  t= f ->time({ts:te})

  res=True
  res@gsnDraw=False
  res@gsnFrame=False

  resP = True
  resP@txString="time series"

  res@vpWidthF=4

  plot=new(40,graphic)                 ; create a plot array
  plot_pdf=new(40,graphic)             ; create a plot array
  plot_ps=new(40,graphic)              ; create a plot array
  plot_eps=new(40,graphic)             ; create a plot array

  n = 0

  do while ( n .lt. n2)
    data = f ->$vNam(n)$({ts:te})
    plot(n) = gsn_csm_xy(wks_x11,t,data,res)
    plot_pdf(n) = gsn_csm_xy(wks_pdf,t,data,res)
    plot_eps(n) = gsn_csm_xy(wks_eps,t,data,res)
    plot_ps(n) = gsn_csm_xy(wks_ps,t,data,res)
    n = n + 1
  end do

  pa = 0
  pb = 5

  do while ( pa .lt. dimsizes(plot) )
    gsn_panel(wks_pdf,plot_pdf(pa:pb),(/6,1/),resP) ; Panel the plots, first using rows x columns, then using number of plots per row
    gsn_panel(wks_eps,plot_eps(pa:pb),(/6,1/),resP)
    gsn_panel(wks_ps,plot_ps(pa:pb),(/6,1/),resP)
    gsn_panel(wks_x11,plot(pa:pb),(/6,1/),resP)
    pa = pa + 6
    pb = pb + 6
  end do

end
