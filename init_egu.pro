;Plot of Initial conditions

path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/check_kx052/sinus/kx052/Bx=90deg/r01/'

init_rho=read_ascii(path+'data/initial_MHD_u0.txt')
init_rho=transpose(init_rho.FIELD001)
init_bx=read_ascii(path+'data/initial_MHD_u5.txt')
init_bx=transpose(init_bx.FIELD001)
init_bz=read_ascii(path+'data/initial_MHD_u7.txt')
init_bz=transpose(init_bz.FIELD001)
init_vx=read_ascii(path+'data/initial_MHD_w0.txt')
init_vx=transpose(init_vx.FIELD001)
init_vy=read_ascii(path+'data/initial_MHD_w1.txt')
init_vy=transpose(init_vy.FIELD001)
init_p=read_ascii(path+'data/initial_MHD_w3.txt')
init_p=transpose(init_p.FIELD001)

y=indgen(800)

loadct,3
!p.multi=[0,6,2]
plot,y,init_rho[10,*],psym=1,ytitle='Rho',xtitle='y'
plot,y,init_bx[10,*],psym=1,ytitle='Bx',xtitle='y'
plot,y,init_bz[10,*],psym=1,ytitle='Bz',xtitle='y'
plot,y,init_vx[10,*],psym=1,ytitle='vx',xtitle='y'
plot,y,init_vy[10,*],psym=1,ytitle='vy',xtitle='y'
plot,y,init_p[10,*],psym=1,ytitle='P',xtitle='y'
plot_image,init_rho,charsize=3,title='Rho'
plot_image,init_bx,charsize=3,title='Bx'
plot_image,init_bz,charsize=3,title='Bz'
plot_image,init_vx,charsize=3,title='vx'
plot_image,init_vy,charsize=3,title='vy'
plot_image,init_p,charsize=3,title='p'
!p.multi=0

set_plot, 'ps'
loadct,33

device, file=path+'plot/init_rho.ps', /color, bits=8
DEVICE, DECOMPOSED = 0
plot_image,init_rho,charsize=1.0
COLOR_BAR,init_rho,0.07,0.12,0.11,0.95,/normal	;stehend
; COLOR_BAR,data_r,0.12,0.96,0.1,0.15,/normal	;liegend
xyouts, 0.4, 0.97, 'KHI Test initial Rho',/normal
device, /close
device, file=path+'plot/init_Bx.ps', /color, bits=8
DEVICE, DECOMPOSED = 0
plot_image,init_bx,charsize=1.0
COLOR_BAR,init_bx,0.07,0.12,0.11,0.95,/normal	;stehend
; COLOR_BAR,data_r,0.12,0.96,0.1,0.15,/normal	;liegend
xyouts, 0.4, 0.97, 'KHI Test initial Bx',/normal
device, /close
device, file=path+'plot/init_Bz.ps', /color, bits=8
DEVICE, DECOMPOSED = 0
plot_image,init_bz,charsize=1.0
COLOR_BAR,init_bz,0.07,0.12,0.11,0.95,/normal	;stehend
; COLOR_BAR,data_r,0.12,0.96,0.1,0.15,/normal	;liegend
xyouts, 0.4, 0.97, 'KHI Test initial Bz',/normal
device, /close
device, file=path+'plot/init_vx.ps', /color, bits=8
DEVICE, DECOMPOSED = 0
plot_image,init_vx,charsize=1.0
COLOR_BAR,init_vx,0.07,0.12,0.11,0.95,/normal	;stehend
; COLOR_BAR,data_r,0.12,0.96,0.1,0.15,/normal	;liegend
xyouts, 0.4, 0.97, 'KHI Test initial vx',/normal
device, /close
device, file=path+'plot/init_vy.ps', /color, bits=8
DEVICE, DECOMPOSED = 0
plot_image,init_vy,charsize=1.0
COLOR_BAR,init_vy,0.07,0.12,0.11,0.95,/normal	;stehend
; COLOR_BAR,data_r,0.12,0.96,0.1,0.15,/normal	;liegend
xyouts, 0.4, 0.97, 'KHI Test initial vy',/normal
device, /close
device, file=path+'plot/init_p.ps', /color, bits=8
DEVICE, DECOMPOSED = 0
plot_image,init_p,charsize=1.0
COLOR_BAR,init_p,0.07,0.12,0.11,0.95,/normal	;stehend
; COLOR_BAR,data_r,0.12,0.96,0.1,0.15,/normal	;liegend
xyouts, 0.4, 0.97, 'KHI Test initial p',/normal
device, /close

set_plot, 'x' 	;to switch back

print, '##### successfully printed #####'
;end