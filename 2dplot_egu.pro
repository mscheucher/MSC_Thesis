;Plot the output data from "MHD_2d"

; path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/study/r05v025/'
; path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/check_kx052/random/r002/'
path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/study-egu/sinus/kx052/Bx=0deg/r01/v2/'


; time=5
time=strcompress(time,/remove)

data_i=read_ascii(path+'data/data_MHD_u0_t'+String(time)+'.*')
; data_i=read_ascii(path+'data/data_MHD_w1_t'+String(time)+'.*')
; data_i=read_ascii(path+'data/data_MHD_u1_t'+String(time)+'.*')
; data_i=read_ascii(path+'data/initial_MHD_w1*')
data_i=transpose(data_i.FIELD001)
data_i=rebin(data_i,n_elements(data_i[*,0])*2,n_elements(data_i[0,*]))

set_plot, 'ps'

device, file=path+'plot/plot_KHI_test_t'+String(time)+'.ps', /color, bits=8
; device, file=path+'plot/plot_KHI_test_t'+String(time)+'_vy.ps', /color, bits=8
; device, file=path+'plot/plot_KHI_test_t'+String(time)+'_rvx.ps', /color, bits=8
; device, file=path+'plot/plot_KHI_test_init_vy.ps', /color, bits=8
DEVICE, DECOMPOSED = 0
loadct,33

plot_image,data_i,charsize=1.5,title=String('t='+String(time)),xticks=2,xtickv=[0,200,400],xtitle='x',ytitle='density'
; contour,data_i,LEVELS = 0.5 + FINDGEN(15) * 0.4,c_color=indgen(15)*5*3.4+1;,/overplot
COLOR_BAR,data_i,0.20,0.25,0.2,0.9,/normal	;stehend
; COLOR_BAR,data_i,0.43,0.76,0.07,0.11,/normal	;liegend
; xyouts,0.45, 0.97, String('t='+String(time)),/normal ; 
; xyouts, 0.25, 0.03, String('15 density contours (0.5-6.5)'),/normal
!p.multi=0
device, /close

; device,file='MHD_Simulations/MHD/2d/KHI_Tests/EGU/plot/plot_KHI_test_profile.ps' 
; plot,y,data_i[10,*],psym=0,ytitle='Density',title='Rho',xtitle='y = 10'
; device, /close
set_plot, 'x' 	;to switch back


print, '##### successfully compiled #####'
; end
