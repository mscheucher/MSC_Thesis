;draw growth rate tables

path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/'

;v05:
run1=[0.0436,0.0422,0.0353,0.0228]
run2=[0.0355,0.0321,0.0235,0.0082]
run3=[0.0623,0.0597,0.0567,0.0504]
run4=[0.0417,0.0396,0.0307,0.0143]

k=[0,5,10,15]
T = FINDGEN(100)/99*25


set_plot, 'ps'
; device, filename=path+'table_v05.eps',/color, bits=8
device, filename=path+'table_bx.eps',/color, bits=8
DEVICE, DECOMPOSED = 0,/encapsulated
loadct,39
plot, k, run1,psym=symcat(16), xtit='Inclination ' + string(greek("phi"))+' of B-field', ytit=string(greek("gamma"))+' a/v!Df!n', yrange=[0,0.09],xrange=[0,20],charsize=1.5, color=0
oplot,k,run2,psym=symcat(16), color=250
oplot,k,run3,psym=symcat(16), color=60
oplot,k,run4,psym=symcat(16), color=150
oplot, T, spline(k,run1,T,1), linestyle=0
oplot, T, spline(k,run2,T,1), linestyle=0, color=250
oplot, T, spline(k,run3,T,1), linestyle=0, color=60
oplot, T, spline(k,run4,T,1), linestyle=0, color=150

oplot,[8.5,9.5],[0.09,0.09],linestyle=0
oplot,[8.5,9.5],[0.085,0.085],linestyle=0, color=250
oplot,[8.5,9.5],[0.08,0.08],linestyle=0, color=60
oplot,[8.5,9.5],[0.075,0.075],linestyle=0, color=150

xyouts,10,0.09,string(greek("rho"))+'!D0!n = 0.1, k!Dx!na = 0.52,  v!D0!n = 0.5', color=0, charsize=1
xyouts,10,0.085,string(greek("rho"))+'!D0!n = 0.1, k!Dx!na = 0.82, v!D0!n = 0.5', color=0, charsize=1
xyouts,10,0.08,string(greek("rho"))+'!D0!n = 0.1, k!Dx!na = 0.52, v!D0!n = 0.25', color=0, charsize=1
xyouts,10,0.075,string(greek("rho"))+'!D0!n = 0.02, k!Dx!na = 0.52, v!D0!n = 0.5', color=0, charsize=1

; oplot,[200,210],[-10.5,-10.5],linestyle=0, color=30
; ; oplot,[200,210],[-11,-11],linestyle=0, color=80
; oplot,[200,210],[-11.5,-11.5],linestyle=0, color=130
; xyouts,220,-10.5," gamma="+string(fit(1)), color=30, charsize=0.5
; ; xyouts,220,-11,"only BL: gamma="+string(fit_bl(1)), color=80, charsize=0.5
; xyouts,220,-11.5," gamma_e="+string(fit_e(1)*0.5), color=130, charsize=0.5

device, /close
set_plot, 'x'

; end
