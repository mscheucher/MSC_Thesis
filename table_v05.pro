;draw growth rate tables

path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/'

;v05:
r002=[0.0238,0.0356,0.0417,0.0437,0.0397]
r005=[0.0272,0.0384,0.0425,0.0436,0.0384]
r01=[0.0311,0.0408,0.0445,0.0432,0.0355]
r05=[0.0349,0.0409,0.0393,0.0317,0.0186]

;v025:
; r002=[0.0346,0.0517,0.0612,0.0615,0.0549]
; r005=[0.0400,0.0552,0.0625,0.0608,0.0504]
; r01=[0.0446,0.0578,0.0623,0.0580,0.0456]
; r05=[0.0415,0.0484,0.0437,0.0277,0.0]

k=[0.22,0.37,0.52,0.67,0.82]
T = FINDGEN(100)/99*0.85 + 0.1


set_plot, 'ps'
device, filename=path+'table_v05.eps',/color, bits=8
; device, filename=path+'table_v025.eps',/color, bits=8
DEVICE, DECOMPOSED = 0,/encapsulated
loadct,39
plot, k, r002,psym=symcat(16), xtit='k!Dx!n a', ytit=string(greek("gamma"))+' a/v!Df!n', yrange=[0,0.05],xrange=[0.,1.],charsize=1.5, color=0
oplot,k,r005,psym=symcat(16), color=250
oplot,k,r01,psym=symcat(16), color=60
oplot,k,r05,psym=symcat(16), color=150
oplot, T, spline(k,r002,T,1), linestyle=0
oplot, T, spline(k,r005,T,1), linestyle=0, color=250
oplot, T, spline(k,r01,T,1), linestyle=0, color=60
oplot, T, spline(k,r05,T,1), linestyle=0, color=150

oplot,[0.3,0.35],[0.015,0.015],linestyle=0
oplot,[0.3,0.35],[0.0125,0.0125],linestyle=0, color=250
oplot,[0.3,0.35],[0.01,0.01],linestyle=0, color=60
oplot,[0.3,0.35],[0.0075,0.0075],linestyle=0, color=150

xyouts,0.4,0.015,string(greek("rho"))+'!D0!n = 0.02', color=0, charsize=1
xyouts,0.4,0.0125,string(greek("rho"))+'!D0!n = 0.05', color=0, charsize=1
xyouts,0.4,0.01,string(greek("rho"))+'!D0!n = 0.1', color=0, charsize=1
xyouts,0.4,0.0075,string(greek("rho"))+'!D0!n = 0.5', color=0, charsize=1

; oplot,[200,210],[-10.5,-10.5],linestyle=0, color=30
; ; oplot,[200,210],[-11,-11],linestyle=0, color=80
; oplot,[200,210],[-11.5,-11.5],linestyle=0, color=130
; xyouts,220,-10.5," gamma="+string(fit(1)), color=30, charsize=0.5
; ; xyouts,220,-11,"only BL: gamma="+string(fit_bl(1)), color=80, charsize=0.5
; xyouts,220,-11.5," gamma_e="+string(fit_e(1)*0.5), color=130, charsize=0.5

device, /close
set_plot, 'x'

; end
