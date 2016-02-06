;get growth rates

; path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/study/r05v025/'
; path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/check_kx052/random/r05/'
path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/study-egu/sinus/kx052/Bx=0deg/r01/v2/'
; path='MHD_Simulations/MHD/2d/KHI_Tests/Test02/kx052/r05/'

t0=20
t1=40

e=0.0001
vy=0.01
;vy_bl=0.01

for i=5,75,5 DO BEGIN
time=strcompress(i,/remove)
data=read_ascii(path+'data/data_MHD_w1_t'+String(time)+'.*')
data=transpose(data.FIELD001)
data2=read_ascii(path+'data/data_MHD_u0_t'+String(time)+'.*')
data2=transpose(data2.FIELD001)

emap=data*data*data2/2
e=[e,max(abs(emap[*,*]))]

vy=[vy,max(abs(data[*,*]))]
;vy_bl=[vy_bl,max(abs(data[*,115:135]))]

ENDFOR

save,e,vy,vy_bl, FILENAME = path+'e+vy.sav'

t=indgen(n_elements(vy))*5
vy=alog(vy)
fit = poly_fit(t(0:(t1-t0)*0.2),vy(t0*0.2:t1*0.2),1)
;vy_bl=alog(vy_bl)
;fit_bl = poly_fit(t(0:(t1-t0)*0.2),vy_bl(t0*0.2:t1*0.2),1)
e=alog(e)
fit_e = poly_fit(t(0:(t1-t0)*0.2),e(t0*0.2:t1*0.2),1)


set_plot, 'ps'
device, filename=path+'growth_rates_egu.eps',/color, bits=8
DEVICE, DECOMPOSED = 0,/encapsulated
loadct,39
plot, t, e,linestyle=0, xtit='t v!Df!n/a', ytit='ln (E!Dy!n)', yrange=[-14,0],charsize=1.5;color=30, xrange=[0,300]
; oplot, t, vy_bl,linestyle=0, color=80
; oplot, t, vy,linestyle=0, color=130
; oplot, t+t0, fit(0)+fit(1)*t, linestyle=1, color=130
; oplot, t+t0, fit_bl(0)+fit_bl(1)*t, linestyle=1, color=80
oplot, t+t0, fit_e(0)+fit_e(1)*t, linestyle=2, color=60,thick=2

;oplot,[75,75],[-15,e(15)],linestyle=2, color=250
;oplot,[130,130],[-15,e(26)],linestyle=2, color=250
;oplot,[75,75],[e(15),e(15)],psym=3, color=250,thick=3
;oplot,[130,130],[e(26),e(26)],psym=3, color=250,thick=3

; oplot,[200,210],[-10.5,-10.5],linestyle=0, color=30
; ; oplot,[200,210],[-11,-11],linestyle=0, color=80
; oplot,[200,210],[-11.5,-11.5],linestyle=0, color=130
; xyouts,220,-10.5," gamma="+string(fit(1)), color=30, charsize=0.5
; ; xyouts,220,-11,"only BL: gamma="+string(fit_bl(1)), color=80, charsize=0.5
xyouts,220,-11.5," gamma_e="+string(fit_e(1)*0.5), color=130, charsize=0.5

device, /close
set_plot, 'x'

print, 'gamma= ', fit(1)
; print, 'gamma_bl= ', fit_bl(1)
print, 'gamma_e= ', fit_e(1)*0.5
end
