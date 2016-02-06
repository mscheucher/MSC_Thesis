;Plot of Initial conditions

path='MHD_Simulations/MHD/2d/KHI_Tests/EGU/study-egu/sinus/kx052/Bx=0deg/r01/v2/'

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

set_plot, 'ps'
loadct,39
device, file=path+'init_profiles.eps', /color, bits=8
DEVICE, DECOMPOSED = 0,/encapsulated

plot,y,init_rho[10,*],linestyle=0,xtitle='y/a',color=0,yrange=[0.0,2.0];, xtickn=[2,4,6,8,10]
oplot,y,init_bx[10,*],linestyle=0,color=60
oplot,y,init_bz[10,*],linestyle=0,color=150
oplot,y,init_vx[10,*],linestyle=0,color=30
oplot,y,init_p[10,*],linestyle=0,color=250

oplot,[500,550],[0.95,0.95],linestyle=0, color=0
oplot,[500,550],[0.9,0.9],linestyle=0, color=60
oplot,[500,550],[0.85,0.85],linestyle=0, color=150
oplot,[500,550],[0.8,0.8],linestyle=0, color=30
oplot,[500,550],[0.75,0.75],linestyle=0, color=250
xyouts,570,0.95,'density', color=0;, charsize=0.5
xyouts,570,0.9,'B!Dx!n', color=0;, charsize=0.5
xyouts,570,0.85,'B!Dz!n', color=0;, charsize=0.5
xyouts,570,0.8,'v!Dx!n', color=0;, charsize=0.5
xyouts,570,0.75,'thermal pressure', color=0;, charsize=0.5


device, /close

set_plot, 'x' 	;to switch back

print, '##### successfully printed #####'
;end
