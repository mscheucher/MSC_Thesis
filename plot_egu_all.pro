;plot all data
pro plot_egu_all

for i=30,70,5 DO BEGIN
time=i
@2dplot_egu
ENDFOR
end
