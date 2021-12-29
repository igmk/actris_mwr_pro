;#######################
PRO MAKE_BGYR_CT,$
;....input
colmin,$
colmax,$
;...output
r,$
g,$ 
b
;generates r,g,b arrays so that lower boundary is blue, upper boundary is red, and zero value of original data is always white

r=INTARR(256)
g=INTARR(256)
b=INTARR(256)

;set colors for lower half of colorbar
r(0)=0
g(0)=0
b(0)=0
r(255)=255
g(255)=255
b(255)=255
r(1)=0
g(1)=0
b(1)=255
r(64)=0
g(64)=255
b(64)=255
r(127)=0
g(127)=255
b(127)=0

colgrad_r=(r(64)-r(1))/63.
colgrad_g=(g(64)-g(1))/63.
colgrad_b=(b(64)-b(1))/63.

FOR i=1,63 DO BEGIN
r(i+1)=r(1)+i*colgrad_r
g(i+1)=g(1)+i*colgrad_g
b(i+1)=b(1)+i*colgrad_b
ENDFOR

colgrad_r=(r(127)-r(64))/63.
colgrad_g=(g(127)-g(64))/63.
colgrad_b=(b(127)-b(64))/63.

FOR i=1,63 DO BEGIN
r(i+64)=r(64)+i*colgrad_r
g(i+64)=g(64)+i*colgrad_g
b(i+64)=b(64)+i*colgrad_b
ENDFOR

;yellow
r(128)=255
g(128)=255
b(128)=0
;orange
r(191)=255
g(191)=128
b(191)=0
;red
r(254)=255
g(254)=0
b(254)=0

colgrad_r=(r(191)-r(128))/63.
colgrad_g=(g(191)-g(128))/63.
colgrad_b=(b(191)-b(128))/63.

FOR i=1,63 DO BEGIN
r(i+128)=r(128)+i*colgrad_r
g(i+128)=g(128)+i*colgrad_g
b(i+128)=b(128)+i*colgrad_b
ENDFOR

colgrad_r=(r(254)-r(191))/63.
colgrad_g=(g(254)-g(191))/63.
colgrad_b=(b(254)-b(191))/63.

FOR i=1,63 DO BEGIN
r(i+191)=r(191)+i*colgrad_r
g(i+191)=g(191)+i*colgrad_g
b(i+191)=b(191)+i*colgrad_b
ENDFOR

r(254)=211
g(254)=211
b(254)=211
END