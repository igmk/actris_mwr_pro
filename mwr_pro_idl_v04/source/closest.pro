;************
PRO CLOSEST,$
;************
;INPUT:
x,$             ;scalar
y,$             ;vector
;OUPUT:
dmin,$          ;smallest distance of x to y
ind             ;index of y corresponding to smallest distance
; $Id: closest.pro,v 1.1 2008/08/25 13:02:23 hatpro Exp $
; Abstract:
; * find closest absolute distance of vector y to scalar x
; * returns index of closest value
; Authors:
; U. Loehnert
; Date:
; 2008-08-11
; Dependencies:
; -
; Changes:
; -
;-

dmin = MIN(ABS(y-x))
ind = WHERE(dmin EQ ABS(y-x))
ind = ind(0) 

RETURN

END