FUNCTION DBZ_2_BYTE, x, xmax, xmin, fac
 ncolor = 254.
 m = (ncolor/(xmax - xmin))*fac; Steigung
; b = -1.*(m*xmin)
 b = (ncolor/2.) - m*(xmax+xmin)/2.
 y = m*x + b
 RETURN, y
END
