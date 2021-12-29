FUNCTION DBZ_2_BYTE_T, z, zmax, zmin, offset=offset

 ncolor = 254
 a = ncolor/(zmax - zmin)
 b = a * zmin
 x = a*z - b + 1
 RETURN, x
END