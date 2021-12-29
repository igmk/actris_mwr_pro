FUNCTION AVERAGE, ARRAY, DIMENSION, MISSING=MISSING
;+
; Project     : SOHO - CDS
;
; Name        : 
;	AVERAGE()
; Purpose     : 
;	Averages an array over one or all of its dimensions.
; Explanation : 
;	Calculates the average value of an array, or calculates the average
;	value over one dimension of an array as a function of all the other
;	dimensions.
; Use         : 
;	Result = AVERAGE( ARRAY )
;	Result = AVERAGE( ARRAY, DIMENSION )
; Inputs      : 
;	ARRAY	  = Input array.  May be any type except string or structure.
; Opt. Inputs : 
;	DIMENSION = Optional dimension to do average over.  Valid inputs are 1
;		    through the total number of dimensions of ARRAY.
; Outputs     : 
;	The average value of the array when called with one parameter.
;
;	If DIMENSION is passed, then the result is an array with all the
;	dimensions of the input array except for the dimension specified,
;	each element of which is the average of the corresponding vector
;	in the input array.
;
;	For example, if A is an array with dimensions of (3,4,5), then the
;	command B = AVERAGE(A,2) is equivalent to
;
;			B = FLTARR(3,5)
;			FOR J = 0,4 DO BEGIN
;				FOR I = 0,2 DO BEGIN
;					B(I,J) = TOTAL( A(I,*,J) ) / 4.
;				ENDFOR
;			ENDFOR
;
; Opt. Outputs: 
;	None.
; Keywords    : 
;	MISSING	= Value signifying missing pixels.  Any pixels with this value
;		  are not included in the average.  If there are no non-missing
;		  pixels, then MISSING is returned.
; Calls       : 
;	None.
; Common      : 
;	None.
; Restrictions: 
;	The dimension specified must be valid for the array passed.
; Side effects: 
;	None.
; Category    : 
;	Utilities, Arrays.
; Prev. Hist. : 
;	Taken from an earlier routine by W. Thompson called AVG, but the
;	definition of the DIMENSION parameter is different to be consistent
;	with current usage in IDL.
; Written     : 
;	William Thompson, GSFC, 9 April 1993.
; Modified    : 
;	Version 1, William Thompson, GSFC, 9 April 1993.
;	Version 2, William Thompson, GSFC, 3 February 1996
;		Added missing keyword.
; Version     : 
;	Version 2, 3 February 1996
;-
;
	ON_ERROR,2
;
;  Check the input parameters.
;
	IF N_PARAMS() LT 1 THEN MESSAGE,	$
		'Syntax:  Result = AVERAGE( ARRAY  [, DIMENSION ] )'
	IF N_ELEMENTS(ARRAY) EQ 0 THEN MESSAGE,	$
		'ARRAY not defined'
;
;  Dimension not passed.  Return a simple average.  If the keyword MISSING was
;  sent, then only average together the points not equal to the missing value.
;  If there are no non-missing pixels, then return the missing value.
;
	IF N_PARAMS(0) EQ 1 THEN BEGIN
	    IF N_ELEMENTS(MISSING) EQ 1 THEN BEGIN
		W = WHERE(ARRAY NE MISSING, COUNT)
		IF COUNT GT 0 THEN AVER = TOTAL(ARRAY(W)) / COUNT ELSE	$
			AVER = MISSING
	    END ELSE AVER = TOTAL(ARRAY) / N_ELEMENTS(ARRAY)
;
;  Dimension passed.  Check DIMENSION, and make sure that ARRAY is an array.
;
	END ELSE BEGIN
		IF N_ELEMENTS(DIMENSION) EQ 0 THEN BEGIN
			MESSAGE,'DIMENSION not defined'
		END ELSE IF N_ELEMENTS(DIMENSION) NE 1 THEN BEGIN
			MESSAGE,'DIMENSION must be a scalar'
		ENDIF
		S = SIZE(ARRAY)
		IF S(0) EQ 0 THEN MESSAGE,'ARRAY must be an array'
;
;  Return an array collapsed along one of the dimensions.  If the keyword
;  MISSING was passed, then 
;
		IF (DIMENSION GE 1) AND (DIMENSION LE S(0)) THEN BEGIN
		    IF N_ELEMENTS(MISSING) EQ 1 THEN BEGIN
;
;  Start by calculating the numerator, substituting 0 where-ever the missing
;  pixel flag is seen.
;
			AVER = ARRAY
			W = WHERE(ARRAY EQ MISSING, COUNT)
			IF COUNT GT 0 THEN AVER(W) = 0
			AVER  = TOTAL(AVER, DIMENSION)
;
;  Next calculate the denominator as the total number of points which are good.
;  Substitute the MISSING pixel value where-ever there are no good pixels to
;  average together.
;
			DENOM = TOTAL(ARRAY NE MISSING, DIMENSION)
			AVER = TEMPORARY(AVER) / (DENOM > 1)
			W = WHERE(DENOM EQ 0, COUNT)
			IF COUNT GT 0 THEN AVER(W) = MISSING
;
;  Otherwise, simply divide the total by the number of pixels along that
;  dimension.
;
		    END ELSE AVER = TOTAL(ARRAY,DIMENSION) / S(DIMENSION)
		END ELSE BEGIN
			MESSAGE,'Dimension out of range'
		ENDELSE
	ENDELSE
;
	RETURN, AVER
	END
