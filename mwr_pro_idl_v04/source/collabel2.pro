PRO COLLABEL2, BarColor, BarLabel,                                             $
    BARTEXT=BarText, CENTERLAB=CenterLab, CHARSIZE=CharSize, COLOR=Color,      $
    DIRECTION=Direction, FONT=Font, LABFORMAT=LabFormat, LABOFFSET=LabOffset,  $
    POSITION=Position

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;   procedure : COLLABEL
;
;   purpose   : This procedure plots a color bar on the current window.
;
;   author    : Rolf Fuhrhop, IfM Kiel
;
;   modified  : Dirk Meetschen, MIUB Bonn
;
;   date      : June 1993
;
;   changed   : April 1998
;
;   parameter : ( input )
;
;             : BarColor      : BYTARR(*)   ; Byte array with color indices
;                                           ; to plot
;
;               BarLabel      : Array(2)    ; First and last label of bar.
;
;
;   keyword   : ( input )
;
;               BarText       : STRING      ; Text label for color bar
;                                           ; def.: ''
;
;               CenterLab     : BYTE        ; center bar label for color box
;
;               CharSize      : FLOAT       ; Character size of text label
;                                           ; ( "BarText" ).
;                                           ; def.: !P.CharSize
;
;               Color         : FLTARR(2)   ; Defines text and line color.
;                                           ; Color(0) -> text color
;                                           ; Color(1) -> line color
;                                           ; def.: [ !P.COLOR, !P.COLOR ]
;
;               Font          : BYTE        ; set font to hardware or software
;                                           ; -1 -> software
;                                           ; 0  -> hardware
;                                           ; def.: !P.FONT
;
;               Direction     : BYTE        ; Defines direction of color bar.
;                                           ; Direction = 0  -> horizontal
;                                           ; Direction = 1  -> vertical
;                                           ; def.: 0
;
;               LabFormat     : STRING      ; String format for color labels
;                                           ; def.: '(I3)'
;
;               LabOffset     : Integer     ; Offset for color bar labels
;                                           ; def.: 1
;               Position      : FLTARR(4)   ; normalized coordinates defining
;                                           ; color bar area.
;                                           ; def.: below current plot region,
;                                           ;       centered.
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
; set error exit
;

  ON_ERROR, 1

;
; check parameter
;

  IF N_PARAMS() LT 1 THEN MESSAGE, 'usage: COLLABEL, BarColor'

;
; get size of color array
;

  ColSize = SIZE( BarColor )
  NCol = ColSize(1)

;
; check keywords
;

  CASE N_ELEMENTS( Color ) of
     0 : Color = [ !P.COLOR, !P.COLOR ]
     1 : Color = [ Color, !P.COLOR ]
  ENDCASE

  IF N_ELEMENTS( CharSize ) EQ 0 THEN CharSize=!P.CharSize
  IF N_ELEMENTS( Font ) GT 0 THEN BEGIN
     FontSav = !P.FONT
     !P.FONT = Font
  ENDIF 
  IF N_ELEMENTS( LabFormat ) EQ 0 THEN LabFormat='(I3)'
  IF N_ELEMENTS( LabOffset ) EQ 0 THEN LabOffset=1

  IF NOT KEYWORD_SET( Direction ) THEN Direction = 0
  IF NOT KEYWORD_SET( Position ) THEN BEGIN
     IF Direction EQ 0 THEN BEGIN
        X = ( !X.WINDOW(1) - !X.WINDOW(0) ) *0.5 + !X.WINDOW(0)
        POSITION = [ X-0.3, !Y.WINDOW(0)-0.08, X + 0.3, !Y.WINDOW(0) - 0.05 ]
     ENDIF ELSE BEGIN
        Y = ( !Y.WINDOW(1) - !Y.WINDOW(0) ) *0.5 + !Y.WINDOW(0)
        POSITION = [ !X.WINDOW(0)- 0.07, Y-0.3, !X.WINDOW(0) - 0.04, Y + 0.3 ]
     ENDELSE
  ENDIF

;
; set data
;

  PSave = !P
  XSave = !X
  YSave = !Y

;
; draw box
;

  IF DIRECTION EQ 0 THEN BEGIN
     X = FINDGEN(NCol)
     Y = FLTARR(NCol)
     XR = [0, NCol]
     YR = [ 0, 1 ]
     XOff = 0
     If KEYWORD_SET( CenterLab ) THEN XOff = XOff + 0.5
     YOff = -1.0 
     Ali = 0.5
     IF KEYWORD_SET(LabFormat) THEN !X.Tickformat=LabFormat
     PLOT, XR, YR, /NOERASE, /NODATA, XStyle=9, YStyle=5,$
        POSITION=Position, XRange=BarLabel, YRange=YR, $
        COLOR=Color(1),XTicks=6,XTicklen=-0.08
  ENDIF ELSE BEGIN
     X = FLTARR(NCol)
     Y = FINDGEN(NCol)
     XR = [ 0, 1 ]
     YR = [0, NCol]
     XOff = -0.2
     YOff = -0.05 
     If KEYWORD_SET( CenterLab ) THEN YOff = YOff + 0.5
     Ali  = 1.0
     IF KEYWORD_SET(LabFormat) THEN !Y.Tickformat=LabFormat
     PLOT, XR, YR, /NOERASE, /NODATA, XStyle=5, YStyle=9,$
        POSITION=Position, XRange=XR, YRange=BarLabel, $
        COLOR=Color(1),YTicks=6,YTicklen=-0.08
  ENDELSE

;
; draw color boxes
;

      xpoly=[Position(0),Position(2),Position(2),Position(0)]
      ypoly=[Position(1),Position(1),Position(3),Position(3)]

      koord_norm=FLTARR(3,4)
      FOR i = 0,3 DO koord_norm(*,i)=[xpoly(i),ypoly(i),0]
;wave:
;      koord_dev=POLY_DEV(koord_norm)
;idl
      koord_dev=CONVERT_COORD(koord_norm,/NORM,/TO_DEVICE)
      POLYFILL,koord_dev(0,*)+[-1,1,1,-1],$
               koord_dev(1,*)+[-1,-1,1,1],Color=Color(1),/Device
      
      IF Direction EQ 0 THEN BEGIN
        xdiff=(Position(2)-Position(0))/FLOAT(NCol+1)
        ydiff=0.
        xpoly=[Position(0),Position(0)+xdiff,Position(0)+xdiff,Position(0)]
       ENDIF ELSE BEGIN
        xdiff=0.
        ydiff=(Position(3)-Position(1))/FLOAT(NCol+1)
        ypoly=[Position(1),Position(1),Position(1)+ydiff,Position(1)+ydiff]
       ENDELSE
      FOR i= 0, NCol-1 DO BEGIN
         xpoly=xpoly+[xdiff,xdiff,xdiff,xdiff]
         ypoly=ypoly+[ydiff,ydiff,ydiff,ydiff]
         POLYFILL,xpoly,ypoly,$
                  Color=BarColor(i),/Normal
       ENDFOR

;
; Title
;

  IF KEYWORD_SET( BarText ) THEN BEGIN
     IF DIRECTION EQ 0 THEN BEGIN
        X = ( Position(2) - Position(0) ) * 0.5 + Position(0)
        Y = Position(1) - 0.25 
        Ori=0.0
     ENDIF ELSE BEGIN
        X = Position(0) - 0.02
        Y = ( Position(3) - Position(1) ) * 0.5 + Position(1)
        Ori=90.0
     ENDELSE
     XYOUTS, X, Y, BarText, ALI=0.5, /NORM, ORIENTATION=Ori, SIZE=CharSize,    $
             COLOR=Color(0)
  ENDIF

;
; finished
;

  IF N_ELEMENTS( Font ) GT 0 THEN !P.FONT = FontSav 
  !P = PSave
  !X = XSave
  !Y = YSave 
  RETURN

END
