PRO LoadCol, File
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;  Diese Routine liest die RGB Farbtabelle aus der Datei 'File' und
;  und setzt die neue Farbtabelle mit TVLCT.
;  Die Anzahl der verfuegbaren Farben ist abhaengig vom Window -
;  System ( DECWindows, MOTIV ). In der Datei werden die RGB Werte
;  fuer 256 Farben gelesen. Ist die tatsaechliche Farbtiefe geringer
;  dann werden entsprechend weiniger Farben gesetzt.
;  Die neu gesetzten Farben werden in der Farbtabelle von PW - Wave
;  ( COMMON colors ) gesetzt.
;
;  Parameter :
;
;     File     : String   : I ; Dateiname mit Farbtabelle
;
;  Common Bloecke
;                                PW - Wave Farbtabelle ( Seite 14-14 )
;     COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ON_ERROR,1                              ; return to main

; definiere Commonblock 

  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

; pruefe Parameter

  IF N_PARAMS() EQ 0 THEN BEGIN
     MESSAGE, 'filename for colortable not defined'
     PRINT, 'use: LOADCOL, filename'
     RETURN
  ENDIF

; pruefe ob Datei existiert

  FStr = FINDFILE(File,COUNT=NF)
  If NF EQ 0 THEN BEGIN
     MESSAGE, 'filename for colortable not found'
     RETURN
  ENDIF

; belege Speicher fuer Farbtabelle aus der Datei

  RR = BYTARR(256) 
  GR = BYTARR(256)
  BR = BYTARR(256)

; lese Farbtabelle aus der Datei

  GET_LUN, Unit                           ; unit fuer Datei
  OPENR, Unit, File                       ; oeffne Datei 
  READU, Unit, RR, GR, BR                 ; lese Farbtabelle
  FREE_LUN, Unit                          ; Freigabe der unit

; speichern der Farbtabelle

  TVLCT, RR, GR, BR 
; NCol = !D.N_COLORS - 1                  ; Farbtiefe - 1

; PW - Wave Farbtabelle definiert ?

  IF N_ELEMENTS(R_ORIG) EQ 0 THEN BEGIN
     R_ORIG = RR                          ; nein ->
     G_ORIG = GR                          ; setze original Werte
     B_ORIG = BR
  ENDIF 

; setze aktuelle Farbwerte

  R_CURR = RR
  G_CURR = GR
  B_CURR = BR
  
  RETURN                                  ; fertig
END; LOADCOL
