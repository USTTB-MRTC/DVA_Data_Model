$PROBLEM
; Model: DVA for DHA  
; Based on: run000
; Description: DVA and PRR model as PD endpoint.
; Author: Mohamed MAIGA and  Sebastian Wicha
; Year: 2025



$INPUT      
	ID          ;conc scneario (full time course)          
	TIME        ;Time after experiment start [h]           
	AMT         ;concentration tested and dosed into CMT           
	CMT         ;1=PD DVA,2=PD PRR         
	EVID        ; 		
	DV          ;DVA: Viable parasite count, PRR: Log (Viable parasite + 1)          
	FLAG 	      ;1= DVA ;2= PRR
		       
	MDV         ;		
	REPLICATE   ;Replicate
	CA          ;Conc DHA 
	CB          ;Conc CQ	
	CC          ;Conc ATO
	CD          ;Conc PYRI	


			
$DATA   DVA_PRR.csv IGNORE=@ 
 ACCEPT=(ID.EQ.1)
; The data are structured as follows:
;; 3D7:
      ; ID.EQ.1 -------> DVA for DHA 
      ; ID.EQ.2 -------> DVA for CQ
      ; ID.EQ.3 -------> DVA for ATO
      ; ID.EQ.4 -------> DVA for PYRI  

      ; ID.EQ.5 -------> PRR for DHA 
      ; ID.EQ.6 -------> PRR for CQ
      ; ID.EQ.7 -------> PRR for ATO
      ; ID.EQ.8 -------> PRR for PYRI 

;; Field Isolates:
      ; ID.EQ.9 ------->  DVA for DHA 
      ; ID.EQ.10 -------> DVA for CQ
      ; ID.EQ.11 -------> DVA for ATO
      ; ID.EQ.12 -------> DVA for PYRI  




$SUBROUTINE ADVAN13 TOL=9

$MODEL      
 COMP=(PDMT)
 COMP=(PDPRR)


         
$PK
  N0 = THETA(1) * EXP(ETA(1))
  KG = THETA(2) * EXP(ETA(2))
  KK1 = THETA(3) * EXP(ETA(3))
  LAGTIME = THETA(4) * EXP(ETA(4))
  BSL = THETA(5) * EXP(ETA(5))
   


  A_0(1) = N0
  



$DES
   IF(T.LE.LAGTIME) THEN ; LAGTIME as ON-OFF
    ONOFF = 0  
   ELSE
    ONOFF = 1
   ENDIF 
   
   IF(A(1).GT.BSL) THEN 
    E = KK1 * ONOFF
   ELSE
	E = 0
   ENDIF
   
   DADT(1) = KG*A(1) - E * A(1) 
   DADT(2) = 0
   
   A1 = A(1)
 
  
  
  
$ERROR
  IPRED = A(1) 
  Y = IPRED * (1+EPS(1)) + EPS(2)


$THETA
(0, 304) ;1_N0 Geomean
(0) FIX ;2_KG
(0, 0.045) ;3_KK1
(0) FIX;4_LAGTIME
(0, 20)  ;5_BSL




$OMEGA 
 0 FIX  ;1_N0
 0 FIX  ;2_KG_
 0 FIX  ;3_KK1
 0 FIX  ;4_LAGTIME
 0 FIX  ;5_BSL




$SIGMA 
 0 FIX 
 0.1



;$SIMULATION (123) SUBPROBLEMS=1 



$ESTIMATION METHOD=1 LAPLACIAN INTER MAXEVAL=9999 NOABORT SIG=2 PRINT=1



$COVARIANCE MATRIX=R PRINT=E



$TABLE ID TIME DV CMT MDV EVID IPRED A1 KG N0 E KK1 CA CB CC CD LAGTIME
NPDE CWRES ONEHEADER NOPRINT FILE=sdtab001



$TABLE KG N0 KK1 FIRSTONLY
ONEHEADER NOPRINT FILE=patab001

