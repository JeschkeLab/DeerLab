#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	355
XMIN	-0.000000
XWID	2832.000000
*
*	Documentational Text:
*
TITL	'DEER_Ringtest2_sampleA_3us_20sc'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'unknown'
IRUNI	''
IIUNI	''
XUNI	'ns'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    epruser
DATE    11/02/20
TIME    12:10:05
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    PLS
OXS1    TADC
AXS1    B0VL
AXS2    NONE
AXS3    
MWPW    0.000
A1CT    1.22807
B0VL    1.22807
A1SW    0.02
MWFQ    0.0
AVGS    20
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              0.000 mW
PowerAtten         60 dB

.DVC     endor, 1.0

EIEENDORFreq       14.9021800 MHz/3.5 kG
EIEIsotope         H1
EIERFSweepDir      Same
EIEStaticField     12140.000 G
EIEStaticRF        38.000 MHz
ENDORType          EIF
RF1Atten           60.000 dB
RF1FreqPos         38.000 MHz
RF1StartFreq       38.000 MHz
RF1SweepWidth      30.000 MHz
RF2Atten           60.000 dB
RF2FreqPos         38.000 MHz
RF2StartFreq       38.000 MHz
RF2SweepWidth      30.000 MHz
RFSrcMixing        Add
SumAtten           10.000 dB
SumAttenStart      10.000 dB
SumAttenWidth      70.000 dB

.DVC     fieldCtrl, 1.0

AllegroMode        False
CenterField        12280.70 G
Delay              0.0 s
FieldFlyback       On
FieldPosition      12280.700 G
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     12280.700 G
SweepDirection     Up
SweepWidth         200.0 G
WidthTM            200.0 G

.DVC     freqCounter, 1.0

FrequencyMon       0.0 GHz
QMonitBridge       Off

.DVC     ftBridge, 1.0


.DVC     ftEpr, 1.0

AutoTimeOut        True
AveragesPerScan    5
ELDORFreqStart     9.300000 GHz
ELDORFreqWidth     0.700000 GHz
FTAcqModeSlct      Run from PulseSPEL
FieldIsStatic      False
GradIntPulse       False
GrdEnable          False
LastXAxis          ?
LastYAxis          ?
PCycleAllowed      True
PCycleOn           False
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     DEER
PlsSPELGlbTxt      \
;     PulseSPEL variables for DEER experiments\n\
;\n\
begin defs\n\
;                                               \n\
;\n\
;optimized for home-made Q band F238, May2010, attention! pulse lengths can be different\n\
;Gunnar Jeschke 2002, Polyhach-Bordignon 2010\n\
; Repetition time\n\
;\n\
srt  = 5000 * srtu ; shot repetition time, srtu= 1 µs\n\
;                 ; usual values for nitroxides: 1000 * srtu at 80 K\n\
                   ; 4000 * srtu at 50 K\n\
;                  ; 10000...40000 * srtu at 15 K\n\
;\n\
; Pulse lengths\n\
;                        \n\
p0   = 24      ; observer pi pulse length (frq1), 12 ns for hard pulses\n\
p1   = 12      ; observer pi/2 pulse length (frq1), 12 ns for hard pulses\n\
p2   = 12      ; pump pi pulse length (frq2), 12 ns for hard pulses\n\
;\n\
; Variable delays\n\
; \n\
d0   = 1220    ; 1344 field-sweep Q band home-made delay, 1312 for DEER with all40 ns\n\
\n\
d1   = 400     ; tau1 (initial value)\n\
d2   = 3000    ; tau2\n\
d3   = 280     ; deadtime delay during which no data are acquired\n\
               ; d3 >= 360, if TWT window must be split (single gate time exceeded)\n\
;\n\
; Time increments\n\
;\n\
d30  = 8       ; t time scale increment of dipolar time evolution\n\
d31  = 16      ; time scale increment for nuclear modulation averaging\n\
d29  = 4       ; time increment for pulse length in nutation\n\
;\n\
; Counters\n\
;\n\
h    = 5       ; number of shots/loop                     ( counter: I )\n\
n    = 20      ; number of sweeps to accumulate           ( counter: K )\n\
m    = 8       ; number of tau1 averages for nuclear modulation\n\
;                suppression\n\
;\n\
; ### Block of additional variables (Bruker definitions) ###\n\
p3   = 40      ;\n\
p7   = 80      ; LeCroy trigger\n\
p9   = 1000    \n\
pg   = 24      ; Integrator Gate Width (is superseeded by DEER pulse program)\n\
;\n\
;\n\
d4   = 2000    ;\n\
d5   = 400     ;\n\
                        ;\n\
dx   = 0       ; t time scale starting value\n\
dy   = 0       ; t1 time scale starting value\n\
               ;\n\
\n\
d20  = 0       ; Initial value of dx\n\
d21  = 0       ; Initial value of dy\n\
                        ;\n\
s    = 300     ; Sweep length (n. of data really taken)   ( counter: X )\n\
t    = 1       ; second time axsis sweep length           ( counter: Y )\n\
              ;\n\
r    = 1       ;\n\
c    = 1       ;\n\
b    = 1       ;\n\
             ;\n\
w    = 180     ; Split point of sweeps (ns) ,i.e. pulse separation where\n\
             ; TWT gate pulse can be split.\n\
;\n\
d6   = 100 us  \n\
end defs\n\
\n
PlsSPELLISTSlct    DEER run
PlsSPELPhPrgEx     Normal
PlsSPELPrg         mayu/deer4p_Q_NO_NO_2010_nobra.exp
PlsSPELPrgTxt      \
;\n\
; DEER setup and run \n\
; G. Jeschke, 2000-2009\n\
;\n\
; ### This assumes that +y (pi), +<x> (+pi/2), -x (-pi/2) are observer frequency channels while\n\
;     -y is a pump frequency (second frequency channel home-made Q band F238)\n\
     \n\
begin defs\n\
 dim s[512]               ;     dimension of data trace [sx] for setup\n\
 dim1 s[512]              ;     dimension [sx] for setup\n\
 dim2 s[512]              ;     dimension [sx] for set-up\n\
 dim3 s[1024]             ;     dimension [sx] for field sweep\n\
 dim4 s[355]              ;     dimension [sx] for DEER\n\
;                               ### must be adapted to delays!\n\
;                               sx =(d1+d2-2*d3)/d30\n\
 dim5 s[193]              ;     dimension for standing DEER\n\
 dim6 s[660,14]           ;     dimension [sx] for field-averaged DEER\n\
 dim7 s[24]               ;     dimension [sx] for nutation\n\
 dim8 s[280]              ;     dimension [sx] for -y eldor nutation\n\
   \n\
end defs\n\
;\n\
; Pulse lengths (all custom at home-made Q-band):\n\
;\n\
; p0    observer pi pulse (frq1), recommendation for nitroxides:        12 ns\n\
; p1    observer pi/2 pulse (frq1), recommendation for nitroxides:      12 ns\n\
; p2    pump pi pulse (frq2), recommendation for nitroxides:            12 ns\n\
;       at a frequency difference of frq1-frq2= 60... 70 MHz \n\
;\n\
; Variable delays:\n\
;\n\
; d0    acquisition trigger,    observation of defence pulses:  0 ns\n\
;                               pulse setup:                    440 ns\n\
;                               field sweep:                    440 ns\n\
;                               DEER runs:                      484 ns\n\
;                               ### depends on video bandwidth!\n\
;                               ### 50 MHz video bandwidth recommended\n\
;                               ### may vary between spectrometers!\n\
; d1    tau1,                   absolute minimum:               80 ns\n\
;                               recommended minimum:            120 ns\n\
;                               recommended standard value:     200 ns\n\
;                               ### short T2 -> minimum d1\n\
;                               ### long deadtime -> longer d1\n\
; d2    tau2,                   typical values:                 1-6 ï¿½s\n\
;                               ### short T2 -> short d2\n\
;                               ### high concentration -> short d2\n\
;                               ### long distances -> long d2\n\
;                               recommendation: increase d2 until echo\n\
;                               amplitude is half of maximum value in\n\
;                               SpecJet window at full receiver gain (66 dB)\n\
;                               for overnight measurements, at 60 dB\n\
;                               receiver gain for measurement times of 1 h,\n\
;                               at 54 dB for measurement times of 15 min\n\
; d3    deadtime,               window at the beginning and end of the\n\
;                               second observer interpulse delay where no\n\
;                               data are acquired, set to 80 ns\n\
;                               ### do not change without need\n\
;\n\
; Variable increments:\n\
;\n\
; d30                           increment for dipolar time evolution, typically 8 ns\n\
; d31                           increment for nuclear modulation averaging,\n\
;                               typically 8 ns\n\
;\n\
;\n\
; Fixed times:\n\
;\n\
; pg    integration gate,       field sweeps:   200 ns\n\
;                               DEER:           p0 (pi pulse length)\n\
; d9    DAF,                    (delay after flash), delay of experiment with\n\
;                               respect to external trigger\n\
;\n\
; Other variables used (definition file settings are superseeded)\n\
;\n\
; p10, p11, p12, d4, d10, d13, d21, d23, d24\n\
; d25, d26 for sensDEER\n\
;\n\
;\n\
;\n\
begin lists "+x (pi/2)"\n\
 ph1 +x\n\
 ph2 +x\n\
 ph3 +x\n\
 asg1 +a\n\
 bsg1 +b\n\
end lists \n\
\n\
\n\
begin lists1 "-x (pi/2)"\n\
 ph1 -x   \n\
 ph2 +x\n\
 ph3 +x                  \n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists1\n\
\n\
begin lists2 "+y (pi frq1)"\n\
 ph1 +x   \n\
 ph2 +y\n\
 ph3 +<x>                  \n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists2\n\
 \n\
begin lists3 "-<y> not used(pi frq1)"\n\
 ph1 -<y>   \n\
 ph2 +x\n\
 ph3 +x                  \n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists3\n\
 \n\
begin lists4 "Field sweep +x/-x (pi frq1)"\n\
 ph1 +x -x  \n\
 ph2 +x\n\
 ph3 +x                  \n\
 asg1 +a -a                   \n\
 bsg1 +b -b \n\
end lists4 \n\
\n\
begin lists5 "-y eldor(pi frq2)"\n\
 ph1 -y  \n\
 ph2 +x\n\
 ph3 +x                  \n\
 asg1 +a                  \n\
 bsg1 +b \n\
end lists5\n\
 \n\
begin lists6 "DEER setup not used (only observer pulses)"\n\
 ph1 -x   \n\
 ph2 +x\n\
 ph3 *                 \n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists6\n\
 \n\
begin lists7 "DEER run"\n\
 ph1 +x -x\n\
 ph2 +y\n\
 ph3 -y                 \n\
 asg1 +a -a                  \n\
 bsg1 +b -b  \n\
end lists7                  \n\
 \n\
begin lists8 "-y eldor nutation"\n\
 ph1 -y \n\
 ph2 +x -x\n\
 ph3 +x        \n\
 asg1 +a -a                  \n\
 bsg1 +b -b\n\
end lists8\n\
\n\
begin lists9 "+y nutation"\n\
 ph1 +y \n\
 ph2 +x\n\
 ph3 +x        \n\
 asg1 +a                  \n\
 bsg1 +b\n\
end lists9\n\
\n\
begin lists10 "DEER run Q"\n\
 ph1 +x -x   \n\
 ph2 +x\n\
 ph3 -y                 \n\
 asg1 +a -a                  \n\
 bsg1 +b -b  \n\
end lists10         \n\
\n\
;\n\
;\n\
; Standing 2-pulse echo for setup of observer pi pulses\n\
;\n\
\n\
begin exp "pi tp=p0 +y" [TRANS QUAD]         ; QUAD detection\n\
\n\
p11=p0\n\
p10=p0\n\
   \n\
    shot i=1 to h              ; accumulation loop\n\
      d9                       ; DAF\n\
      p10 [ph1]                ; 1st pulse and phase program\n\
      d1                       ; tau                       \n\
      p11 [ph2]                ; 2nd pulse and phase program\n\
      d1                       ; tau\n\
      d0                       ; constant acquisition delay\n\
      dig [sg1]                ; acquisition\n\
    next i                     ; end of accumulation loop \n\
   \n\
end exp   \n\
\n\
;\n\
;\n\
; Standing 2-pulse echo for setup of observer pi/2 pulse \n\
;\n\
\n\
begin exp1 "pi/2 tp=p1 +x/-x" [TRANS QUAD]         ; QUAD detection\n\
\n\
p11=2*p1\n\
p10=p1\n\
   \n\
    shot i=1 to h              ; accumulation loop\n\
      d9                       ; DAF\n\
      p10 [ph1]                ; 1st pulse and phase program\n\
      d1                       ; tau                       \n\
      p11 [ph2]                ; 2nd pulse and phase program\n\
      d1                       ; tau\n\
      d0                       ; constant acquisition delay\n\
      dig [sg1]                ; acquisition\n\
    next i                     ; end of accumulation loop \n\
   \n\
end exp1   \n\
\n\
;\n\
;\n\
; Standing inversion recovery for setup of pump pi pulse \n\
\n\
\n\
begin exp2 "pi tp=p2 -y eldor frq2" [TRANS QUAD]         ; QUAD detection\n\
\n\
p11=2*p1\n\
p10=p1\n\
p12=p2\n\
d13=400\n\
d10=d0\n\
dx=0\n\
    shot i=1 to h              ; accumulation loop\n\
      d9                       ; DAF\n\
      p12 [ph1]\n\
      d13\n\
      p10 [ph2]                ; 1st pulse and phase program\n\
      d1                       ; tau                       \n\
      p11 [ph3]                ; 2nd pulse and phase program\n\
      d1                       ; tau\n\
      d10                      ; constant acquisition delay\n\
      dig [sg1]                ; acquisition\n\
    next i                     ; end of accumulation loop\n\
   \n\
end exp2   \n\
\n\
\n\
; 2-Pulse Field Sweep\n\
\n\
\n\
begin exp3 "Field Sweep +x/-x pg=200ns" [INTG QUAD]        ; QUAD detection\n\
\n\
pg=200                         ; integration gate set to 200\n\
p11=2*p1\n\
p10=p1\n\
\n\
 for k=1 to n                  ; averaging loop\n\
\n\
  totscans(n)                  ; output of total number of scans\n\
\n\
  bsweep x=1 to sx             ; sweep loop\n\
    shot i=1 to h              ; accumulation loop\n\
      d9                       ; DAF\n\
      p10 [ph1]                 ; 1st pulse and phase program\n\
      d1                       ; tau                       \n\
      p11 [ph2]                 ; 2nd pulse and phase program\n\
      d1                       ; tau\n\
      d0                       ; constant acquisition delay\n\
      acq [sg1]                ; acquisition\n\
    next i                     ; end of accumulation loop \n\
   next x                      ; end of sweep loop\n\
\n\
  scansdone(k)                 ; output of scans done\n\
\n\
 next k\n\
\n\
end exp3 \n\
\n\
\n\
; Four-pulse DEER\n\
\n\
\n\
begin exp4 "DEER" [INTG QUAD]\n\
\n\
 totscans(n)                   ; output of total number of scans\n\
 for k=1 to n                  ; averaging loop\n\
\n\
 pg=p0;                        ; integration gate, standard 32 ns\n\
d25=200-pg                    ; trigger correction for shorter gate, not is use, remember to change d0 by hand \n\
d25=d25/2\n\
d25=d25-p0\n\
d25=d25-4\n\
 d4=d1+d2                      ; second interpulse delay of observer\n\
;                                sequence, available observer window\n\
 d4=d4-d3                      ; subtract deadtime\n\
 d13=d3                        ; initial delay between second observer\n\
                               ; pulse and pump pulse\n\
 d21=d1                        ; initial tau1 value\n\
\n\
   for j=1 to m                ; tau1 averaging loop for suppression of\n\
;                                nuclear modulations\n\
      d23=d13                  ; reset initial pump pulse delay\n\
      d24=d4                   ; reset initial delay of last observer\n\
;                                pulse\n\
      sweep x=1 to sx          ; sweep loop\n\
       shot i=1 to h           ; accumulation loop\n\
         d9                    ; DAF\n\
         p1 [ph1]              ; 1st observer pulse (pi/2) and phase program\n\
         d21                   ; tau1\n\
         p0 [ph2]              ; 2nd observer pulse (pi) and phase program\n\
         d23                   ; pump pulse delay, initially deadtime\n\
         p2 [ph3]              ; pump pulse (pi)\n\
         d24                   ; tau1+tau2-d13\n\
         p0 [ph2]             ; 3rd observer pulse (pi) and phase program\n\
         d2                    ; tau2\n\
         d0                    ; constant acquisition delay\n\
         d25                   ; trigger correction\n\
         acq [sg1]             ; acquisition\n\
       next i                  ; end of accumulation loop\n\
       dx=dx+d30               ; define time axis increment\n\
       d23=d23+d30             ; increment pump pulse delay\n\
       d24=d24-d30             ; decrement delay after pump pulse   \n\
      next x                   ; end of sweep loop  \n\
      d21=d21+d31              ; increment tau1 for nucl. mod. averaging   \n\
      d13=d13+d31              ; pump pulse delay must be incremented, too\n\
   next j                      ; end of nucl. mod. averaging loop\n\
   scansdone(k)                ; output of scans done  \n\
 next k                        ; end of scan accumulation loop \n\
end exp4  \n\
\n\
\n\
begin exp5 "standing DEER" [INTG QUAD]\n\
\n\
 for k=1 to n                  ; averaging loop\n\
\n\
 pg=2*p0;                        ; integration gate, for homemade Q band is 12 ns\n\
d13=d1\n\
d4=d2\n\
d25=200-pg                    ; trigger correction for shorter gate, not is use, remember to change d0 by hand \n\
d25=d25/2\n\
d25=d25-p0\n\
d25=d25-4\n\
 totscans(n)                   ; output of total number of scans\n\
 \n\
   sweep x=1 to sx             ; sweep loop\n\
    shot i=1 to h              ; accumulation loop\n\
      d9                       ; DAF\n\
      p1 [ph1]                 ; 1st pulse and phase program\n\
      d1                       ; tau1\n\
      p0 [ph2]                 ; 2nd pulse and phase program\n\
      d13                       \n\
      p2 [ph3]                 ; pi pulse on second frequency\n\
      d4                       ; tau1+tau2-d13\n\
      p0 [ph2]\n\
      d2\n\
      d0                       ; constant acquisition delay\n\
      d25                      ; trigger correction shorter gate\n\
      acq [sg1]                ; acquisition\n\
    next i                     ; end of accumulation loop\n\
    dx=dx+d30                  ; define trigger increment   \n\
   next x                      ; end of sweep loop\n\
   dx=0\n\
  \n\
  scansdone(k)                 ; output of scans done\n\
  \n\
 next k\n\
  \n\
end exp5  \n\
\n\
\n\
begin exp6 "DEER field avg." [INTG QUAD]\n\
\n\
 totscans(n)                   ; output of total number of scans\n\
 for k=1 to n                  ; averaging loop\n\
\n\
 bsweep y=1 to sy\n\
 \n\
 pg=2*p0;                        ; integration gate, standard 32 ns\n\
 d25=200-pg                    ; trigger correction for shorter gate, not is use, remember to change d0 by hand \n\
 d25=d25/2\n\
 d25=d25-p0\n\
 d25=d25-4\n\
\n\
 d4=d1+d2                      ; second interpulse delay of observer\n\
;                                sequence, available observer window\n\
 d4=d4-d3                      ; subtract deadtime\n\
 d13=d3                        ; initial delay between second observer\n\
                               ; pulse and pump pulse\n\
 d21=d1                        ; initial tau1 value\n\
 \n\
   for j=1 to m                ; tau1 averaging loop for suppression of\n\
;                                nuclear modulations\n\
      d23=d13                  ; reset initial pump pulse delay\n\
      d24=d4                   ; reset initial delay of last observer\n\
;                                pulse\n\
      sweep x=1 to sx          ; sweep loop\n\
       shot i=1 to h           ; accumulation loop\n\
         d9                    ; DAF\n\
         p1 [ph1]              ; 1st observer pulse (pi/2) and phase program\n\
         d21                   ; tau1\n\
         p0 [ph2]              ; 2nd observer pulse (pi) and phase program\n\
         d23                   ; pump pulse delay, initially deadtime\n\
         p2 [ph3]              ; pump pulse (pi)\n\
         d24                   ; tau1+tau2-d13\n\
         p0 [ph2]             ; 3rd observer pulse (pi) and phase program\n\
         d2                    ; tau2\n\
         d0                    ; constant acquisition delay\n\
         d25                   ; trigger correction shorter gate\n\
         acq [sg1]             ; acquisition\n\
       next i                  ; end of accumulation loop\n\
       dx=dx+d30               ; define time axis increment\n\
       d23=d23+d30             ; increment pump pulse delay\n\
       d24=d24-d30             ; decrement delay after pump pulse   \n\
      next x                   ; end of sweep loop  \n\
      d21=d21+d31              ; increment tau1 for nucl. mod. averaging   \n\
      d13=d13+d31              ; pump pulse delay must be incremented, too\n\
   next j                      ; end of nucl. mod. averaging loop\n\
   next y\n\
   scansdone(k)                ; output of scans done  \n\
 next k                        ; end of scan accumulation loop \n\
end exp6  \n\
\n\
;\n\
;\n\
; -y eldor nutation for testing inversion efficiency and pump field \n\
\n\
\n\
begin exp7 "-y eldor  +y nutation" [INTG QUAD]         ; QUAD detection\n\
\n\
\n\
p10=2*p1\n\
pg=p10\n\
d30=d29*sx\n\
d30=d30+200\n\
d25=200-pg                    ; trigger correction for shorter gate\n\
d25=d25/2\n\
for k=1 to n\n\
  dx=0\n\
  p20=0\n\
  totscans(n)\n\
  sweep x=1 to sx\n\
    shot i=1 to h              ; accumulation loop\n\
      d9                       ; DAF\n\
      p20 [ph1]\n\
      d30\n\
      p1 [ph2]                 ; 1st pulse and phase program\n\
      d1                       ; tau                       \n\
      p10 [ph3]                 ; 2nd pulse and phase program\n\
      d1                       ; tau\n\
      d0                       ; constant acquisition delay\n\
      d25                      ; trigger correction shorter gate\n\
      acq [sg1]                ; acquisition\n\
    next i                     ; end of accumulation loop\n\
    dx=dx+d29\n\
    p20=p20+d29\n\
   next x\n\
scansdone(k)\n\
next k\n\
end exp7   \n\
\n\
\n\
begin exp8 "DEER-II" [INTG QUAD]\n\
\n\
 for k=1 to n                  ; averaging loop\n\
\n\
 pg=2*p0;                        ; integration gate, for homemade Q band is 12 ns\n\
d13=d1\n\
d4=d2\n\
; d25=200-pg                    ; trigger correction for shorter gate, not is use, remember to change d0 by hand \n\
; d25=d25-p0\n\
; d25=d25/2\n\
 totscans(n)                   ; output of total number of scans\n\
   d23=d13\n\
   d24=d4\n\
 \n\
   sweep x=1 to sx             ; sweep loop\n\
    shot i=1 to h              ; accumulation loop\n\
      d9                       ; DAF\n\
      p1 [ph1]                 ; 1st pulse and phase program\n\
      d1                       ; tau1\n\
      p0 [ph2]                 ; 2nd pulse and phase program\n\
      d23                       \n\
      p2 [ph3]                 ; pi pulse on second frequency\n\
      d24                       ; tau1+tau2-d13\n\
      p0 [ph2]\n\
      d2\n\
      d0                       ; constant acquisition delay\n\
;      d25                      ; trigger correction shorter gate\n\
      acq [sg1]                ; acquisition\n\
    next i                     ; end of accumulation loop\n\
    dx=dx+d30                  ; define trigger increment   \n\
    d23=d23+d30             ; increment pump pulse delay\n\
    d24=d24-d30             ; decrement delay after pump pulse   \n\
   next x                      ; end of sweep loop\n\
   dx=0\n\
  \n\
  scansdone(k)                 ; output of scans done\n\
  \n\
 next k\n\
  \n\
end exp8  \n\
\n\
\n
Psd1               {2;33,4;0} 0,0,1908,2308,2588,5708,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,24,12,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd10              {2;33,4;0} 0,1900,2300,5700,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd11              {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd12              {2;33,4;0} 0,1900,2580,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd13              {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd14              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd15              {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd16              {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd17              {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd18              {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd19              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd2               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd20              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd21              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd22              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd23              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd24              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd25              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd26              {2;33,4;0} 0,0,1780,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4236,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd27              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd28              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd29              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd3               {2;33,4;0} 0,0,1780,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4236,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd4               {2;33,4;0} 0,0,1732,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4104,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd5               {2;33,4;0} 0,700,8780,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd6               {2;33,4;0} 0,1900,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd7               {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd8               {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd9               {2;33,4;0} 0,1900,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadDetect         True
RF1Prg             {2;33,3;0.000[MHz]} 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
RF2Prg             {2;33,3;0.000[MHz]} 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
ReplaceMode        Off
ShotRepTime        5100.00 us
ShotsPLoop         5
SmoothAllowed      False
SmoothPoints       1
SptProgress        99 %
StochMode          False
SweepsPExp         20
TriggerTimeOut     19 s
XAxisQuant         Time
XSpecRes           355
YAxisQuant         Magnetic Field
YSpecRes           1

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         20
NbScansDone        20
NbScansToDo        20
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       1

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
