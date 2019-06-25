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
XPTS	380
XMIN	-0.000000
XWID	3032.000000
*
*	Documentational Text:
*
TITL	'DEER Mixture Birad. 1.94/2.8/3.8 nm 9.3331/9.2681 GHz, 80 K'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'Time'
IRUNI	''
IIUNI	''
XUNI	'ns'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    gunnar
DATE    07/23/05
TIME    15:47:02
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    PLS
OXS1    TADC
AXS1    B0VL
AXS2    NONE
AXS3    
MWPW    2.002e-07
A1CT    0.33
B0VL    0.33
A1SW    0.015
MWFQ    9.312e+09
AVGS    5
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
Power              0.0002002 mW
PowerAtten         60.0 dB

.DVC     endor, 1.0

EIEENDORFreq       14.9021800 MHz/3.5 kG
EIEIsotope         H1
EIERFSweepDir      Same
EIEStaticField     3480.000 G
EIEStaticRF        1.000 MHz
RF1Atten           60.000 dB
RF1FreqPos         1.000 MHz
RF1StartFreq       1.000 MHz
RF1SweepWidth      30.000 MHz
RF2Atten           60.000 dB
RF2FreqPos         1.000 MHz
RF2StartFreq       1.000 MHz
RF2SweepWidth      30.000 MHz
RFSrcMixing        Add
SumAtten           10.000 dB
SumAttenStart      3.000 dB
SumAttenWidth      77.000 dB

.DVC     fieldCtrl, 1.0

CenterField        3300.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
SweepDirection     Up
SweepWidth         150.0 G

.DVC     freqCounter, 1.0

FrequencyMon       9.312 GHz
QMonitBridge       On

.DVC     ftBridge, 1.0

Attenuation        0.00 dB
ELDORAtt           30 dB
FrequencyA         9.300000 GHz
VideoBW            50 MHz
VideoGain          57 dB

.DVC     ftEpr, 1.0

AveragesPerScan    30
ELDORFreqStart     9.300000 GHz
ELDORFreqWidth     0.700000 GHz
FTAcqModeSlct      Run from PulseSPEL
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     DEER
PlsSPELGlbTxt      \
;     PulseSPEL variables for DEER experiments\n\
;\n\
begin defs\n\
;                                               \n\
;\n\
;                                                         Gunnar Jeschke, 2002\n\
; ### Block of used variables ###\n\
;\n\
; Repetition time\n\
;\n\
srt  = 1000 * srtu ; shot repetition time, srtu= 1 µs\n\
;                  ; usual values for nitroxides: 1000 * srtu at 80 K\n\
;                  ; 10000...40000 * srtu at 15 K\n\
;\n\
; Pulse lengths\n\
;                        \n\
p0   = 32      ; observer pi pulse length (frq1)\n\
p1   = 32      ; observer pi/2 pulse length (frq1) \n\
p2   = 12      ; pump pi pulse length (frq2)\n\
;\n\
; Variable delays\n\
;\n\
d0   = 484     ; trigger delay\n\
d1   = 200     ; tau1 (initial value)\n\
d2   = 3000    ; tau2\n\
d3   = 80      ; deadtime delay during which no data are acquired\n\
;\n\
; Time increments\n\
;\n\
d30  = 8       ; t time scale increment of dipolar time evolution\n\
d31  = 8       ; time scale increment for nuclear modulation averaging\n\
;\n\
; Counters\n\
;\n\
h    = 30      ; number of shots/loop                     ( counter: I )\n\
n    = 5       ; number of sweeps to accumulate           ( counter: K )\n\
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
PlsSPELPrg         PulseSPEL/deer4p_2003.exp
PlsSPELPrgTxt      \
;\n\
; DEER setup and run \n\
; G. Jeschke, 2000-2003\n\
;\n\
; ### This assumes that +x, -x, +<x> are observer frequency channels while\n\
;     -<y> is a pump frequency (second frequency channel)\n\
;     The pulse sequence must be adjusted for the commercial Bruker ELDOR\n\
;     setup!\n\
;     pulse lengths, delays etc. are optimzed for X band\n\
;     for best results, a 3 mm split-ring resonator is recommended\n\
;  \n\
begin defs\n\
 dim s[512]               ;     dimension of data trace [sx] for setup\n\
 dim1 s[512]              ;     dimension [sx] for setup\n\
 dim2 s[512]              ;     dimension [sx] for set-up\n\
 dim3 s[1024]             ;     dimension [sx] for field sweep\n\
 dim4 s[380]              ;     dimension [sx] for DEER\n\
;                               ### must be adapted to delays!\n\
;                               sx =(d1+d2-2*d3)/d30\n\
 dim5 s[193]              ;     dimension for standing DEER\n\
 dim6 s[380,2]            ;     dimension [sx,sy] for sensitive DEER\n\
   \n\
end defs\n\
;\n\
; Pulse lengths:\n\
;\n\
; p0    observer pi pulse (frq1), recommendation for nitroxides:        32 ns\n\
; p1    observer pi/2 pulse (frq1), recommendation for nitroxides:      32 ns\n\
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
; d2    tau2,                   typical values:                 1-6 µs\n\
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
begin lists2 "+<x> (pi frq1)"\n\
 ph1 +<x>   \n\
 ph2 +<x>\n\
 ph3 +x                  \n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists2\n\
 \n\
begin lists3 "-<x> (pi frq1)"\n\
 ph1 -<x>   \n\
 ph2 -<x>\n\
 ph3 +x                  \n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists3\n\
 \n\
begin lists4 "-<y> (pi frq2)"\n\
 ph1 -<y>  \n\
 ph2 +<x>\n\
 ph3 +<x>                  \n\
 asg1 +a                  \n\
 bsg1 +b \n\
end lists4\n\
 \n\
begin lists5 "DEER setup (only observer pulses)"\n\
 ph1 +x   \n\
 ph2 +<x>\n\
 ph3 *                 \n\
 asg1 +a                   \n\
 bsg1 +b  \n\
end lists5\n\
 \n\
begin lists6 "DEER run"\n\
 ph1 +x -x   \n\
 ph2 +<x>\n\
 ph3 -<y>                 \n\
 asg1 +a -a                  \n\
 bsg1 +b -b  \n\
end lists6                  \n\
 \n\
begin lists7 "sensDEER sym"\n\
 ph1 +x -x   \n\
 ph2 +<x> +<x> -<x> -<x>\n\
 ph3 -<y>                 \n\
 asg1 +a -a +a -a                  \n\
 bsg1 +b -b +b -b  \n\
end lists7     \n\
\n\
begin lists8 "+<y> (pi frq2)"\n\
 ph1 +<y>  \n\
 ph2 +<x>\n\
 ph3 +<x>                  \n\
 asg1 +a                  \n\
 bsg1 +b \n\
end lists8\n\
 \n\
;\n\
;\n\
; Standing 2-pulse echo for setup of observer pi pulses\n\
;\n\
\n\
begin exp "pi tp=p0 +<x>/-<x>" [TRANS QUAD]         ; QUAD detection\n\
\n\
p11=p0\n\
p10=p0/2\n\
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
begin exp2 "pi tp=p2 -<y> frq2" [TRANS QUAD]         ; QUAD detection\n\
\n\
p11=p0\n\
p10=p0/2\n\
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
begin exp3 "Field Sweep +<x> pg=200ns" [INTG QUAD]        ; QUAD detection\n\
\n\
pg=200                         ; integration gate\n\
p10=p0/2\n\
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
      p0 [ph2]                 ; 2nd pulse and phase program\n\
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
 pg=p0;                        ; integration gate, standard 32 ns\n\
 d13=d1\n\
 d4=d2\n\
\n\
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
; Sensitive four-pulse DEER\n\
\n\
\n\
begin exp6 "sensDEER" [INTG QUAD]\n\
\n\
 for k=1 to n                  ; averaging loop\n\
\n\
 pg=p0;                        ; integration gate, standard 32 ns\n\
 d4=d1+d2                      ; second interpulse delay of observer\n\
;                                sequence, available observer window\n\
 d21=d1+p1                     ; initial tau1 value\n\
 d13=d2+p0                     ; initial tau2\n\
\n\
 totscans(n)                   ; output of total number of scans\n\
   for j=1 to m                ; tau1 averaging loop for suppression of\n\
;                                nuclear modulations\n\
      d24=d4                   ; reset initial delay of last observer\n\
;                                pulse\n\
      d26=0                    ; set pump pulse increment to zero\n\
      for y=1 to sy            ; loop DEER & reference signal\n\
         d23=d21                  ; reset initial pump pulse delay\n\
         d24=d13                  ; reset tau2 to tau2_0+mod. avg. delay\n\
         d25=d24                  ; reset initial delay of last observer pulse \n\
         sweep x=1 to sx          ; sweep loop\n\
          shot i=1 to h           ; accumulation loop\n\
            d9                    ; DAF\n\
            p1 [ph1]              ; 1st observer pulse (pi/2) and phase program\n\
            d21                   ; tau1\n\
            p0 [ph2]              ; 2nd observer pulse (pi) and phase program\n\
            d23                   ; pump pulse delay, initially tau1\n\
            p2 [ph3]              ; pump pulse (pi)\n\
            d24                   ; tau2-t\n\
            p0 [+<x>]             ; 3rd observer pulse (pi) and phase program\n\
            d25                   ; tau2\n\
            d0                    ; constant acquisition delay\n\
            acq [sg1]             ; acquisition\n\
          next i                  ; end of accumulation loop\n\
          dx=dx+d30               ; define time axis increment\n\
          d23=d23+d26             ; increment pump pulse delay\n\
          d24=d24-d26             ; decrement delay after pump pulse\n\
          d24=d24+d30             ; increment tau2 (delay after pump pulse)\n\
          d25=d25+d30             ; increment tau2 (delay after last observer pulse)\n\
         next x                   ; end of sweep loop\n\
         d26=d26+d30              ; set pump pulse increment to DEER increment\n\
      next y                      ; end of reference & DEER loop  \n\
      d13=d13+d31              ; pump pulse delay must be incremented, too\n\
   next j                      ; end of nucl. mod. averaging loop\n\
   scansdone(k)                ; output of scans done  \n\
 next k                        ; end of scan accumulation loop \n\
end exp6  \n\
\n\
\n\
\n\
\n
Psd1               0,0,292,548,684,3804,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,36,16,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd10              0,300,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd11              0,300,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd12              0,304,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd13              0,300,692,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd14              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd15              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd16              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd17              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd18              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd19              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd2               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd20              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd21              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd22              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd23              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd24              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd25              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd26              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd3               0,0,132,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3804,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd4               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3728,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd5               0,8,7004,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd6               0,300,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd7               0,300,556,3812,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd8               0,300,300,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd9               0,300,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadDetect         True
RF1Prg             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
RF2Prg             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
ReplaceMode        Off
ShotRepTime        1020.00 us
ShotsPLoop         30
SweepsPExp         5
TriggerTimeOut     10 s
XAxisQuant         Time
XSpecRes           380
YAxisQuant         Magnetic Field
YSpecRes           1

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         5
NbScansDone        5
NbScansToDo        5
ReplaceMode        Off

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
