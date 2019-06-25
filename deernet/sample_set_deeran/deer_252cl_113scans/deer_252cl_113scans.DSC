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
XPTS	189
XMIN	-100.000000
XWID	2256.000000
*
*	Documentational Text:
*
TITL	'DEER'
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
OPER    xuser
DATE    12/22/10
TIME    00:03:48
CMNT    Deer
SAMP    WT1
SFOR    
STAG    L
EXPT    
OXS1    
AXS1    
AXS2    
AXS3    
MWPW    2e-07
A1CT    0.33822
B0VL    0.3355589
A1SW    0.01
MWFQ    9.503071e+09
AVGS    0
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
Power              0.0002000 mW
PowerAtten         60 dB

.DVC     fieldCtrl, 1.0

CenterField        3382.20 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
SweepDirection     Up
SweepWidth         100.0 G

.DVC     freqCounter, 1.0

FrequencyMon       9.503071 GHz
QMonitBridge       On

.DVC     ftBridge, 1.0

Attenuation        22.00 dB
ELDORAtt           20 dB
FrequencyA         9.427790 GHz
VideoBW            200 MHz
VideoGain          42 dB

.DVC     ftEpr, 1.0

AveragesPerScan    600
FTEzB0SettlDel     0 s
FTEzB0Step         0 G
FTEzCurrELDOR      9.427790 GHz
FTEzD11deer        200 ns
FTEzD20deer        1200 ns
FTEzD3deer         100 ns
FTEzDXdeer         12 ns
FTEzDelay0         0 ns
FTEzDelay1         0 ns
FTEzDelay2         0 ns
FTEzDeltaX         0 ns
FTEzDeltaY         0 ns
FTEzELDORStart     9.100000 GHz
FTEzELDORWidth     0.000000 GHz
FTEzENDORDetDel    0 ns
FTEzExtTrgDel      0 ns
FTEzIntgBase       4 ns
FTEzIntgGate       52 ns
FTEzMWPiHalf       16 ns
FTEzP20deer        32 ns
FTEzP3             0 ns
FTEzRF1Pi          0 ns
FTEzRF1PlsSep      4 ns
FTEzRF1Pump        0.500 MHz
FTEzRF1Start       0.500 MHz
FTEzRF1Width       0.000 MHz
FTEzSRELDORExp     False
FTEzSXdeer         105
FTEzSXdeerExp      True
FTEzSleep          0 s
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     4P ELDOR Setup
PlsSPELLISTSlct    8-step
PlsSPELPhPrgEx     Normal
QuadDetect         True
ReplaceMode        Off
ShotRepTime        1499.40 us
ShotsPLoop         50
SweepsPExp         1
TriggerTimeOut     10 s
XSpecRes           512
YSpecRes           1

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         0
NbScansDone        0
NbScansToDo        0
ReplaceMode        Off

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
