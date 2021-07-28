

//============================================ Basic info===========================================\\
There are three directorys in repository:

1. CalculationResults contains R and Z (xx, xy, yx,yy) values obtained for DTM model
2. CalculationExample contains the example of files prepared for calculation running
3. Code directory contains the source code for calculation modules that are needed to run calculation. 
All programs except calcimpedance are written with c++. 
MKL is needed for OutputNu assembling.


//========================================= Configuration files ====================================\\
To run calculation one shold prepare 5 text files:
- frequencies
- xyzVectorB
- z_sig_2d
- objects
- settings.cfg

//============================= frequencies ==========================\\
'frequencies' file contains the list of frequencies to be calculated. 
Format is following:
<freq. count (FN)>
<freq1>
<freq2>
...
<freqFN>

//============================= xyzVectorB ==========================\\
'xyzVectorB' file contains the receivers positions
Format is following:
<receivers count (RN)>
<x1> <z1> <z1>
<x2> <z2> <z2>
...
<xRN> <zRN> <zRN>

//============================= z_sig_2d ==========================\\
'z_sig_2d' file contains the environment parameters (conductivity layers)
Format is following:
<layers count (LN)>
<z1> <cond1 (under z1)>
<z2> <cond2 (under z2)>
...
<zLN> <condLN (under zLN)>

Note that layers must be ordered from lowest to highest. Air is not included, last z should be 0.

//============================= objects ==========================\\
'objects' file contains objects parameters (3D heterogeneties)
Format is following:
<objects count (ON)>
for each object:
<x1> <x2> <y1> <y2> <z1> <z2> <conductivity>

x1, x2, y1, y2, z1, z2 are objects bounds. Note that x1 < x2, y1 < y2, z1 < z2.

//============================= settings.cfg ==========================\\
'settings.cfg' file contains mesh and SLAE solver parameters
Format is following:
<latheral mesh step (HXY)> <vertical mesh step (HZ)>
<latheral mesh sparse (kXY)> <vertical mesh sparse (kZ)>
<latheral regular step domain (RXY)> <regular step domain (RZ)>
<latheral calculation domain (DXY)> <calculation domain (DZ)>
<required SLAE relative residual (usually 1e-4 is enough)>
<max. SLAE solving iterations>


//============================================ mesh structure ===========================================\\
Mesh is built as follows:

↑                     
| Z   - - - - - - - - ------------------------------------------------
|     ↑                |          |                       |          |
|     |                |        ↑ |                       |          |
|     |                |        | |                       |          |
|     |                |     kZ | |                       |    kXY   |
|     |                |        | |                       | ------>  |
|     |DZ     - - - - -|----------|----------------------------------|
|     |       ↑        |          |         HX            |          |
|     |       |        |          |                       |          |
|     |       |RZ      |          |  Regular step domain  |          |
|     |       |        |          |                       |          |
|     |       |        |          |                     HZ|          |
|     ↓       ↓        |          |                       |          |
|     - - - - - - - - -|- - - - - | - - - -.......- - - - |          |
|     ↑       ↑        |          |       RECEIVERS       |          |
|     |       |        |          |           |           |          |
|     |       |RZ      |          |                       |          |
|     |       |        |          |           |           |          |
|     |       |        |          |                       |          |
|     |       ↓        |          |           |           |          |
|     |DZ     - - - - -|----------|----------------------------------|
|     |                |  <------ |           |           | |        |
|     |                |    kXY   |                       | |        |
|     |                |          |           |           | | kZ     |
|     |                |          |                       | ↓        |
|     ↓                |          |           |           |          |
|     - - - - - - - - ------------------------------------------------
|                      |          |           |           |          |
|                                      RXY         RXY                
|                      |          |<--------->|<--------->|          |
|                                                                     
|                      |                      |                      |
|                                 DXY                    DXY             
|                      |<-------------------->|<-------------------->|
|                                            
|--------------------------------------------------------------------------------------------->
                                                                                          X(Y)


//============================================ How to run calculation ===================================\\											  
To run calculation follow the steps:
1. Build all programs
2. Create empty directory for calculation running (let's denote it as Calculation directory)
3. Create 'Programs' directory inside of Calculation directory
4. Put all calculation programs (BuildMatrix.exe, calcimpedance.exe, COCR_FP.exe, Harm1D.exe, 
ImpedanceZToEdsAll.exe, OutputNu.exe, RegularMeshBuilder.exe) to 'Programs'
5. (If needed) Put required dll files to 'Programs' directory
5. Put CalcStarter.exe to Calculation directory
5. Put all prepared text files (frequencies, xyzVectorB, z_sig_2d, objects, settings.cfg) to Calculation directory
6. Run CalcStarter.exe and wait for it to complete
7. Look for result file (edsall0_0_imp) in Results directory


//============================================ Results file ==============================================\\
'edsall0_0_imp' file contains all signals in receivers.
File consists of small tables for each signal and positions.
The table format is followng:
<title line 1 (ignore it)>
<title line 2 (ignore it)>
<title line 3 (ignore it)>
<title line 4 (ignore it)>
<title line 5 (ignore it)>
<title line 6 (ignore it)>
<freq1> <some floating point value (ignore it)> <some floating point value (ignore it)> <signalValue1>
<freq2> <some floating point value (ignore it)> <some floating point value (ignore it)> <signalValue2>
...
<freqFN> <some floating point value (ignore it)> <some floating point value (ignore it)> <signalValueFN>

Curves are grouped by signals. Groups are ordered as follows:
Zxx
Zxy
Zyx
Zyy
Fxx
Fxy
Fyx
Fyy
Rxx
Rxy
Ryx
Ryy
Fixx
Fixy
Fiyx
Fiyy

Inside each group curves are ordered by receivers. For example for Zxx there will be a curves for
rec1, rec2, ... , recRN
Then curves for Zxy go (rec1, rec2, ... , recRN), then Zyx and so on.