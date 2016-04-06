@echo off
@SET Compare=FC /T
@SET BeamDyn=..\bin\BeamDyn_Driver_Win32.exe
@set CompareFile=CertTest.out
@SET Editor=NotePad.EXE


@IF EXIST %CompareFile%  DEL %CompareFile%

%BeamDyn% Dvr_5MW_Static.inp
%Compare% Dvr_5MW_Static.out                     Results\Dvr_5MW_Static.out               >> %CompareFile%
%Compare% Static_BeamDyn_Input_5MW.inp.sum       Results\Static_BeamDyn_Input_5MW.inp.sum >> %CompareFile%


%BeamDyn% Dvr_5MW_Dynamic.inp
%Compare% Dvr_5MW_Dynamic.out                   Results\Dvr_5MW_Dynamic.out               >> %CompareFile%
%Compare% Dynamic_BeamDyn_Input_5MW.inp.sum     Results\Dynamic_BeamDyn_Input_5MW.inp.sum >> %CompareFile%

%Editor% %CompareFile%