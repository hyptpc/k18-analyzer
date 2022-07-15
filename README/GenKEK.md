---
stylesheet: https://cdnjs.cloudflare.com/ajax/libs/github-markdown-css/2.10.0/github-markdown.min.css
body_class: markdown-body
pdf_options:
  format: A4
  margin: 24mm 16mm
  displayHeaderFooter: true
  headerTemplate: |-
    <style>
      section {
        margin: 0 auto;
        font-family: system-ui;
        font-size: 11px;
      }
    </style>
    <section>
    </section>
  footerTemplate: |-
    <section>
      <div>
        <span class="pageNumber"></span>
        <!-- /<span class="totalPages"></span> -->
      </div>
    </section>
---


GenFit implement into the J-PARC E42 K1.8 Analyzer
====================

<div style="text-align: right;">
 2022.07.15
 </div><br>

## Environment setting

   Nothing to do. Same as the K1.8-analyzer

## Features

   - GenFit package
   - GenKEK codes
   - GenKEK dsts
   - HypTPC geometry GDML

## Install

   Install the GenFit package and the GenKEK follows.

   How-to Compile:
   ```sh
   $> cp Makefile.genfit Makefile
   $> make
   $> make pcms
   ```
   You can complie seperately 1. K1.8-Analyer(k18ana) or 2. GenKEK(genfit)

   ```yml
   all: k18ana genfit
   k18ana: lib usr dst
   genfit: genkek genfit_dst
   ```

## Library

   lib/libK18Analyzer.a : for K1.8-analyzer compiling \
   lib/ligGenKEK.a : for K1.8-analyzer & GenFit package & GenKEK compiling


## Paramters

   Add "TPCGDML" in the conf param file
   ```yml
   TPCGDML:  param/geometry/hypypcGeo.gdml
   ```

   And add "Fitter" and "Iteration" in the USER param file
   ```yml
   #GenFit
   Fitter 0 #KalmanFitterRefTrack
   #Fitter 1 #KalmanFitter
   #Fitter 2 #DAF w/ RefTrack
   #Fitter 3 #DAF w/o RefTrack
   nIteration 5 20
   ```

## dst rule

   GenKEK dsts have name starting with "Genfit" and should be placed at "dst/" \
   e.g.) dst/GenfitHelixGeant4.cc, dst/GenfitSkeleton.cc

GenKEK development guide
====================

## GenKEK directory

   All files are placed in "genfit/genkek/". \
   e.g.) genfit/genkek/include/HypTPCTask.hh, genfit/genkek/src/HypTPCTask.cc

## Units

   GenFit : GeV/c, ns, cm, kGauss \
   K1.8Ana : GeV/c, ns, mm, T

## Features

   Environment setting and Translation

   - HypTPC geometry file (param/geometry/hyptpcGeo.gdml)
   - HypTPCFieldMan : HS field management
   - HypTPCSpacePointMeasurement & HypTPCHit : Translating HitPos&Resolution information into the GenFit format
   - HypTPCFitter : Handling GenFit Fitting algorithms

   Main Parts for the track fitting \
   Inheritance ( HypTPCTrack -> HypTPCFitProcess -> HypTPCTask )

   - HypTPCTrack : Track container \
   GenFit tracks can be provided several TrackReps to describe the same track in order to fit different particle hypotheses(Pion, Proton ...) -> Find best reslt (defult setting is finding minChi2) \
   TPCLocalTrack/TPCLocalTrackHelix should provide "PDGcode" and initial "Position seed" & "Momentum seed" from pre-fitting.
   - HypTPCFitProcess : Handling the fitting process
   - HypTPCTask : Useful functions (Handling the fit results)
   e.g.) get fitting results (chi2, ndf, tof, mom, length, residuals...) or useful functions (extrapolation...)

## Development guide

   Please add more functions in the HypTPCTask class. \
   Get a fitted-track from the container and use GenFit functions to work the way you want. \
   1. Directly use the genfit::Track or 2. Get genfit::FitStatus from the track and use it or 3. Get genfit::AbsTrackRep to use hypothese and track parameterization

   You can find most useful GenFit funtions in the follows.
   ```yml
   genfit/core/include/track.h .. : useful track functions
   genfit/core/include/MeasuredStateOnPlane.h, StateOnPlane.h ... : functions for State vector
   genfit/core/include/AbsTrackRep : extrapolation and others
   ```
   Expecially core directory has most useful functions. \
   All header files have discription.

   e.g. 1) HypTPC::GetTrackLength()
   ```yml
   genfit::Track* fittedTrack = GetTrack(trackid);
   length = 10*fittedTrack -> getTrackLen(nullptr,start,end); //cm -> mm
   return length;
   ```

   e.g. 2) HypTPC::GetChi2()
   ```yml
   genfit::FitStatus *fitStatus = GetFitStatus(trackid);
   if(fitStatus) chi2 = fitStatus -> getChi2();
   ```