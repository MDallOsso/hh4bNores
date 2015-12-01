# hh4bNores
macros and code for non resonant hh->4b analysis in CMS Run2.

-- work in progress... update for V15 vhbb ntuples --
-- V12 -> Apr15 -- freezed and used with Phys14 samples.
-- V13 -> Sept15
-- V15 -> Nov15

BASELINE:

- step0  -> Step0 ntuples production from vhbb official ntuples.
	    no trigger request or kin cuts. Jets level info only.

- kinSel -> Step1 ntuples production: trigger, kinematical cuts, jets coupling
            OUTPUT: Step1 ntuples in 'data' folder + Histos (.root) in 'plots' folder + log in 'plots' folder 
            hh4b_13TeV_kinSel.cc -> produce Step1 ntuples
            Plot_kinSel_new.cc -> produce final plots from Step1 histos.

ADDITIONAL CODE:

- tmva -> hh4b_13TeV_tmva.cc , code to apply a tmva (BDT, Likelihood) to selected variables from Step1 ntuples.
          see also https://github.com/cgottard/hh4b_8TeV_CarloThesis.     

- fit -> to do a 2D binned likelihood fit + to use combine to extract final limit

- jetMatch -> preliminary study on different method to select 4th jet - to be updated

- studies -> random stuff


