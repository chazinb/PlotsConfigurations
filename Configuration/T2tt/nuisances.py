# lnN

nuisances['lumiMoriond17']  = {
    'name'  : 'lumiMoriond17', 
    'samples'  : {
        #https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
        'T2tt'                             : '1.026',
        '03_VZ'                            : '1.026',
        '04_TTTo2L2Nu'                     : '1.026',
        '06_WW'                            : '1.026',
        '05_ST'                            : '1.026',
        '09_TTW'                           : '1.026',
        '10_TTZ'                           : '1.026',
        '07_ZJets'                         : '1.026',
        '02_WZTo3LNu'                      : '1.026',
        },
    'type'  : 'lnN',
    }

nuisances['normVZ']  = {
    'name'  : 'normVZ', 
    'samples'  : {
        '03_VZ'                            : '1.5',
        },
    'type'  : 'lnN',
    }

nuisances['normTTW']  = {
    'name'  : 'normTTW', 
    'samples'  : {
        '09_TTW'                           : '1.5',
        },
    'type'  : 'lnN',
    }

nuisances['normWZ']  = {
    'name'  : 'normTTW', 
    'samples'  : {
        '02_WZTo3LNu'                      : '1.5',
        },
    'type'  : 'lnN',
    }

nuisances['normTTZ']  = {
    'name'  : 'normTTZ', 
    'samples'  : {
        '10_TTZ'                           : '1.5',
        },
    'type'  : 'lnN',
    }

nuisances['normDY']  = {
    'name'  : 'normDY', 
    'samples'  : {
        '07_ZJets'                         : '1.5',
        },
    'type'  : 'lnN',
    }

nuisances['normTTbar']  = {
    'name'  : 'normTTbar', 
    'samples'  : {
        '04_TTTo2L2Nu'                     : '1.1',
        },
    'type'  : 'lnN',
    }

nuisances['normWW']  = {
    'name'  : 'normWW', 
    'samples'  : {
        '06_WW'                            : '1.1',
        },
    'type'  : 'lnN',
    }

nuisances['normTW']  = {
    'name'  : 'normTW', 
    'samples'  : {
        '05_ST'                            : '1.1',
        },
    'type'  : 'lnN',
    }

# Shape

nuisances['Reco']  = {
    'name'  : 'Reco', 
    'samples'  : {
        'T2tt'                             : '1.',
        '03_VZ'                            : '1.',
        '04_TTTo2L2Nu'                     : '1.',
        '06_WW'                            : '1.',
        '05_ST'                            : '1.',
        '09_TTW'                           : '1.',
        '10_TTZ'                          : '1.',
        '07_ZJets'                         : '1.',
        '02_WZTo3LNu'                      : '1.',
        },
    'type'  : 'shape',
    }

nuisances['Idiso']  = {
    'name'  : 'Idiso', 
    'samples'  : {
        'T2tt'                             : '1.',
        '03_VZ'                            : '1.',
        '04_TTTo2L2Nu'                     : '1.',
        '06_WW'                            : '1.',
        '05_ST'                            : '1.',
        '09_TTW'                           : '1.',
        '10_TTZ'                          : '1.',
        '07_ZJets'                         : '1.',
        '02_WZTo3LNu'                      : '1.',
        },
    'type'  : 'shape',
    }

nuisances['Fastsim']  = {
    'name'  : 'Fastsim', 
    'samples'  : {
        'T2tt'                             : '1.',
        },
    'type'  : 'shape',
    }

#nuisances['Trigger']  = {
#    'name'  : 'Trigger', 
#    'samples'  : {
#        'T2tt'                             : '1.',
#        '03_VZ'                            : '1.',
#        '04_TTTo2L2Nu'                     : '1.',
#        '06_WW'                            : '1.',
#        '05_ST'                            : '1.',
#        '09_TTW'                           : '1.',
#        '10_TTZ'                          : '1.',
#        '07_ZJets'                         : '1.',
#        '02_WZTo3LNu'                      : '1.',
#        },
#    'type'  : 'shape',
#    }

nuisances['Trigger']  = {
    'name'  : 'Trigger', 
    'samples'  : {
        'T2tt'                             : '1.02',
        '03_VZ'                            : '1.02',
        '04_TTTo2L2Nu'                     : '1.02',
        '06_WW'                            : '1.02',
        '05_ST'                            : '1.02',
        '09_TTW'                           : '1.02',
        '10_TTZ'                           : '1.02',
        '07_ZJets'                         : '1.02',
        '02_WZTo3LNu'                      : '1.02',
        },
    'type'  : 'lnN',
    }

nuisances['Btag']  = {
    'name'  : 'Btag', 
    'samples'  : {
        'T2tt'                             : '1.',
        '03_VZ'                            : '1.',
        '04_TTTo2L2Nu'                     : '1.',
        '06_WW'                            : '1.',
        '05_ST'                            : '1.',
        '09_TTW'                           : '1.',
        '10_TTZ'                          : '1.',
        '07_ZJets'                         : '1.',
        '02_WZTo3LNu'                      : '1.',
        },
    'type'  : 'shape',
    }

nuisances['BtagFS']  = {
    'name'  : 'BtagFS', 
    'samples'  : {
        'T2tt'                             : '1.',
        },
    'type'  : 'shape',
    }

nuisances['JES']  = {
    'name'  : 'JES', 
    'samples'  : {
        'T2tt'                             : '1.',
        '03_VZ'                            : '1.',
        '04_TTTo2L2Nu'                     : '1.',
        '06_WW'                            : '1.',
        '05_ST'                            : '1.',
        '09_TTW'                           : '1.',
        '10_TTZ'                           : '1.',
        '07_ZJets'                         : '1.',
        '02_WZTo3LNu'                      : '1.',
        },
    'type'  : 'shape',
    }

# Stat

nuisances['stat']  = {
    # apply to the following samples: name of samples here must match keys in samples.py
    'samples'  : {
        
        'T2tt'                           : {
            'typeStat' : 'bbb',
            },
        
        '04_TTTo2L2Nu': {
            'typeStat' : 'bbb',
                         },  
        '05_ST': {
            'typeStat' : 'bbb',
                         },  
        '06_WW': {
            'typeStat' : 'bbb',
                         }, 
        '03_VZ': {
            'typeStat' : 'bbb',
                         },
        '09_TTW': {
            'typeStat' : 'bbb',
                         }, 
        '10_TTZ': {
            'typeStat' : 'bbb',
                         },  
        '07_ZJets': {
            'typeStat' : 'bbb',
                         }, 
        '02_WZTo3LNu': {
            'typeStat' : 'bbb',
                         },                                
        },
    'type'  : 'shape'
    }

# MT2ll correlate regions
"""
nuisances['MT2llBin4Top']  = {
    'name'  : 'MT2llBin4Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5Top']  = {
    'name'  : 'MT2llBin5Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6Top']  = {
    'name'  : 'MT2llBin6Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7Top']  = {
    'name'  : 'MT2llBin7Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4WW']  = {
    'name'  : 'MT2llBin4WW',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5WW']  = {
    'name'  : 'MT2llBin5WW',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6WW']  = {
    'name'  : 'MT2llBin6WW',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7WW']  = {
    'name'  : 'MT2llBin7WW',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }
"""
# MT2ll correlate bins
"""
nuisances['MT2llBin4Bin5Bin6Bin7SR1Top']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7SR1Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4Bin5Bin6Bin7SR2Top']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7SR2Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4Bin5Bin6Bin7SR3Top']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7SR3Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4Bin5Bin6Bin7SR1WW']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7SR1WW',  
    'samples'  : {
        '06_WW'        : '0.05',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4Bin5Bin6Bin7SR2WW']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7SR2WW',  
    'samples'  : {
        '06_WW'        : '0.05',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4Bin5Bin6Bin7SR3WW']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7SR3WW',  
    'samples'  : {
        '06_WW'        : '0.05',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }
"""
# MT2ll correlate regions and bins
"""
nuisances['MT2llBin4Bin5Bin6Bin7Top']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7Top',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4Bin5Bin6Bin7WW']  = {
    'name'  : 'MT2llBin4Bin5Bin6Bin7WW',  
    'samples'  : {
        '06_WW' : '0.05',
        },
    'type'  : 'shape',
    }
"""
# MT2ll correlate nothing

nuisances['MT2llBin4TopSR1']  = {
    'name'  : 'MT2llBin4TopSR1',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4TopSR2']  = {
    'name'  : 'MT2llBin4TopSR2',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4TopSR3']  = {
    'name'  : 'MT2llBin4TopSR3',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5TopSR1']  = {
    'name'  : 'MT2llBin5TopSR1',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6TopSR1']  = {
    'name'  : 'MT2llBin6TopSR1',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.15',
        '05_ST'        : '0.15',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7TopSR1']  = {
    'name'  : 'MT2llBin7TopSR1',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.30',
        '05_ST'        : '0.30',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5TopSR2']  = {
    'name'  : 'MT2llBin5TopSR2',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6TopSR2']  = {
    'name'  : 'MT2llBin6TopSR2',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.15',
        '05_ST'        : '0.15',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7TopSR2']  = {
    'name'  : 'MT2llBin7TopSR2',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.30',
        '05_ST'        : '0.30',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5TopSR3']  = {
    'name'  : 'MT2llBin5TopSR3',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.05',
        '05_ST'        : '0.05',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6TopSR3']  = {
    'name'  : 'MT2llBin6TopSR3',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.15',
        '05_ST'        : '0.15',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7TopSR3']  = {
    'name'  : 'MT2llBin7TopSR3',  
    'samples'  : {
        '04_TTTo2L2Nu' : '0.30',
        '05_ST'        : '0.30',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4WWSR1']  = {
    'name'  : 'MT2llBin4WWSR1',  
    'samples'  : {
        '06_WW'        : '0.05',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4WWSR2']  = {
    'name'  : 'MT2llBin4WWSR2',  
    'samples'  : {
        '06_WW'        : '0.05',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin4WWSR3']  = {
    'name'  : 'MT2llBin4WWSR3',  
    'samples'  : {
        '06_WW'        : '0.05',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5WWSR1']  = {
    'name'  : 'MT2llBin5WWSR1',  
    'samples'  : {
        '06_WW' : '0.05',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6WWSR1']  = {
    'name'  : 'MT2llBin6WWSR1',  
    'samples'  : {
        '06_WW' : '0.15',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7WWSR1']  = {
    'name'  : 'MT2llBin7WWSR1',  
    'samples'  : {
        '06_WW' : '0.30',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5WWSR2']  = {
    'name'  : 'MT2llBin5WWSR2',  
    'samples'  : {
        '06_WW' : '0.05',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6WWSR2']  = {
    'name'  : 'MT2llBin6WWSR2',  
    'samples'  : {
        '06_WW' : '0.15',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7WWSR2']  = {
    'name'  : 'MT2llBin7WWSR2',  
    'samples'  : {
        '06_WW' : '0.30',
        },
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin5WWSR3']  = {
    'name'  : 'MT2llBin5WWSR3',  
    'samples'  : {
        '06_WW' : '0.05',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin6WWSR3']  = {
    'name'  : 'MT2llBin6WWSR3',  
    'samples'  : {
        '06_WW' : '0.15',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

nuisances['MT2llBin7WWSR3']  = {
    'name'  : 'MT2llBin7WWSR3',  
    'samples'  : {
        '06_WW' : '0.30',
        },
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shape',
    }

# VR --> From fit (not used)
"""
nuisances['TopMT2llShape']  = {
    'name'  : 'TopMT2llShape', 
    'rootfile'  : 'VR1TagFit', 
    'samples'  : {
        '04_TTTo2L2Nu' : '1.',
        '05_ST'        : '1.',
        },
    'cuts'  : {
        'VR1_Tag_ee' : '1.00',
        'VR1_Tag_em' : '1.00',
        'VR1_Tag_mm' : '1.00',
        'VR1_NoTag_ee' : '1.00',
        'VR1_NoTag_em' : '1.00',
        'VR1_NoTag_mm' : '1.00',
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shapefit',
    }

nuisances['WWMT2llShape']  = {
    'name'  : 'WWMT2llShape', 
    'rootfile'  : 'VR1NoTagFit', 
    'samples'  : {
        '06_WW' : '1.',
        },
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_Tag_em' : '1.00',
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        'SR1_NoTag_em' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        'SR2_Tag_ee' : '1.00',
        'SR2_Tag_em' : '1.00',
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        'SR2_NoTag_em' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        'SR3_Tag_ee' : '1.00',
        'SR3_Tag_em' : '1.00',
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        'SR3_NoTag_em' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    'type'  : 'shapefit',
    }
"""

# rateParam ttbar

nuisances['Topnorm_SR1_em']  = {
    'name'  : 'Topnorm_SR1_em',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_em' : '1.00',
        'SR1_NoTag_em' : '1.00',
        },
    }

nuisances['Topnorm_SR1_ee']  = {
    'name'  : 'Topnorm_SR1_ee',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        },
    }

nuisances['Topnorm_SR1_mm']  = {
    'name'  : 'Topnorm_SR1_mm',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    }

nuisances['Topnorm_SR2_em']  = {
    'name'  : 'Topnorm_SR2_em',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_em' : '1.00',
        'SR2_NoTag_em' : '1.00',
        },
    }

nuisances['Topnorm_SR2_ee']  = {
    'name'  : 'Topnorm_SR2_ee',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        },
    }

nuisances['Topnorm_SR2_mm']  = {
    'name'  : 'Topnorm_SR2_mm',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    }

nuisances['Topnorm_SR3_em']  = {
    'name'  : 'Topnorm_SR3_em',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_em' : '1.00',
        'SR3_NoTag_em' : '1.00',
        },
    }

nuisances['Topnorm_SR3_ee']  = {
    'name'  : 'Topnorm_SR3_ee',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        },
    }

nuisances['Topnorm_SR3_mm']  = {
    'name'  : 'Topnorm_SR3_mm',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    }

nuisances['Topnorm_VR1_em']  = {
    'name'  : 'Topnorm_VR1_em',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'VR1_Tag_em' : '1.00',
        'VR1_NoTag_em' : '1.00',
        },
    }

nuisances['Topnorm_VR1_ee']  = {
    'name'  : 'Topnorm_VR1_ee',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'VR1_Tag_ee' : '1.00',
        'VR1_NoTag_ee' : '1.00',
        },
    }

nuisances['Topnorm_VR1_mm']  = {
    'name'  : 'Topnorm_VR1_mm',
    'samples'  : {
        '04_TTTo2L2Nu' : '1.00',
        #'05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'VR1_Tag_mm' : '1.00',
        'VR1_NoTag_mm' : '1.00',
        },
    }

# rateParam tW

nuisances['tWnorm_SR1_em']  = {
    'name'  : 'Topnorm_SR1_em',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_em' : '1.00',
        'SR1_NoTag_em' : '1.00',
        },
    }

nuisances['tWnorm_SR1_ee']  = {
    'name'  : 'Topnorm_SR1_ee',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        },
    }

nuisances['tWnorm_SR1_mm']  = {
    'name'  : 'Topnorm_SR1_mm',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    }

nuisances['tWnorm_SR2_em']  = {
    'name'  : 'Topnorm_SR2_em',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_em' : '1.00',
        'SR2_NoTag_em' : '1.00',
        },
    }

nuisances['tWnorm_SR2_ee']  = {
    'name'  : 'Topnorm_SR2_ee',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        },
    }

nuisances['tWnorm_SR2_mm']  = {
    'name'  : 'Topnorm_SR2_mm',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    }

nuisances['tWnorm_SR3_em']  = {
    'name'  : 'Topnorm_SR3_em',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_em' : '1.00',
        'SR3_NoTag_em' : '1.00',
        },
    }

nuisances['tWnorm_SR3_ee']  = {
    'name'  : 'Topnorm_SR3_ee',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        },
    }

nuisances['tWnorm_SR3_mm']  = {
    'name'  : 'Topnorm_SR3_mm',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    }

nuisances['tWnorm_VR1_em']  = {
    'name'  : 'Topnorm_VR1_em',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'VR1_Tag_em' : '1.00',
        'VR1_NoTag_em' : '1.00',
        },
    }

nuisances['tWnorm_VR1_ee']  = {
    'name'  : 'Topnorm_VR1_ee',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'VR1_Tag_ee' : '1.00',
        'VR1_NoTag_ee' : '1.00',
        },
    }


nuisances['tWnorm_VR1_mm']  = {
    'name'  : 'Topnorm_VR1_mm',
    'samples'  : {
        #'04_TTTo2L2Nu' : '1.00',
        '05_ST' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'VR1_Tag_mm' : '1.00',
        'VR1_NoTag_mm' : '1.00',
        },
    }

# rateParam WW

nuisances['WWnorm_SR1_em']  = {
    'name'  : 'WWnorm_SR1_em',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_em' : '1.00',
        'SR1_NoTag_em' : '1.00',
        },
    }

nuisances['WWnorm_SR1_ee']  = {
    'name'  : 'WWnorm_SR1_ee',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_ee' : '1.00',
        'SR1_NoTag_ee' : '1.00',
        },
    }

nuisances['WWnorm_SR1_mm']  = {
    'name'  : 'WWnorm_SR1_mm',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR1_Tag_mm' : '1.00',
        'SR1_NoTag_mm' : '1.00',
        },
    }

nuisances['WWnorm_SR2_em']  = {
    'name'  : 'WWnorm_SR2_em',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_em' : '1.00',
        'SR2_NoTag_em' : '1.00',
        },
    }

nuisances['WWnorm_SR2_ee']  = {
    'name'  : 'WWnorm_SR2_ee',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_ee' : '1.00',
        'SR2_NoTag_ee' : '1.00',
        },
    }

nuisances['WWnorm_SR2_mm']  = {
    'name'  : 'WWnorm_SR2_mm',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR2_Tag_mm' : '1.00',
        'SR2_NoTag_mm' : '1.00',
        },
    }

nuisances['WWnorm_SR3_em']  = {
    'name'  : 'WWnorm_SR3_em',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_em' : '1.00',
        'SR3_NoTag_em' : '1.00',
        },
    }

nuisances['WWnorm_SR3_ee']  = {
    'name'  : 'WWnorm_SR3_ee',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_ee' : '1.00',
        'SR3_NoTag_ee' : '1.00',
        },
    }

nuisances['WWnorm_SR3_mm']  = {
    'name'  : 'WWnorm_SR3_mm',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        'SR3_Tag_mm' : '1.00',
        'SR3_NoTag_mm' : '1.00',
        },
    }

nuisances['WWnorm_VR1_em']  = {
    'name'  : 'WWnorm_VR1_em',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        #'VR1_Tag_em' : '1.00',
        'VR1_NoTag_em' : '1.00',
        },
    }

nuisances['WWnorm_VR1_ee']  = {
    'name'  : 'WWnorm_VR1_ee',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        #'VR1_Tag_ee' : '1.00',
        'VR1_NoTag_ee' : '1.00',
        },
    }

nuisances['WWnorm_VR1_mm']  = {
    'name'  : 'WWnorm_VR1_mm',
    'samples'  : {
        '06_WW' : '1.00',
        },
    'type'  : 'rateParam',
    'cuts'  : {
        #'VR1_Tag_mm' : '1.00',
        'VR1_NoTag_mm' : '1.00',
        },
    }
