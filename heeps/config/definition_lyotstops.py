
# oversize for the spiders is calculated with respect to the enveloppe 
# but HEEPS compute it with respect to the spider width. 
# Hence an additional factor is needed to account for that: (enveloppe_spider / diameter_nominal)
# 0.26% is an average for the two type of spiders
#    (308 - 202) / 38542mm and (475-386)/38542
dspi = 0.0025 

STOPS_LM = {
            # 'PPS-LM': {#'par' : (-0.0322, 0, 0), 
            #            'force_sym': False, 
            #            'asym_wfs': None
            #            },
            'SPM-LM': {'par' : (0.01, 0.01, 0.0126+dspi), 
                       'force_sym': False, 
                       'asym_wfs': None,
                       'f_lyot_stop':''
                       },
            'RLS-LM': {'par' : (0.0492, 0.04, 0.0249+dspi), 
                       'force_sym': False, 
                       'asym_wfs': 'one_spider',
                       'asym_widths': [0.137, 0.0, 0, 0, 0, 0], # 5% transmission loss
                       'f_lyot_stop':'pupil/ls/ls_RLS-LM_10000.fits.tar.gz'
                       },
            'CLS-LM': {'par' : (0.0159, 0.09, 0.0126+dspi), 
                       'force_sym': False, 
                       'asym_wfs': 'mickey',
                       'asym_widths': [0.42, 0.0, 0.42, 0, 0, 0], # 3% transmission loss
                       'f_lyot_stop':''
                       },
            'CLS-LM-b': {'par' : (0.0224, 0.09, 0.0245+dspi), 
                         'force_sym': True, 
                         'asym_wfs': 'mickey',
                         'asym_widths': [0.5275, 0.0, 0.5275, 0, 0, 0], # 8% transmission loss
                         'f_lyot_stop':''
                         },
            'ULS-LM': {'par' : (0.1084, 0.0791, 0.0245+dspi),
                       'force_sym': False,
                       'asym_wfs': 'one_spider',
                       'asym_widths': [0.118, 0.0, 0, 0, 0, 0], # 5% transmission loss
                       'f_lyot_stop':''
                       }
            }

STOPS_N = {
            # 'PPS-N': {#'par' : (-0.0322, 0, 0),
            #              'force_sym': False,
            #              'asym_wfs': None
            #              },
           'SPM-N': {'par' : (0.01, 0.01, 0.0126+dspi),
                     'force_sym': False,
                     'asym_wfs': None,
                     'f_lyot_stop':''
                     },
           'SPM-N-b': {'par' : (0.0283, 0.0268, 0.0357+dspi),
                       'force_sym': True,
                       'asym_wfs': 'one_spider',
                       'asym_widths': [0.125, 0.0, 0, 0, 0, 0], # 5% transmission loss
                       'f_lyot_stop':''
                       },
           'CLS-N': {'par' : (0.0176, 0.09, 0.0126+dspi),
                     'force_sym': False,
                     'asym_wfs': 'mickey',
                     'asym_widths': [0.42, 0.0, 0.42, 0, 0, 0], # 3% transmission loss
                     'f_lyot_stop':''
                     },
           'CLS-N-b': {'par' : (0.0283, 0.09, 0.0357+dspi),
                       'force_sym': True,
                       'asym_wfs': 'mickey',
                       'asym_widths': [0.529, 0.0, 0.529, 0, 0, 0], # 8% transmission loss
                       'f_lyot_stop':''
                       },
           'ULS-N': {'par' : (0.1036, 0.0810, 0.0357+dspi),
                     'force_sym': False,
                     'asym_wfs': 'one_spider',
                     'asym_widths': [0.126, 0.0, 0, 0, 0, 0], # 5% transmission loss
                     'f_lyot_stop':''
                     }
            }

STOPS_LMS = {
            # 'PPS-LMS': {#'par' : (-0.0322, 0, 0),
            #              'force_sym': True,
            #              'asym_wfs': None},
             'SPM-LMS': {'par' : (0.01, 0.01, 0.0126+dspi),
                         'force_sym': True,
                         'asym_wfs': None,
                         'f_lyot_stop':''},
             'SPM-LMS-b': {'par' : (0.0256, 0.0241, 0.0308+dspi),
                         'force_sym': True,
                         'asym_wfs': None,
                         'f_lyot_stop':''},
             'RLS-LMS': {'par' : (0.0483, 0.04, 0.028+dspi),
                         'force_sym': True,
                         'asym_wfs': None,
                         'f_lyot_stop':''},
             'CLS-LMS': {'par' : (0.0168, 0.09, 0.0126+dspi), #'par' : (0.0233, 0.09, 0.0288+dspi),
                         'force_sym': True,
                         'asym_wfs': None,
                         'f_lyot_stop':''},
             'CLS-LMS-b': {'par' : (0.0233, 0.09, 0.0288+dspi),
                           'force_sym': True,
                           'asym_wfs': 'mickey',
                           'asym_widths': [0.53, 0.0, 0.53, 0, 0, 0], # 8% transmission loss
                           'f_lyot_stop':''
                           },
             'ULS-LMS': {'par' : (0.1084, 0.0791, 0.0308+dspi),
                         'force_sym': True,
                         'asym_wfs': None,
                         'f_lyot_stop':''}
    }