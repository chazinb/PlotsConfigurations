# variables

#variables = {}
    
variables['nvtx']  = {   'name': 'nvtx',      
                        'range' : (40,0,40),  
                        'xaxis' : 'nvtx', 
                         'fold' : 3
                        }
                        
variables['mll']  = {   'name': 'mll',            #   variable name    
                        'range' : (40,0,200),    #   variable range
                        'xaxis' : 'm_{ll} [GeV]',  #   x axis name
                         'fold' : 3
                        }

variables['pt1']  = {   'name': 'std_vector_lepton_pt[0]',     
                        'range' : (100,0,200),   
                        'xaxis' : 'p_{T} 1st lep',
                        'fold'  : 3                         
                        }

variables['pt2']  = {   'name': 'std_vector_lepton_pt[1]',     
                        'range' : (100,0,200),   
                        'xaxis' : 'p_{T} 2nd lep',
                        'fold'  : 3                         
                        }


variables['eta1']  = {   'name': 'std_vector_lepton_eta[0]',     
                        'range' : (100,-3,3),   
                        'xaxis' : '#eta 1st lep',
                        'fold'  : 3                         
                        }

variables['eta2']  = {   'name': 'std_vector_lepton_eta[1]',     
                        'range' : (100,-3,3),   
                        'xaxis' : '#eta 2nd lep',
                        'fold'  : 3                         
                        }

variables['csvv2ivf_1']  = { 
                        'name': 'std_vector_jet_csvv2ivf[0]',     
                        'range' : (100,0,1),   
                        'xaxis' : 'csvv2ivf jet 1st',
                        'fold'  : 3                         
                        }

variables['pfmet']  = { 
                        'name': 'pfType1Met',     
                        'range' : (100,0,300),   
                        'xaxis' : 'pfmet [GeV]',
                        'fold'  : 3                         
                        }

#variables['jetpt1']  = {
                        #'name': 'std_vector_jet_pt[0]',     
                        #'range' : (100,0,200),   
                        #'xaxis' : 'p_{T} 1st jet'
                        #}


variables['njet']  = {
                        'name': 'njet',     
                        'range' : (20,0,20),   
                        'xaxis' : 'Number of jets',
                        'fold' : 2   # 0 = not fold (default), 1 = fold underflowbin, 2 = fold overflow bin, 3 = fold underflow and overflow
                        }




variables['jetpt1']  = {
                        'name': 'std_vector_jet_pt[0]',     
                        'range' : (20,0,200),   
                        'xaxis' : 'p_{T} 1st jet',
                        'fold' : 2   # 0 = not fold (default), 1 = fold underflowbin, 2 = fold overflow bin, 3 = fold underflow and overflow
                        }

variables['jetpt2']  = {
                        'name': 'std_vector_jet_pt[1]',     
                        'range' : (20,0,200),   
                        'xaxis' : 'p_{T} 2nd jet',
                        'fold' : 2   # 0 = not fold (default), 1 = fold underflowbin, 2 = fold overflow bin, 3 = fold underflow and overflow
                        }
