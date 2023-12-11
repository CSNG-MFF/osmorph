paramsets = []
for neuron in ['L5','L23']:
    for ChRdist in ['uniform','shemesh_supfig9b_exp_yoff','shemesh_supfig9b_exp_lin_yoff']:
        for diam in [25,50,100,200,400]:
            for NA in [0.1,0.22,0.39,0.5]:
                if neuron == 'L5':
                    l_min = 1e-6
                    l_max = 1e-2
                else:
                    l_min = 1e-7
                    l_max = 7e-1
                paramsets.append(
                    ''.join(
                        [
                    'hoc_file-', neuron, '--light_model-foutz_et_al2012--chanrhod_distribution-', 
                    ChRdist, '--chanrhod_expression-130000000000--fiber_diameter-',str(diam),
                    '--fiber_NA-',str(NA), '--stim_duration-200.0--l_min-',str(lmin),
                    '--l_max-',str(l_max),'--APC_desired-10'
                        ]
                    )
                )
