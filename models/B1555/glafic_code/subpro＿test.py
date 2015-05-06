
    
    # B1555 constraints
    x= np.arrange(12)
    y= np.array([0.0,0.0,17.0,0.0726,0.0480,10.54,0.4117,-0.0280,8.619,0.1619,-0.3680,1.462])
    err= np.array([1.0e-6,1.0e-6,1.70,0.01,0.01,0.97,0.01,0.01,0.83,0.06,0.06,0.13])

        
        params = [float(v) for v in line.split()]
        mody = params[0] + params[1]*x
        chi2 = np.sum( ((y-mody) / err)**2 )
        lnprob = -0.5*chi2
        sys.stdout.write(str(lnprob)+'\n')
        sys.stdout.flush()

