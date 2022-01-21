import numpy as np
import pandas as pd
import itertools
import sys
import copy

#Come constant parameters
def module10(filedata,filebase,fileout):
    with open(fileout, 'w') as f: 
        f.write('%-8s,%-30s,%-8s,%11s,    %s\n' % ("COUNTRY","SECTOR.ACTIVITY.TECHNO","GNFR","COST(\N{euro sign})","Emissions abated(kg)"))
    
    ITEMAX = 1
    VALMAX = 1.e20
    VALMIN = -VALMAX
    NEARZERO = 0.0e-10
    EPS = 0.
    ERROR = 1.
    # First, some useful procedure 
    # Define a sequence 
    def sort_uniq(sequence):
        return (x[0] for x in itertools.groupby(sorted(sequence)))
    #
    # Cut string when '.' is found
    def searchpol(ch):
        i = 0
        polname=''
        while ch[i] != '.':
            polname += ch[i]
            i+=1
        return(polname)
    #
    # Key Generator in the dataset, e.g. NOx.TRA_RD_LD4T.GAS.LFEUVI
    def keygen(iden,ind,matrix,liste):
        out=matrix[iden,ind[liste[0]]]
        for i in range(1,len(liste)):
            out = out + '.' + matrix[iden,ind[liste[i]]]
        return out
    #
    # Read cost database and input file
    df = pd.read_csv(filebase, sep=r',', engine='python') 
    #dt = pd.read_csv(filedata, sep=r',', header=None, engine='python') 
    dt = pd.read_csv(filedata, sep=r',', engine='python') 
    costbase = np.array(df)
    ncoldf = np.shape(costbase)[1] 
    nrowdf = np.shape(costbase)[0] 
    
    #Manage header, create a dictionary of column number
    header_costbase = np.array(list(df.columns))
    origdata = np.array(dt)
    counlis = origdata[:,0] 
    print(list(sort_uniq(counlis)))
    for country in list(sort_uniq(counlis)):
        inputdata=origdata[(origdata[:,0] == country)]
        ncase = np.shape(inputdata)[0] 
        #print(inputdata)
        #First screening to check inputdata on maximum abatement
        #
        maxabat = {}   #Max abatement per GNFR
        teatech = {}   #Abatement per techno as sorted in the database
        keycost = {}   #Techno key: GAINS/ACTIVITY/SECTOR
        indexdf = {}
        keygnfr = {}
        for idx in range (ncoldf):
            result = np.where(header_costbase == header_costbase[idx])
            indexdf[header_costbase[idx]]=result[0][0]
        colcoun = costbase[:,indexdf['EU_COUNTRY_CODE']] # Country
        colpolu = costbase[:,indexdf['POL']]             # Pollutant
        colgnfr = costbase[:,indexdf['GNFR']]            # GNFR
        #  
        rankprior = np.empty((ncase,1),int)
        before = ''
        for idx in range(ncase):
            if (inputdata[idx,1] == before):
                rankprior[idx,0] = 0
            else:
                rankprior[idx,0] = 1
                before = inputdata[idx,1]
            
            subset = costbase[(colcoun == inputdata[idx,0]) & (colpolu == inputdata[idx,2]) & (colgnfr == inputdata[idx,1])]
            nrowsubset = np.shape(subset)[0]
            #Initialisation of delta emission abated
            keyg=inputdata[idx,2]+'.'+inputdata[idx,1]
            maxabat[keyg] = VALMIN
            if nrowsubset > 0:
                for idy in range(nrowsubset):
                    key  = keygen(idy,indexdf,subset,('POL','GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                    keyc = keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                    keyr = keygen(idy,indexdf,subset,('POL','GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                    teatech[key] =  subset[idy,indexdf['tea']]
                    keycost[key] =  keyc
                    keygnfr[key] =  keyg
        for idx in range(len(np.array(list(keycost.items())))):
            keyp = list(keycost.items())[idx][0]
            keys = list(keycost.values())[idx]
            if teatech[keyp] > maxabat[keygnfr[keyp]]:
                maxabat[keygnfr[keyp]] = teatech[keyp]
        ninputdata = np.empty((0,5))
        nrankprior = np.empty((0,1),int)
        for idx in range(ncase):
            keyg=inputdata[idx,2]+'.'+inputdata[idx,1]
            #if maxabat[keyg] <= 0.:
            #   print("Skip",keyg,".....")
            #   continue 
            ninputdata = np.append(ninputdata, [np.array(inputdata[idx,:])] , axis=0)
            nrankprior = np.append(nrankprior, [np.array(rankprior[idx,:])] , axis=0)
        #
        #I make sure I have a priority pollutant for a GNFR
        ncase = np.shape(ninputdata)[0] 
        before = ''
        keyprio = {}
        for idx in range(ncase):
           keyg=ninputdata[idx,2]+'.'+ninputdata[idx,1]
           if ninputdata[idx,1] != before:
              nrankprior[idx,0] = 1
           before = ninputdata[idx,1]
           keyprio[keyg] = nrankprior[idx,0]
        #
        #print(keyprio)
        #Here I load the inputdata after screening with strickly positive max abatement factor
        #ninputdata is the inputdata
        #
        ini_inputdata = copy.deepcopy(ninputdata)
        polist = list(sort_uniq(ninputdata[:,2]))
        ptotcrit = -1.
        for ite in range(ITEMAX):
            var4ranking=np.empty((ncase))
            rdea = {}      #delta emission abated per pollutant, gains and activity sector and abatement technology
            rdeac = {}      #delta emission abated per pollutant, gains and activity sector and abatement technology
            targetk = {}
            margcost = {}
            maxabat = {}   #Max abatement per GNFR
            teatech = {}   #Abatement per techno as sorted in the database
            keycost = {}   #Techno key: GAINS/ACTIVITY/SECTOR
            keygnfr = {}
            remgain = {}   #remaining GAINS emissions after abatement
            remtech = {}   #remaining GAINSup to technology emissions after abatement
            totperc = {}
            topreme = {}
            remgnfr = {}   #remaining GNFR emissions after abatement
            savgnfr = {}   #Initial emissions
            costech = {}   #cost per abated emissions (Euros/kg)
            isacost = {}    #Check if a cost is available for the global key
            keygnfr = {}   #GNFR key: POL/GNFR
            keyrlis = {}   #GNFR key: POL/GNFR
            kkkgnfr = {}   #GNFR key: POL/GNFR
            priorpo = {}   #Main pollutant
            otherpo = {}   #Polutant concerned by the techno
            isabated = {}   #Polutant concerned by the techno
            otherre = {}   #Removal efficiency of polutant concerned by the techno
            remetec = {}   #removal factor for main pollutant
            minrtec = {}   #min removal factor for the techno
            gnslist = {}   #list of keys per GNFR
            gnfcost = {}   #GNFR cost
            gnspols = {}   #List of all potential pollutants affected by GNFR
            reme = {}      #removal efficiency of the key
            efac = {}      #Emission factor
            ifte = {}      #Is the techno used or not 
            perc = {}      #Percentage avec techno in the activity, could be used or not
            critkey = {}
            deltacr = {}
            backuptech = {}
            rdealist = {}
            #
            #Here I sort the inputs by Total cost by Total Emission Abated
            #and some initializations
            colcoun = costbase[:,indexdf['EU_COUNTRY_CODE']] # Country
            colpolu = costbase[:,indexdf['POL']]             # Pollutant
            colgnfr = costbase[:,indexdf['GNFR']]            # GNFR
            
            for idx in range(ncase): 
                subset = costbase[(colcoun == ninputdata[idx,0]) & (colpolu == ninputdata[idx,2]) & (colgnfr == ninputdata[idx,1])]
                nrowsubset = np.shape(subset)[0]
                #Initialisation of delta emission abated
                keyg=ninputdata[idx,2]+'.'+ninputdata[idx,1]
                critkey[keyg] = ninputdata[idx,4]
                targetk[keyg] = (1. - ninputdata[idx,4]/100.)*ninputdata[idx,3]
                deltacr[keyg] = 0.
                remgnfr[keyg] = ninputdata[idx,3]
                maxabat[keyg] = VALMIN
                savgnfr[keyg] = ninputdata[idx,3]
                gnslist[keyg] = ([]) 
                gnfcost[keyg] = 0.
                gnspols[keyg] = ([]) 
                if nrowsubset > 0:
                    for idy in range(nrowsubset):
                        key  = keygen(idy,indexdf,subset,('POL','GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                        keyc = keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                        keyr = keygen(idy,indexdf,subset,('POL','GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                        margcost[keyc] = VALMIN
                        rdea[key] = 0.
                        rdeac[key] = 0.
                        isacost[keyc] = False  #Only for the priority pollutant and if the techno is used
                        ifte[key] = True      #False if the techno is used
                        reme[key] = subset[idy,indexdf['rem_eff']] 
                        #perc[key] = 1. - subset[idy,indexdf['perc']] 
                        perc[key] = 0.
                        remetec[keyc] = VALMIN
                        minrtec[keyc] = VALMAX
                        otherpo[keyc] = ([]) 
                        isabated[keyc] = False
                        otherre[keyc] = ([]) 
                        rdealist[keyc] = ([]) 
                        efac[key] = subset[idy,indexdf['ef_by_act']] 
                        remgain[keyr] = subset[idy,indexdf['fraction']]*ninputdata[idx,3]
                        remtech[key] = subset[idy,indexdf['perc']]*remgain[keyr]/100.
                        if (keyr in totperc):
                            totperc[keyr] = subset[idy,indexdf['perc']]/100. + totperc[keyr]
                        else:
                            totperc[keyr] = subset[idy,indexdf['perc']]/100.
                        if (keyr in topreme):
                            if (reme[key] >= topreme[keyr][0]):
                                topreme[keyr] = [reme[key],key]
                        else:
                            #totperc[keyr] = NEARZERO
                            topreme[keyr] = [reme[key],key]
                        #print(keyr,topreme[keyr])
                        costech[key] = subset[idy,indexdf['cost_per_dea']]
                        teatech[key] =  subset[idy,indexdf['tea']]
                        keycost[key] =  keyc
                        keygnfr[key] =  keyg
                        keyrlis[key] =  keyr
                        kkkgnfr[keyc] =  keyg
                        gnslist[keyg].append(keyc)
                        teatech[key] =  subset[idy,indexdf['tea']]
                    #idmax = subset[:,indexdf['tea']].argmax()
                    #var4ranking[idx] = subset[idmax,indexdf['tc']]/subset[idmax,indexdf['tea']]
                    idmax = subset[:,indexdf['cost_per_dea']].argmax()
                    var4ranking[idx] = subset[idmax,indexdf['cost_per_dea']]
                else:
                    print("Check your inputs: GNFR or Pollutant or Country is not present in the database: %s" % (keyg))
                # Herebelow I perform the last calculation of the emission by technologies remtech
                if nrowsubset > 0:
                    for idy in range(nrowsubset):
                        key  = keygen(idy,indexdf,subset,('POL','GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                        keyc = keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                        keyr = keygen(idy,indexdf,subset,('POL','GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                        #print(keyr,totperc[keyr])
                        #if totperc[keyr] > 0.:
                        #    remtech[key] = remtech[key] /totperc[keyr]
            #Add column for ranking
            listofkeys = list(reme.keys())
            dummy = np.c_[ninputdata,var4ranking]
            #Flip the rows
            input_sorted = np.flipud(dummy[dummy[:,5].argsort()][:,:])[:,:5]
            #input_sorted = np.copy(ninputdata)
            sortedlist = input_sorted[:,1]+'.'+input_sorted[:,2]
            nosortedlist = ninputdata[:,1]+'.'+ninputdata[:,2]
            #print(sortedlist)
            conv = ([])
            for idx in range(len(sortedlist)):
                conv.append(list(sortedlist).index(nosortedlist[idx]))
            #
            #Here I search the priority pollutant based on removal efficiency
            #
            for idx in range(len(np.array(list(keycost.items())))):
                keyp = list(keycost.items())[idx][0]
                keys = list(keycost.values())[idx]
                if teatech[keyp] > maxabat[keygnfr[keyp]]:
                    maxabat[keygnfr[keyp]] = teatech[keyp]
                if reme[keyp] > remetec[keys]:
                    remetec[keys] = reme[keyp]
                    priorpo[keys] = keyp[0:keyp.index('.')]
                dum = otherpo[keys]
                dur = otherre[keys]
                ppo = keyp[0:keyp.index('.')]
                dum.append(ppo)
                dur.append(reme[keyp])
                otherpo[keys] = dum
                otherre[keys] = dur
                if reme[keyp] < minrtec[keys]:
                    minrtec[keys] = reme[keyp]
            #
            #Here I search if a I am a key pollutant with possible negative impact on other pollutants
            #for idx in range(len(np.array(list(gnslist.items())))):
            for idx in range(len(list(gnslist.items()))):
                keyg = list(gnslist.items())[idx][0]
                count = {}
                for p in polist: 
                    count[p] = 0
                for idy in range(len(list(gnslist.values())[idx])):
                        count[priorpo[list(gnslist.values())[idx][idy]]] += 1
                for idy in range(len(list(gnslist.values())[idx])):
                    for idz in range(len(otherpo[list(gnslist.values())[idx][idy]])):
                        if otherpo[list(gnslist.values())[idx][idy]][idz] not in gnspols[keyg]:
                            gnspols[keyg].append(otherpo[list(gnslist.values())[idx][idy]][idz])
            for idx in range(len(np.array(list(keycost.items())))):
                keyp = list(keycost.items())[idx][0]
                keys = list(keycost.values())[idx]
                if (reme[keyp] > 0.) & (minrtec[keys] < 0.):
                    #backuptech[keyp] = True
                    backuptech[keyp] = False
                else:
                    backuptech[keyp] = False
            #Print Maximal reachable abatament
            #print("********Max abatement (%)*********")
            #for idx in range(len(np.array(list(maxabat.items())))):
            #    print("%-13s : %8.1f " % (list(maxabat.items())[idx][0],list(maxabat.values())[idx]))
            #print("**********************************")
            input_sorted = copy.deepcopy(ninputdata)
            #input_sorted = dummy[dummy[:,5].argsort()][:,:5]
            #ncase = np.shape(ninputdata)[0] 
            #print("Start Optimization...")
            #Now Let's start the optimization
            
            for idx in range(ncase): #Loop over cases in input
                subset = costbase[(colcoun == input_sorted[idx,0]) & (colpolu == input_sorted[idx,2]) & (colgnfr == input_sorted[idx,1])]
                nrowsubset = np.shape(subset)[0]
                #
                #First round I follow the ranking in the database
                #
                idy = 0
                keyg = input_sorted[idx,2]+'.'+input_sorted[idx,1]
                target = (1. - input_sorted[idx,4]/100.)*input_sorted[idx,3]
                #findpol = key[0:key.index('.')]
                #findgfr = key[key.index('.')+1:]
                #findact = findgfr[findgfr.index('.')+1:]
                #findgfr = findgfr[0:findgfr.index('.')]
                #findact = findact[0:findact.index('.')]
                #print(key,findpol,findgfr,findact)
                while (remgnfr[keyg] > target) & (idy < nrowsubset):
                    key  = keygen(idy,indexdf,subset,('POL','GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                    keyr = keygen(idy,indexdf,subset,('POL','GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                    idk = listofkeys.index(key)
                    #print(idk)
                    if ifte[key]: # & ( reme[key] > 0.):  #I do not account for negative removal efficiency
                    #if ifte[key] & ( reme[key] > 0.):  #I do not account for negative removal efficiency
                        #print("------------YYYYY",input_sorted[idx,2],remgnfr[keyg],target,ifte[key],reme[key],backuptech[key])
                        if not backuptech[key]:
                            ifte[key] = False
                            xfrac = 1.0 #Eventually we can cut the reduction
                            # From Here I calculate the impact on other pollutants if Positive impact only
                            for pol in polist:
                                if pol != input_sorted[idx,2]: #But I skip the impact on the current pollutants of course !!
                                    keyg_pol = pol+'.'+input_sorted[idx,1]
                                    key_pol = pol+'.'+keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO')) 
                                    keyr_pol = pol+'.'+keygen(idy,indexdf,subset,('GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                                    if (key_pol in reme):
                                        if ifte[key_pol]:
                                            ifte[key_pol] = False
                                            if reme[key_pol] != 100.:
                                                #print("++++++++aaaa",pol,key,key_pol)
                                                rdea[key_pol] = ( ( topreme[keyr_pol][0] - reme[key_pol] ) / (100. - reme[key_pol]) ) * remtech[key_pol]
                                                remgnfr[keyg_pol] -= rdea[key_pol]
                                                remgain[keyr_pol] -= rdea[key_pol]
                                                remtech[key_pol] -= rdea[key_pol]
                                            #else:
                                                #print("++++++++bbbb",pol,key,key_pol)
                            if reme[key] != 100.:
                                rdea[key] = ( ( topreme[keyr][0] - reme[key] ) / (100. - reme[key]) ) * remtech[key]
                                remgnfr[keyg] -= rdea[key]
                                remgain[keyr] -= rdea[key]
                                remtech[key] -= rdea[key]
                                rdeac[key] = rdea[key]
                                #print("GGGGG",key,remgnfr[keyg],topreme[keyr][0],reme[key],rdea[key],target)
                    idy += 1
                if (idy >= 1) & (remgnfr[keyg] <= target) & (idy <= nrowsubset): # & (keyprio[keyg] == 1): # Only if I am a priority pollutant
                    print("Hello")
                    idy -= 1
                    key  = keygen(idy,indexdf,subset,('POL','GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                    keyr = keygen(idy,indexdf,subset,('POL','GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                    prdea = rdea[key]
                    ifte[key] = False
                    #print("Target1",key,remgnfr[keyg],target,rdea[key],prdea,xfrac)
                    rdea[key] = rdea[key] + remgnfr[keyg] - target
                    rdeac[key] = rdea[key]
                    remgnfr[keyg] = remgnfr[keyg] - rdea[key] + prdea
                    remgain[keyr] = remgain[keyr] - rdea[key] + prdea
                    xfrac = 1.0
                    if prdea != 0.:
                        xfrac = rdea[key]/prdea
                    #print("Target1",key,remgnfr[keyg],target,rdea[key],prdea,xfrac)
                    for pol in polist:
                        if pol != input_sorted[idx,2]: #But I skip the impact on the current pollutants of course !!
                            keyg_pol = pol+'.'+input_sorted[idx,1]
                            key_pol = pol+'.'+keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                            keyr_pol = pol+'.'+keygen(idy,indexdf,subset,('GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                            #print("....",keyg_pol,remgnfr[keyg_pol],xfrac)
                            if (key_pol in reme):
                                ifte[key_pol] = False
                                #prdea = xfrac * (reme[key_pol] / 100. ) * remgain[keyr_pol] # Here possibly we can have a negative removal efficiency
                                prdea = xfrac * rdea[key_pol] # Here possibly we can have a negative removal efficiency
                                remgnfr[keyg_pol] = remgnfr[keyg_pol] + rdea[key_pol] - prdea
                                remgain[keyr_pol] = remgain[keyr_pol] + rdea[key_pol] - prdea
                                remtech[key_pol]  = remtech[key_pol]  + rdea[key_pol] - prdea
                                rdea[key_pol] = prdea
                #
            #print("END ROUND 1",remgnfr)
                #Second round, Check is backup is available 
                #Backup means positive removal for the pollutant of interest but potential negative effect on other except proprity pollutant
            for idx in range(ncase): #Loop over cases in input
                target = (1. - input_sorted[idx,4]/100.)*input_sorted[idx,3]
                subset = costbase[(colcoun == input_sorted[idx,0]) & (colpolu == input_sorted[idx,2]) & (colgnfr == input_sorted[idx,1])]
                keyg = input_sorted[idx,2]+'.'+input_sorted[idx,1]
                nrowsubset = np.shape(subset)[0]
                for idy in range(nrowsubset):
                    key  = keygen(idy,indexdf,subset,('POL','GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                    if backuptech[key] & (remgnfr[keyg] >= target):
                        print(">>>>COUCOU",key,backuptech[key],remgnfr[keyg],target)
                        keyr = keygen(idy,indexdf,subset,('POL','GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                        if ifte[key] & ( reme[key] > 0.): #I do not account for negative removal efficiency
                            zfrac = 1. #default value 
                            #Here I check for the main polltant selected if the backup techno does not affect the removal 
                            for pol in polist:
                                if pol != input_sorted[idx,2]: #But I skip the impact on the current pollutants of course !!
                                    keyg_pol = pol+'.'+input_sorted[idx,1]
                                    key_pol = pol+'.'+keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                                    keyr_pol = pol+'.'+keygen(idy,indexdf,subset,('GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                                    if (key_pol in reme):
                                        if (keyprio[keyg_pol] == 1) & (reme[key_pol] <= 0.): #in that case I remove the impact of the techno
                                            zfrac = 0.
                            for pol in polist: #Here I check all pollutant
                                if pol != input_sorted[idx,2]: #But I skip the impact on the current pollutants of course !!
                                    keyg_pol = pol+'.'+input_sorted[idx,1]
                                    key_pol = pol+'.'+keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO')) 
                                    keyr_pol = pol+'.'+keygen(idy,indexdf,subset,('GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                                    if  (key_pol in reme):
                                        #if keyprio[keyg] == 1: #
                                        ifte[key_pol] = False
                                        rdea[key_pol] = zfrac*( ( topreme[keyr_pol][0] - reme[key_pol] ) / (100. - reme[key_pol]) ) * remtech[key_pol]
                                        remgnfr[keyg_pol] -= rdea[key_pol]
                                        remgain[keyr_pol] -= rdea[key_pol]
                                        #print(">>>>COUCOU",key,key_pol,reme[key_pol])
                            # From Here I calculate the impact on other pollutants if any
                            rdea[key] = zfrac*( ( topreme[keyr][0] - reme[key] ) / (100. - reme[key]) ) * remtech[key]
                            ifte[key] = not bool(round(zfrac))
                            #rdea[key] = (perc[key] / 100.)*( subset[idy,indexdf['rem_eff']] / 100. ) * remgain[keyr] # Fraction of GNFR * Total Emission in kg
                            remgnfr[keyg] -= rdea[key]
                            remgain[keyr] -= rdea[key]
                            rdeac[key] = rdea[key]
                            #print("END ROUND 2a",remgnfr[keyg],remgain[keyr])
                            if (remgnfr[keyg] < target) & (idy <= nrowsubset):# & (keyprio[keyg] == 1): # Only if I am a priority pollutant
                                prdea = rdea[key]
                                ifte[key] = False
                                rdea[key] = rdea[key] + remgnfr[keyg] - target
                                #print("TARGET1",key,remgnfr[keyg],target,rdea[key],prdea,xfrac)
                                remgnfr[keyg] = remgnfr[keyg] - rdea[key] + prdea
                                remgain[keyr] = remgain[keyr] - rdea[key] + prdea
                                xfrac = 1.0
                                #print("TARGET2",key,remgnfr[keyg],target,rdea[key],prdea,xfrac)
                                rdeac[key] = rdea[key]
                                if prdea != 0.:
                                    xfrac = rdea[key]/prdea
                                for pol in polist:
                                    if pol != input_sorted[idx,2]: #But I skip the impact on the current pollutants of course !!
                                        keyg_pol = pol+'.'+input_sorted[idx,1]
                                        key_pol = pol+'.'+keygen(idy,indexdf,subset,('GAINS_SECTOR','GAINS_ACTIVITY','TECHNO'))
                                        keyr_pol = pol+'.'+keygen(idy,indexdf,subset,('GNFR','GAINS_SECTOR','GAINS_ACTIVITY'))
                                        if (key_pol in reme):
                                            #print(">>>>CCCC",key,reme[key],xfrac)
                                            ifte[key_pol] = False
                                            #prdea = xfrac * (reme[key_pol] / 100. ) * remgain[keyr_pol] # Here possibly we can have a negative removal efficiency
                                            prdea = xfrac * rdea[key_pol] # Here possibly we can have a negative removal efficiency
                                            remgnfr[keyg_pol] = remgnfr[keyg_pol] + rdea[key_pol] - prdea
                                            remgain[keyr_pol] = remgain[keyr_pol] + rdea[key_pol] - prdea
                                            rdea[key_pol] = prdea
                        
            #print("END ROUND 2b",remgnfr)
            #Here I start to calculate the cost by sector/activity/techniques and by pollutants
            #but I take the maximum cost in case of muti pollutant case
            #
            totcrit = 0.
            for idx in range(len(np.array(list(savgnfr.items())))):
                key = list(savgnfr.items())[idx][0]
                pol = key[0:key.index('.')]
                gnf = key[key.index('.')+1:]
                ntar = list(remgnfr.values())[idx]-list(savgnfr.values())[idx]
                ntar = -100. * ntar / list(savgnfr.values())[idx] #Effective positive emission reduction
                if (ntar - critkey[key] < -EPS) & (ntar < maxabat[key]):
                #if (ntar - ini_inputdata[idx,4] < -EPS) & (ntar < maxabat[key]):
                    totcrit = totcrit - ntar + critkey[key] #I sum over all pollutants that not reach the initial criteria
                    #totcrit = totcrit - ntar + ini_inputdata[idx,4] #I sum over all pollutants that not reach the initial criteria
                    for idy in range(len(np.array(list(savgnfr.items())))): #From here I test if another pollutant for 
                        rkey = list(savgnfr.items())[idy][0]                #a the given sector
                        rpol = rkey[0:rkey.index('.')]
                        rgnf = rkey[rkey.index('.')+1:]
                        rtar = list(remgnfr.values())[idy]-list(savgnfr.values())[idy]
                        rtar = -100. * rtar / list(savgnfr.values())[idy]
                        #Test on same sector, different pollutant, effective reduction above target, not already calculated
                        if (gnf == rgnf) & (rpol != pol) & (rtar - critkey[rkey] >= -EPS) & (deltacr[rkey] == 0.) & (maxabat[rkey] > rtar):
                        #print(pol,rpol,rtar,ini_inputdata[idy,4])
                        #if (gnf == rgnf) & (rpol != pol) & (rtar - ini_inputdata[idy,4] >= -EPS) & (deltacr[rkey] == 0.) & (maxabat[rkey] > rtar):
                            #deltacr[rkey] = (maxabat[rkey] - rtar)/50.
                            deltacr[rkey] = rtar/100.
            
            #print("CRITERIA=",totcrit)            
            #if (totcrit == 0.) | (ptotcrit == totcrit):
            if (totcrit < ERROR) | (np.sum(np.array(list(deltacr.values()))) == 0.):
                break  
            ptotcrit = totcrit
            for idx in range(len(np.array(list(savgnfr.items())))):
                key = list(savgnfr.items())[idx][0]
                ninputdata[idx,4] = ninputdata[idx,4] - deltacr[key]
        #
        #End of iteration loop
        #I keep the maximal cost for eacj technologie and avoid double counting
        #
        #for idx in range(len(np.array(list(otherpo.items())))):
        for idx in range(len(list(otherpo.items()))):
            keyc = list(otherpo.items())[idx][0]
            for idy in range(len(otherpo[keyc])):
                keyp = otherpo[keyc][idy]+'.'+keyc
                tmp = VALMIN
                if ( topreme[keyrlis[keyp]][0] != reme[keyp] ):
                    tmp = rdeac[keyp]*( costech[topreme[keyrlis[keyp]][1]]*topreme[keyrlis[keyp]][0] - costech[keyp] * reme[keyp] ) / \
                      ( topreme[keyrlis[keyp]][0] - reme[keyp] )
                if (tmp > margcost[keyc]) & (ifte[keyp] == False):
                    margcost[keyc] = tmp
                    isacost[keyc] = True
        #
        for idx in range(len(np.array(list(margcost.items())))):
            keys = list(margcost.items())[idx][0]
            if isacost[keys]:
                #print(keys,margcost[keys])
                gnfcost[kkkgnfr[keys]] = gnfcost[kkkgnfr[keys]] + margcost[keys] 
        #
        totabated = np.sum(np.array(list(rdea.values())))
        if (totabated != 0.0) :
            currency_string = "{:s} Total cost: \N{euro sign} {:,.0f} - [ \N{euro sign}/kg abated on average: {:,.2f}] - Final Abat.: {:,.2f} %". \
                           format(str(country),np.sum(np.array(list(gnfcost.values()))), \
                           np.sum(np.array(list(gnfcost.values())))/totabated,100.*totabated/np.sum(ini_inputdata[:,3]))
        else:
            currency_string = "{:s} Total cost: \N{euro sign} {:,.0f} - [ \N{euro sign}/kg abated on average: {:s}] - Final Abat.: {:,.2f} %". \
                           format(str(country),np.sum(np.array(list(gnfcost.values()))), \
                           "ND",100.*totabated/np.sum(ini_inputdata[:,3]))
        #
        print("*********************")
        print(currency_string)
        print("*********************")
        print("******Summary********")
        print("STOP after %s/%s iterations" % (ite+1,ITEMAX))
        print("%-12s %19s %9s %5s" \
            % ('POL.GNFR','EMISabat.','EFFabat','INIabat'))
        for idx in range(len(np.array(list(savgnfr.items())))):
            key = list(savgnfr.items())[idx][0]
            diff = list(remgnfr.values())[idx]-list(savgnfr.values())[idx]
            print("%-13s : %11.1f Tons (%5.1f %1s [%5.1f]" \
            % (key,diff/1000,100*diff/list(savgnfr.values())[idx],"%)",-ini_inputdata[idx,4]))
        #
        for key, value in margcost.items(): 
            for idx in range(len(otherpo[key])):
                keyp=otherpo[key][idx]+'.'+key
                dum = (str(otherpo[key][idx])+':'+str(round(rdea[keyp],9)))
                rdealist[key].append(dum)
                if rdea[keyp] > 0.: # I check if an emission was available to abate
                    isabated[key] = True
        #
        with open(fileout, 'a+') as f: 
            #f.write('%-8s,%-30s,%-8s,%11s,    %s\n' % ("COUNTRY","SECTOR.ACTIVITY.TECHNO","GNFR","COST(\N{euro sign})","Emissions abated(kg)"))
            for key, value in margcost.items(): 
                if isacost[key] & isabated[key]:
                    f.write('%-3s,%-30s,%-8s,%11.0f,  %s\n' % (country,key,kkkgnfr[key][kkkgnfr[key].index('.')+1:], value, str(rdealist[key]).replace('\'',' ').replace('[',' ').replace(']',' ')))
