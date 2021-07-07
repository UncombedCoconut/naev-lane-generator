# This file generates safe lanes

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lgs
#from scipy.linalg import null_space
import scipy.linalg as slg
import math
import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import time

def createFactions():
    '''Creates the dico of lane-making factions'''
    factions = [
                "Empire",
                "Soromid",
                "Dvaered",
                "Za'lek",
                "Collective", # TODO : see if this one is right
                "Sirius",
                "Frontier",
                "Goddard",
                "Proteron",
                "Thurion",
               ]

    return {name: i for (i, name) in enumerate(factions)}

FACTION_COLORS = ["g","orange","brown","darkred","silver","aqua","y","b","purple","grey"]

# Creates anchors to prevent the matrix from being singular
# Anchors are jumpoints, there is 1 per connected set of systems
# A Robin condition will be added to these points and we will check the rhs 
# thanks to them because in u (not utilde) the flux on these anchors should
# be 0, otherwise, it means that the 2 non-null terms on the rhs are on
# separate sets of systems.
def createAnchors(): # TODO : understand what's going on with anchors
    anchorSys = [
                 "Alteris",
                 "Flow",
                 "Zied",        # TODO : be sure this one is necessary
                 "Qorel",
                ]
    
    anchorJps = [
                 "Delta Pavonis",
                 "Aesria",
                 "Pudas",
                 "Firk",
                ]
    
    anchorAst = [
                 "Darkshed",
                 "Sevlow",
                 "Bon Sebb",
                 "Qorellia",
                ]
    
    return (anchorSys, anchorJps, anchorAst)

# Reads all the assets
def readAssets( path ):
    assets = {}
    
    for fileName in os.listdir(path):
        tree = ET.parse((path+fileName))
        root = tree.getroot()
        
        name = root.attrib['name']
        
        pos = root.find('pos')
        if pos == None: # It's a virtual asset
            x = None
            y = None
        else:
            x = float(pos.find('x').text)
            y = float(pos.find('y').text)
        
        presence = root.find('presence')
        if presence == None: # Inhabited
            faction = 'nobody'
            population = 0
            ran = 0
        else:
            faction = presence.find('faction').text
            population = float(presence.find('value').text)
            ran = int(presence.find('range').text)

        assets[name] = ( x, y, faction, population, ran )
    
    return assets
    

# Compute insystem paths. TODO maybe : use Delauney triangulation instead ?
def inSysStiff( nodess, factass, g2ass, loc2globNs ):
    #stiff = sp.csr_matrix() # Basic stiffness matrix (without lanes)
    si  = []  # Element to build the sparse default matrix
    sj  = []
    sv  = []
    sil = [] # Lanes in local numerotation
    sjl = []
    sr  = [] # Lanes reserved for a faction
    loc2globs = []
    system = [] # Tells which system the lane is in
    
    minangle = 10 # Minimal angle between 2 lanes (deg)
    crit = math.cos(math.pi*minangle/180)
    
    i  = 0
    ii = 0

    for k in range(len(nodess)):
        nodes = nodess[k]
        loc2glob = []
        loc2globN = loc2globNs[k]
        #print(nodes)
        for n in range(len(nodes)):
            xn = nodes[n][0]
            yn = nodes[n][1]
            
            na = g2ass[loc2globN[n]]  # Find global asset numerotation of local node
            #print(na)
            if na>=0:
                fn = factass[na]
            else: # It's not an asset
                fn = -1
            
            for m in range(n): # Because symmetry
                #if n >= m: 
                #    continue
                xm = nodes[m][0]
                ym = nodes[m][1]
                
                ma = g2ass[loc2globN[m]]
                if ma>=0:
                    fm = factass[ma]
                else: # It's not an asset
                    fm = -1
                
                lmn = math.hypot( xn-xm, yn-ym )
                
                # Check if very close path exist already.
                forbitten = False
                for p in range(len(nodes)):
                    if n==p or m==p:
                        continue
                    xi = nodes[p][0]
                    yi = nodes[p][1]
                    
                    # The scalar product should not be too close to 1
                    # That would mean we've got a flat triangle
                    lni = math.hypot(xn-xi, yn-yi)
                    lmi = math.hypot(xm-xi, ym-yi)
                    dpn = ((xn-xm)*(xn-xi) + (yn-ym)*(yn-yi)) / ( lmn * lni )
                    dpm = ((xm-xn)*(xm-xi) + (ym-yn)*(ym-yi)) / ( lmn * lmi )
                    if (dpn > crit and lni < lmn) or (dpm > crit and lmi < lmn):
                        forbitten = True
                        break
                    
                if forbitten:
                    continue
                
                si.append( i+n )
                sj.append( i+m )
                sv.append( 1/lmn )
                sil.append( n )
                sjl.append( m )
                sr.append((fn,fm)) # This tells who is allowed to build along the lane
                system.append(k)
                
#                if fn==0 or fm==0:
#                   print(fm,fn)
                
                loc2glob.append(ii)
                ii += 1
                
        i += len(nodes)
        loc2globs.append(loc2glob)

    return ( si, sj, sv, sil, sjl, loc2globs, sr, system )

# Reads all the systems
class Systems:
  def __init__( self, path, factions ):
    assets  = readAssets( '../naev/dat/assets/' )
    anchorSys, anchorJps, anchorAst = createAnchors()
    
    sysdict = {} # This dico will give index of systems
    self.sysnames = [] # List giving the invert of sysdict
    self.nodess = [] # This is a list of nodes in systems
    
    jpnames = [] # List of jump points
    jpdicts = [] # List of dict (invert of jpnames)
    
    autoposs = []
    radius = [] # List of radius of systems
    xlist = []
    ylist = []
    self.presences = [] # List of presences in systems
    
    nodespres = [] # Presence per node
    nodesfact = [] # Faction per node
    
    self.loc2globs = [] # Gives the global numerotation from local one for nodes
    jp2locs = [] # Gives the local numerotation from jump number
    self.ass2g   = [] # Gives the global numerotation from asset one
    g2ass   = [] # Gives asset numerotation from global one
    g2sys   = [] # Gives the system from global numerotation
    
    sysass = [] # Assets names per system
    
    nsys = len(os.listdir(path))
    connect = np.zeros((nsys,nsys)) # Connectivity matrix for systems. TODO : use sparse
    
    self.anchors = []
    presass = [] # List of presences in assets
    factass = [] # Factions in assets. (I'm sorry for all these asses)
    
    i = 0 # No of system
    nglob = 0 # Global nb of nodes
    nasg = 0 # Global nb of assest
    for fileName in os.listdir(path):
        tree = ET.parse((path+fileName))
        root = tree.getroot()
        
        name = root.attrib['name']
        #print(name)
        sysdict[name] = i
        self.sysnames.append(name)
        
        pos = root.find('pos')
        xlist.append( float(pos.find('x').text) )
        ylist.append( float(pos.find('y').text) )
        
        general = root.find('general')
        radius.append( float(general.find('radius').text) )
        
        # Store list of nodes (without jumps)
        nodes = []
        presence = [0] * len(factions)
        assts = root.find('assets')
        loc2glob = []
        sysas = []
        nass = 0
        
        if assts != None:
            aslist = assts.findall('asset')
            for pnt in aslist :
                asname = pnt.text
                sysas.append(asname)
                info = assets[asname]
                if info[3] > 0: # Check the asset is actually inhabited
                    #if info[2] in factions.keys(): # Increment presence
                     #   presence[ factions[info[2]] ] += info[3]
                    if info[0] != None: # Check it's not a virual asset
                        nodes.append( (info[0], info[1]) )
                        loc2glob.append(nglob)
                        self.ass2g.append(nglob)
                        g2ass.append(nasg)
                        g2sys.append(i)
                        if info[2] in factions.keys(): # Store presence
                            presass.append(info[3])
                            factass.append(factions[info[2]])
                                                    
                            #if asname == "Raelid Outpost":
                                #print(len(factass))
                                #print(factions[info[2]])
                        else: # Someone that does not build
                            presass.append(info[3])
                            factass.append(-1)
#                        if (name in anchorSys and (asname in anchorAst: # Add to anchors
#                            self.anchors.append( nglob )
#                        nodespres[nglob] = info[3]
#                        nodesfact[nglob] = info[2]
                        nglob += 1
                        nass += 1
                        nasg += 1
        
        presence = [max(0,j) for j in presence] # Ensure presences are >= 0
        self.presences.append(presence)
        sysass.append(sysas)
        
        
        # Store jump points.
        jpname = []
        jpdict = {}
        autopos = []
        jp2loc = []
        jumps = root.find('jumps')
        jplist = jumps.findall('jump')
        jjj = 0
#        for jj in range(len(jplist)) :
#            jpt = jplist[jj]
        for jpt in jplist :
            
            hid = jpt.find('hidden')
            if hid != None:  # Jump is hidden : don't consider it
                continue

            xit = jpt.find('exitonly')
            if xit != None:  # Jump is exit only : don't consider it
                continue
            
            pos = jpt.find('pos')
            if pos == None: # Autopos is activated: need a second loop
                nodes.append( (None, None) ) # Unknown x and y for now
                autopos.append(True)
            else: # Position is given : create a node
                x = float(pos.attrib['x'])
                y = float(pos.attrib['y'])
                nodes.append( (x, y) )
                autopos.append(False)
                
            jp2loc.append(nass+jjj)
            loc2glob.append(nglob)
            g2ass.append(-1)
            g2sys.append(i)
            nglob += 1
            jpname.append(jpt.attrib['target'])
            jpdict[jpt.attrib['target']] = jjj
            
            jjj += 1
          
        autoposs.append(autopos)
        jpnames.append(jpname)
        jpdicts.append(jpdict)
        jp2locs.append(jp2loc)
        self.nodess.append(nodes)
        self.loc2globs.append(loc2glob)
        i += 1
     
    
    # Second loop over systems in order to create nodes associated to jumps
    # And to populate the connectivity matrix (for ranged presence)
    si0 = [] #
    sj0 = [] #
    sv0 = [] # For the construction of the sparse weighted connectivity matrix
    
    for i in range(nsys):
        
        jpname = jpnames[i]
        autopos = autoposs[i]
        
        loc2globi = self.loc2globs[i]
        jp2loci = jp2locs[i]
        namei = self.sysnames[i]
        #print(namei)
        for j in range(len(jpname)):
            k = sysdict[jpname[j]] # Get the index of target
            connect[i,k] = 1 # Systems connectivity
            connect[k,i] = 1
            
            jpnamek = jpdicts[k]
            loc2globk = self.loc2globs[k]
            jp2lock = jp2locs[k]

            if (namei in anchorSys) and (jpname[j] in anchorJps): # Add to anchors
                self.anchors.append( loc2globi[jp2loci[j]] )
            
            if autopos[j]: # Compute autopos stuff
                theta = math.atan2( ylist[k]-ylist[i], xlist[k]-xlist[i] )
                x = radius[i] * math.cos(theta)
                y = radius[i] * math.sin(theta)
                self.nodess[i][jp2loci[j]] = (x,y) # Now we have the position

            if not ( namei in jpnamek.keys() ):
                continue # It's an exit-only : we don't count this one as a link
            
            m = jpnamek[namei] # Index of system i in system k numerotation
            
            si0.append(loc2globi[jp2loci[j]]) # Assets connectivity (jp only)
            sj0.append(loc2globk[jp2lock[m]]) # 
            sv0.append( 1/1000 ) # TODO : better value
                
    # Remove the redundant info because right now, we have i->j and j->i
    while k<len(si0):
        if (si0[k] in sj0):
            si0.pop(k)
            sj0.pop(k)
            sv0.pop(k)
            k -= 1
        k += 1
    
    # Compute distances. 
    distances = sp.csgraph.dijkstra(connect)

    #print(distances)
    # Use distances to compute ranged presences
    for i in range(nsys):
        for j in range(nsys):
            sysas = sysass[j]
            for k in range(len(sysas)):
                info = assets[sysas[k]] # not really optimized, but should be OK
                if info[2] in factions.keys():
                    fact = factions[ info[2] ]
                    pres = info[3]
                    ran = info[4]
                    d = distances[i,j]
                    if d <= ran:
                        #self.presences[i][fact] += pres/(2**d)
                        self.presences[i][fact] += pres / (1+d)
                        
#        if self.sysnames[i] == "Raelid":#"Alteris":
#            print(i)
#            print(self.presences[i])
        # TODO maybe : ensure positive presence
        

    # Get the stiffness inside systems
    self.internal_lanes = inSysStiff( self.nodess, factass, g2ass, self.loc2globs )

    # Merge both and generate matrix
    si = si0 + self.internal_lanes[0]
    sj = sj0 + self.internal_lanes[1]
    sv = sv0 + self.internal_lanes[2]
    sz = len(si)
    
    # Build the sparse matrix
    sii = [0]*4*sz
    sjj = [0]*4*sz
    svv = [0.]*4*sz
    for k in range(sz):
        sii[4*k]   = si[k]
        sjj[4*k]   = si[k]
        svv[4*k]   = sv[k]
        
        sii[4*k+1] = sj[k]
        sjj[4*k+1] = sj[k]
        svv[4*k+1] = sv[k]
        
        sii[4*k+2] = si[k]
        sjj[4*k+2] = sj[k]
        svv[4*k+2] = - sv[k]
        
        sii[4*k+3] = sj[k]
        sjj[4*k+3] = si[k]
        svv[4*k+3] = - sv[k]

    self.stiff = sp.csr_matrix( ( svv, (sii, sjj) ) )

    self.default_lanes = (si0,sj0,sv0)
    self.assts = (presass,factass)
    self.sysdist = (distances,g2sys)


# Computes the conductivity matrix
def buildStiffness( default_lanes, internal_lanes, activated, alpha, anchors ):
    
    def_si  = default_lanes[0]
    def_sj  = default_lanes[1]
    def_sv  = default_lanes[2]
    
    int_si  = internal_lanes[0]
    int_sj  = internal_lanes[1]
    int_sv0 = internal_lanes[2]
    #print(np.max(np.c_[internal_lanes[2]]))
    # Take activated lanes into account
    int_sv = int_sv0.copy()
    for i in range(len(int_sv0)):
        if activated[i]:
            int_sv[i] = (alpha+1)*int_sv0[i]
    
    # Assemble both lists
    si = def_si + int_si
    sj = def_sj + int_sj
    sv = def_sv + int_sv
    sz = len(si)
    
    # Build the sparse matrix
    sii = [0]*4*sz
    sjj = [0]*4*sz
    svv = [0.]*4*sz
    for k in range(sz):
        sik = si[k]
        sjk = sj[k]
        svk = sv[k]
        
        sii[4*k]   = sik
        sjj[4*k]   = sik
        svv[4*k]   = svk
        
        sii[4*k+1] = sjk
        sjj[4*k+1] = sjk
        svv[4*k+1] = svk
        
        sii[4*k+2] = sik
        sjj[4*k+2] = sjk
        svv[4*k+2] = - svk
        
        sii[4*k+3] = sjk
        sjj[4*k+3] = sik
        svv[4*k+3] = - svk

    # impose Robin condition at anchors (because of singularity)
    mu = max(sv)  # Just a parameter to make a Robin condition that does not spoil the spectrum of the matrix
    for i in anchors:
        sii.append(i)
        sjj.append(i)
        svv.append(mu)
        #stiff[i,i] = stiff[i,i] + mu

    stiff = sp.csr_matrix( ( svv, (sii, sjj) ) )
    # Rem : it may become mandatory at some point to impose nb of dofs

    return stiff

# Gives the matrix that computes penibility form potentials
# By chance, this does not depend on the presence of lane
def PenMat( nass, ndof, internal_lanes, utilde, ass2g, assts, sysdist ):
    presass = assts[0]
    factass = assts[1]
    nfact = max(factass)+1
    
    # First : compute the (sparse) matrix that transforms utilde into u
    pi = []
    pj = []
    pv = []
    
    di = [[] for i in range(nfact)]
    dj = [[] for i in range(nfact)]
    dv = [[] for i in range(nfact)]
    
    distances = sysdist[0]
    g2sys = sysdist[1]
    
    for i in range(nass):
        ai = ass2g[i] # Find corresponding dof
        facti = factass[i]
        #print(facti)
        presi = presass[i]        
        for j in range(i):# Stuff is symmetric, we use that
            if utilde[ai,j] <= 1e-12: # i and j are disconnected. Dont consider this couple
                continue
            ij = i*nass + j # Multiindex
            
            # Just a substraction
            pi.append(i)
            pj.append(ij)
            pv.append(1)
            pi.append(j)
            pj.append(ij)
            pv.append(-1)
            
            factj = factass[j]
            presj = presass[j]
            
            # More the assets are far from each other less interesting it is
            # This avoids strange shapes for peripheric systems
#            aj = ass2g[j]
#            dis = distances[g2sys[ai],g2sys[aj]]
#            if dis==0:
#                dis = 1
            dis = 1
                    
            # Build the diagonal ponderators
            if (facti == factj):
                di[facti].append(ij)
                dj[facti].append(ij)
                dv[facti].append( (presi+presj) / dis )
                
            else: # Foreign assets are not so interesting
                di[facti].append(ij)
                dj[facti].append(ij)
                dv[facti].append( presi / dis )
                
                di[factj].append(ij)
                dj[factj].append(ij)
                dv[factj].append( presj / dis ) 
            
    P = sp.csr_matrix( ( pv, (pi, pj) ) )
    # P*utilde^T = u^T
    
    #print(di[0][100])
    #print(di[1][100])
    
    D = [None]*nfact
    for k in range(nfact):
        D[k] = sp.csr_matrix( ( dv[k], (di[k], dj[k]) ), (P.shape[1], P.shape[1]) )
    #print(lgs.norm(D[0]-D[1], 'fro'))
    
    si = internal_lanes[0]
    sj = internal_lanes[1]
    
    # Then : the matrix that gives penibility on each internal lane from u
    qi = []
    qj = []
    qv = []
    for i in range(len(si)):
        # Just a substraction (again)
        qi.append(i)
        qj.append(si[i]) # TODO : check the transpose
        qv.append(1)
        qi.append(i)
        qj.append(sj[i])
        qv.append(-1)
            
    Q = sp.csr_matrix( ( qv, (qi, qj) ), (len(si),ndof) )
    # Q*u = p
    
    return (P,Q,D)


# Get the gradient from state and adjoint state
# Luckly, this does not depend on the stiffness itself (by linearity wrt. stiffness)
def getGradient( internal_lanes, u, lamt, alpha, PP, PPl, pres_0 ):
    si = internal_lanes[0]
    sj = internal_lanes[1]
    sv = internal_lanes[2]
    sr = internal_lanes[6] # Tells who has right to build on each lane
    sy = internal_lanes[7] # Tells the system
    
    sz = len(si)
#    print(u.shape, sz)
    sh = u.shape
    sh = sh[0]
    
    #lam = lamt.dot(PP) # 0.01 s
    #LUT = lam.dot(u.transpose())
    g = np.zeros((sz,1)) # Just a vector

    nfact = len(PPl)
    gl = []
    #t0 = time.perf_counter()
    #ut = u.transpose()
    for k in range(nfact): # .2
        lal = lamt.dot(PPl[k])
        #LUTl = lal.dot(ut)

        glk = np.zeros((sz,1))

    #for k in range(nfact): 
        for i in range(sz):
            
            if pres_0[sy[i]][k] <= 0: # Does this faction have presence here ?
                continue
            
            # Does this faction have the right to build here ?
            cond = (sr[i][0] == k) or (sr[i][1] == k) or ((sr[i][0] == -1) and (sr[i][1] == -1))
            if not cond: # No need to compute stuff if we've dont have the right to build here !
                continue
            
            sis = si[i]
            sjs = sj[i]
            
            LUTll = np.dot( lal[[sis,sjs],:] , u[[sis,sjs],:].T )
            glk[i] = alpha*sv[i] * ( LUTll[0,0] + LUTll[1,1] - LUTll[0,1] - LUTll[1,0] )
            
            #glk[i] = alpha*sv[i] * ( LUTl[sis, sis] + LUTl[sjs, sjs] - LUTl[sis, sjs] - LUTl[sjs, sis] )
            
        gl.append(glk)
    #t2 = time.perf_counter()
    #print(t2-t0)
    return (g,gl)
        

# Activates the best lane in each system
def activateBest( internal_lanes, g, activated, Lfaction, nodess ):
    # Find lanes to keep in each system
    sv  = internal_lanes[2]
    #sil = internal_lanes[3]
    #sjl = internal_lanes[4]
    lanesLoc2globs = internal_lanes[5]
    nsys = len(lanesLoc2globs)
    
    g1 = g * np.c_[sv] # Short lanes are more interesting
    
    for i in range(nsys):
        lanesLoc2glob = lanesLoc2globs[i]
        #nodes = nodess[i]
        #aloc = activated[lanesLoc2glob]
        aloc = [activated[k] for k in lanesLoc2glob]
        
        if not (False in aloc): # There should be something to actually activate
            continue
        
        gloc = g1[lanesLoc2glob]
        ind1 = np.argsort(gloc.transpose()) # For some reason, it does not work without transpose :/
        ind = [lanesLoc2glob[k] for k in ind1[0]]

        # Find a lane to activate
        for k in ind:#range(len(ind)):
            if activated[k] == False: # One connot activate something that is already active
                break
                    
        #if admissible: # Because it's possible noone was admissible
        activated[k] = True
        Lfaction[k] = 0
        
    return 1


# Activates the best lane in each system for each faction
def activateBestFact( internal_lanes, g, gl, activated, Lfaction, nodess, pres_c, pres_0 ):
    # Find lanes to keep in each system
    sv  = internal_lanes[2]
    #sil = internal_lanes[3]
    #sjl = internal_lanes[4]
    lanesLoc2globs = internal_lanes[5]
    sr = internal_lanes[6]
    nsys = len(lanesLoc2globs)
    
    nfact = len(pres_c[0])
    price = .007  # TODO : tune metric price
    
    g1 = g * np.c_[sv] # Short lanes are more interesting
    
    g1l = [None]*nfact
    for k in range(nfact):
        g1l[k] = gl[k] * np.c_[sv]
        #print(np.linalg.norm(gl[k]))
    
    #print(sv)
    
    for i in range(nsys):
        lanesLoc2glob = lanesLoc2globs[i]
        #nodes = nodess[i]
        #aloc = activated[lanesLoc2glob]
        aloc = [activated[k] for k in lanesLoc2glob]
        
        if not (False in aloc): # There should be something to actually activate
            continue

        # Faction with more presence choose first
        ploc = pres_0[i]
        sind = np.argsort(ploc)
        sind = np.flip(sind)
        
#        if i==320: #116=Alteris 320=Raelid
#            print(pres_c[i])
#            print(sind)
        
        for ff in range(nfact):
            f = sind[ff]

            ntimes = 1 # Nb of lanes per iteration for this faction
            
            for lane in range(ntimes):
                if pres_c[i][f] <= 0.: # This faction has no presence
                    continue
                
                gloc = g1l[f][lanesLoc2glob]
                #gloc = g1[lanesLoc2glob]
                ind1 = np.argsort(gloc.transpose()) # For some reason, it does not work without transpose :/
                ind = [lanesLoc2glob[k] for k in ind1[0]]
        
                # Find a lane to activate
                for k in ind:#range(len(ind)):
                    cond = (sr[k][0] == f) or (sr[k][1] == f) or ((sr[k][0] == -1) and (sr[k][1] == -1))
                    #if (activated[k] == False) and (pres_c[i][f] >= 1/sv[k] * price):     # One connot activate something that is already active           
                    if (activated[k] == False) and (pres_c[i][f] >= 1/sv[k] * price) \
                      and cond:                   
                        pres_c[i][f] -= 1/sv[k] * price
                        activated[k] = True
                        Lfaction[k] = f
                        break
        
    return 1


# Display the active lanes
def printLanes( internal_lanes, activated, Lfaction, nodess, sysnames ):
    
#    lanes2print = ["Arcturus", "Delta Pavonis", "Gamma Polaris", "Goddard",\
#                   "Alteris", "Cygnus", "Za'lek", "Raelid", "Armorhead",\
#                   "Pudas", "Aesir", "Eye of Night", "Pilatis", "Pisces Prime",\
#                   "Fidelis"]
    #lanes2print = ["Arcturus", "Gamma Polaris", "Goddard", "Za'lek", "Alteris", "Pilatis"]
    #lanes2print = ["Eye of Night", "Armorhead", "Gamma Polaris"]
    #lanes2print = ["Raelid"]
    lanes2print = ["Fidelis"]
    
    nsys = len(nodess)
    sil = internal_lanes[3]
    sjl = internal_lanes[4]
    lanesLoc2globs = internal_lanes[5]
    
    for i in range(nsys):
        if not (sysnames[i] in lanes2print):
            continue
        
        print(sysnames[i])
        nodes = nodess[i]
        lanesLoc2glob = lanesLoc2globs[i]
        aloc = [activated[k] for k in lanesLoc2glob]
        #floc = [Lfaction[k] for k in lanesLoc2glob]
        
        plt.figure()
        xlist = [nodes[no][0] for no in range(len(nodes))]
        ylist = [nodes[no][1] for no in range(len(nodes))]
        plt.scatter(xlist, ylist, color='b')
        
        for j in range(len(aloc)):
            if aloc[j]:
                jj = lanesLoc2glob[j]

                no1 = sil[jj]
                no2 = sjl[jj]
            
                x1 = nodes[no1][0]
                x2 = nodes[no2][0]
                y1 = nodes[no1][1]
                y2 = nodes[no2][1]
            
                col = FACTION_COLORS[ Lfaction[jj] ]
                plt. plot([x1,x2], [y1,y2], color=col)
        
        # Pass to global
        plt.show()


# Optimize the lanes
def optimizeLanes( systems, alpha, ndof ):
    sz = len(systems.internal_lanes[0])
    activated = [False] * sz # Initialization : no lane is activated
    Lfaction = [-1] * sz;
    pres_c = systems.presences.copy()
    nfact = len(systems.presences[0])
    
    nass = len(systems.ass2g)
    
    # TODO : It could be interesting to use sparse format for the RHS as well (see)
    ftilde = np.eye( ndof )
    ftilde = ftilde[:,systems.ass2g] # Keep only lines corresponding to assets
    
    niter = 20
    for i in range(niter):
        stiff = buildStiffness( systems.default_lanes, systems.internal_lanes, activated, alpha, systems.anchors ) # 0.02 s

        # Compute direct and adjoint state
        if i >= 1:
            utildp = utilde
        #start = time.perf_counter()
        #print(stiff.count_nonzero())
        utilde = lgs.spsolve( stiff, ftilde ) # 0.11 s
        #end = time.perf_counter()

        # Check stopping condition
        nu = np.linalg.norm(utilde,'fro')
        dr = nu # Just for i==0
        if i >=1:
            dr = np.linalg.norm(utildp-utilde,'fro')
            
        #print(dr/nu)
        if dr/nu <= 1e-12:
            break

        # Compute QQ and PP, and use utilde to detect connected parts of the mesh
        # It's in the loop, but only happends once
        if i == 0: # .5 s
            Pi = PenMat( nass, ndof, systems.internal_lanes, utilde, systems.ass2g, systems.assts, systems.sysdist )
            P = Pi[0]
            Q = Pi[1] 
            D = Pi[2]
            
            QQ = (Q.transpose()).dot(Q)
            #QQ = QQ.todense()
            #start = time.perf_counter()
            PP0 = P.dot(P.transpose())
            #print(PP0.shape)
            PP = PP0.todense() # PP is actually not sparse
            #end = time.perf_counter()
            #print(end-start)
            PPl = [None]*nfact
            for k in range(nfact): # Assemble per-faction ponderators
                Pp = P.dot(D[k])
                PP0 = Pp.dot(Pp.transpose())
                #print(PP0.count_nonzero())
                PPl[k] = PP0.todense() # Actually, for some factions, sparse is better

        rhs = - QQ.dot(utilde)  
        lamt = lgs.spsolve( stiff, rhs ) #.11 s # TODO if possible : reuse cholesky factorization
        
        # Compute the gradient.
        gNgl = getGradient( systems.internal_lanes, utilde, lamt, alpha, PP, PPl, pres_c ) # 0.2 s
        g = gNgl[0]
        gl = gNgl[1]

        # Activate one lane per system
        #activateBest( systems.internal_lanes, g, activated, Lfaction, nodess ) # 0.01 s
        activateBestFact( systems.internal_lanes, g, gl, activated, Lfaction, systems.nodess, pres_c, systems.presences ) # 0.01 s

    #print(np.max(np.c_[systems.internal_lanes[2]]))
        #print(end - start)

    # And print the lanes
    printLanes( systems.internal_lanes, activated, Lfaction, systems.nodess, systems.sysnames )
    print(np.linalg.norm(utilde,'fro'))
    print(i+1)
    #print(np.linalg.norm(utilde-utilde.transpose(),'fro'))
    
    return (activated, Lfaction)

if __name__ == "__main__":
    factions = createFactions()
    systems = Systems( '../naev/dat/ssys/', factions )
    
    alpha = 9. # Efficiency parameter for lanes    
    
    ndof = systems.stiff.shape[0]
    
    # Run the optim algo
    a = time.perf_counter()
    act = optimizeLanes( systems, alpha, ndof )
    b = time.perf_counter()
    print(b-a," s")
    
    activated = act[0]
    Lfaction = act[1]
    # TODO: put that into xml
    
    
def uselessStuff():
    Pi = PenMat( systems.stiff.shape[0], systems.internal_lanes ) # Get the ponderator matrix. TODO : ti should not be stiff.shape, but the nb of assets
    P = Pi[0]
    Q = Pi[1]
    
    mu = systems.stiff.max()  # Just a parameter to make a Robin condition that does not spoil the spectrum of the matrix
    stiffk = systems.stiff.copy()
    for i in systems.anchors:
        stiffk[i,i] = stiffk[i,i] + mu
        
    #A = stiffk.todense()
    
    #theta, U = np.linalg.eig(A)
    #np.eigvals     #null_space(A)
    
    ftilde = np.eye( stiffk.shape[0] ) # TODO : only those corresponding to assets
    utilde = lgs.spsolve( stiffk, ftilde )
    # TODO maybe : at first iterate, identify non-connex parts
    #utt    = np.linalg.solve( A, ftilde[:,1] )
    #lam = np.matmult( np.matmul(np.transpose(Pi),Pi) , utilde )
    QQ = (Q.transpose()).dot(Q)
    #QQ = np.transpose(Q).dot(Q)
    PP = P.dot(P.transpose())
    u1t = (PP.transpose()).dot(utilde)
    u1 = u1t.transpose()    # Manoeuvers in order to have always sparse matrices first
    rhs = - QQ.dot(u1)
    lam = lgs.spsolve( stiffk, rhs ) # TODO : reuse cholesky factorization
    
    g = getGradient( systems.internal_lanes, utilde, lam, alpha )
    
    # Find lanes to keep in each system
    sv, sil, sjl, lanesLoc2globs = systems.internal_lanes[2:6]
    nsys = len(systems.nodess)
    tokeeps = []
    
    g1 = g * np.c_[sv] # Short lanes are more interesting
    
    #n = 0
    lanes2print = ["Arcturus", "Delta Pavonis", "Gamma Polaris", "Goddard", "Alteris", "Cygnus"]
    for i in range(nsys): 
        if not (systems.sysnames[i] in lanes2print):
            continue
        
#        loc2glob = loc2globs[i]
        lanesLoc2glob = lanesLoc2globs[i]
        nodes = systems.nodess[i]
        
        print(systems.sysnames[i])
        
        gloc = g1[lanesLoc2glob]
        ind1 = np.argsort(gloc.transpose()) # For some reason, it does not work without transpose :/
        ind = [lanesLoc2glob[k] for k in ind1[0]]
        
        # Find which lanes to activate. TODO : better way
        ng = gloc.shape[0]
        #nk = int(ng/5) # Keep half of lanes
        nk = min(1,ng)
        tokeep = [ind[k] for k in range(nk)]
        
        plt.figure()
        for j in range(nk):
            no1 = sil[tokeep[j]]
            no2 = sjl[tokeep[j]]
            
            x1 = nodes[no1][0]
            x2 = nodes[no2][0]
            y1 = nodes[no1][1]
            y2 = nodes[no2][1]
            
            plt. plot([x1,x2], [y1,y2])
            
        xlist = [nodes[no][0] for no in range(len(nodes))]
        ylist = [nodes[no][1] for no in range(len(nodes))]
        plt.scatter(xlist,ylist)
        #print(tokeep)
        # Pass to global
        plt.show()
        tokeeps.append(tokeep)
    
    #computeCostFunction( systems.stiff, factions )
    #isti = lg.inv(systems.stiff)
    #e = np.eye( systems.stiff.shape[0] )
    # TODO : find the kernel
    #isti = lgs.spsolve( systems.stiff, e )
    #u, s, vh = lgs.svds(systems.stiff)
    
#    tree = ET.parse('../naev/dat/ssys/alteris.xml')
#    root = tree.getroot()
#    assets = root.find('assets')
