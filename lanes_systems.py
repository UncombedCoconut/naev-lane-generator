import os
import xml.etree.ElementTree as ET


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

class Systems:
    '''Readable representation of the systems.'''
    def __init__( self, path, factions ):
        self.assets  = readAssets( '../naev/dat/assets/' )

        self.sysdict = {} # This dico will give index of systems
        self.sysnames = [] # List giving the invert of self.sysdict
        self.nodess = [] # This is a list of nodes in systems

        self.jpnames = [] # List of jump points
        self.jpdicts = [] # List of dict (invert of self.jpnames)

        self.autoposs = []
        self.radius = [] # List of radius of systems
        self.xlist = []
        self.ylist = []
        self.presences = [] # List of presences in systems

        self.loc2globs = [] # Gives the global numerotation from local one for nodes
        self.jp2locs = [] # Gives the local numerotation from jump number
        self.ass2g   = [] # Gives the global numerotation from asset one
        self.g2ass   = [] # Gives asset numerotation from global one
        self.g2sys   = [] # Gives the system from global numerotation

        self.sysass = [] # Assets names per system

        nsys = len(os.listdir(path))

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
            self.sysdict[name] = i
            self.sysnames.append(name)

            pos = root.find('pos')
            self.xlist.append( float(pos.find('x').text) )
            self.ylist.append( float(pos.find('y').text) )

            general = root.find('general')
            self.radius.append( float(general.find('radius').text) )

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
                    info = self.assets[asname]
                    if info[3] > 0: # Check the asset is actually inhabited
                        #if info[2] in factions.keys(): # Increment presence
                         #   presence[ factions[info[2]] ] += info[3]
                        if info[0] != None: # Check it's not a virual asset
                            nodes.append( (info[0], info[1]) )
                            loc2glob.append(nglob)
                            self.ass2g.append(nglob)
                            self.g2ass.append(nasg)
                            self.g2sys.append(i)
                            if info[2] in factions.keys(): # Store presence
                                presass.append(info[3])
                                factass.append(factions[info[2]])

                                #if asname == "Raelid Outpost":
                                    #print(len(factass))
                                    #print(factions[info[2]])
                            else: # Someone that does not build
                                presass.append(info[3])
                                factass.append(-1)
                            nglob += 1
                            nass += 1
                            nasg += 1

            presence = [max(0,j) for j in presence] # Ensure presences are >= 0
            self.presences.append(presence)
            self.sysass.append(sysas)


            # Store jump points.
            jpname = []
            jpdict = {}
            autopos = []
            jp2loc = []
            jumps = root.find('jumps')
            jplist = jumps.findall('jump')
            jjj = 0
#            for jj in range(len(jplist)) :
#                jpt = jplist[jj]
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
                self.g2ass.append(-1)
                self.g2sys.append(i)
                nglob += 1
                jpname.append(jpt.attrib['target'])
                jpdict[jpt.attrib['target']] = jjj

                jjj += 1

            self.autoposs.append(autopos)
            self.jpnames.append(jpname)
            self.jpdicts.append(jpdict)
            self.jp2locs.append(jp2loc)
            self.nodess.append(nodes)
            self.loc2globs.append(loc2glob)
            i += 1

        self.assts = (presass,factass)
