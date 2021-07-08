import matplotlib.pyplot as plt

FACTION_COLORS = ["g","orange","brown","darkred","silver","aqua","y","b","purple","grey"]


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
