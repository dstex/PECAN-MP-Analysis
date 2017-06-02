def getPECANparams(flight, probe):
    
    if flight == '20150617':
        # Start and end time of spirals (in seconds)
        startT = [14203, 16750, 18576, 20353, 21846, 23588, 24995]
        endT = [14917, 17948, 19577, 21440, 22905, 24570, 25895]

        # Flight-level data index pertaining to start and end of spirals
        startFLix = [18317, 20864, 22690, 24467, 25960, 27702, 29109]
        endFLix = [19031, 22062, 23691, 25554, 27019, 28684, 30009]

        # Time/[temp at same time] when first evidence of ice was observed
        mlBotTemp = []

        mlBotTime = []

        # MCS stage of evolution. 
        # F = Formative, M = Mature, W = Weakening
        mcsStg = []

        # Location of spiral relative to MCS structure. 
        # T = Transition Zone, S = Enhanced Stratiform Region, A = Rear Anvil
        sprlZone = []
        
        if probe == 'CIP':
            intar_threshold = [1.9e-5, 1.2e-5, 7.8e-5, 3.7e-6, 7e-6, 8.7e-6, 3e-6]
        
        elif probe == 'PIP':
            intar_threshold = []
        
    elif flight == '20150620':
        # Start and end time of spirals (in seconds)
        startT = [17757, 18230, 19284, 20845, 21877, 23890, 27947]
        endT = [18229, 18783, 20227, 21826, 22968, 24451, 28702]

        # Flight-level data index pertaining to start and end of spirals
        startFLix = [19948, 20421, 21475, 23036, 24068, 26081, 30138]
        endFLix = [20420, 20974, 22418, 24017, 25159, 26642, 30893]

        # Time/[temp at same time] when first evidence of ice was observed
        mlBotTemp = []

        mlBotTime = []

        # MCS stage of evolution. 
        # F = Formative, M = Mature, W = Weakening
        mcsStg = []

        # Location of spiral relative to MCS structure. 
        # T = Transition Zone, S = Enhanced Stratiform Region, A = Rear Anvil
        sprlZone = []
        
        if probe == 'CIP':
            intar_threshold = [-999, 7e-6, 3e-5, 6e-6, 6e-6, 1.5e-5, 3e-6] # Need value for new first spiral
        
        elif probe == 'PIP':
            intar_threshold = []
    
    
    elif flight == '20150701':
        # Start and end time of spirals (in seconds)
        startT = [21882]
        endT = [22838]

        # Flight-level data index pertaining to start and end of spirals
        startFLix = [19653]
        endFLix = [20609]

        # Time/[temp at same time] when first evidence of ice was observed
        mlBotTemp = []

        mlBotTime = []

        # MCS stage of evolution. 
        # F = Formative, M = Mature, W = Weakening
        mcsStg = []

        # Location of spiral relative to MCS structure. 
        # T = Transition Zone, S = Enhanced Stratiform Region, A = Rear Anvil
        sprlZone = []
        
        if probe == 'CIP':
            intar_threshold = [3e-5]
        
        elif probe == 'PIP':
            intar_threshold = []
    
    
    elif flight == '20150702':
        # Start and end time of spirals (in seconds)
        startT = [13336, 15350, 16892]
        endT = [14551, 16501, 18196]

        # Flight-level data index pertaining to start and end of spirals
        startFLix = [13383, 15397, 16939]
        endFLix = [14598, 16548, 18243]

        # Time/[temp at same time] when first evidence of ice was observed
        mlBotTemp = []

        mlBotTime = []

        # MCS stage of evolution. 
        # F = Formative, M = Mature, W = Weakening
        mcsStg = []

        # Location of spiral relative to MCS structure. 
        # T = Transition Zone, S = Enhanced Stratiform Region, A = Rear Anvil
        sprlZone = []
        
        if probe == 'CIP':
            intar_threshold = [1.7e-4, 3.9e-5, 3.9e-5]
        
        elif probe == 'PIP':
            intar_threshold = []
    
    
    elif flight == '20150706':
        # Start and end time of spirals (in seconds)
        startT = [11987, 12653, 15824, 16805, 20419, 21384, 22950, 23879]
        endT = [12652, 13534, 16804, 17615, 21349, 22393, 23878, 24875]

        # Flight-level data index pertaining to start and end of spirals
        startFLix = [11789, 12455, 15626, 16607, 20221, 21186, 22752, 23681]
        endFLix = [12454, 13336, 16606, 17417, 21151, 22195, 23680, 24677]

        # Time/[temp at same time] when first evidence of ice was observed
        mlBotTemp = [-999, 1.9967, 3.4365, 4.9645, 1.9219, 3.5011, 2.8945, 3.6202] # Need value for new first spiral

        mlBotTime = [-999, 13098, 16234, 17252, 20851, 21872, 23424, 24383] # Need value for new first spiral
        
        # Time/[temp at same time] when first evidence of liquid water was observed
        mlTopTemp = [-999, 0.8366, 0.03, 0.2736, 0.329, 0.5315, 0.122, 0.8] # Need value for new first spiral

        mlTopTime = [-999, 13053, 16337, 17157, 20914, 21796, 23492, 24293] # Need value for new first spiral

        # MCS stage of evolution. 
        # F = Formative, M = Mature, W = Weakening
        mcsStg = ['F','F', 'F', 'F', 'M', 'M', 'M', 'M']

        # Location of spiral relative to MCS structure. 
        # T = Transition Zone, S = Enhanced Stratiform Region, A = Rear Anvil
        sprlZone = ['T','T', 'T', 'T', 'S', 'S', 'S', 'S']
        
        if probe == 'CIP':
            intar_threshold = [-999, 2.7e-5, 2.7e-5, 3.4e-5, 3.3e-6, 9e-6, 2.5e-5, 2.5e-5] # Need value for new first spiral
        
        elif probe == 'PIP':
            intar_threshold = []
    
    
    elif flight == '20150709':
        # Start and end time of spirals (in seconds)
        startT = [8857, 9584, 10643, 11371, 12769, 13665, 14517, 15379, 16585, 17482,
                  19045, 19960, 20891, 21834, 22729, 23546]
        endT = [9582, 10417, 11369, 12160, 13532, 14516, 15377, 16150, 17481, 18226,
                19959, 20687, 21833, 22363, 23544, 24074]

        # Flight-level data index pertaining to start and end of spirals
        startFLix = [12024, 12751, 13810, 14538, 15936 ,16832 ,17684 ,18546 ,19752, 20649,
                     22212, 23127, 24058, 25001, 25896, 26713]
        endFLix = [12749, 13584, 14536, 15327, 16699, 17683, 18544, 19317, 20648, 21393,
                   23126, 23854, 25000, 25530, 26711, 27241]

        # Time/[temp at same time] when first evidence of ice was observed
        mlBotTemp = [0.827, 2.76, 1.584, 2.2527, 2.8054, 4.7628, 3.823, 5.7988, 3.7277,
                     5.6487, 2.3316, 4.7067, 3.6857, 4.9586, 2.5909, 5.7844]

        mlBotTime = [9166, 10102, 10997, 11769, 13072, 14202, 14867, 15861, 16956, 17932,
                     19429, 20445, 21312, 22169, 23075, 23942]
        
        # Time/[temp at same time] when first evidence of liquid water was observed
        mlTopTemp = [0.0104, 0.0034, 0.0423, 0.1119, 0.4757, 3.4426, 1.3057, 2.5534, 0.4398,
                     2.0643, 0.0656, 1.5315, 0.6095, 2.740, 0.6454, 2.5692]

        mlTopTime = [9185, 10028, 11036, 11711, 13146, 14163, 14903, 15820, 17038,
                     17866, 19494, 20391, 21391, 22136, 23121, 23882]

        # MCS stage of evolution. 
        # F = Formative, M = Mature, W = Weakening
        mcsStg = ['F', 'F', 'F', 'F', 'F', 'M', 'M', 'M', 'M',
                  'M', 'M', 'M', 'M', 'M', 'M', 'M']

        # Location of spiral relative to MCS structure. 
        # T = Transition Zone, S = Enhanced Stratiform Region, A = Rear Anvil
        sprlZone = ['T', 'T', 'T', 'T', 'S', 'S', 'S', 'S', 'S',
                    'S', 'S', 'S', 'S', 'S', 'S', 'S']
        
        if probe == 'CIP':
            intar_threshold = [6.9e-5, 5.4e-5 ,6.9e-5, 3.4e-5, 3.4e-5, 2.7e-5, 1e-4, 8.7e-5, 4.3e-5,
                               1e-4, 3.4e-5, 3.4e-5, 3.4e-5, 2.7e-5, 3.4e-5, 3.4e-5]
        
        elif probe == 'PIP':
            intar_threshold = []

    
    return startT, endT, startFLix, endFLix, intar_threshold, mlBotTime, mcsStg, sprlZone