__author__ = 'Helena Mendes-Soares'

import os as os
import csv as csv
import cobra as cobra
import itertools as itertools
import numpy as np
import pandas as pd
import copy as copy
from math import exp

####### Prepare models that will be used in the analysis #######

def getEXrxns(folderLocation = './models/', speciesID = '85962.8', modelType = 'sbml', EXrxnsTag = 'EX_'):

    '''
    Function to extract the IDs of the exchange reactions of a particular model and adds them to a previous list in the file EX_rxns.txt.
    This is needed to make sure that all the diet files contain all the exchange reactions that may be found in the community models.
    This allows us to make sure that every single compound that may be exchanged is available to the models.

    :param folder_location: specification of where the single species models can be be found. Default is singleSpeciesModels folder.
    :param speciesID: genomeID in the NBCI format. Default is 83333.111 (E. coli)
    :param modelType: format the model is in. Can be sbml or xml, mat or json. Default is sbml.
    :param EXrxnsTag: Exchange reactions in models tend to have a specific tag. Takes by default the tag 'EX_'

    '''

    # Determine the full path to the models to be used.
    modelPath = folderLocation + speciesID + '.' + modelType

    # Import the model using cobrapy.
    if modelType == 'sbml' or modelType == 'xml':
        model = cobra.io.read_sbml_model(modelPath)
    elif modelType == 'mat':
        model = cobra.io.load_matlab_model(modelPath)
    else:
        model = cobra.io.load_json_model(modelPath)

    print model.id
    # Get the list of all reactions in the model
    rxns = model.reactions

    # Create an empty list to store only the exchange reactions
    exRxnsNew = []

    # Store only the exchange reactions in ex_rxns from list of all reactions
    for i in rxns:
        if i.id.startswith(EXrxnsTag):
            exRxnsNew.append(i.id)



    # Check to see if an Outputs folder exists, if not, create one.
    if not os.path.exists('./Outputs/'):
        os.makedirs('./Outputs/')


    # Define which file will contain the list of exchange reactions
    exFile = './Outputs/EX_rxns.txt'

    # Create a list with the exchange reactions already contained in the file
    exRxnsCurrent = []

    try:
        with open(exFile, 'r') as file:
            for line in file:
                line = line.rstrip()
                exRxnsCurrent.append(line)

    except:
        print "The file EX_rxns.txt did not exist and will be created in this iteration."
        pass

    # Append the exchange reactions of the current model to the previous reactions
    for rxn in exRxnsNew:
        exRxnsCurrent.append(rxn)

    # Clean list to remove duplicates
    exRxnsCurrent = list(set(exRxnsCurrent))

    # Print the updated list of exchange reactions into the file with the list of exchange reactions.
    # This file is re-created in each run of the function.
    exs = open(exFile, 'w')

    for line in exRxnsCurrent:
        exs.write(line)
        exs.write('\n')

    # Close the file
    exs.close()


###### Create the Diet file for time point 0 #########

def createNewDiet(oldDietFile = './support/DietsOriginal.txt', exRxnsFile = './Outputs/EX_rxns.txt', DietFile = './Outputs/Calculations/Diets0.txt'):

    oldDiet = open(oldDietFile,'r')
    exRxns = open(exRxnsFile,'r')
    Diet = open(DietFile,'w')

    lines = []
    for i in oldDiet:
        i = i.rstrip().split()
        lines.append(i)


    #find the reactions in the exRxns list that match reactions existing in the diet file

    rxns = []

    for j in exRxns:
        j = j.rstrip()
        rxns.append(j)

    matches = []

    for i in lines:
        for j in rxns:
            if j == i[0]:
                matches.append(j)

    #get the reactions from the exRxns that found no match in the diets file. These will have to be added to the diet file with a flux of 0
    no_matches = list(set(rxns) - set(matches))

    for i in no_matches:
        new_line = [i,'NA',0]
        lines.append(new_line)



    oldDiet.close()
    exRxns.close()

    for i in lines:
        print>>Diet, i[0], '\t', i[1], '\t', i[2]

    Diet.close()


###### Tag the reaction and metabolites in each model with the model ID ##################

def tagModel(folderLocation = './models/', speciesID = '85962.8', modelType = 'sbml', outputFolder = './Outputs/taggedModels/'):

    # Open original species model from ModelSEED:
    locationSpeciesModel = folderLocation + '%s.%s' %(speciesID,modelType)

    if modelType == 'sbml' or modelType == 'xml':
        speciesModel = cobra.io.read_sbml_model(locationSpeciesModel)
    elif modelType == 'mat':
        speciesModel = cobra.io.load_matlab_model(locationSpeciesModel)
    else:
        speciesModel = cobra.io.load_json_model(locationSpeciesModel)

    print speciesModel.id

    # Tag all metabolites with speciesID
    for met in range(len(speciesModel.metabolites)):
        old_met = speciesModel.metabolites[met].id
        new_met = speciesID + '_' + old_met

        speciesModel.metabolites[met].id = new_met


    # Tag all reactions with speciesID
    for rxn in range(len(speciesModel.reactions)):
        old_rxn = speciesModel.reactions[rxn].id
        new_rxn = speciesID + '_' + old_rxn

        speciesModel.reactions[rxn].id = new_rxn


    # Some small hack that is required in cobrapy for id of Reaction object to be part of Model object as well.
    speciesModel.repair()

    # Time to export the tagged models. Check if a directory where to put these models exists already. If not, create one.
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)

    # Finally, export the tagged models.
    outputFile = outputFolder + '%s_tagged.sbml' %speciesID

    try:
        cobra.io.write_sbml_model(speciesModel, outputFile)
        print "We finished exporting the file %s_tagged.sbml" % speciesID
    except:
        print "Something is wrong with the model"


def createExternalModels():
    # First, get the exchange reactions
    exRxns = []

    with open('./Outputs/EX_rxns.txt', 'r') as file:
        for line in file:
            line = line.rstrip()
            exRxns.append(line)

    # Then create a metabolic model just with the exchange reactions
    exModel = createEXmodel(exRxns)

    # Reverse the ???????? in the model with just the exchange reactions
    revModel = createReverseEXmodel(exRxns)



####### Create all 2-species community models ######


def allPairComModels(folderLocation = './Outputs/taggedModels/',outputFolder = './Outputs/communityModels/'):
    '''
    This function creates all possible two-species communities for the species in the community.

    :param folderLocation: path to the folder containing the tagged metabolic models of individual species in a SBML format
    :param outputFolder: path to the folder that will store the two-species community metabolic models.
    :return set of two-species community metabolic models
    '''


    # Create a list with all possible combinations of two species
    pairsSpecies = createAllPairs(folderLocation)

    # Get the exchange reactions and reverse exchange reactions models and put them in a Model object
    exModel = cobra.io.read_sbml_model('./Outputs/Exchange_Reactions_Model.sbml')
    revExModel = cobra.io.read_sbml_model('./Outputs/Reverse_Exchange_Reactions_Model.sbml')


    # Go down the list containing pairs of species and create a two-species community model.
    for i in range(len(pairsSpecies)):
        modelA = folderLocation + '%s.sbml' %pairsSpecies[i][0]
        modelA = str(modelA)

        modelB = folderLocation + '%s.sbml' %pairsSpecies[i][1]
        modelB =str(modelB)

        try:
            createCommunityModel(modelA,modelB, exModel, revExModel, outputFolder)
        except Exception as e:
            print e





####### Estimate Biomasses for all species in each 2-species community ######


def runMMinte(comFolder = '../Outputs/comModels/', timePoint = 0, outInter = '../Outputs/Calculations/effects', EXFluxesFile = '../Outputs/Calculations/exchangeFluxes'):
    calculateGR(timePoint, comFolder, EXFluxesFile)
    evaluateInteractions(timePoint,outInter)




####### Turn the output with the growth rates/biomass information into a matrix #######


def createGRmatrix(timePoint, pathToFile = './Outputs/Calculations/outputGRs'):
    genomeAdata = []
    genomeBdata = []
    grAmixdata = []
    grBmixdata = []
    grAsolodata = []
    grBsolodata = []
    pathToFile = pathToFile + str(timePoint) +'.txt'
    with open(pathToFile,'r') as iTableFile:
        next(iTableFile)
        for line in iTableFile:
            line = line.rstrip()
            line = line.split()
            genomeAdata.append(str(line[1]))
            genomeBdata.append(str(line[2]))
            grAmixdata.append(float(line[3]))
            grBmixdata.append(float(line[4]))
            grAsolodata.append(float(line[5]))
            grBsolodata.append(float(line[6]))
        raw_data_df1 = {'Growth_of': genomeAdata, 'In_presence_of': genomeBdata, 'gr': grAmixdata}
        raw_data_df2 = {'Growth_of': genomeBdata, 'In_presence_of': genomeAdata, 'gr': grBmixdata}
        raw_data_df3 = {'Growth_of': genomeAdata, 'In_presence_of': genomeAdata, 'gr': grAsolodata}
        raw_data_df4 = {'Growth_of': genomeBdata, 'In_presence_of': genomeBdata, 'gr': grBsolodata}
        df1 = pd.DataFrame(raw_data_df1, columns=['Growth_of', 'In_presence_of', 'gr'])
        df2 = pd.DataFrame(raw_data_df2, columns=['Growth_of', 'In_presence_of', 'gr'])
        df3 = pd.DataFrame(raw_data_df3, columns=['Growth_of', 'In_presence_of', 'gr'])
        df4 = pd.DataFrame(raw_data_df4, columns=['Growth_of', 'In_presence_of', 'gr'])
        new_df = df1
        new_df = new_df.append(df2, ignore_index=True)
        new_df = new_df.append(df3, ignore_index=True)
        new_df = new_df.append(df4, ignore_index=True)
        new_df = new_df.pivot_table(index = 'Growth_of', columns = 'In_presence_of', values = 'gr')
        new_df.to_csv('./Outputs/Calculations/GRMatrix%s.csv'%timePoint, header=True)
        return new_df





####### Turn the output with the information about the effect of each species on the others into a matrix #######


def createEffectsmatrix(timePoint, pathToFile = './Outputs/Calculations/effects'):
    genomeAdata = []
    genomeBdata = []
    percentChangeA = []
    percentChangeB = []
    fullPath = pathToFile + str(timePoint) + '.txt'
    with open(fullPath, 'r') as iTableFile:
        next(iTableFile)
        for line in iTableFile:
            line = line.rstrip()
            line = line.split()
            genomeAdata.append(str(line[1]))
            genomeBdata.append(str(line[2]))
            percentChangeA.append(float(line[7]))
            percentChangeB.append(float(line[8]))
        raw_data_df1 = {'Percent_change_in_growth_of': genomeAdata, 'Because_of': genomeBdata, 'change': percentChangeA}
        raw_data_df2 = {'Percent_change_in_growth_of': genomeBdata, 'Because_of': genomeAdata, 'change': percentChangeB}
        df1 = pd.DataFrame(raw_data_df1, columns=['Percent_change_in_growth_of', 'Because_of', 'change'])
        df2 = pd.DataFrame(raw_data_df2, columns=['Percent_change_in_growth_of', 'Because_of', 'change'])
        new_df = df1
        new_df = new_df.append(df2, ignore_index=True)
        new_df = new_df.pivot_table(index='Percent_change_in_growth_of', columns='Because_of', values='change')

        new_df = new_df.replace(np.nan,'1', regex = True)

        new_df.to_csv('./Outputs/Calculations/effectsMatrix%s.csv'%timePoint, header=True)
        return new_df





####### Run one iteration of the Leslie-Gower model to determine population size changes #############

def LeslieGowerModel(timePoint, initFile = './Outputs/Calculations/popDensities.txt' , effectsFile = './Outputs/Calculations/effectsMatrix', grMatrixFile = './Outputs/Calculations/GRMatrix', k = 1, timeStep = 0.5):
    '''
    :param t: timepoint we want to start the simulation from (so Nt and we're going to get Nt+1)
    :param initFile: population densities file
    :param effectsFile: matrix of effects of each species on the growth of every other species
    :param grMatrixFile: matrix of growth rates of all the species in the presence of every other species
    :param k: carrying capacity
    :return:
    '''
    effectsFile = effectsFile + str(timePoint) + '.csv'
    grMatrixFile = grMatrixFile +str(timePoint) + '.csv'

    #first create a vector with the sum of the effects. It multiplies the matrix of effects by initial populations conditions times. So do matrix multiplication

    #get the population densities at the beginning of this time point
    init = []
    with open(initFile,'r') as file:
        for line in file:
            line = line.rstrip()
            line = line.split('\t')
            if line[1] == str(timePoint):
                init.append(float(line[2]))
    #print len(init)

    # get matrix of effects
    effects = []
    with open(effectsFile,'r') as file:
        file.readline()
        for line in file:
            line = line.rstrip()
            line = line.split(',')
            #print(line)
            effects.append(line[1:])

    print len(effects)

    #convert init and effects into float objects
    init = np.array(init, dtype = float)
    effects = np.array(effects, dtype = float)

    #actual Effects (21 January 2017). This is the decrease in growth of the focal species due to the others. Previously, we had the total growth of the focal species in the presence of the other, which i don't think made much sense for being in the numeratos.
    #Actually, I'm still not completely sure about this.
    effects = 1 - effects



    #calculate the vector of the total effects of other species on each of our focal species
    sumEffects = np.dot(effects,init) #Ok, this is the right way.

    #get the information relative to how much biomass is created per species under the present conditions (Bt)
    speciesBiomasses = extractBiomass(grMatrixFile) #remember that column 1 is the speciesIDs and column 2 is biomasses

    #get just the biomasses in a vector
    speciesIDs = []
    Bt = []

    for line in speciesBiomasses:
        speciesIDs.append(line[0])
        Bt.append(line[1])

    speciesIDs = speciesIDs[1:]
    Bt = Bt[1:]
    Bt = np.array(Bt, dtype= float)

    #reduce the size of the time step
    Bt = Bt * timeStep

    #create a vector of lambdas
    #FIXME Attention: lbds should be equal to 1 + (Bt over the initial population size used for the calculation of Bt, so 1). This way, Bt will also be independent of populations size
    #The initial calculations were using the below and that may be why everything always seemed so strange
    #lbds = 1 + Bt/init

    #I think it is supposed to be
    lbds = 1 + Bt




    #create a vector of alphas
    alphas = (lbds-1)/k


    #Or, normalize the effects. This has no bases in the model whatsoever, but makes the effects all positive
    #mini = np.amin(effects)
    #maxi = np.amax(effects)
    #effects_norm = (effects - mini)/(maxi-mini)
    #sumEffects = np.dot(effects_norm,init)


    #create a vector with the values of the population sizes at Nt+1
    Nafter = (lbds * init) / (1 + alphas * init + sumEffects)


    #output the informationabout the pop densities after this time point into the popDensitites file
    finalPopFile = initFile

    finalPop = open(finalPopFile,'a')


    for i in range(len(Nafter)):
        #if Nafter[i] < 0 or str(Nafter[i]) == 'inf' or str(Nafter[i]) == 'nan': @FIXME changed to less than 1e-12 because that's the biomass of one cell
        if Nafter[i] < 1e-12 or str(Nafter[i]) == 'inf' or str(Nafter[i]) == 'nan':
            seq = (speciesIDs[i],str(timePoint+1),str(0))
            new_line = '\t'.join(seq)
        else:
            seq = (speciesIDs[i], str(timePoint + 1), str(Nafter[i]))
            new_line = '\t'.join(seq)
        print >>finalPop, new_line


    finalPop.close()




####### Update the metabolite availabilities in the environment given the predicted population changes #########

def findFluxesPerBiomass(timePoint, exchangeFluxesFile = './Outputs/Calculations/exchangeFluxes', exchangePerBiomassFile = './Outputs/Calculations/exchangePerBiomass'):
    import numpy as np


    #open exchageFluxesPerBiomass file where we will dump the new exchange flux values
    exchangePerBiomassFile = exchangePerBiomassFile + str(timePoint) + '.txt'
    exchangePerBiomass = open(exchangePerBiomassFile, 'w')


    #create a dataframe with the uniqueExchangeFluxes for the specific time point
    exchangeFluxes = []
    uniqueExchangeFluxes = []
    speciesIDs = []
    biomassFluxes = []


    exchangeFluxesFile = exchangeFluxesFile + str(timePoint) + '.txt'

    with open(exchangeFluxesFile, 'r') as file:
        for line in file:
            line = line.rstrip()
            line = line.split()
            exchangeFluxes.append(line)


    for line in exchangeFluxes:
        if int(line[0]) == timePoint:
            if line[1] in speciesIDs:
                continue
            else:
                speciesIDs.append(line[1])

    for line in exchangeFluxes:
        if 'biomass' in line[2] or 'Biomass' in line[2]:
                biomassFluxes.append(line[2:])

    for line in exchangeFluxes:
        if int(line[0]) == timePoint:
            if line in uniqueExchangeFluxes:
                continue
            else:
                uniqueExchangeFluxes.append(line)




    new_exchangeFluxes = []
    for i in speciesIDs:
        fluxes = []
        rxns = []
        biomass = []
        for line in uniqueExchangeFluxes:
            if line[1] == i:
                fluxes.append(float(line[3]))
                rxns.append(line[2])
        for item in biomassFluxes:
            if str(i) in item[0]:
                if float(item[1]) <= 1e-12: #@FIXME issue here
                    biomass = float(1e-12)
                    print "Warning, changed biomass to 1e-12 even though it was 0 or less than 0 for species  %s" %str(i)
                else:
                    biomass = float(item[1])

        fluxes = np.array(fluxes)
        new_fluxes = list(fluxes / biomass)

        for j in range(len(rxns)):
            new_line = [str(i), rxns[j], new_fluxes[j]]
            new_exchangeFluxes.append(new_line)


    for i in range(len(new_exchangeFluxes)):
        print >>exchangePerBiomass, timePoint, '\t', new_exchangeFluxes[i][0], '\t', new_exchangeFluxes[i][1], '\t', new_exchangeFluxes[i][2]

    exchangePerBiomass.close()


def newDiet(timePoint = 0, dietFile = './Outputs/Calculations/Diets',popChangesFile = None, rxnChangesFile = None):

    import copy

    if rxnChangesFile == None:
        rxn_changes = totalRxnChanges(timePoint = timePoint)
    else:
        rxn_changes = []
        with open(rxnChangesFile, 'r') as file:
            for line in file:
                line = line.rstrip()
                line = line.split()
                rxn_changes.append(line)


    dietOld = []
    old_dietFile = dietFile + str(timePoint) + '.txt'
    with open(old_dietFile, 'r') as file:
        for line in file:
            line = line.rstrip()
            line = line.split()
            dietOld.append(line)


    dietNew = copy.deepcopy(dietOld)

    rxns_in_change = []
    rxns_in_diet = []

    for i in rxn_changes:
        new_i = i[0]

        #if timePoint == 0:
        #    new_i = i[0]
        #else:
        #    new_i = i[0][2:-2]
        rxns_in_change.append(new_i)

    for j in dietNew:
        rxns_in_diet.append(j[0])

    rxns_not_in_change = list(set(rxns_in_diet) - set(rxns_in_change))

    #@todo, the output is still weird. Fix this from here.
    for rxn in rxns_not_in_change:
        new_line = [rxn, float(0)]
        rxn_changes.append(new_line)

    new_lines = []
    for change in rxn_changes:
        for line in dietNew:
            if change[0] in line[0]:
                try:
                    new_value = float(line[-1]) - float(change[-1]) #replaced + by -. Ok, this makes sense. Ok, had change and line in the worng order so everything was getting positive?
                    if new_value > 0:
                        new_line = [line[0],line[1],0]
                        new_lines.append(new_line)
                    else:
                        new_line = [line[0],line[1], new_value]
                        new_lines.append(new_line)
                    #line.append(new_value)
                except:
                    continue


    new_dietFile = dietFile + str(timePoint+1) + '.txt'


    updateDietFile = open(new_dietFile, 'w')

    print>>updateDietFile, 'RxnID','\t','Name','\t','Flux'

    for line in new_lines:
        print>> updateDietFile, line[0], '\t', line[1], '\t', line[2]

    updateDietFile.close()

    return dietNew

    pass



#########################
### Support functions ###
#########################


def getListOfModels(comFolder):
    '''
    Required by calculateGR function.
    This function creates a list with all the community models that will be used in the analysis. It creates this list by listing the metabolic models in SBML format present in the user specified folder that contains the community models.

    :param comFolder: path to the folder that contains the community metabolic models.
    :return listOfModelPaths: list object containing the filenames for all the models that will be analysed in the function calculateGR
    '''

    import os

    #Create a list object with the paths for all the community model files that are contained in the community models folder
    path = comFolder
    listOfFiles = os.listdir(path)

    listOfModelPaths = []

    for file in listOfFiles:
        if file.endswith('.sbml') or file.endswith('.xml'):
            pathToFile = path + file
            listOfModelPaths.append(pathToFile)

    return listOfModelPaths










############### See if the next two can be replaced with Mike's code ######
#@TODO replace this code with Mike's

def calculateGR(timePoint, comFolder = './Outputs/communityModels/', EXFluxesFileName = './Outputs/Calculations/exchangeFluxesTemp'):
    '''
    In this function we use cobrapy to calculate the growth rates of the two species that make up the
    two species community metabolic models under particular metabolite availability conditions. We start
    by loading the community model (full model) into 3 distinct Model objects with cobrapy. We then change
    the fluxes of the exchange reactions of the external model so they have lower bounds corresponding to
    whichever 'Diet' condition the user specifies. We then run a flux balance analysis on the full model,
    optimizing the biomass reactions of the two species that make up the community at the same time. It
    then remove all reactions whose IDs start with modelA from the model in modelMinusA, thus leaving only
    the reactions from modelB and from the external compartment. It then runs a FBA on it, maximizing the
    biomass reaction for the model with the tag modelB. It then does the same thing but for reactions tagged
    with modelB on modelMinusB. The optimal flux values for the biomass reactions of each species resulting
    from optimization in the full model and in each model containing only one species, which correspond to
    predicted growth rates, are then exported to a table in the tab-delimited text formal to a folder
    chosen by the user.

    :param timePoint
    :param comFolder
    :param EXFluxesFileName
    :return
    '''

    import cobra

    #print 'We will now calculate the growth rates of the two species in a community model in the presence and absence of the other species'
    print 'This is time point %.02f' %timePoint

    outputGRs = './Outputs/Calculations/outputGRs%s.txt' %timePoint
    growthRatesFile = open(outputGRs, 'w')

    EXFluxesFile = EXFluxesFileName + str(timePoint) + '.txt'

    EXFluxesFinal = open(EXFluxesFile, 'w')

    print>> growthRatesFile, 'ModelName', '\t', 'ObjFunctionSpeciesA', '\t', 'ObjFunctionSpeciesB', '\t', 'GRSpeciesAFull', '\t', 'GRSpeciesBFull', '\t', 'GRASolo', '\t', 'GRBSolo'

    # Create a list of all the models that will be analysed
    allModels = getListOfModels(comFolder)

    for item in range(len(allModels)):

        '''
        @summary: load the model, all versions that will be manipulated in the analysis.
        '''

        # Import the models with cobrapy
        modelFull = cobra.io.read_sbml_model(allModels[item])
        modelMinusA = cobra.io.read_sbml_model(allModels[item])
        modelMinusB = cobra.io.read_sbml_model(allModels[item])

        #print 'We successfully loaded the file %s into three different Model objects with ids, %s, %s, and %s. They should all have the same id.' % (allModels[item], modelFull.id, modelMinusA.id, modelMinusB.id)

        # Determine what the objective function is. It should be composed of two reactions, the biomass reactions for each of the species that compose the model. Store the biomass reactions in a new variable to be used later.
        ObjKeys = modelFull.objective.keys()
        idObjKeys = ObjKeys[0].id, ObjKeys[1].id #I don't think this is needed. Is it? Yes it is [line 820].

        #get the tags used in the model to differentiate the reactions and metabolites from each individual model
        modelA_ID = ObjKeys[0].id
        modelB_ID = ObjKeys[1].id


        #get the species ID for modelA and modelB. This was fixed to take the AGORA models
        modelA_temp = modelA_ID.split('_')[:-1]
        modelA = '_'.join(modelA_temp)

        modelB_temp = modelB_ID.split('_')[:-1]
        modelB = '_'.join(modelB_temp)


        # Open the metabolite conditions file ('Diet')
        dietFile = './Outputs/Calculations/Diets'+ str(timePoint) + '.txt' #for bioreactor conditions comment this line and uncomment the next one
        #dietFile = '../SupportFiles/Diets'+ str(0) + '.txt'
        dietValues = open(dietFile, 'r')
        dietValues.readline()


        # Default timepoint is 0. After choosing the timepoint change the lower bounds for the exchange reactions of the external compartment.

        #print 'The diet chosen for this particular run of the functions was %s' %timePoint
        #print 'The diet file open for this particular run of the function is %s' %dietFile
        print 'We are analysing the model %s' %modelFull.id

        for line in dietValues:
            try:
                new_line = line.rstrip('\n').split()
                #new_line[0] = new_line[0][:-3] This is for modelSEED models
                new_line[0] = new_line[0] #this is for agora models
                #print new_line[0]
                newTimePoint = 2
                #print modelFull.reactions.get_by_id(new_line[0]).lower_bound
                modelFull.reactions.get_by_id(new_line[0]).lower_bound = float(new_line[newTimePoint])
                #print modelFull.reactions.get_by_id(new_line[0]).lower_bound
                modelMinusA.reactions.get_by_id(new_line[0]).lower_bound = float(new_line[newTimePoint])
                modelMinusB.reactions.get_by_id(new_line[0]).lower_bound = float(new_line[newTimePoint])
            except:
                continue

        dietValues.close()

        #print 'We finished changing the lower bounds for the fluxes of the exchange reactions in the models to better fit the availability of metabolites for the microbial communities we are simulating the growth of. '

        # Run FBA on Full model
        modelFull.optimize()


        if modelFull.solution.x is not None:
            print 'We finished running the FBA on the full model. The status of the solution is %s and the value of the solution found is %f' % (modelFull.solution.status, modelFull.solution.f)
        else:
            print 'Full model is infeasible'
        # Find which reactions are tagged as being from modelA, store them in a list, then create the modelMinusA, that is, remove all reactions that are part of one of the species in the model.
        # Run FBA on that reduced model.

        #print 'We are now going to remove all reactions that are tagged as being of species A from the model.'



        listSilentRxnsA = []


        for item in modelMinusA.reactions:
            item = str(item)
            if item.startswith(modelA):
                listSilentRxnsA.append(item)

        #print 'We found %d reactions tagged as being of species A' % len(listSilentRxnsA)

        for j in listSilentRxnsA:
            rxn = j.strip()
            deadRxnA = modelMinusA.reactions.get_by_id(rxn)
            deadRxnA.remove_from_model()

        #print 'We are now going to remove all metabolites that are tagged as being of species A from the model.'

        listSilentMetsA = []

        for item in modelMinusA.metabolites:
            item = str(item)
            if item.startswith(modelA):
                listSilentMetsA.append(item)

        #print 'We found %d metabolites tagged as being of species A' % len(listSilentMetsA)

        for j in listSilentMetsA:
            met = j.strip()
            deadMetA = modelMinusA.metabolites.get_by_id(met)
            deadMetA.remove_from_model()

        modelMinusA.repair()

        #print 'We finished removing all reactions and metabolites from the model. The model now has %s reactions, %s metabolites and the reaction on the objective function in %s' % (len(modelMinusA.reactions), len(modelMinusA.metabolites), modelMinusA.objective)


        modelMinusA.optimize()



        if modelMinusA.solution.x is not None:
            print 'We finished running the FBA on the model without reactions of species A. The status of the solution is %s and the value of the solution found is %f' % (modelMinusA.solution.status, modelMinusA.solution.f)
        else:
            print 'Model minus A is infeasible'



        # Find which reactions are tagged as being from modelB, store them in a list, then create the modelMinusB, that is, remove all reactions that are part of one of the species in the model.
        # Run FBA on that reduced model.

        #print 'We are now going to remove all reactions that are tagged as being of species B from the model.'

        listSilentItemsB = []

        for item in modelMinusB.reactions:
            item = str(item)
            if item.startswith(modelB):
                listSilentItemsB.append(item)

        #print 'We found %d reactions tagged as being of species B' % len(listSilentItemsB)

        for j in listSilentItemsB:
            rxn = j.strip()
            deadRxnB = modelMinusB.reactions.get_by_id(rxn)
            deadRxnB.remove_from_model()

        #print 'We are now going to remove all metabolites that are tagged as being of species B from the model.'

        listSilentMetsB = []

        for item in modelMinusB.metabolites:
            item = str(item)
            if item.startswith(modelB):
                listSilentMetsB.append(item)

        #print 'We found %d metabolites tagged as being of species B' % len(listSilentMetsB)

        for j in listSilentMetsB:
            met = j.strip()
            deadMetB = modelMinusB.metabolites.get_by_id(met)
            deadMetB.remove_from_model()

        modelMinusB.repair()

        #print 'We finished removing all reactions and metabolites from the model. The model now has %s reactions, %s metabolites and the reaction on the objective function in %s' % (len(modelMinusB.reactions), len(modelMinusB.metabolites), modelMinusB.objective)


        modelMinusB.optimize()



        if modelMinusB.solution.x is not None:
            print 'We finished running the FBA on the model without reactions of species A. The status of the solution is %s and the value of the solution found is %f' % (modelMinusB.solution.status, modelMinusB.solution.f)
        else:
            print 'Model minus B is infeasible'


        # Get the x_dict values (fluxes) for the reactions listed under idObjKeys for all three models.
        # Output them to a file with a table that has the information about the model, the species ID in the model, and the growth rates of each species in the full model and in isolation.


        '''new part
        Export the fluxes of the exchange reactions and the biomass reaction
        '''


        for i in range(len(modelMinusA.reactions)):
            if modelMinusA.reactions[i].id.startswith('EX_') and 'biomass' not in modelMinusA.reactions[i].id:
            #if modelMinusA.reactions[i].id.startswith('EX_'):
                if modelMinusA.solution.x is None:
                    print >>EXFluxesFinal, timePoint, '\t', modelMinusA.id, '\t', modelB,'\t', modelMinusA.reactions[i].id,'\t', int(0)
                else:
                    print >>EXFluxesFinal, timePoint, '\t', modelMinusA.id, '\t', modelB,'\t', modelMinusA.reactions[i].id,'\t', modelMinusA.solution.x[i]

            elif str(modelB_ID) in modelMinusA.reactions[i].id:
            #elif 'biomass_' in modelMinusA.reactions[i].id: #maybe replace 'biomass' for modelB_ID? and tehn just not output anything that says biomass
                if modelMinusA.solution.x is None or modelMinusA.solution.x[i] < 1e-12:
                    print >> EXFluxesFinal, timePoint, '\t', modelMinusA.id, '\t', modelB, '\t', modelMinusA.reactions[i].id, '\t', int(0)
                else:
                    print >> EXFluxesFinal, timePoint, '\t', modelMinusA.id, '\t', modelB, '\t', modelMinusA.reactions[i].id, '\t', modelMinusA.solution.x[i]
            else:
                continue

        for j in range(len(modelMinusB.reactions)):
            if modelMinusB.reactions[j].id.startswith('EX_') and 'biomass' not in modelMinusB.reactions[j].id:
            #if modelMinusB.reactions[j].id.startswith('EX_'):
                if modelMinusB.solution.x is None:
                    print >>EXFluxesFinal, timePoint, '\t', modelMinusB.id, '\t', modelA,'\t', modelMinusB.reactions[j].id,'\t', int(0)
                else:
                    print >>EXFluxesFinal, timePoint, '\t', modelMinusB.id, '\t', modelA,'\t', modelMinusB.reactions[j].id,'\t', modelMinusB.solution.x[j]

            elif str(modelA_ID) in modelMinusB.reactions[j].id:
            #elif 'biomass_' in modelMinusB.reactions[j].id:
                if modelMinusB.solution.x is None or modelMinusB.solution.x[j] < 1e-12:
                    print >> EXFluxesFinal, timePoint, '\t', modelMinusB.id, '\t', modelA, '\t', modelMinusB.reactions[j].id, '\t', int(0)
                else:
                    print >> EXFluxesFinal, timePoint, '\t', modelMinusB.id, '\t', modelA, '\t', modelMinusB.reactions[j].id, '\t', modelMinusB.solution.x[j]
            else:
                continue


        #@TODO output this next part to a different file
        #print>>EXFluxesFinal, timePoint, '\t', modelFull.id,'\t', modelFull.id, '\t', '%s_biomass' %modelFull.id, '\t', modelFull.solution.f
        #for k in range(len(modelFull.reactions)):
        #    if modelFull.reactions[k].id.startswith('EX_'):
        #        print >> EXFluxesFinal, timePoint, '\t', modelFull.id, '\t', modelFull.id, '\t', modelFull.reactions[k].id, '\t', modelFull.solution.x[k]
            #elif 'biomass_' in modelFull.reactions[k].id:
            #    if modelFull.solution.f < 0.001:
            #        print >> EXFluxesFinal, timePoint, '\t', modelFull.id, '\t', modelFull.id, '\t', 'TotalBiomass', '\t', int(0)
            #    else:
            #        print >> EXFluxesFinal, timePoint, '\t', modelFull.id, '\t', modelFull.id, '\t', 'TotalBiomass', '\t', modelFull.solution.f
        #    else:
        #        continue






        #prepare to export the growth rate values
        ObjA = []
        ObjB = []

        if idObjKeys[0].startswith(modelA):
            ObjA = idObjKeys[0]
        else:
            ObjB = idObjKeys[0]

        if idObjKeys[1].startswith(modelB):
            ObjB = idObjKeys[1]
        else:
            ObjA = idObjKeys[1]

        #print 'We matched the reactions in the objective function to the model they came from. %s was originally from modelA and %s was originally from modelB' % (ObjA, ObjB)

        if modelFull.solution.x is None or modelFull.solution.x_dict[ObjA] < 1e-12:
            grAfull = float(0)
        else:
            grAfull = modelFull.solution.x_dict[ObjA]

        if modelFull.solution.x is None or modelFull.solution.x_dict[ObjB] < 1e-12:
            grBfull = float(0)
        else:
            grBfull = modelFull.solution.x_dict[ObjB]

        #print "We are going to create a table with the information about the growth of A and B alone and in the presence of B and A respectively."



        if modelMinusA.solution.x is not None:
            if ObjA in modelMinusA.solution.x_dict: #this shouldn't happen ever
                if modelMinusA.solution.x is None or modelMinusA.solution.x_dict[ObjA] < 1e-12:
                    grAMinusA = float(0)
                else:
                    grAMinusA = modelMinusA.solution.x_dict[ObjA]
            else:
                grAMinusA = 'Solo'
        else:
            grAMinusA = float(0)

        if modelMinusA.solution.x is not None:
            if ObjB in modelMinusA.solution.x_dict:
                if modelMinusA.solution.x is None or modelMinusA.solution.x_dict[ObjB] < 1e-12:
                    grBMinusA = float(0)
                else:
                    grBMinusA = modelMinusA.solution.x_dict[ObjB]
            else:
                grBMinusA = 'Solo'
        else:
            grBMinusA = float(0)


        if modelMinusB.solution.x is not None:
            if ObjA in modelMinusB.solution.x_dict: #this shouldn't happen ever
                if modelMinusB.solution.x is None or modelMinusB.solution.x_dict[ObjA] < 1e-12:
                    grAMinusB = float(0)
                else:
                    grAMinusB = modelMinusB.solution.x_dict[ObjA]
            else:
                grAMinusB = 'Solo'
        else:
            grAMinusB = float(0)

        if modelMinusB.solution.x is not None:
            if ObjB in modelMinusB.solution.x_dict:
                if modelMinusB.solution.x is None or modelMinusB.solution.x_dict[ObjB] < 1e-12:
                    grBMinusB = float(0)
                else:
                    grBMinusB = modelMinusB.solution.x_dict[ObjB]
            else:
                grBMinusB = 'Solo'
        else:
            grBMinusB = float(0)


        if grAMinusA != 'Solo' or grBMinusB != 'Solo':
            print 'There is a problem with the attribution of growth rate values to their respective species in model {0} {1} {2}.'.format(modelFull.id, grAMinusA, grBMinusB)

        modelID = modelFull.id
        #organisms = modelID.split('X')

        print>> growthRatesFile, modelID, '\t', ObjA, '\t', ObjB, '\t', grAfull, '\t', grBfull, '\t', grAMinusB, '\t', grBMinusA

    #print 'We finished calculating the growth rates of the species in isolation and when in the presence of another species and dumped the information to the file: %s' % growthRatesFile

    growthRatesFile.close()
    EXFluxesFinal.close()


def evaluateInteractions(timePoint, outInter ='./Outputs/Calculations/effects'):
    #print "will now determine the interactions"
    '''
    This function goes over the file with the growth rates of the species that make up a two-species
    community model and determines the kind of interaction occurring in between the two species. The
    types interactions are determined according to the paper by Heinken and Thiele 2015 AEM. These
    are determined based on the amplitude of change in growth rate of species in the presence and
    absence of another species in the community (>10% of change in growth of the particular species
    when in the presence of another species relative to the absence of another species indicates
    significant interaction), and the sign of the change (positive or negative). The information about
    the calculations of change and the type of interaction predicted in each community is added to the
    original table with the growth rates.

    :param inGRs: path to the file with the table listing the growth rates of the two species in a two-species community metabolic model in the presence and absence of another species in the community.
    :param outInter: path to the file that will contain the information contained in the file with growth rates, plus information regarding the the types of interactions predicted to be occurring in the community
    :return outInter: file with the interactions that are predicted to be occurring between species in a two-species community.
    '''

    #logging.info("We will use the information on the growth rates of the species in file %s to determine what kind of interaction is occurring between the organisms. We will output the table of interactions to %s. We will also count how many instances of each type of interaction are found" %outInter)

    inGRs = './Outputs/Calculations/outputGRs%s.txt'%timePoint

    grFile = open(inGRs, 'r')

    interactionsTableFile = open(outInter+ str(timePoint)+'.txt', 'a')

    print>> interactionsTableFile, 'Model', '\t', 'GenomeIDSpeciesA', '\t', 'GenomeIDSpeciesB', '\t', 'GRSpeciesAFull', '\t', 'GRSpeciesBFull', '\t', 'GRASolo', '\t', 'GRBSolo', '\t', 'PercentChangeRawA', '\t', 'PercentChangeRawB'

    next(grFile)


    for item in grFile:

        item = item.replace("_", ".")
        item = item.replace("A.", "")
        item = item.replace(".model", "")

        try:
            itemNew = item.split('\t')
            item = item.rstrip()

            '''
            # Calculation of the effect of the presence of a competing species in the growth rate of species A
            if float(itemNew[5]) != 0:
                percentChangeRawA = (float(itemNew[3]) - float(itemNew[5])) / float(itemNew[5])
            else:
                percentChangeRawA = (float(itemNew[3]) - float(itemNew[5])) / float(1e-25)

            # Calculation of the effect of the presence of a competing species in the growth rate of species B
            if float(itemNew[6]) != 0:
                percentChangeRawB = (float(itemNew[4]) - float(itemNew[6])) / float(itemNew[6])
            else:
                percentChangeRawB = (float(itemNew[4]) - float(itemNew[6])) / float(1e-25)
            '''

            #Alternative calculations. The effect is calculated as how much of the solo growth is the full growth.
            #This makes the effect vary from 0 to infinity

            #Calculation of the effect on growth of A
            if float(itemNew[5]) !=0:
                percentChangeRawA = (float(itemNew[3]))/(float(itemNew[5]))
            else:
                percentChangeRawA = (float(itemNew[3])) / (float(1e25))

            # Calculation of the effect on growth of B
            if float(itemNew[6]) != 0:
                percentChangeRawB = (float(itemNew[4])) / (float(itemNew[6]))
            else:
                percentChangeRawB = (float(itemNew[4])) / (float(1e25))


        except Exception as e:
            print e

        # Create the interactions table. See what the growth rates file looks like, and get the appropriate columns from there, and merge the information with the information in the file with the growth rates.

        print>> interactionsTableFile, item, '\t', percentChangeRawA, '\t', percentChangeRawB

    # Report the counts for each interaction type.


    interactionsTableFile.close()

######################### Finish of replacing with Mike's code










def extractBiomass(grMatrixFile):
    '''
    :param grMatrixFile: file with matrix created using createMatrix file that shows the growth rates of each species in the presence of all others.
    :return: merged: data frame with two columns. the first column has the species IDs and the second column has the biomass created.
    '''

    #open csv file
    matrixData = open(grMatrixFile,'r')

    #get the data from the file
    mtx = []
    for row in matrixData:
        new_row = row.rstrip()
        new_row = new_row.split(',')
        mtx.append(new_row)

    #create an array jsut with the speciesIDs
    speciesIDs = []
    for i in range(len(mtx)):
        speciesIDs.append(mtx[i][0])
    speciesIDs[0] = 'SpeciesIDs'

    #create an array with the biomass created for each species in isolation
    col2 = []
    for i in range(len(mtx)):
        for j in range(len(mtx)):
            if i == j:
                col2.append(mtx[i][j])
    col2[0] = 'Biomass'

    #merge the information matching the biomass values to the speciesID
    merged = []
    for i in range(len(speciesIDs)):
        new_line = [speciesIDs[i],col2[i]]
        merged.append(new_line)


    return merged


def popChanges(timePoint = 0, popDensities = './Outputs/Calculations/popDensities.txt'):

    import pandas as pd
    dataFrame = []

    with open(popDensities, 'r') as file:
        for line in file:
            line = line.rstrip()
            line = line.split()
            dataFrame.append(line)

    #get just the data for calculation of population size changes
    reduced_dataFrame = []

    for line in dataFrame:
        if line[1] == str(timePoint+1) or line[1] == str(timePoint):
            reduced_dataFrame.append(line)
        else:
            continue

    #turn data into a data frame
    reduced_dataFrame = pd.DataFrame(reduced_dataFrame, columns = ['SpeciesID', 'TimePoint', 'PopSize'])

    #turn data frame from long format to wide format
    reduced_dataFrame = reduced_dataFrame.pivot(index = 'SpeciesID', columns = 'TimePoint', values = 'PopSize')

    #create the final data table just with the changes in population size for the time points of interest
    diff = []
    for i in range(len(reduced_dataFrame)):
        delta = float(reduced_dataFrame.iloc[i][1]) - float(reduced_dataFrame.iloc[i][0])
        new_line = [reduced_dataFrame.index[i], delta]
        diff.append(new_line)

    #change the genome IDs to match NCBI format
    final_diff = []

    for i in range(len(diff)):
        item = diff[i][0]
        #item = item.split('_')
        #new_line = [item[0], diff[i][1]]
        new_line = [item,diff[i][1]]
        final_diff.append(new_line)


    return final_diff


def totalRxnChanges(timePoint= 0, exchangePerBiomassFile = './Outputs/Calculations/exchangePerBiomass', popChangesFile = None, timeStep = 0.5):

    exchangePerBio = []
    exchangePerBiomassFile = exchangePerBiomassFile + str(timePoint) + '.txt'

    with open(exchangePerBiomassFile,'r') as file:
        for line in file:
            line = line.rstrip()
            line = line.split()
            exchangePerBio.append(line)


    if popChangesFile == None:
        pop_changes = popChanges()
    else:
        pop_changes = []
        with open(popChangesFile, 'r') as file:
            for line in file:
                line = line.rstrip()
                line = line.split()
                pop_changes.append(line)


    #@TODO add timeStep here (I think)
    new_table = []
    for line in exchangePerBio:
        for species in pop_changes:
            if line[1] in species[0]:
                total_changeSpecies = float(line[3]) * float(species[1]) * timeStep
                new_line = [line[0], line[1], line[2], line[3], total_changeSpecies]
                new_table.append(new_line)


    list_rxns = []

    for i in new_table:
        if i[2] in list_rxns:
            continue
        else:
            list_rxns.append(i[2])

    total_change = []

    for i in list_rxns:
        fluxes = []
        for j in new_table:
            if j[2] == i:
                fluxes.append(j[4])
        fluxes_vector = np.array(fluxes)
        total_flux = float(np.sum(fluxes_vector))
        if total_flux > 0:
            #new_line = [i, 0] #@todo, here, this needs to be changed. Putting 0 here is not updating the diets file to take in consideration what is produced by the bugs
            new_line = [i, total_flux]
            total_change.append(new_line)
        else:
            new_line = [i, total_flux]
            total_change.append(new_line)

    totalChangeFile = './Outputs/Calculations/totalFluxChanges%s.txt' %timePoint
    totalChangesFluxes = open(totalChangeFile,'w')
    for line in total_change:
        print>>totalChangesFluxes, line

    totalChangesFluxes.close()


    return total_change


def compare_models(model1, model2):


    print('Comparing first model {0} to second model {1}'.format(model1.id, model2.id))

    # Compare the compartments.
    for id in model1.compartments:
        if id in model2.compartments:
            if model1.compartments[id] != model2.compartments[id]:
                print('Compartment ID {0} has different names: {1} != {2}'
                      .format(id, model1.compartments[id], model2.compartments[id]))
        else:
            print('Compartment ID {0} is not in second model'.format(id))
    for id in model2.compartments:
        if id not in model1.compartments:
            print('Compartment ID {0} is not in first model'.format(id))

    # Compare the metabolites.
    for metabolite1 in model1.metabolites:
        if model2.metabolites.has_id(metabolite1.id):
            metabolite2 = model2.metabolites.get_by_id(metabolite1.id)
            # if metabolite1.name != metabolite2.name:
            #     print('Metabolite ID {0} has different names: {1} != {2}'
            #           .format(metabolite1.id, metabolite1.name, metabolite2.name))
            if metabolite1.compartment != metabolite2.compartment:
                print('Metabolite ID {0} has different compartments: {1} != {2}'
                      .format(metabolite1.id, metabolite1.compartment, metabolite2.compartment))
            # if metabolite1.formula != metabolite2.formula:
            #     print('Metabolite ID {0} has different formulas: {1} != {2}'
            #           .format(metabolite1.id, metabolite1.formula, metabolite2.formula))
        else:
            print('Metabolite ID {0} is not in second model'.format(metabolite1.id))
    for metabolite2 in model2.metabolites:
        if not model1.metabolites.has_id(metabolite2.id):
            print('Metabolite ID {0} is not in first model'.format(metabolite2.id))

    # Compare the reactions.
    # Special consideration for exchange reactions
    for reaction1 in model1.reactions:
        if model2.reactions.has_id(reaction1.id):
            reaction2 = model2.reactions.get_by_id(reaction1.id)
            # if reaction1.name != reaction2.name:
            #     print('Reaction ID {0} has different names: {1} != {2}'
            #           .format(reaction1.id, reaction1.name, reaction2.name))
            if reaction1.reaction != reaction2.reaction:
                print('Reaction ID {0} has different definitions: {1} != {2}'
                      .format(reaction1.id, reaction1.reaction, reaction2.reaction))
        else:
            print('Reaction ID {0} is not in second model'.format(reaction1.id))
    for reaction2 in model2.reactions:
        if not model1.reactions.has_id(reaction2.id):
            print('Reaction ID {0} is not in first model'.format(reaction2.id))

    #Compare the reaction bounds
    #Lower bounds
    for reaction1 in model1.reactions:
        if model2.reactions.has_id(reaction1.id):
            reaction2 = model2.reactions.get_by_id(reaction1.id)
            if reaction1.lower_bound != reaction2.lower_bound:
                print('Lower bound of reaction {0} is different between the two reactions: {1} != {2}'.format(reaction1.id, reaction1.lower_bound, reaction2.lower_bound))

    #Upper bounds
    for reaction1 in model1.reactions:
        if model2.reactions.has_id(reaction1.id):
            reaction2 = model2.reactions.get_by_id(reaction1.id)
            if reaction1.upper_bound != reaction2.upper_bound:
                print('Lower bound of reaction {0} is different between the two reactions: {1} != {2}'.format(reaction1.id, reaction1.upper_bound, reaction2.upper_bound))

    #Compare fluxes
    for reaction1 in model1.reactions:
        if model2.reactions.has_id(reaction1.id):
            reaction2 = model2.reactions.get_by_id(reaction1.id)
            if model1.solution.x_dict[reaction1.id] != model2.solution.x_dict[reaction2.id]:
                print('Solutions for reaction {0} are different: In model1 it has solution {1} and in model2 it has solution {2}.'.format(reaction1.id, model1.solution.x_dict[reaction1.id], model2.solution.x_dict[reaction2.id]))

    return


def createAllPairs(folderLocation = './models/'):

    '''
    This function lists the files for species metabolic models that exist in the folder specified in folderLocation and creates a list with all possible pairwise combinations of species models..

    :param folderLocation: path to the folder containing all the individual species metabolic models.
    :return pairs: list with all possible pairwise combinations of the species, that all combinations without replacement.

    '''


    # List all the IDs of all the species that have a model in the original models folder.

    listSpecies = []

    for file in os.listdir(folderLocation):
        if file.endswith('.mat') or file.endswith('.xml'):
            listSpecies.append(file[:-4])
        elif file.endswith('.sbml') or file.endswith('.json'):
            listSpecies.append(file[:-5])

    # Create all pairwise combinations of species IDs (filenames)
    pairs = []

    for c in itertools.combinations(listSpecies,2):
        pairs.append(c)

    return pairs


def createEXmodel(exRxns, outputFolder = './Outputs/'):
    '''
    This function takes the list of exchange reactions in the file EX_rxns.txt and creates a Model object using cobrapy composed of those reactions with the upper bound flux values of 1000, lower bound flux values of -1000, and objective coefficient of 0, and one metabolite as being uptaken by the reaction (stoichiometric coefficient of -1). This is a model composed solely of exchange reactions and it's the model for the extra compartment created for the full community model

    :param exRxns: list of exchange reactions of all the species in the community
    :return SBML model with only the exchange reactions. It's the metabolic model of the extra compartment.
    '''


    # Create the new model name
    exchangeModel = cobra.Model('Model with the exchange reactions only')
    exchangeModel.id = 'Exchange_Reactions_Model'

    # Add the reactions and metabolites to the model
    for i in exRxns:
        new_i = str(i)
        new_i = new_i[3:]
        new_met = cobra.Metabolite(new_i)

        rxn = cobra.Reaction(i)
        rxn.lower_bound = -1000.000
        rxn.upper_bound = 1000.000
        rxn.objective_coefficient = 0.000
        rxn.add_metabolites({new_met:-1.0})

        exchangeModel.add_reaction(rxn)

    # Export the model to a SBML file
    outputFile = outputFolder + str(exchangeModel.id) + '.sbml'
    cobra.io.write_sbml_model(exchangeModel, outputFile)


def createReverseEXmodel(exRxns, outputFolder = './Outputs/'):
    '''
    This function takes the list of exchange reactions in the file EX_rxns.txt and creates a Model object using cobrapy composed of those reactions with the upper bound flux values of 1000, lower bound flux values of -1000, and objective coefficient of 0, and one metabolite as being produced by the reaction (stoichiometric coefficient of 1). This is a model composed solely of exchange reactions and will be used to add the metabolite information to the excahnge reactions of species models in the 2-species community models.

    :param exRxns: list of exchange reactions of all the species in the community
    :return SBML model with only the exchange reactions. It's the metabolic model of the extra compartment.
        '''

    # Create the new model name
    revExchangeModel = cobra.Model('Model with the exchange reactions only')
    revExchangeModel.id = 'Reverse_Exchange_Reactions_Model'

    # Add the reactions and metabolites to the model
    for i in exRxns:
        new_i = str(i)
        new_i = new_i[3:]
        new_met = cobra.Metabolite(new_i)

        rxn = cobra.Reaction(i)
        rxn.lower_bound = -1000.000
        rxn.upper_bound = 1000.000
        rxn.objective_coefficient = 0.000
        rxn.add_metabolites({new_met: 1.0})

        revExchangeModel.add_reaction(rxn)

    # Export the model to a SBML file
    outputFile = outputFolder + str(revExchangeModel.id) + '.sbml'
    cobra.io.write_sbml_model(revExchangeModel, outputFile)


def addEXMets2SpeciesEX(reverseEXmodel, speciesModel):
    '''
    This function takes the model with exchange reactions where the metabolite is produced (output from function createReverseEXmodel) and a species model, and adds the metabolite from the reverse model to the exchange reactions of the species model. For instance:

    Reaction :  modelB_EX_cpd11588_e0 got the cpd11588_e0[u] added.
                 'model_B_cpd11588_e0 <=> cpd11588_e0[u]'

    This way, when a compound is exported to the extracellular environment, it is automatically transformed into a form that is common to all members in the community.

    :param reverseEXmodel: cobrapy Model object containing only exchange reactions with the production of their respective metabolites
    :param speciesModel: Model object of a particular species.
    :return speciesModel: Model object of a particular species with updated exchange reactions.
    '''

    for j in range(len(reverseEXmodel.reactions)):
        exRxn = str(reverseEXmodel.reactions[j])

        for i in range(len(speciesModel.reactions)):
            rxn = str(speciesModel.reactions[i])
            if exRxn in rxn:
                new_met = reverseEXmodel.reactions[j].metabolites
                speciesModel.reactions[i].add_metabolites(new_met)
                speciesModel.reactions[i].lower_bound = -1000.000
                speciesModel.reactions[i].upper_bound = 1000.000

    return speciesModel


def createCommunityModel(modelFileA, modelFileB, exModelCobra, revModelCobra, comFolder = './Outputs/communityModels'):
    '''
    This function takes advantage of the outputs of all the functions defined previously to actually piece together the individual species models and the extra compartment model.


    :param modelFileA: path to the metabolic model of species A in SBML format
    :param modelFileB: path to the metabolic model of species B in SBML format
    :param comFolder:  path to the folder where the metabolic models of the two-species communities will be stored.

    :return two-species community model: in SBML format exported to the folder designated by the user (comFolder) to store these models
    '''

    exModel = exModelCobra
    revExModel = revModelCobra

    #import the model file into the Model object model1 using cobrapy
    try:
        if modelFileA.endswith('.mat'):
            model1 = cobra.io.load_matlab_model(modelFileA)
        elif modelFileA.endswith('.xml') or modelFileA.endswith('.sbml'):
            model1 = cobra.io.read_sbml_model(modelFileA)
        elif modelFileA.endswith('.json'):
            model1 = cobra.io.load_json_model(modelFileA)
        else:
            print "not able to find model %s" %modelFileA
    except Exception as e:
        print e

    #import the model file into the Model object model2 using cobrapy
    if modelFileB.endswith('.mat'):
        model2 = cobra.io.load_matlab_model(modelFileB)
    elif modelFileB.endswith('.xml') or modelFileB.endswith('.sbml'):
        model2 = cobra.io.read_sbml_model(modelFileB)
    elif modelFileB.endswith('.json'):
        model2 = cobra.io.load_json_model(modelFileB)
    else:
        print "not able to find model %s" %modelFileB



    # Create a communityID to identify the output files belonging to each 2-species community created
    communityID = model1.id+ 'X' + model2.id


    # Add the metabolites of the external model to the exchange reactions of each species.

    new_m1 = addEXMets2SpeciesEX(revExModel,model1)
    new_m2 = addEXMets2SpeciesEX(revExModel,model2)


    mix = new_m1
    mix.id = communityID
    mix.add_reactions(new_m2.reactions)
    mix.add_metabolites(new_m2.metabolites)
    mix.add_reactions(exModel.reactions)
    mix.add_metabolites(exModel.metabolites)

    # If the community folder doesn't exist, create one
    if not os.path.exists(comFolder):
        os.makedirs(comFolder)


    # Export the newly created community model to its folder. The models should then be ready to be further analyzed on Widget 5
    cobra.io.write_sbml_model(mix, "%s/community%s.sbml" %(comFolder,communityID))

def cleanExFluxes(timePoint, ExFluxesFile = './Outputs/Calculations/exchangeFluxesTemp', cleanExFluxesFile = './Outputs/Calculations/exchangeFluxes'):

    # If I recall, this is selecting the first exchange flux for a metabolite.

    ExFluxesPath = ExFluxesFile + str(timePoint) + '.txt'
    cleanExFluxesPath = cleanExFluxesFile + str(timePoint) + '.txt'

    cleanEx = []
    ExFluxes = open(ExFluxesPath,'r')
    cleanExFluxes = open(cleanExFluxesPath,'w')

    tempList = []

    for line in ExFluxes:
        line = line.rstrip().split()
        identifier = [line[0],line[2],line[3]]
        new_line = 'XXX'.join(identifier), line[4]
        tempList.append(new_line)

    identifierList = []

    for i in tempList:
        identifierList.append(i[0])

    identifierList = list(set(identifierList))

    for i in identifierList:
        for j in tempList:
            if i == j[0]:
                i = i.split('XXX')
                new_ex = i[0],i[1],i[2],j[1]
                cleanEx.append(new_ex)
                break

    for i in cleanEx:
        print>>cleanExFluxes, i[0],'\t',i[1],'\t',i[2],'\t',i[3]

    cleanExFluxes.close()







