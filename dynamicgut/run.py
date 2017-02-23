'''
After:
getEXrxns
createNewDiet
tagModel
createExternalModels
allPairComModels
'''

from .helpers import calculateGR, evaluateInteractions, createGRmatrix, createEffectsmatrix, LeslieGowerModel, \
    cleanExFluxes, findFluxesPerBiomass, newDiet, getEXrxns, tagModel, createExternalModels, allPairComModels


# how about prepare(), run(), visualize()
# can you run in the same folder again?

def prepare(analysis_folder):
    models = [
        'Bacteroides_caccae_ATCC_43185',
        'Bacteroides_ovatus_ATCC_8483',
        'Bacteroides_thetaiotaomicron_VPI_5482',
        'Clostridium_symbiosum_ATCC_14940',
        'Collinsella_aerofaciens_ATCC_25986',
        'Desulfovibrio_piger_ATCC_29098',
        'Escherichia_coli_str_K_12_substr_MG1655',
        'Eubacterium_rectale_ATCC_33656',
        'Marvinbryantia_formatexigens_I_52_DSM_14469'
    ]

    for m in models:
        getEXrxns(analysis_folder, m, 'xml')

    # Think I can use Diets0.txt

    for m in models:
        tagModel(analysis_folder, m, 'xml')

    createExternalModels()

    allPairComModels()

timeInterval = range(0, 10)

# Can replace tagModel, allPairComModels, calculateGR, evaluateInteractions, with mminte code.
# Where to put make_media_from_models()?

# Would be best to pass around DataFrames with results.

def runSim(timeInterval):
    for i in timeInterval:
        calculateGR(i)
        evaluateInteractions(i)
        createGRmatrix(i)
        createEffectsmatrix(i)
        LeslieGowerModel(i)
        cleanExFluxes(i)  # fast fix. Not correct.
        findFluxesPerBiomass(i)
        newDiet(i)
