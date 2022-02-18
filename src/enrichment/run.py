from enrich import expressionLoader, convertLogFC, interactionLoader, graphCreate, addWeightsFromMap, edgeWeightAdd, \
    nodeWeightAdd, recalNodeWeights, printing, concordance
from src.enrichment.constants import EXPRESSION, INTERACTION, PATHWAY

result = expressionLoader(EXPRESSION)

result = convertLogFC(result)
gene_interaction_map = interactionLoader(INTERACTION)

G = graphCreate(PATHWAY)
G = addWeightsFromMap(G)
G = edgeWeightAdd(G)
G = nodeWeightAdd(G, result)
G = recalNodeWeights(G)
printing(G)
concordance(G, result)


