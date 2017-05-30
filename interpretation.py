import json

#IN fact the vast majority of this file should be built from the metamodel JSON


ALLELE_REGISTRY_ID_KEY = '@id'
#TODO: Build this off the context somehow.
DMWG_INTERPRETATION_TYPE = 'VariantInterpretation'
DMWG_AGENT_TYPE = 'Agent'
DMWG_CONTRIBUTION_TYPE = 'Contribution'
DMWG_CANONICAL_ALLELE_TYPE = 'CanonicalAllele'
DMWG_CONTRIBUTION_KEY = 'contribution'
DMWG_CONTRIBUTION_ONDATE_KEY = 'onDate'
DMWG_INTERPRETER_ROLE = 'CG-contributory-role:interpreter'
DMWG_ASSESSOR_ROLE = 'CG-contributory-role:assessor'
DMWG_CURATOR_ROLE = 'CG-contributory-role:curator'
DMWG_ROLE_KEY = 'role'
DMWG_AGENT_KEY = 'agent'
DMWG_AGENT_NAME_KEY = 'name'
DMWG_CLINICAL_SIGNIFICANCE_KEY = 'clinicalSignificance'
DMWG_VARIANT_KEY = 'variant'
DMWG_TYPE_KEY = 'type'
DMWG_ID_KEY = 'id'
DMWG_CONDITION_KEY = 'condition'
DMWG_CONDITION_KEY = 'condition'
DMWG_CONDITION_TYPE = 'MendelianCondition'
DMWG_DISEASE_TYPE = 'Disease'
DMWG_CONDITION_EXPLANATION_KEY = 'explanation'
DMWG_CONDITION_DISEASE_KEY = 'disease'
DMWG_DISEASE_NAME_KEY = 'name'
DMWG_HAS_EVIDENCE_KEY = 'hasEvidence'
DMWG_HAS_SUPPORTING_EVIDENCE_KEY = 'hasSupportingEvidence'
DMWG_EVIDENCE_LINE_TYPE = 'EvidenceLine'
DMWG_EVIDENCE_LINE_STRENGTH_KEY = 'strength'
DMWG_CRITERION_ASSESSMENT_TYPE = 'CriterionAssessment'
DMWG_CRITERIA_ASSESSMENT_EXPLANATION_KEY = 'explanation'
#Do we need this?
#DMWG_CRITERIA_ASSESSMENT_VARIANT_KEY = 'variant'
DMWG_SUPPORTING_DATA_KEY = 'supportingData'
DMWG_CRITERION_TYPE = 'Criterion' 
DMWG_CRITERION_DESCRIPTION_KEY = 'description'
DMWG_CRITERION_SHORT_DESCRIPTION_KEY = 'shortDescription'
DMWG_CRITERION_DEFAULT_STRENGTH_KEY = 'defaultStrength'
DMWG_CRITERION_TARGET_INTERPRETATION_KEY = 'targetInterpretation'
DMWG_CRITERION_KEY = 'criterion' 
DMWG_OUTCOME_KEY = 'outcome'
DMWG_MET = 'CG-criterion-outcome:met'
DMWG_NOT_MET = 'CG-criterion-outcome:not-met'
DMWG_ALLELE_FREQUENCY_TYPE = 'AlleleFrequency'
DMWG_ASCERTAINMENT_KEY = 'ascertainment'
DMWG_ESP = 'ESP'
DMWG_EXAC = 'ExAC'
DMWG_1000_GENOMES = 'ThousandGenomes'
DMWG_1000_GENOMES_ESP = 'ThousandGenomesESP6500'
DMWG_POPULATION_KEY = 'population'
DMWG_ALLELE_COUNT_KEY = 'alleleCount'
DMWG_ALLELE_NUMBER_KEY = 'alleleNumber'
DMWG_ALLELE_FREQUENCY_KEY = 'alleleFrequency'
DMWG_COMBINED_POP = 'CG-population-type:combined'
#review population types
DMWG_EUROPEAN_AMERICAN_POP = 'CG-population-type:EuropeanAmerican'
DMWG_AFRICAN_AMERICAN_POP = 'CG-population-type:AfricanAmerican'
DMWG_LATINO_POP = 'CG-population-type:Latino'
DMWG_AFRICAN_POP = 'CG-population-type:African'
DMWG_EUROPEAN_POP = 'CG-population-type:European'
DMWG_EAST_ASIAN_POP = 'CG-population-type:EastAsian'
DMWG_FINNISH_POP = 'CG-population-type:Finnish'
DMWG_NON_FINNISH_EURO_POP = 'CG-population-type:NonFinnishEuropean'
DMWG_SOUTH_ASIAN_POP = 'CG-population-type:SouthAsian'
DMWG_OTHER_POP = 'CG-population-type:Other'
DMWG_HOMOZYGOUS_INDIVIDUAL_COUNT_KEY = 'homozygousAlleleIndiviualCount'
DMWG_HETEROZYGOUS_INDIVIDUAL_COUNT_KEY = 'heterozygousAlleleIndiviualCount'
DMWG_CONSERVATION_TYPE = 'Conservation'
DMWG_CONSERVATION_METHOD_KEY = 'method'
DMWG_CONSERVATION_SCORE_KEY = 'score'
DMWG_INSILICO_TYPE = 'InSilicoPrediction'
DMWG_INSILICO_TOOL_KEY = 'tool'
DMWG_INSILICO_VALUE_KEY = 'value'
DMWG_INSILICO_PREDICTION_CLASS_KEY = 'prediction_type'
DMWG_MISSENSE_EFFECT_PREDICTOR = 'CG:me'


class Node:
    def __init__(self):
        pass
    def add_contribution(self,contribution):
        self.data[DMWG_CONTRIBUTION_KEY] = contribution
    #move id tag into entity? Or some kind of flag?
    def get_id(self): 
        try:
            return self.data[DMWG_ID_KEY]
        except:
            return self.data[ALLELE_REGISTRY_ID_KEY]

class Interpretation(Node):
    def __init__(self,iri):
        self.data = {}
        self.data[DMWG_ID_KEY] = iri
        self.data[DMWG_TYPE_KEY] = DMWG_INTERPRETATION_TYPE
        self.data[DMWG_HAS_EVIDENCE_KEY] = []
    #TODO: VALIDATE
    def add_clinical_significance(self, significance):
        self.data[DMWG_CLINICAL_SIGNIFICANCE_KEY] = significance
    def add_variant(self,variant):
        self.data[DMWG_VARIANT_KEY] = variant
    def add_condition(self,condition):
        self.data[DMWG_CONDITION_KEY] = condition
    def add_criterion_assessment(self, assessment, strength):
        evidence_line = EvidenceLine( assessment )
        evidence_line.add_strength(strength)
        self.data[DMWG_HAS_EVIDENCE_KEY].append(evidence_line)

class CriterionAssessment(Node):
    def __init__(self,iri):
        self.data = {}
        self.data[DMWG_ID_KEY] = iri
        self.data[DMWG_TYPE_KEY] = DMWG_CRITERION_ASSESSMENT_TYPE
        self.data[ DMWG_HAS_SUPPORTING_EVIDENCE_KEY ] = []
    def add_outcome(self,outcome):
        self.data[DMWG_OUTCOME_KEY] = outcome
    #TODO: Review with how this is really supposed to look
    def add_criterion(self,criterion):
        self.data[DMWG_CRITERION_KEY] = criterion
    def add_explanation(self, explanation):
        self.data[DMWG_CRITERIA_ASSESSMENT_EXPLANATION_KEY] = explanation
    def add_variant(self,variant):
        self.data[DMWG_VARIANT_KEY] = variant
    def add_information_nodes(self, infolist):
        for info in infolist:
            evidence_line = EvidenceLine( info )
            self.data[DMWG_HAS_SUPPORTING_EVIDENCE_KEY].append(evidence_line)

class EvidenceLine(Node):
    def __init__(self,assessment):
        self.data = {}
        self.data[DMWG_TYPE_KEY] = DMWG_EVIDENCE_LINE_TYPE
        self.data[ DMWG_SUPPORTING_DATA_KEY ] = [ assessment ]
    def add_strength(self, strength):
        self.data[DMWG_EVIDENCE_LINE_STRENGTH_KEY] = strength

class Contribution(Node):
    def __init__(self, agent, when, role):
        self.data = {}
        self.data[DMWG_TYPE_KEY] = DMWG_CONTRIBUTION_TYPE
        self.data[DMWG_AGENT_KEY] = agent
        self.data[DMWG_CONTRIBUTION_ONDATE_KEY] = when
        self.data[DMWG_ROLE_KEY] = role
    def get_agent(self): return self.data[DMWG_AGENT_KEY]
    def get_ondate(self): return self.data[DMWG_CONTRIBUTION_ONDATE_KEY]
    def get_role(self): return self.data[DMWG_ROLE_KEY]

class Agent(Node):
    def __init__(self, iri):
        self.data={}
        self.data[DMWG_ID_KEY] = iri
        self.data[DMWG_TYPE_KEY] = DMWG_AGENT_TYPE
    def set_name(self,name):
        self.data[DMWG_AGENT_NAME_KEY] = name
    def get_name(self): return self.data[DMWG_AGENT_NAME_KEY]

class Disease(Node):
    def __init__(self, iri, name):
        self.data={}
        self.data[DMWG_ID_KEY] = iri
        self.data[DMWG_TYPE_KEY] = DMWG_DISEASE_TYPE
        self.data[DMWG_DISEASE_NAME_KEY] = name

#TODO: Obv need to fix with CodableConcept
#We're using a blank node ID here because we're not putting a deref id, but we have this
#node many places.
class MendelianCondition(Node):
    count = 1
    def __init__(self):
        self.data={}
        self.data[DMWG_ID_KEY] = '_:MendelianCondition_%d' % MendelianCondition.count
        self.data[DMWG_TYPE_KEY] = DMWG_CONDITION_TYPE
        self.data[DMWG_CONDITION_DISEASE_KEY] = []
        MendelianCondition.count += 1
    def add_disease(self, disease):
        self.data[DMWG_CONDITION_DISEASE_KEY].append(disease) 

#What's the point here?  We want variant to be a subclass of Node so that our serializer can know whether
# to rewrite it or not.  We could and perhaps should add all the methods to build the variant correctly, but
# not right now.  
#TODO Pull out strings.
class Variant(Node):
    def __init__(self, ar_rep):
        self.data={}
        for k,v in ar_rep.items():
            self.data[k] = v
    def get_allele( self, reference ):
        allele_list = self.data['genomicAlleles']
        for allele in allele_list:
            if allele['referenceGenome'] == reference:
                return allele['coordinates'][0]['allele'] #why is coords a list?
        raise Exception
    def get_ref_allele( self, reference ):
        allele_list = self.data['genomicAlleles']
        for allele in allele_list:
            if allele['referenceGenome'] == reference:
                return allele['coordinates'][0]['referenceAllele'] #why a list?
        raise Exception


class AlleleFrequency(Node):
    def __init__(self, ascertainment):
        self.data = {}
        self.data[DMWG_TYPE_KEY] = DMWG_ALLELE_FREQUENCY_TYPE
        self.data[DMWG_ASCERTAINMENT_KEY] = ascertainment
    def set_variant(self,variant):
        self.data[DMWG_VARIANT_KEY] = variant
    def set_population(self, population):
        self.data[DMWG_POPULATION_KEY] = population
    def set_allele_count(self, allele_count):
        self.data[DMWG_ALLELE_COUNT_KEY] = allele_count
    def set_allele_number(self, allele_number):
        self.data[DMWG_ALLELE_NUMBER_KEY] = allele_number
    def set_allele_frequency(self,allele_frequency):
        self.data[DMWG_ALLELE_FREQUENCY_KEY] = allele_frequency
    def set_heterozygous_count(self, count):
        self.data[DMWG_HETEROZYGOUS_INDIVIDUAL_COUNT_KEY] = count
    def set_homozygous_count(self, count):
        self.data[DMWG_HOMOZYGOUS_INDIVIDUAL_COUNT_KEY] = count
    def get_allele_count(self): return self.data[DMWG_ALLELE_COUNT_KEY]
    def get_allele_number(self):
        if DMWG_ALLELE_NUMBER_KEY in self.data:
            return self.data[DMWG_ALLELE_NUMBER_KEY] 
        return 0

class Conservation(Node):
    def __init__ (self):
        self.data = {}
        self.data[DMWG_TYPE_KEY] = DMWG_CONSERVATION_TYPE
    def set_variant(self,variant):
        self.data[DMWG_VARIANT_KEY] = variant
    def set_conservation_method(self, method):
        self.data[DMWG_CONSERVATION_METHOD_KEY] = method
    def set_conservation_score(self, score):
        self.data[DMWG_CONSERVATION_SCORE_KEY] = float(score)
    def set_conserved(self, isconserved):
        self.data[DMWG_ISCONSERVED_SCORE_KEY] = isconserved

class InSilicoPrediction(Node):
    def __init__(self):
        self.data = {}
        self.data[DMWG_TYPE_KEY] = DMWG_INSILICO_TYPE
    def set_variant(self,variant):
        self.data[DMWG_VARIANT_KEY] = variant
    def set_prediction_method(self, method):
        self.data[DMWG_INSILICO_TOOL_KEY] = method
    def set_prediction_value(self, score):
        self.data[DMWG_INSILICO_VALUE_KEY] = score
    def set_prediction_type(self, ptype):
        self.data[DMWG_INSILICO_PREDICTION_CLASS_KEY] = ptype

class CodableConcept(Node):
    def __init__(self,ccid):
        self.data = {}
        self.data[DMWG_TYPE_KEY] = DMWG_CODABLE_CONCEPT_TYPE
        self.data[DMWG_ID_KEY] = ccid
    def add_coding(self,coding):
        self.data[DMWG_CODING] = coding

def Coding(Node):
    def __init__(self,codingid,system,code,display):
        self.data = {}
        self.data[DMWG_TYPE_KEY] = DMWG_CODING_TYPE
        self.data[DMWG_ID_KEY] = codingid
        self.data[DMWG_SYSTEM_KEY] = system
        self.data[DMWG_CODE_KEY] = code
        self.data[DMWG_DISPLAY_KEY] = display


#This declutters a bit by only printing the full node the first time we encounter it.  
#But that's not necessarily in the highest node, because of the ordering of the keys.
#Should probably be smarter, but it works for the moment
class InterpretationEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        super(InterpretationEncoder, self).__init__(*args, **kwargs)
        self.written_nodes = set()
    def default(self,obj):
        if isinstance(obj,Node):
            if obj not in self.written_nodes:
                self.written_nodes.add(obj)
                return obj.data
            else:
                try:
                    return obj.data[DMWG_ID_KEY]
                except:
                    return obj.data[ALLELE_REGISTRY_ID_KEY]
        else:
            return json.JSONEncoder.default(self,obj)
