from interpretation_generated import *
from coding_generated import *
import json
from CodingFactory import the_factory

#UtilityMethods for wrapping things in Evidence Lines
def add_criterion_assessment(interpretation, assessment, strength):
    evidence_line = EvidenceLine()
    evidence_line.add_information( assessment )
    evidence_line.set_evidenceStrength( strength )
    interpretation.add_evidence( evidence_line )

def add_informations( assessment, information_list):
    for information in information_list:
        evidence_line = EvidenceLine()
        evidence_line.add_information( information )
        assessment.add_evidence( evidence_line )

DMWG_CURATOR_ROLE = 'curator'
DMWG_INTERPRETER_ROLE = 'interpreter'
DMWG_ASSESSOR_ROLE = 'assessor'
#Utility methods for creating contributions
def create_contribution(agent, ondate, role):
    contribution = Contribution()
    contribution.set_agent(agent)
    contribution.set_onDate(ondate)
    #roleconcept = CodeableConcept()
    # Populate from somewhere?
    #coding = create_coding('http://clinicalgenome.org/datamodel/contributory-role',role,role)
    #roleconcept.add_coding(coding)
    contribution.set_role(role)
    return contribution

#Utility method for creating diseases
def create_orphanet_disease(code, name):
    cc = CodeableConcept()
    coding = create_coding('http://www.orpha.net/', name, code )
    cc.add_coding(coding)
    return cc

#Utility method to make sure that coding gets the ID set correctly
def create_coding(system,display,code):
    iri = system+code
    coding = Coding(iri)
    coding.set_display(display)
    coding.set_code(code)
    coding.set_system(system)

def read_criteria():
    inf = file('ValueSets/Criterion.json','r')
    jcrit = json.load(inf)
    inf.close()
    criteria = {}
    for rulenum in jcrit:
        crit = jcrit[rulenum]
        criterion = Criterion(crit['id'])
        criterion.set_description (crit['description'] )
        criterion.set_shortDescription ( crit['shortDescription'] )
        criterion.set_defaultStrength( crit['defaultStrength'] )
        criteria[rulenum] = criterion
    return criteria

#What's the point here?  We want variant to be a subclass of Node so that our
#serializer can know whether to rewrite it or not.  We could and perhaps
#should add all the methods to build the variant correctly, but not right now.  
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



#This declutters a bit by only printing the full node the first time we
#encounter it.  But that's not necessarily in the highest node, because of the
#ordering of the keys.  Should probably be smarter, but it works for the moment
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
                    #try to write it just by id
                    return obj.get_id()
                except:
                    #but if it doesn't have an id (like a codeable concept, just put all the data)
                    return obj.data
                #try:
                #    return obj.data[DMWG_ID_KEY]
                #except:
                #    return obj.data[ALLELE_REGISTRY_ID_KEY]
        else:
            return json.JSONEncoder.default(self,obj)

