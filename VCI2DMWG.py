#!/usr/bin/env python

import sys
from collections import defaultdict
import requests
import hashlib
import csv
from clingen_interpretation.interpretation_generated import *
from clingen_interpretation.interpretation_extras import *
from clingen_interpretation.interpretation_constants import *
from clingen_interpretation.Allele import Variant
import argparse
import logging

IRI_BASE='https://vci.clinicalgenome.org'
VCI_ID_KEY = '@id'
VCI_TYPE_KEY = '@type'
VCI_CONTRIBUTION_KEY = 'submitted_by'
VCI_LAST_MODIFIED_KEY = 'last_modified'
VCI_AGENT_NAME_KEY = 'title'
VCI_AGENT_AFFILIATION_KEY = 'affiliation'
VCI_AUTOCLASSIFICATION_KEY = 'autoClassification'
VCI_ALTEREDCLASSIFICATION_KEY = 'alteredClassification'
VCI_EVIDENCE_SUMMARY_KEY = 'evidenceSummary'
VCI_APPROVAL_REVIEW_DATE_KEY = 'approvalReviewDate'
VCI_CLASSIFICATION_APPROVER_KEY = 'classificationApprover'
VCI_APPROVAL_DATE_KEY = 'approvalDate'
VCI_APPROVAL_SUBMITTER_KEY = 'approvalSubmitter'
VCI_AFFILIATION_KEY = 'affiliation'
VCI_VARIANT_KEY = 'variant'
VCI_VARIANT_TYPE = 'variant'
VCI_DISEASE_TYPE = 'disease'
VCI_AGENT_TYPE = 'user'
VCI_CANONICAL_ID_KEY = 'carId'
VCI_HGVS_NAMES_KEY = 'hgvsNames'
VCI_GENOMIC_HGVS_38_KEY = 'GRCh38'
VCI_CONDITION_KEY = 'disease'
VCI_DISEASE_ONTOLOGY_KEY = 'ontology'
VCI_DISEASE_TERM_KEY = 'term'
VCI_DISEASE_ID_KEY = 'diseaseId'
VCI_EVALUATION_KEY = 'evaluations'
VCI_EVALUATION_VARIANT_KEY = 'variant'
VCI_CRITERIA_KEY = 'criteria'
VCI_CRITERIA_STATUS_KEY = 'criteriaStatus'
VCI_CRITERIA_NOT_EVALUATED = 'not-evaluated'
VCI_MET = 'met'
VCI_NOT_MET = 'not-met'
VCI_EVALUATION_EXPLANATION_KEY = 'explanation'
VCI_EVIDENCE_DESCRIPTION_KEY = 'evidenceDescription'
VCI_CRITERIA_MODIFIER_KEY = 'criteriaModifier'
VCI_MODIFIER_KEY = 'modifier'
VCI_FREQUENCY_KEY = 'population'
VCI_POPULATION_DATA_KEY = 'populationData'
VCI_COMPUTATIONAL_KEY = 'computational'
VCI_COMPUTATIONAL_DATA_KEY = 'computationalData'
VCI_COMPUTATIONAL_ALLELE_KEY = 'variant'
VCI_CLINGEN_COMPUTATION_KEY = 'clingen'
VCI_OTHER_COMPUTATION_KEY = 'other_predictors'
VCI_CONSERVATION_KEY = 'conservation'
VCI_CONSERVATION_DATA_KEY = 'conservationData'
VCI_ESP_KEY = 'esp'
VCI_EXAC_KEY = 'exac'
VCI_1000_GENOMES_KEY = 'tGenomes'
VCI_FREQUENCY_ALLELE_KEY='_extra'
VCI_COMBINED_POP = '_tot'
VCI_EUROPEAN_AMERICAN_POP = 'ea'
VCI_AFRICAN_AMERICAN_POP = 'aa'
VCI_CHROMOSOME_KEY = 'chrom'
VCI_REF_ALLELE_KEY = 'ref'
VCI_ALT_ALLELE_KEY = 'alt'
VCI_HG19_START_KEY = 'hg19_start'
VCI_EXAC_START_KEY = 'pos'
VCI_ALLELE_COUNT_KEY = 'ac'
VCI_ALLELE_NUMBER_KEY = 'an'
VCI_ALLELE_FREQUENCY_KEY = 'af'
VCI_GENOTYPE_COUNT_KEY = 'gc'
VCI_EXAC_AFRICAN = 'afr'
VCI_EXAC_LATINO_POP = 'amr'
VCI_EXAC_EAST_ASIAN_POP = 'eas'
VCI_EXAC_FINNISH_POP = 'fin'
VCI_EXAC_NON_FINNISH_EURO_POP = 'nfe'
VCI_EXAC_SOUTH_ASIAN_POP = 'sas'
VCI_EXAC_OTHER_POP = 'oth'
VCI_HOMOZYGOUS_GENOTYPE_COUNT_KEY = 'hom'
VCI_1000_GENOMES_DBSNP_KEY = 'name'
VCI_1000_GENOMES_AFRICAN_POP = 'afr'
VCI_1000_GENOMES_LATINO_POP = 'amr'
VCI_1000_GENOMES_EAST_ASIAN_POP = 'eas'
VCI_1000_GENOMES_SOUTH_ASIAN_POP = 'sas'
VCI_1000_GENOMES_EURO_POP = 'eur'
VCI_1000_GENOMES_ESP_AA_POP = 'espaa'
VCI_1000_GENOMES_ESP_EA_POP = 'espea'
VCI_SCORE_KEY = 'score'
VCI_PREDICTION_KEY = 'prediction'
VCI_EXTRA_EVIDENCE_KEY = 'extra_evidence_list'
VCI_ARTICLES_KEY = 'articles'
VCI_PMID_KEY = 'pmid'
VCI_CATEGORY_KEY = 'category'
VCI_SUBCATEGORY_KEY = 'subcategory'
VCI_MODEINHERITANCE_KEY = 'modeInheritance'
VCI_MODEINHERITANCE_ADJECTIVE_KEY = 'modeInheritanceAdjective'
VCI_PROVISIONAL_VARIANT_KEY = 'provisional_variant'

VCI_MISSENSE_EFFECT_PREDICTOR = 'missense_predictor'
VCI_SPLICE_EFFECT_PREDICTOR = 'splice'

DMWG_ESP='ESP ascertainment method'
DMWG_1KGESP='1000 Genomes ascertainment method'
DMWG_1KG='1000 Genomes ascertainment method'
DMWG_EXAC='ExAC ascertainment method'

# should probably be in a separate file...
systems = {
    'MONDO': 'http://monarchinitiative.org/disease/',  # FIXME - the VCI puts something like MONDO_0009861, but the url should be like MONDO:0009861 (but html escaped...)
    'DOID': 'http://www.disease-ontology.org/term/', # FIXME- the VCI puts something like DOID_1321, but the url should be like DOID:1321 (but html escaped...)
    'OMIM': 'http://omim.org/entry/',
    'Orphanet': ' http://www.orpha.net' # FIXME: is there a stable IRI per term?
}


term_map = { VCI_MET: 'Met', \
             VCI_NOT_MET: 'Not Met', \
             VCI_MISSENSE_EFFECT_PREDICTOR: 'missense effect', \
             VCI_SPLICE_EFFECT_PREDICTOR: 'splicing prediction' \
             }

extra_evidence_map = { ('population','population'): ['BA1','PM2','BS1'],\
                       ('predictors','functional-conservation-splicing-predictors'): ['PP3','BP4','BP1','PP2'], \
                       ('predictors','other-variants-in-codon'): ['PM5','PS1'], \
                       ('predictors','null-variant-analysis'): ['PVS1'], \
                       ('predictors','molecular-consequence-silent-intron'): ['BP7'], \
                       ('predictors','molecular-consequence-inframe-indel'): ['BP3','PM4'], \
                       ('experimental','hotspot-functional-domain'): ['PM1'], \
                       ('experimental','experimental-studies'): ['BS3','PS3'], \
                       ('case-segregation','observed-in-healthy'): ['BS2'], \
                       ('case-segregation','case-control'): ['PS4'], \
                       ('case-segregation','segregation-data'): ['BS4','PP1'], \
                       #The following is a typo that occurs in the VCI data:
                       ('case-segregation','segreagtion-data'): ['BS4','PP1'], \
                       ('case-segregation','de-novo'): ['PM6','PS2'], \
                       ('case-segregation','allele-data'): ['BP2','PM3'], \
                       ('case-segregation','alternate-mechanism'): ['BP5'], \
                       ('case-segregation','specificity-of-phenotype'): ['PP4'], \
                       ('case-segregation','reputable-source'): ['BP6','PP5'] }


def read_affiliations():
    with file('Affiliation_id_name_lookup.js','r') as inf:
        affs = json.load(inf)
    return { a['affiliation_id']: a['affiliation_fullname'] for a in affs }
    
#Just going to read this thing at global scope.  It's bad practice, but 
# really this file should be read from a service somewhere anyway
affiliations_id_to_name = read_affiliations()

def get_chromosome_name(chromosome, version):
    v37chromosomes={'1':'NC_000001.10',\
                    '2':'NC_000002.11',\
                    '3':'NC_000003.11',\
                    '4':'NC_000004.11',\
                    '5':'NC_000005.9',\
                    '6':'NC_000006.11',\
                    '7':'NC_000007.13',\
                    '8':'NC_000008.10',\
                    '9':'NC_000009.11',\
                    '10':'NC_000010.10',\
                    '11':'NC_000011.9',\
                    '12':'NC_000012.11',\
                    '13':'NC_000013.10',\
                    '14':'NC_000014.8',\
                    '15':'NC_000015.9',\
                    '16':'NC_000016.9',\
                    '17':'NC_000017.10',\
                    '18':'NC_000018.9',\
                    '19':'NC_000019.9',\
                    '20':'NC_000020.10',\
                    '21':'NC_000021.8',\
                    '22':'NC_000022.10',\
                    'X':'NC_000023.10',\
                    'Y':'NC_000024.9',\
                    'M':'NC_012920.1'}
    if version in ( 'hg19', 'GRCh37' ):
        return v37chromosomes[ str(chromosome) ]
    else:
        raise Exception

def get_canonical_id(hgvs):
    # send a GET request with parameter
    url = 'http://reg.genome.network/allele?hgvs='
    # convert symbol > to special code %3E
    url += requests.utils.quote(hgvs)
    res = requests.get(url)
    txt = res.text
    cardata = json.loads(txt)
    return cardata

#We don't need to entity map everything, just some of the data nodes.
# For instance, not the interpretation.  When it occurs multiple times,
# it causes problems because the evaluations at each level
# are represented differently.
class EntityMap:
    EMtypes = set([VCI_VARIANT_TYPE, VCI_AGENT_TYPE, VCI_DISEASE_TYPE])
    def __init__(self, source, idtag='@id'):
        self.entities = defaultdict(dict)
        self.transformed={}
        self.idtag = idtag
        self.walk(source)
    def walk(self,source):
        if isinstance(source,list):
            for e in source:
                self.walk(e)
        elif isinstance(source,dict):
            self.register(source)
            for k,v in source.items():
                self.walk(v)
        else:
            pass
    def register(self,node):
        if self.idtag in node:
            if len(self.EMtypes.intersection(set(node[VCI_TYPE_KEY]))) != 0:
                atid = fully_qualify(node[self.idtag])
                entity = self.entities[atid]
                for key in node:
                    if key in entity:
                        if entity[key] != node[key]:
                            print key, '\n--------\n',entity[key], '\n---------------\n',node[key]
                            raise Exception('Incoherent Nodes')
                    else:
                        entity[key] = node[key]
                #print atid,entity
    def get_entity(self,eid):
        return self.entities[eid]
    def get_transformed(self,eid):
        if eid in self.transformed:
            return self.transformed[eid]
        return None
    def add_transformed(self,eid,entity):
        self.transformed[eid] = entity

def canonicalizeVariant(rep):
    orig_carid = rep[VCI_CANONICAL_ID_KEY]
    #We want to get the id/representation from the Baylor Allele Registry.
    hgvs38 = rep[VCI_HGVS_NAMES_KEY][VCI_GENOMIC_HGVS_38_KEY]
    baylor_car_rep = get_canonical_id(hgvs38)
    baylor_carid = baylor_car_rep['@id']
    #Make sure it's the same id
    if (orig_carid) != '':
        #If the original is de-curied, then we need to compare to that...
        if not orig_carid.startswith('http://'):
            compare_id = baylor_carid.split('/')[-1]
        else:
            compare_id = baylor_carid
        if orig_carid != compare_id:
            logging.error('Original ID: %s.    Final ID: %s' % (orig_carid, baylor_carid) )
            raise Exception
    #compact CAR id:
    baylor_car_rep['@id'] = 'CAR:{}'.format(baylor_carid.split('/')[-1])
    return baylor_car_rep

def fully_qualify(iri):
    fqiri = iri
    if fqiri.startswith('/'):
        fqiri = IRI_BASE + fqiri
    return fqiri

def get_id(source):
    if isinstance(source, dict):
        sid = fully_qualify(source[VCI_ID_KEY])
    else:
        sid = fully_qualify(source)
    return sid


def add_contribution( user_input, target, ondate, entities, role ):
    userid = get_id(user_input)
    user = entities.get_entity(userid)
    #We're going to make some assumptions: 
    # 1. if we want an affiliation associated with a particular contribution
    #    then the affiliation tag will be local. i.e. we are not going to grab
    #    the affiliation from some other place that the user occurs
    # 2. The user_input here might be a dict, could be a string
    # 3. affiliation can sometimes be a list and sometimes be a single value :(
    # 4. The agent_for is not formally added to the DMWG model, so that there
    #    is not a set_affiliations method on contributions already.  We
    #    will probably need to modify this down the road.
    # 5. That information about the affiliation will be found in a local data 
    #    file (which needs to be made accessible in another way) 
    affiliations = None
    if isinstance(user_input,dict) and VCI_AGENT_AFFILIATION_KEY in user_input:
        affiliations = user_input[VCI_AGENT_AFFILIATION_KEY]
        if isinstance(affiliations,str):
            affiliations = [affiliations]
        affiliations = [ affiliations_id_to_name[a] for a in affiliations ]
    #print 'User:',user
    if user == {}:
        agent = Agent()
        agent.set_label( user_input )
    else:
        agent = entities.get_transformed(userid)
        if agent is None:
            agent = Agent(userid)
            try:
                agent.set_label( user[VCI_AGENT_NAME_KEY] )
            except KeyError:
                #sometimes we might not have the name.. oh well.
                pass
            entities.add_transformed(userid, agent)
    contribution = create_contribution(agent, ondate, role)
    if affiliations is not None:
        contribution.data['agent_for'] = affiliations
    target.add_contribution(contribution)

def add_contributions( source, target, entities, ondate, role ):
    if isinstance( source, list ):
        for single_source in source:
            add_contribution( single_source, target, ondate, entities, role)
    else:
        add_contribution( source, target,  ondate, entities, role)

#Ignore these keys:
# interpretation_genes: unused in the VCI
# actions: tracking user behavior
# markAsProvisional: UI/process related
# submitted_by: this is who started the interpretion.  we care who ended it. (in prov_variant)
# uuid: redundant with @id
# schema_version: not relevant to transformed interp
# evaluation_count: we can figure out from the list of evaluations
# status: not tracking status
# interpretation_status
# last_modified: using the information from the provisional_variant
# audit: Interal information
# interpretation disease: redundant with disease
# date created: Not tracking user behavior
# provisional_count: should be 1:1
#Still need to handle:
# modeInheritance:
# modeInheritanceAdjective: for things like "with maternal imprinting"
def transform_root(vci):
    vi = VariantPathogenicityInterpretation( fully_qualify(vci[VCI_ID_KEY]) )
    modestring = vci[VCI_MODEINHERITANCE_KEY]
    modeadjstring = vci[VCI_MODEINHERITANCE_ADJECTIVE_KEY]
    if modeadjstring != '':
        raise Exception('Mode Inheritance Adjective not handled')
    if '(HP:' in modestring:
        p = modestring.split()
        for pi in p:
            if pi.startswith('(HP:'):
                mode = pi[1:-1]
    else:
        mode = modestring
    return vi, mode

#Ignore:
# status: not tracking status
# uuid: redudant with id
# @id: we don't need an id for this entity which we are merging with the parent
#  Confirmed with Karen on 4/19/2017 that the id for tracking submission is the interpretation id.
# interpretation_associated: just a reverse link to the parent
# schema_version: not relevant to the transformed interpretation
# date_create: just tracking end dates
#May need to handle
# alteredClassification_present
# reason_present
def transform_provisional_variant(vci_pv , interpretation, entities ):
    if isinstance(vci_pv, list):
        if len(vci_pv) > 1:
            print '???'
            exit()
        vci_pv = vci_pv[0]
    interpreter = ''
    if VCI_PROVISIONAL_VARIANT_KEY in vci_pv:
        interpreter = vci_pv[VCI_CLASSIFICATION_APPROVER_KEY]
    elif VCI_APPROVAL_SUBMITTER_KEY in vci_pv:
        interpreter = vci_pv[VCI_APPROVAL_SUBMITTER_KEY]
    else:
        interpreter = vci_pv[VCI_CONTRIBUTION_KEY]

    interpreted_date = ''
    if VCI_APPROVAL_REVIEW_DATE_KEY in vci_pv:
        interpreted_date = vci_pv[VCI_APPROVAL_REVIEW_DATE_KEY]
    elif VCI_APPROVAL_DATE_KEY in vci_pv:
        interpreted_date = vci_pv[VCI_APPROVAL_DATE_KEY]
    else:
        interpreted_date = vci_pv[VCI_LAST_MODIFIED_KEY]

    add_contributions( interpreter, interpretation, entities, interpreted_date, DMWG_INTERPRETER_ROLE)

    if VCI_EVIDENCE_SUMMARY_KEY in vci_pv:
        interpretation.set_description( vci_pv[VCI_EVIDENCE_SUMMARY_KEY])

    interpretation.set_statementOutcome( convert_significance(vci_pv) )

def convert_significance(vci_provisional_variant):
    significance = ''
    if VCI_ALTEREDCLASSIFICATION_KEY in vci_provisional_variant:
        significance = vci_provisional_variant[VCI_ALTEREDCLASSIFICATION_KEY]
    else:
        significance = vci_provisional_variant[VCI_AUTOCLASSIFICATION_KEY]
    if significance == 'Uncertain significance - conflicting evidence':
        significance = 'LOINC:LA26333-7'
    if significance == 'Uncertain significance - insufficient evidence':
        significance = 'LOINC:LA26333-7'
    return significance

def transform_variant(variant,entities):
    vci_variant_id = get_id(variant)
    dmwg_variant = entities.get_transformed(vci_variant_id)
    if dmwg_variant is None:
        vci_variant = entities.get_entity(vci_variant_id)
        ar_variant = canonicalizeVariant( vci_variant )
        dmwg_variant = Variant(ar_variant)
        entities.add_transformed(vci_variant_id, dmwg_variant)
    return dmwg_variant

#Evaluation keys:
## Ones we won't use
# date_created
# schema_version
# interpretation_associated
# uuid
# @type
# status
# evidence_type: going to be clear from the data
## Go into contribution
# last_modified X
# submitted_by X
## Identifier X
# @id
## transform
# criteria X
# criteriaStatus X
# explanation X
# modifier: related to criteriaModifier for the UI
# variant X
# population: this is actually allele frequency data (in populations, not a population itself)
def transform_evaluation(vci_evaluation, interpretation, entities, criteria):
    dmwg_assessment = CriterionAssessment( vci_evaluation[VCI_ID_KEY] )
    criterion = criteria[ vci_evaluation[ VCI_CRITERIA_KEY] ]
    dmwg_assessment.set_criterion( criterion )
    dmwg_assessment.set_statementOutcome( term_map[ vci_evaluation[ VCI_CRITERIA_STATUS_KEY ] ] )
    dmwg_assessment.set_description( vci_evaluation[ VCI_EVALUATION_EXPLANATION_KEY] )
    dmwg_assessment.set_variant( transform_variant( vci_evaluation[VCI_EVALUATION_VARIANT_KEY ] , entities) )
    #Have to do a little work to figure out the strength
    #VCI has both a modifier and a criteria modifier.  If these both exist they should be the
    #same, but if one or the other does not exist, then we should use the other one.
    crit_mod = ''
    mod = ''
    if VCI_CRITERIA_MODIFIER_KEY in vci_evaluation:
        crit_mod = vci_evaluation[VCI_CRITERIA_MODIFIER_KEY]
    if VCI_MODIFIER_KEY in vci_evaluation:
        mod = vci_evaluation[VCI_MODIFIER_KEY]
    if crit_mod != mod:
        raise Exception
    defaultStrength = criterion.get_defaultStrength()
    strength = transform_strength( crit_mod, defaultStrength )
    ##TODO: Also add contribution to the evidence line if crit_mod != ''.  A little tricky since I hid the evidence line, but it's gettable
    add_contributions( vci_evaluation[VCI_CONTRIBUTION_KEY], dmwg_assessment, entities,vci_evaluation[VCI_LAST_MODIFIED_KEY], DMWG_ASSESSOR_ROLE)
    #Now the evidence
    if VCI_FREQUENCY_KEY in vci_evaluation:
        frequencies = transform_frequency( vci_evaluation[VCI_FREQUENCY_KEY],  entities)
        add_evidenceItems( dmwg_assessment, frequencies )
    if VCI_COMPUTATIONAL_KEY in vci_evaluation:
        predictions = transform_computational( vci_evaluation[VCI_COMPUTATIONAL_KEY], entities )
        add_evidenceItems( dmwg_assessment, predictions )
    add_criterion_assessment(interpretation, dmwg_assessment, strength)
    return vci_evaluation[VCI_CRITERIA_KEY], dmwg_assessment

def transform_computational(source, entities):
    predictions = []
    vci_variant = source[VCI_VARIANT_KEY]
    #the form of this call suggests that the transform* functions should be member functions of entitymap
    dmwg_variant = transform_variant(vci_variant,entities)
    compdata = source[VCI_COMPUTATIONAL_DATA_KEY]
    if VCI_CONSERVATION_KEY in compdata:
        predictions += transform_conservation_data( compdata[VCI_CONSERVATION_KEY] , dmwg_variant )
    if VCI_CLINGEN_COMPUTATION_KEY in compdata:
        predictions += transform_clingen_comp_data( compdata[VCI_CLINGEN_COMPUTATION_KEY], dmwg_variant )
    if VCI_OTHER_COMPUTATION_KEY in compdata:
        predictions += transform_other_comp_data( compdata[VCI_OTHER_COMPUTATION_KEY], dmwg_variant )
    add_contributions_to_data( source , predictions, entities )
    return predictions

#The VCI insilico predictions have a score attribute as well as a prediction.
# It looks like the prediction is always generic text that "higher score = higher pathogenicity"
# So we will just keep thte score piece.
#In the event that this changes, and the prediction becomes a qualitative prediction,
# we'll make a prediction (with the qualitative score)->evidence line -> ISPScore(score)
#See transform_other_comp_data for an example
def transform_clingen_comp_data( source,variant):
    predictions = []
    for pred in source:
        score = source[pred][VCI_SCORE_KEY]
        if score is not None:
            if source[pred]['prediction'] != "higher score = higher pathogenicity":
                print '!',source[pred]['prediction']
            #Have to sort out the scores....
            dmwg_prediction = InSilicoPredictionScoreStatement()
            #We really also want to set the transcript, but the VCI is not returning that informaiton
            #dmwg_prediction.set_transcript()
            dmwg_prediction.set_canonicalAllele(variant)
            dmwg_prediction.set_algorithm(pred)
            dmwg_prediction.set_prediction(score)
            dmwg_prediction.set_predictionType( term_map[VCI_MISSENSE_EFFECT_PREDICTOR] )
            predictions.append(dmwg_prediction)
    return predictions

def transform_other_comp_data( source, variant ):
    predictions = []
    for pred in source:
        scores = source[pred][ VCI_SCORE_KEY ]
        predsv  = source[pred][ VCI_PREDICTION_KEY ]
        if (scores is not None) and (not isinstance(scores, list)):
            scores = [ scores ]
        if scores is not None and predsv is not None:
            preds  = predsv.split(',')
        elif scores is None and predsv is not None:
            preds = predsv.split(',')
            scores = [ None for p in preds ]
        elif scores is not None and predsv is None:
            preds = [None for s in scores]
        else:
            #both none, ignore
            continue
        for (s,p) in zip(scores,preds):
            dmwg_prediction = dmwg_prediction_score = None
            if s is not None:
                dmwg_prediction_score = InSilicoPredictionScoreStatement()
                dmwg_prediction_score.set_predictionType( term_map[ VCI_MISSENSE_EFFECT_PREDICTOR]  )
                dmwg_prediction_score.set_canonicalAllele( variant )
                dmwg_prediction_score.set_algorithm( pred )
                dmwg_prediction_score.set_prediction(s)
            if p is not None:
                dmwg_prediction = InSilicoPredictionStatement()
                dmwg_prediction.set_predictionType( term_map[ VCI_MISSENSE_EFFECT_PREDICTOR]  )
                dmwg_prediction.set_canonicalAllele( variant )
                dmwg_prediction.set_algorithm( pred )
                dmwg_prediction.set_statementOutcome(p)
            if dmwg_prediction is None:
                predictions.append(dmwg_prediction_score)
            else:
                if dmwg_prediction_score is not None:
                    evidence_line = EvidenceLine()
                    evidence_line.add_evidenceItem( dmwg_prediction_score )
                    dmwg_prediction.add_evidenceLine( evidence_line )
                predictions.append(dmwg_prediction)
    return predictions


#The values from VCI are not including a true/false on whether the thing is conserved. But we'll want that.
#When we get that, this will have to change.  What it will become is an AlleleConservation object
# with an evidence line to an AlleleConservationScore, which will be defined as below
def transform_conservation_data( source, variant ):
    results = []
    for constool in source:
        dmwg_conservation = AlleleConservationScoreStatement()
        dmwg_conservation.set_allele(variant)
        dmwg_conservation.set_algorithm(constool)
        dmwg_conservation.set_score(source[constool])
        results.append(dmwg_conservation)
    return results

def transform_frequency( source, entities):
    popdata    = source[VCI_POPULATION_DATA_KEY]
    vci_variant = source[VCI_VARIANT_KEY]
    #the form of this call suggests that the transform* functions should be member functions of entitymap
    dmwg_variant = transform_variant(vci_variant,entities)
    frequencies = []
    if VCI_ESP_KEY in popdata:
        esp_frequencies = transform_esp_data( popdata[VCI_ESP_KEY],  dmwg_variant )
        frequencies += esp_frequencies
    if VCI_EXAC_KEY in popdata:
        exac_frequencies = transform_exac_data( popdata[VCI_EXAC_KEY], dmwg_variant )
        frequencies += exac_frequencies
    if VCI_1000_GENOMES_KEY in popdata:
        tg_frequencies = transform_1000_genomes_data(popdata[VCI_1000_GENOMES_KEY], dmwg_variant)
        frequencies += tg_frequencies
    #Add contributors to data nodes
    add_contributions_to_data( source , frequencies, entities )
    return frequencies

def add_contributions_to_data( data_source, data_targets, entities):
    submitters = data_source[VCI_CONTRIBUTION_KEY]
    modtime    = data_source[VCI_LAST_MODIFIED_KEY]
    for data in data_targets:
        add_contributions( submitters, data, entities, modtime, DMWG_CURATOR_ROLE)


#Clean up the VCI keys and decide if they are per experiment or not.
def convert_esp_pop(pop):
    if pop == VCI_COMBINED_POP:
        return 'combined'
    return 'EVS:%s' % pop.upper()

def convert_exac_pop(pop):
    if pop == VCI_COMBINED_POP: return 'combined'
    if pop == VCI_EXAC_OTHER_POP: return 'other'
    return 'GNOMAD:%s' % pop

def convert_1000_genomes_pop(pop):
    if pop == VCI_COMBINED_POP: return 'combined'
    if pop in [ VCI_1000_GENOMES_ESP_AA_POP, VCI_1000_GENOMES_ESP_EA_POP ]:
        return 'EVS:%s' % pop[-2:].upper()
    return 'IGSR:%s' % pop

def transform_1000_genomes_data( source, dmwg_variant ):
    alt_allele = dmwg_variant.get_allele('GRCh37') #right?
    frequencies = []
    for pop in source:
        if pop != VCI_FREQUENCY_ALLELE_KEY:
            if pop in [VCI_1000_GENOMES_ESP_EA_POP, VCI_1000_GENOMES_ESP_AA_POP]:
                dmwg_af = PopulationAlleleFrequencyStatement()
                dmwg_af.set_ascertainment(DMWG_1KGESP)
            else:
                dmwg_af = PopulationAlleleFrequencyStatement()
                dmwg_af.set_ascertainment(DMWG_1KG)
            frequencies.append(dmwg_af)
            dmwg_af.set_allele( dmwg_variant )
            dmwg_af.set_population(  convert_1000_genomes_pop(pop) )
            if len( source[pop][VCI_ALLELE_COUNT_KEY] ) == 0:
                dmwg_af.set_alleleNumber( 0 )
            else:
                dmwg_af.set_alleleNumber( sum( source[pop][VCI_ALLELE_COUNT_KEY].values() ) )
            #TODO: Should the element be not here, or should it be 0 or should it be null?
            if dmwg_af.get_alleleNumber() > 0:
                if alt_allele in source[pop][VCI_ALLELE_FREQUENCY_KEY]:
                    dmwg_af.set_alleleFrequency( source[pop][VCI_ALLELE_FREQUENCY_KEY][alt_allele] )
                else:
                    dmwg_af.set_alleleFrequency( 0.  )
            if alt_allele in source[pop][VCI_ALLELE_COUNT_KEY]:
                dmwg_af.set_alleleCount( source[pop][VCI_ALLELE_COUNT_KEY][alt_allele] )
            else:
                dmwg_af.set_alleleCount( 0 )
            hkey = '%s|%s' % (alt_allele,alt_allele)
            if hkey in source[pop][VCI_GENOTYPE_COUNT_KEY]:
                dmwg_af.set_homozygousAlleleIndividualCount(  source[pop][VCI_GENOTYPE_COUNT_KEY][hkey] )
            else:
                dmwg_af.set_homozygousAlleleIndividualCount( 0 )
    return frequencies


def transform_exac_data( source, dmwg_variant ):
    frequencies = []
    for pop in source:
        if pop != VCI_FREQUENCY_ALLELE_KEY:
            dmwg_af = PopulationAlleleFrequencyStatement()
            dmwg_af.set_ascertainment(DMWG_EXAC)
            frequencies.append(dmwg_af)
            dmwg_af.set_population( convert_exac_pop(pop) )
            dmwg_af.set_allele( dmwg_variant )
            if VCI_ALLELE_COUNT_KEY in source[pop]:
                dmwg_af.set_alleleCount(  source[pop][VCI_ALLELE_COUNT_KEY] )
            else:
                dmwg_af.set_alleleCount( 0 )
            if VCI_ALLELE_NUMBER_KEY in source[pop]:
                dmwg_af.set_alleleNumber(source[pop][VCI_ALLELE_NUMBER_KEY])
            else:
                dmwg_af.set_alleleNumber(0)
            if dmwg_af.get_alleleNumber() > 0:
                dmwg_af.set_alleleFrequency(source[pop][VCI_ALLELE_FREQUENCY_KEY])
            if VCI_HOMOZYGOUS_GENOTYPE_COUNT_KEY in source[pop]:
                dmwg_af.set_homozygousAlleleIndividualCount( source[pop][VCI_HOMOZYGOUS_GENOTYPE_COUNT_KEY] )
            else:
                dmwg_af.set_homozygousAlleleIndividualCount( 0 )
    return frequencies


def transform_esp_data(source,dmwg_variant):
    ref_allele = dmwg_variant.get_ref_allele('GRCh37')
    alt_allele = dmwg_variant.get_allele('GRCh37')
    frequencies = []
    for pop in source:
        if pop != VCI_FREQUENCY_ALLELE_KEY:
            dmwg_af = PopulationAlleleFrequencyStatement( )
            dmwg_af.set_ascertainment(DMWG_ESP)
            frequencies.append(dmwg_af)
            dmwg_af.set_population( convert_esp_pop(pop) )
            dmwg_af.set_allele( dmwg_variant )
            #Assumes ESP uses v37
            if (len(source[pop][VCI_ALLELE_COUNT_KEY]) > 0) and (alt_allele in source[pop][VCI_ALLELE_COUNT_KEY]):
                dmwg_af.set_alleleCount( source[pop][VCI_ALLELE_COUNT_KEY][alt_allele] )
            else:
                dmwg_af.set_alleleCount( 0 )
            dmwg_af.set_alleleNumber( sum(source[pop][VCI_ALLELE_COUNT_KEY].values()) )
            if dmwg_af.get_alleleNumber() > 0:
                f = 1.*dmwg_af.get_alleleCount() / dmwg_af.get_alleleNumber()
                dmwg_af.set_alleleFrequency(f)
            if len(source[pop][VCI_GENOTYPE_COUNT_KEY]) > 0:
                homkey = '%s%s' % (alt_allele,alt_allele)
                dmwg_af.set_homozygousAlleleIndividualCount(source[pop][VCI_GENOTYPE_COUNT_KEY][homkey])
                hetkey1 = '%s%s' % (ref_allele,alt_allele)
                hetkey2 = '%s%s' % (alt_allele,ref_allele)
                if hetkey1 in  source[pop][VCI_GENOTYPE_COUNT_KEY]:
                    dmwg_af.set_heterozygousAlleleIndividualCount(source[pop][VCI_GENOTYPE_COUNT_KEY][hetkey1])
                elif hetkey2 in  source[pop][VCI_GENOTYPE_COUNT_KEY]:
                    dmwg_af.set_heterozygousAlleleIndividualCount( source[pop][VCI_GENOTYPE_COUNT_KEY][hetkey2] )
    f = frequencies[0]
    return frequencies


def transform_strength(modifier, defaultStrength):
    """defaultStrength is one of our DomainEntities.  If we need to change the stength,
    we are going to get the display, and change it, then pass back that string, which will
    be turned back into an object in the back end"""
    if modifier == '':
        return defaultStrength
    path = defaultStrength.get_label().split(' ')[0]
    if modifier == 'strong':
        mod = 'Strong'
    elif modifier == 'supporting':
        mod = 'Supporting'
    elif modifier == 'moderate':
        mod = 'Moderate'
    elif modifier == 'very-strong':
        mod = 'Very Strong'
    else:
        print modifier
        exit()
    strength = ' '.join( [path, mod] )
    return strength

#Returns a dictionary that maps from criterion id (e.g. "PS2") to
# a dmwg assessment entity
def transform_evaluations(evaluation_list,interpretation,entities):
    evaluation_map = {}
    criteria = read_criteria()
    if not isinstance(evaluation_list, list):
        raise Exception
    for vci_eval in evaluation_list:
        if vci_eval[ VCI_CRITERIA_STATUS_KEY ] == VCI_CRITERIA_NOT_EVALUATED:
                continue
        crit_id, evaluation = transform_evaluation(vci_eval,interpretation,entities,criteria)
        evaluation_map[crit_id] = evaluation
    return evaluation_map

def transform_articles( article_list, interpretation, entities ):
    sourcelist = []
    for article in article_list:
        pmid = article[ VCI_PMID_KEY ]
        atid = 'https://www.ncbi.nlm.nih.gov/pubmed/%s' % pmid
        #isource = InformationSource(atid)
        sourcelist.append(atid)
    return sourcelist

# The evidence in the VCI is not a child node of the evaluation.  It is linked to the relevant
# rules via the category and subcategory proerties
def transform_evidence(extra_evidence_list, interpretation, entities, evalmap):
    for ee_node in extra_evidence_list:
        info = Statement()
        add_contributions_to_data( ee_node, [info], entities )
        info.set_description( ee_node[ VCI_EVIDENCE_DESCRIPTION_KEY] )
        sources = transform_articles(ee_node[ VCI_ARTICLES_KEY], interpretation, entities )
        for source in sources:
            info.add_source(source)
        key = ( ee_node[VCI_CATEGORY_KEY], ee_node[ VCI_SUBCATEGORY_KEY] )
        possible_rules = extra_evidence_map[key]
        found = False
        for rule in possible_rules:
            if rule in evalmap:
                found = True
                dmwg_assessment = evalmap[rule]
                add_evidenceItems( dmwg_assessment, [info])
        if not found:
            print  "Did not find any evaluated criteria for this data: %s "% ee_node['uuid']


#Note that we don't necessarily have the full VCI disease node when we get into this function.
# It could be anything from a bare IRI to a full node or anything in between.  The first two lines
# standardize the "local" information to the "global" information node.
def transform_condition(vci_local_disease,interpretation,entities,mode):
    vci_disease_id = get_id(vci_local_disease)
    vci_disease = entities.get_entity(vci_disease_id)
    dmwg_disease = entities.get_transformed(vci_disease_id)
    if dmwg_disease is None:
        disease_ontology = 'MONDO:'
        #disease_ontology = systems.get(disease_ontology, disease_ontology)
        #print "VCI_DISEASE_ONTOLOGY_KEY="+VCI_DISEASE_ONTOLOGY_KEY
        #disease_ontology = vci_disease[VCI_DISEASE_ONTOLOGY_KEY]
        #print disease_ontology
        disease_code = vci_disease[VCI_DISEASE_ID_KEY]
        if not disease_code.startswith('MONDO'):
            raise Exception("Expected a MONDO disease identifier")
        disease_code = disease_code.split('_')[-1]
        disease_name = vci_disease[VCI_DISEASE_TERM_KEY]
        dmwg_disease = create_dmwg_disease(disease_ontology, disease_code, disease_name)
        entities.add_transformed(vci_disease_id, dmwg_disease)
    dmwg_condition = GeneticCondition()
    dmwg_condition.add_disease(dmwg_disease)
    if mode!= '': dmwg_condition.set_inheritancePattern( mode )
    interpretation.add_condition(dmwg_condition)

def transform(jsonf):
    vci = json.load(jsonf)
    if VCI_MODEINHERITANCE_KEY not in vci:
        vci[VCI_MODEINHERITANCE_KEY] = ''
    entities = EntityMap(vci)
    interpretation, inheritance = transform_root(vci)
    if VCI_PROVISIONAL_VARIANT_KEY in vci:
        transform_provisional_variant(vci[VCI_PROVISIONAL_VARIANT_KEY],interpretation,entities)
    else:
        logging.warning('No provisional_variant element in the input.  Proceeding, but clinical significance will not be set.')
    try:
        variant = transform_variant(vci[VCI_VARIANT_KEY],entities)
        interpretation.set_variant(variant)
    except KeyError:
        logging.warning('No variant found. Proceeding.')
    try:
        transform_condition(vci[VCI_CONDITION_KEY], interpretation, entities, inheritance)
    except KeyError:
        logging.warning('No condition found. Proceeding')
    try:
        eval_map = transform_evaluations(vci[VCI_EVALUATION_KEY], interpretation, entities)
    except KeyError:
        logging.warning('No criteria evaluations found.  Proceeding')
        eval_map = {}
    try:
        transform_evidence(vci[VCI_EXTRA_EVIDENCE_KEY], interpretation, entities, eval_map)
    except KeyError:
        logging.warning('No evidence found.  Proceeding')
    return interpretation,entities

def transform_json_file(inf, outf, out_style):
    interp,ents = transform(inf)
    #The idea of flatten would be that the root node would contain fully specified descriptions of all the entities, and then
    # the interpretation, which would be written using only IDs.
    # To implement this, just create a new dict, put the interpretation into it, and write the entites into it from the EntityMap
    # Then dump that envelope, and modify the interpretaion encoder to have a mode that keeps track of depth and writes full nodes
    # for the top level.
    # The other option is not to include this at all, and rely on JSON-LD libraries to do any appropriate flattening.
    # TODO: Decide and implement (if required) or remove option (if not)
    if out_style == 'flat':
        raise Exception('flatten not implemented yet')
    json.dump(interp,outf,sort_keys=True, indent=4, separators=(',', ': '), cls=InterpretationEncoder, out_style=out_style)

def test():
    transform_json_file(open('test_data/test_interp_1.vci.json'), open('test_data/test_interp_1.dmwg.json'), 'first')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input',  nargs='?', type=argparse.FileType('r'), default=sys.stdin,
            help='Path to an input JSON file created by the VCI (defaults to stdin)')
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
            help='Output path for DMWG JSON file to be created (defaults to stdout)')
    parser.add_argument("-s", "--output-style", type=str, choices=['full', 'first', 'flat'], help="full: expand all nodes, first: expand first node, flat: define entities outside the interpretation", default = 'first')
    args = parser.parse_args()
    transform_json_file(args.input, args.output, args.output_style)
