#!/usr/bin/env python3

import sys
from collections import defaultdict
import hashlib
import csv
from clingen_interpretation.interpretation_generated import *
from clingen_interpretation.interpretation_extras import *
from clingen_interpretation.interpretation_constants import *
from clingen_interpretation.Allele import *
import argparse
import logging
import re
import datetime
from dateutil.parser import parse

IRI_BASE = 'https://vci.clinicalgenome.org'
VCI_PK_KEY = 'PK'
VCI_ITEM_TYPE_KEY = 'item_type'
VCI_CLINVAR_VARIANT_TITLE = 'clinvarVariantTitle'
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
# VCI_VARIANT_SOURCE = 'source'
VCI_DISEASE_TYPE = 'disease'
VCI_AGENT_TYPE = 'user'
VCI_CLINVAR_ID_KEY = 'clinvarVariantId'
VCI_CANONICAL_ID_KEY = 'carId'
VCI_HGVS_NAMES_KEY = 'hgvsNames'
VCI_GENOMIC_HGVS_38_KEY = 'GRCh38'
VCI_OTHERS_HGVS_KEY = 'others'
VCI_CONDITION_KEY = 'disease'
VCI_DISEASE_ONTOLOGY_KEY = 'ontology'
VCI_DISEASE_TERM_KEY = 'term'
# VCI_DISEASE_ID_KEY = 'diseaseId'
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
# VCI_MODIFIER_KEY = 'modifier'
VCI_FREQUENCY_KEY = 'population'
VCI_POPULATION_DATA_KEY = 'populationData'
VCI_COMPUTATIONAL_KEY = 'computational'
VCI_COMPUTATIONAL_DATA_KEY = 'computationalData'
VCI_CLINGEN_COMPUTATION_KEY = 'clingen'
VCI_OTHER_COMPUTATION_KEY = 'other_predictors'
VCI_CONSERVATION_KEY = 'conservation'
VCI_ESP_KEY = 'esp'
VCI_EXAC_KEY = 'exac'
VCI_1000_GENOMES_KEY = 'tGenomes'
VCI_COMBINED_POP = '_tot'
VCI_ALLELE_COUNT_KEY = 'ac'
VCI_ALLELE_NUMBER_KEY = 'an'
VCI_ALLELE_FREQUENCY_KEY = 'af'
VCI_GENOTYPE_COUNT_KEY = 'gc'
VCI_EXAC_OTHER_POP = 'oth'
VCI_HOMOZYGOUS_GENOTYPE_COUNT_KEY = 'hom'
VCI_1000_GENOMES_ESP_AA_POP = 'espaa'
VCI_1000_GENOMES_ESP_EA_POP = 'espea'
VCI_SCORE_KEY = 'score'
VCI_PREDICTION_KEY = 'prediction'
VCI_V1_EXTRA_EVIDENCE_KEY = 'extra_evidence_list'
VCI_CURATED_EVIDENCE_KEY = 'curated_evidence_list'
VCI_ARTICLES_KEY = 'articles'
VCI_PMID_KEY = 'pmid'
VCI_CATEGORY_KEY = 'category'
VCI_SUBCATEGORY_KEY = 'subcategory'
VCI_MODEINHERITANCE_KEY = 'modeInheritance'
VCI_MODEINHERITANCE_ADJECTIVE_KEY = 'modeInheritanceAdjective'
VCI_PROVISIONAL_VARIANT_KEY = 'provisionalVariant'
VCI_V1_PROVISIONAL_VARIANT_KEY = 'provisional_variant'

VCI_MISSENSE_EFFECT_PREDICTOR = 'missense_predictor'
VCI_SPLICE_EFFECT_PREDICTOR = 'splice'

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

#This maps categories and sub-categories to rules.
# The value is a list of tuples.  Each tuple contains a set of criteria only one of which can
# actually be true.
extra_evidence_map = { ('population','population'): [('BA1','PM2','BS1')],\
                       ('predictors','functional-conservation-splicing-predictors'): [('PP3','BP4'),('BP1','PP2')], \
                       ('predictors','other-variants-in-codon'): [('PM5',),('PS1',)], \
                       ('predictors','null-variant-analysis'): [('PVS1',)], \
                       ('predictors','molecular-consequence-silent-intron'): [('BP7',)], \
                       ('predictors','molecular-consequence-inframe-indel'): [('BP3','PM4')], \
                       #The following is a typo that occurs in the VCI data:
                       ('experimental','hotspot-functiona-domain'): [('PM1',)], \
                       ('experimental','experimental-studies'): [('BS3','PS3')], \
                       ('case-segregation','observed-in-healthy'): [('BS2',)], \
                       ('case-segregation','case-control'): [('PS4',)], \
                       ('case-segregation','segregation-data'): [('BS4',),('PP1',)], \
                       #The following is a typo that occurs in the VCI data:
                       ('case-segregation','segreagtion-data'): [('BS4',),('PP1',)], \
                       ('case-segregation','de-novo'): [('PM6',),('PS2',)], \
                       ('case-segregation','allele-data'): [('BP2',),('PM3',)], \
                       ('case-segregation','alternate-mechanism'): [('BP5',)], \
                       ('case-segregation','specificity-of-phenotype'): [('PP4',)], \
                       ('case-segregation','reputable-source'): [('BP6','PP5')] }

AFFILIATION_IRI_NAMESPACE = 'http://curation.clinicalgenome.org/affiliation/'
def read_affiliations():
    affiliation_file = os.path.join(os.path.dirname(__file__),"Affiliation_id_name_lookup.js")
    with open(affiliation_file,'r') as inf:
        affs = json.load(inf)

    adict = {}
    adict[None] = None
    for a in affs:
        affiliation_public_id = ''
        guideline = {}
        if 'affiliation_vcep_id' in a:
            affiliation_public_id = a['affiliation_vcep_id']
        elif 'affiliation_id' in a:
            affiliation_public_id = a['affiliation_id']
        if 'guideline_name' in a:
            guideline['name'] = a['guideline_name']
        if 'guideline_url' in a:
            guideline['url'] = a['guideline_url']
        adict[a['affiliation_id']] = { \
            'id': AFFILIATION_IRI_NAMESPACE + affiliation_public_id, \
            'label': a['affiliation_fullname'], \
            'guideline': guideline }
    return adict

#Just going to read this thing at global scope.  It's bad practice, but
# really this file should be read from a service somewhere anyway
affiliations_id_lookup = read_affiliations()

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

#We don't need to entity map everything, just some of the data nodes.
# For instance, not the interpretation.  When it occurs multiple times,
# it causes problems because the evaluations at each level
# are represented differently.
class EntityMap:
    EMtypes = set([VCI_VARIANT_TYPE, VCI_AGENT_TYPE, VCI_DISEASE_TYPE])
    def __init__(self, source, idtag=VCI_PK_KEY):
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
            # type_set = self.EMtypes.intersection(set(node[VCI_TYPE_KEY]))
            # if len(type_set) > 1:
            #     raise Exception("ERROR: More than one unique type match found for a single node.")
            # if len(type_set) == 1:
            #     etype = type_set.pop()
            if VCI_ITEM_TYPE_KEY in node and node[VCI_ITEM_TYPE_KEY] in self.EMtypes:
                etype = node[VCI_ITEM_TYPE_KEY]
                atid = fully_qualify(node[self.idtag])
                entity = self.entities[atid]
                for key in node:
                    if key in entity:
                        if entity[key] != node[key]:
                            print(key, '\n--------\n',entity[key], '\n---------------\n',node[key])
                            raise Exception('Incoherent Nodes')
                    else:
                        entity[key] = node[key]
    def get_entity(self,eid):
        return self.entities[eid]
    def get_transformed(self,eid):
        if eid in self.transformed:
            return self.transformed[eid]
        return None
    def add_transformed(self,eid,entity):
        self.transformed[eid] = entity

def fully_qualify(iri):
    fqiri = iri
    # if fqiri.startswith('/'):
    #     fqiri = IRI_BASE + fqiri

    # Extract PK/UUID from string that's assumed to be in older "@id" format (/[data type]/[data UUID]/)
    if fqiri.startswith('/'):
        fqiri_parts = fqiri.split('/')
        if len(fqiri_parts) > 2:
            if re.match(r'[a-f0-9]{8}\-[a-f0-9]{4}\-[a-f0-9]{4}\-[a-f0-9]{4}\-[a-f0-9]{12}', fqiri_parts[2], re.I):
                fqiri = fqiri_parts[2]

    return fqiri

def get_id(source):
    if isinstance(source, dict):
        sid = fully_qualify(source[VCI_PK_KEY])
    else:
        sid = fully_qualify(source)
    return sid

def get_affiliation(affiliation_id):
    if affiliation_id not in affiliations_id_lookup:
        raise Exception('Affiliation_id '+str(affiliation_id) + ' was not found in affiliation lookup file.')
    return affiliations_id_lookup[affiliation_id]

def add_contribution( user_input, affiliation, target, ondate, entities, role):

    if user_input is None and affiliation is None:
        raise Exception('Required User or Affiliation not provided when adding contributions.')

    # We're going to make some assumptions:
    # 1. at least one of user_input or affiliation_id is required
    # 2. if both exist then the affiliation is the "agent_for" and the "user_input" is the agent.
    # 3. if only the affiliation_id is available it is the agent, otherwise the "user_input" is the agent.
    # 4. affiliation is either a single value or None
    # 5. The agent_for is not formally added to the ClinGen-SEPIO model, so that there
    #    is not a set_affiliations method on contributions already.  We
    #    will probably need to modify this down the road.
    # 6. That information about the affiliation will be found in a local data
    #    file (which needs to be made accessible in another way)

    if user_input:
        # if there's also an affiliation then set it as the agent_for of the agent
        # this does presume that any single 'user' will always have one and
        # only one agent for the processing of the input file
        agent_for = None
        if affiliation is not None:
            agent_for = entities.get_transformed(affiliation['id'])
            if agent_for is None:
                agent_for = create_agent(affiliation['id'], affiliation['label'])
                entities.add_transformed(agent_for.get_id(), agent_for)

        userid = get_id(user_input)
        user = entities.get_entity(userid)

        if user  == {}:
            agent = create_agent(None, user_input, agent_for)
            entities.add_transformed(agent.get_id(), agent)
        else:
            agent = entities.get_transformed(userid)
            if agent is None:
                agent_name = None
                if VCI_AGENT_NAME_KEY in user:
                    agent_name = user[VCI_AGENT_NAME_KEY]
                agent = create_agent(userid, agent_name, agent_for)
                entities.add_transformed(agent.get_id(), agent)
    else:
        agent = entities.get_transformed(affiliation['id'])
        if agent is None:
            agent = create_agent(affiliation['id'], affiliation['label'])
            entities.add_transformed(agent.get_id(), agent)

    contribution = create_contribution( agent, ondate, role )

    target.add_contribution(contribution)

def add_contributions( source, affiliation_id, target, entities, ondate, role ):
    if isinstance( source, list ):
        for single_source in source:
            add_contribution( single_source, affiliation_id, target, ondate, entities, role)
    else:
        add_contribution( source, affiliation_id, target,  ondate, entities, role)

def convert_moi(moi, moiadj):
  m = re.search('^([A-Za-z \-]+)\(HP:[0-9]{7}\)$', moi)
  if m is None:
    if moi == 'Other':
      moi_label = moiadj
    else:
      moi_label = moi
  else:
    if moiadj:
      moi_label = m.group(1) + '('+moiadj+')'
    else:
      moi_label = m.group(1)
  return moi_label.strip()

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
def transform_root(vci):
    vi = VariantPathogenicityInterpretation( fully_qualify(vci[VCI_PK_KEY]) )
    modestring = vci[VCI_MODEINHERITANCE_KEY]
    modeadjstring = vci[VCI_MODEINHERITANCE_ADJECTIVE_KEY]
    mode = convert_moi(modestring, modeadjstring)
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
def transform_provisional_variant(vci_pv , interpretation, entities, publish_datetime):
    if isinstance(vci_pv, list):
        if len(vci_pv) > 1:
            exit()
        vci_pv = vci_pv[0]

    # approval date (aka last evaluated date)
    #  1st if approval review date exists in provisional variant section
    #  2nd if apporval date exists in provisional variant section
    #  last the last modified date of the provisional variant
    approval_date = ''
    if VCI_APPROVAL_REVIEW_DATE_KEY in vci_pv:
        approval_date = vci_pv[VCI_APPROVAL_REVIEW_DATE_KEY]
    elif VCI_APPROVAL_DATE_KEY in vci_pv:
        approval_date = vci_pv[VCI_APPROVAL_DATE_KEY]
    else:
        approval_date = vci_pv[VCI_LAST_MODIFIED_KEY]

    # approver can be either the affiliation or the user
    # if the affiliation exists then only use it for the approver
    # otherwise use the user associated by the following rules
    #  1st if an 'affiliation_id' in the provisional variant section
    #  2nd if there is a 'classificationApprover' in the provisional variant section
    #  3rd if there is a 'approvalSubmitter' in the provisional variant section
    #  4th defaults to the 'submitted_by" in the provisional variant section
    #
    approving_user = None
    affiliation_id = None
    if VCI_AFFILIATION_KEY in vci_pv:
        affiliation_id = vci_pv[VCI_AFFILIATION_KEY]
    elif VCI_CLASSIFICATION_APPROVER_KEY in vci_pv:
        approving_user = vci_pv[VCI_CLASSIFICATION_APPROVER_KEY]
    elif VCI_APPROVAL_SUBMITTER_KEY in vci_pv:
        approving_user = vci_pv[VCI_APPROVAL_SUBMITTER_KEY]
    else:
        approving_user = vci_pv[VCI_CONTRIBUTION_KEY]

    affiliation = None
    if affiliation_id:
        affiliation = get_affiliation(affiliation_id)
        if affiliation:
            guideline = affiliation['guideline']
            url = None
            if 'url' in guideline:
                url = guideline['url']
            if 'name' in guideline:
                assertion_method = create_assertion_method ( guideline['name'], url )
                interpretation.set_assertionMethod(assertion_method)

    add_contributions( approving_user, affiliation, interpretation, entities, approval_date, PROP_APPROVER_ROLE )

    # publisher is same as approver but for current datetime.
    now = datetime.datetime.now()
    add_contributions( approving_user, affiliation, interpretation, entities, publish_datetime.isoformat(), PROP_PUBLISHER_ROLE)

    if VCI_EVIDENCE_SUMMARY_KEY in vci_pv:
        interpretation.set_description( vci_pv[VCI_EVIDENCE_SUMMARY_KEY])

    interpretation.set_statementOutcome( convert_significance(vci_pv) )

# TODO talk to Ronak about these specialized significance terms.
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
    canonical_variant = entities.get_transformed(vci_variant_id)
    if canonical_variant is None:
        vci_variant = entities.get_entity(vci_variant_id)

        if VCI_CLINVAR_VARIANT_TITLE in vci_variant and vci_variant[VCI_CLINVAR_VARIANT_TITLE]:
            preferred_name = vci_variant[VCI_CLINVAR_VARIANT_TITLE]
        else:
            # for now - assume the first "NM_" or "NR_" transcript in the "hgvsNames.others" list is preferred
            preferred_name = None
            for other_hgvs in vci_variant[VCI_HGVS_NAMES_KEY][VCI_OTHERS_HGVS_KEY]:
                if re.search('^N[MR]_.*', other_hgvs):
                    preferred_name = other_hgvs
            if not preferred_name:
                # if no "NM_" or "NR_" is found in the transcripts then use 'GRCh38' as the default.
                preferred_name = vci_variant[VCI_HGVS_NAMES_KEY][VCI_GENOMIC_HGVS_38_KEY]

        # if vci_variant[VCI_VARIANT_SOURCE] == 'ClinVar':
        #     identifier = "%s:%s" % ('ClinVar', vci_variant['clinvarVariantId'])
        if VCI_CLINVAR_ID_KEY in vci_variant and vci_variant[VCI_CLINVAR_ID_KEY]:
            identifier = "%s:%s" % ('ClinVar', vci_variant[VCI_CLINVAR_ID_KEY])
        else:
            identifier = "%s:%s" % ('CAR', vci_variant[VCI_CANONICAL_ID_KEY])

        hgvs_names = vci_variant[VCI_HGVS_NAMES_KEY]
        dbsnp_ids = vci_variant['dbSNPIds']

        # build the following structure to initialize any VCI canonical allele
        # 1. 'identifier' 'CAR:CA999999' (CAid) or 'ClinVar:9999999' (variationId)
        # 2. 'hgvs_names' <array of entries, some genome build-specific + 'others' nested array)
        # 3. 'dbsnp_ids'  <array of dbsnp ids in rs form without rs prefix>
        # 4. 'preferred_name' <clinvar variant title if available, if blank find and pass the b38 hgvs expression>
        canonical_variant = CanonicalAllele( identifier=identifier, \
                                             hgvs_names=hgvs_names, \
                                             dbsnp_ids=dbsnp_ids, \
                                             preferred_name=preferred_name )

        entities.add_transformed(vci_variant_id, canonical_variant)
    return canonical_variant

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
    assessment = CriterionAssessment( vci_evaluation[VCI_PK_KEY] )
    criterion = criteria[ vci_evaluation[ VCI_CRITERIA_KEY] ]
    assessment.set_criterion( criterion )
    assessment.set_statementOutcome( term_map[ vci_evaluation[ VCI_CRITERIA_STATUS_KEY ] ] )
    if VCI_EVALUATION_EXPLANATION_KEY in vci_evaluation:
        assessment.set_description( vci_evaluation[ VCI_EVALUATION_EXPLANATION_KEY] )
    else:
        assessment.set_description( '' )
    assessment.set_variant( transform_variant( vci_evaluation[VCI_EVALUATION_VARIANT_KEY ] , entities) )
    #Have to do a little work to figure out the strength
    #VCI has both a modifier and a criteria modifier.  If these both exist they should be the
    #same, but if one or the other does not exist, then we should use the other one.
    crit_mod = ''
    # mod = ''
    if VCI_CRITERIA_MODIFIER_KEY in vci_evaluation:
        crit_mod = vci_evaluation[VCI_CRITERIA_MODIFIER_KEY]
    # if VCI_MODIFIER_KEY in vci_evaluation:
    #     mod = vci_evaluation[VCI_MODIFIER_KEY]
    # if crit_mod != mod:
    #     raise Exception
    defaultStrength = criterion.get_defaultStrength()
    strength = transform_strength( crit_mod, defaultStrength )

    ##REMOVING contributions from all but the top level interpretation. Only supporting Approver and Publisher contributions at the top-level for now
    ##TODO: Also add contribution to the evidence line if crit_mod != ''.  A little tricky since I hid the evidence line, but it's gettable
    ##affiliation_id = None
    ##if VCI_AFFILIATION_KEY in vci_evaluation:
    ##    affiliation_id = vci_evaluation[VCI_AFFILIATION_KEY]
    ##add_contributions( vci_evaluation[VCI_CONTRIBUTION_KEY], assessment, entities,vci_evaluation[VCI_LAST_MODIFIED_KEY], PROP_ASSESSOR_ROLE,affiliation_id)

    ## REMOVING logic to add frequency and prediction data until we can define the precise association between criterion assessments and specific evidence items.
    # #Now the evidence
    #  if VCI_FREQUENCY_KEY in vci_evaluation:
    #      frequencies = transform_frequency( vci_evaluation[VCI_FREQUENCY_KEY],  entities)
    #      add_evidenceItems( assessment, frequencies )
    #  if VCI_COMPUTATIONAL_KEY in vci_evaluation:
    #      predictions = transform_computational( vci_evaluation[VCI_COMPUTATIONAL_KEY], entities )
    #      add_evidenceItems( assessment, predictions )

    add_criterion_assessment(interpretation, assessment, strength)
    return vci_evaluation[VCI_CRITERIA_KEY], assessment

def transform_computational(source, entities):
    predictions = []
    vci_variant = source[VCI_VARIANT_KEY]
    #the form of this call suggests that the transform* functions should be member functions of entitymap
    variant = transform_variant(vci_variant,entities)
    compdata = source[VCI_COMPUTATIONAL_DATA_KEY]
    if VCI_CONSERVATION_KEY in compdata:
        predictions += transform_conservation_data( compdata[VCI_CONSERVATION_KEY] , variant )
    if VCI_CLINGEN_COMPUTATION_KEY in compdata:
        predictions += transform_clingen_comp_data( compdata[VCI_CLINGEN_COMPUTATION_KEY], variant )
    if VCI_OTHER_COMPUTATION_KEY in compdata:
        predictions += transform_other_comp_data( compdata[VCI_OTHER_COMPUTATION_KEY], variant )
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
                print('!',source[pred]['prediction'])
            #Have to sort out the scores....
            prediction = InSilicoPredictionScoreStatement()
            #We really also want to set the transcript, but the VCI is not returning that informaiton
            #prediction.set_transcript()
            prediction.set_canonicalAllele(variant)
            prediction.set_algorithm(pred)
            prediction.set_prediction(score)
            prediction.set_predictionType( term_map[VCI_MISSENSE_EFFECT_PREDICTOR] )
            predictions.append(prediction)
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
            prediction = prediction_score = None
            if s is not None:
                prediction_score = InSilicoPredictionScoreStatement()
                prediction_score.set_predictionType( term_map[ VCI_MISSENSE_EFFECT_PREDICTOR]  )
                prediction_score.set_canonicalAllele( variant )
                prediction_score.set_algorithm( pred )
                prediction_score.set_prediction(s)
            if p is not None:
                prediction = InSilicoPredictionStatement()
                prediction.set_predictionType( term_map[ VCI_MISSENSE_EFFECT_PREDICTOR]  )
                prediction.set_canonicalAllele( variant )
                prediction.set_algorithm( pred )
                prediction.set_statementOutcome(p)
            if prediction is None:
                predictions.append(prediction_score)
            else:
                if prediction_score is not None:
                    evidence_line = EvidenceLine()
                    evidence_line.add_evidenceItem( prediction_score )
                    prediction.add_evidenceLine( evidence_line )
                predictions.append(prediction)
    return predictions

#The values from VCI are not including a true/false on whether the thing is conserved. But we'll want that.
#When we get that, this will have to change.  What it will become is an AlleleConservation object
# with an evidence line to an AlleleConservationScore, which will be defined as below
def transform_conservation_data( source, variant ):
    results = []
    for constool in source:
        conservation = AlleleConservationScoreStatement()
        conservation.set_allele(variant)
        conservation.set_algorithm(constool)
        conservation.set_score(source[constool])
        results.append(conservation)
    return results

def transform_frequency( source, entities):
    popdata    = source[VCI_POPULATION_DATA_KEY]
    vci_variant = source[VCI_VARIANT_KEY]
    #the form of this call suggests that the transform* functions should be member functions of entitymap
    variant = transform_variant(vci_variant,entities)
    frequencies = []
    if VCI_ESP_KEY in popdata:
        esp_frequencies = transform_esp_data( popdata[VCI_ESP_KEY],  variant )
        frequencies += esp_frequencies
    if VCI_EXAC_KEY in popdata:
        exac_frequencies = transform_exac_data( popdata[VCI_EXAC_KEY], variant )
        frequencies += exac_frequencies
    if VCI_1000_GENOMES_KEY in popdata:
        tg_frequencies = transform_1000_genomes_data(popdata[VCI_1000_GENOMES_KEY], variant)
        frequencies += tg_frequencies
    #Add contributors to data nodes
    add_contributions_to_data( source , frequencies, entities )
    return frequencies

def add_contributions_to_data( data_source, data_targets, entities):
    submitters = data_source[VCI_CONTRIBUTION_KEY]
    modtime    = data_source[VCI_LAST_MODIFIED_KEY]
    for data in data_targets:
        add_contributions( submitters, None, data, entities, modtime, PROP_CURATOR_ROLE)

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

def transform_1000_genomes_data( source, variant ):
    alt_allele = variant.get_allele('GRCh37') #right?
    frequencies = []
    for pop in source:
        if pop != VCI_FREQUENCY_ALLELE_KEY:
            if pop in [VCI_1000_GENOMES_ESP_EA_POP, VCI_1000_GENOMES_ESP_AA_POP]:
                af = PopulationAlleleFrequencyStatement()
                af.set_ascertainment(METHOD_1KGESP)
            else:
                af = PopulationAlleleFrequencyStatement()
                af.set_ascertainment(METHOD_1KG)
            frequencies.append(af)
            af.set_allele( variant )
            af.set_population(  convert_1000_genomes_pop(pop) )
            if len( source[pop][VCI_ALLELE_COUNT_KEY] ) == 0:
                af.set_alleleNumber( 0 )
            else:
                af.set_alleleNumber( sum( source[pop][VCI_ALLELE_COUNT_KEY].values() ) )
            #TODO: Should the element be not here, or should it be 0 or should it be null?
            if af.get_alleleNumber() > 0:
                if alt_allele in source[pop][VCI_ALLELE_FREQUENCY_KEY]:
                    af.set_alleleFrequency( source[pop][VCI_ALLELE_FREQUENCY_KEY][alt_allele] )
                else:
                    af.set_alleleFrequency( 0.  )
            if alt_allele in source[pop][VCI_ALLELE_COUNT_KEY]:
                af.set_alleleCount( source[pop][VCI_ALLELE_COUNT_KEY][alt_allele] )
            else:
                af.set_alleleCount( 0 )
            hkey = '%s|%s' % (alt_allele,alt_allele)
            if hkey in source[pop][VCI_GENOTYPE_COUNT_KEY]:
                af.set_homozygousAlleleIndividualCount(  source[pop][VCI_GENOTYPE_COUNT_KEY][hkey] )
            else:
                af.set_homozygousAlleleIndividualCount( 0 )
    return frequencies

def transform_exac_data( source, variant ):
    frequencies = []
    for pop in source:
        if pop != VCI_FREQUENCY_ALLELE_KEY:
            af = PopulationAlleleFrequencyStatement()
            af.set_ascertainment(METHOD_EXAC)
            frequencies.append(af)
            af.set_population( convert_exac_pop(pop) )
            af.set_allele( variant )
            if VCI_ALLELE_COUNT_KEY in source[pop]:
                af.set_alleleCount(  source[pop][VCI_ALLELE_COUNT_KEY] )
            else:
                af.set_alleleCount( 0 )
            if VCI_ALLELE_NUMBER_KEY in source[pop]:
                af.set_alleleNumber(source[pop][VCI_ALLELE_NUMBER_KEY])
            else:
                af.set_alleleNumber(0)
            if af.get_alleleNumber() > 0:
                af.set_alleleFrequency(source[pop][VCI_ALLELE_FREQUENCY_KEY])
            if VCI_HOMOZYGOUS_GENOTYPE_COUNT_KEY in source[pop]:
                af.set_homozygousAlleleIndividualCount( source[pop][VCI_HOMOZYGOUS_GENOTYPE_COUNT_KEY] )
            else:
                af.set_homozygousAlleleIndividualCount( 0 )
    return frequencies

def transform_esp_data(source,variant):
    ref_allele = variant.get_ref_allele('GRCh37')
    alt_allele = variant.get_allele('GRCh37')
    frequencies = []
    for pop in source:
        if pop != VCI_FREQUENCY_ALLELE_KEY:
            af = PopulationAlleleFrequencyStatement( )
            af.set_ascertainment(METHOD_ESP)
            frequencies.append(af)
            af.set_population( convert_esp_pop(pop) )
            af.set_allele( variant )
            #Assumes ESP uses v37
            if (len(source[pop][VCI_ALLELE_COUNT_KEY]) > 0) and (alt_allele in source[pop][VCI_ALLELE_COUNT_KEY]):
                af.set_alleleCount( source[pop][VCI_ALLELE_COUNT_KEY][alt_allele] )
            else:
                af.set_alleleCount( 0 )
            af.set_alleleNumber( sum(source[pop][VCI_ALLELE_COUNT_KEY].values()) )
            if af.get_alleleNumber() > 0:
                f = 1.*af.get_alleleCount() / af.get_alleleNumber()
                af.set_alleleFrequency(f)
            if len(source[pop][VCI_GENOTYPE_COUNT_KEY]) > 0:
                homkey = '%s%s' % (alt_allele,alt_allele)
                af.set_homozygousAlleleIndividualCount(source[pop][VCI_GENOTYPE_COUNT_KEY][homkey])
                hetkey1 = '%s%s' % (ref_allele,alt_allele)
                hetkey2 = '%s%s' % (alt_allele,ref_allele)
                if hetkey1 in  source[pop][VCI_GENOTYPE_COUNT_KEY]:
                    af.set_heterozygousAlleleIndividualCount(source[pop][VCI_GENOTYPE_COUNT_KEY][hetkey1])
                elif hetkey2 in  source[pop][VCI_GENOTYPE_COUNT_KEY]:
                    af.set_heterozygousAlleleIndividualCount( source[pop][VCI_GENOTYPE_COUNT_KEY][hetkey2] )
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
    elif modifier == 'stand-alone':
        mod = 'Stand Alone'
    else:
        exit()
    strength = ' '.join( [path, mod] )
    return strength

#Returns a dictionary that maps from criterion id (e.g. "PS2") to
# a criterion assessment entity
def transform_evaluations(evaluation_list,interpretation,entities):
    evaluation_map = {}
    criteria = read_criteria()
    if not isinstance(evaluation_list, list):
        raise Exception
    for vci_eval in evaluation_list:
        if VCI_CRITERIA_STATUS_KEY not in vci_eval or vci_eval[ VCI_CRITERIA_STATUS_KEY ] == VCI_CRITERIA_NOT_EVALUATED:
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
        # add_contributions_to_data( ee_node, [info], entities )
        info.set_description( ee_node[ VCI_EVIDENCE_DESCRIPTION_KEY] )
        sources = transform_articles(ee_node[ VCI_ARTICLES_KEY], interpretation, entities )
        for source in sources:
            info.add_source(source)
        key = ( ee_node[VCI_CATEGORY_KEY], ee_node[ VCI_SUBCATEGORY_KEY] )
        #possible_rules will be a list of tuples of rules.  Each tuple is a group that can
        #only have one thing met.
        possible_rules = extra_evidence_map[key]
        #This is a map from possible_rule tuples to assessments (met/not met) to a list of rules
        # So if a given rule tuple has only met rules then it will be [tuple][met]=[met assessments] and [tuple][not met] = []
        # But if there are also not met rules thn it will be [tuple][not met]=[not met assessments]
        rule_groups = defaultdict( lambda: defaultdict (list) )
        found = False
        for rule_set in possible_rules:
            for rule in rule_set:
                if rule in evalmap:
                    found = True
                    assessment = evalmap[rule]
                    rule_groups[rule_set][assessment.get_statementOutcome().get_label()].append(assessment)
                    #add_evidenceItems( assessment, [info])
        if not found:
            print("Did not find any evaluated criteria for this data: %s "% ee_node['uuid'])
        else:
            for rule_set in rule_groups:
                if len(rule_groups[rule_set]['Met']) > 0:
                    #attach evidence to the met rules
                    for assessment in rule_groups[rule_set]['Met']:
                        add_evidenceItems( assessment,[info])
                else:
                    #But if there are no met rule, attach to the not met rules
                    for assessment in rule_groups[rule_set]['Not Met']:
                        add_evidenceItems( assessment,[info])

#Note that we don't necessarily have the full VCI disease node when we get into this function.
# It could be anything from a bare IRI to a full node or anything in between.  The first two lines
# standardize the "local" information to the "global" information node.
def transform_condition(vci_local_disease,interpretation,entities,mode):
    vci_disease_id = get_id(vci_local_disease)
    vci_disease = entities.get_entity(vci_disease_id)
    disease = entities.get_transformed(vci_disease_id)
    if disease is None:
        disease_ontology = 'MONDO:'
        disease_code = vci_disease[VCI_PK_KEY]
        if not disease_code.startswith('MONDO'):
            raise Exception("Expected a MONDO disease identifier")
        disease_code = disease_code.split('_')[-1]
        disease_name = vci_disease[VCI_DISEASE_TERM_KEY]
        disease = create_disease(disease_ontology, disease_code, disease_name)
        entities.add_transformed(vci_disease_id, disease)
    condition = GeneticCondition()
    condition.add_disease(disease)
    if mode!= '': condition.set_inheritancePattern( mode )
    interpretation.add_condition(condition)

def transform(jsonf, payload, publish_datetime):
    vci = None

    if jsonf:
        vci = json.load(jsonf)

    if payload:
        vci = json.loads(payload)
        #vci = json.loads(payload)['interpretation']

    if VCI_MODEINHERITANCE_KEY not in vci:
        vci[VCI_MODEINHERITANCE_KEY] = ''
    if VCI_MODEINHERITANCE_ADJECTIVE_KEY not in vci:
        vci[VCI_MODEINHERITANCE_ADJECTIVE_KEY] = ''
    entities = EntityMap(vci)
    interpretation, inheritance = transform_root(vci)
    if VCI_PROVISIONAL_VARIANT_KEY in vci:
        transform_provisional_variant(vci[VCI_PROVISIONAL_VARIANT_KEY],interpretation,entities, publish_datetime)
    elif VCI_V1_PROVISIONAL_VARIANT_KEY in vci:
        transform_provisional_variant(vci[VCI_V1_PROVISIONAL_VARIANT_KEY],interpretation,entities, publish_datetime)
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
    except (KeyError, IndexError):
        logging.warning('No criteria evaluations found.  Proceeding')
        eval_map = {}
    try:
        if VCI_CURATED_EVIDENCE_KEY in vci:
            transform_evidence(vci[VCI_CURATED_EVIDENCE_KEY], interpretation, entities, eval_map)
        else:
            transform_evidence(vci[VCI_V1_EXTRA_EVIDENCE_KEY], interpretation, entities, eval_map)
    except KeyError:
        logging.warning('No evidence found.  Proceeding')

    return interpretation,entities

def transform_json_file(inf, outf, out_style, publish_datetime):
    interp, ents = transform(inf, None, publish_datetime)
    #The idea of flatten would be that the root node would contain fully specified descriptions of all the entities, and then
    # the interpretation, which would be written using only IDs.
    # To implement this, just create a new dict, put the interpretation into it, and write the entites into it from the EntityMap
    # Then dump that envelope, and modify the interpretaion encoder to have a mode that keeps track of depth and writes full nodes
    # for the top level.
    # The other option is not to include this at all, and rely on JSON-LD libraries to do any appropriate flattening.
    # TODO: Decide and implement (if required) or remove option (if not)
    if out_style == 'flat':
        raise Exception('flatten not implemented yet')
    json.dump(interp, outf, sort_keys=True, indent=4, separators=(',', ': '), cls=InterpretationEncoder, out_style=out_style, ensure_ascii=False)

def transform_json_input(payload, out_style):
    # Temporarily passing in "now" timestamp (probably should use/supply publication timestamp from VCI)
    interp,ents = transform(None, payload, datetime.datetime.now())
    if out_style == 'flat':
        raise Exception('flatten not implemented yet')
    return json.dumps(interp, sort_keys=True, indent=4, separators=(',', ': '), cls=InterpretationEncoder, out_style=out_style, ensure_ascii=False)

def valid_date(s):
    try:
        return parse(s)
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)

def test():
    transform_json_file(open('test_data/test_interp_1.vci.json'), open('test_data/test_interp_1.cg-sepio.json'), 'first', datetime.datetime.now())
#    transform_json_file(open('test-vci.json'), open('test-cgsepio.json', 'w'), 'first', datetime.datetime.now())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input',  nargs='?', type=argparse.FileType('r'), default=sys.stdin,
            help='Path to an input JSON file created by the VCI (defaults to stdin)')
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
            help='Output path for ClinGen-SEPIO JSON file to be created (defaults to stdout)')
    parser.add_argument("-s", "--output-style", type=str, choices=['full', 'first', 'flat'],
            help="full: expand all nodes, first: expand first node, flat: define entities outside the interpretation", default = 'first')
    parser.add_argument("-p",
                    "--publish-datetime",
                    help="The publish date of the records, defaults to now - format YYYY-MM-DD or YYYY-MM-DDThh:mm:ss.fffZ",
                    required=False,
                    type=valid_date,
                    default=datetime.datetime.now())
    args = parser.parse_args()
    transform_json_file(args.input, args.output, args.output_style, args.publish_datetime)
