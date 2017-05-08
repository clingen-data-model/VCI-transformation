#To run these tests, run "pytest" from the command line
from VCI2DMWG import *
from interpretation import *
import json

class MockNode(Node):
    def add_contribution(self,c):
        self.contribution = c

def test_add_contrib():
    ondate = 'DATESTRING'
    role = 'ROLESTRING'
    source = {}
    sourcetxt = '''
    {
    "submitted_by": {
         "@type": [ "user", "item" ],
         "@id": "/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/",
         "title": "Amanda Thomas",
         "first_name": "Amanda",
         "last_name": "Thomas",
         "lab": "/labs/curator/",
         "uuid": "1db40a7e-e6dc-46dd-9e17-44efbf53d508",
         "email": "amanda.thomas@aruplab.com"
         }
    }
    '''
    source=json.loads(sourcetxt)    
    ents = EntityMap(source)
    target = MockNode()
    add_contribution( source['submitted_by'], target, ondate, ents, role)
    contribution = target.contribution
    assert isinstance(contribution, Contribution)
    agent = contribution.get_agent()
    assert agent.get_name() == 'Amanda Thomas'
    assert agent.get_id() == IRI_BASE + '/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/'
    assert contribution.get_ondate() == ondate
    assert contribution.get_role() == role

def test_add_contrib_IRI():
    tlist = []
    ondate = 'DATESTRING'
    role = 'ROLESTRING'
    source = {}
    sourcetxt = '''
    {
    "previous_user":{
        "@id":"/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/",
        "@type": ["user"],
        "title": "Jed Clampett"
    },
    "submitted_by":  "/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/"
    }
    '''
    source=json.loads(sourcetxt)
    ents = EntityMap(source)
    subby=source['submitted_by']
    target = MockNode()
    add_contribution( subby, target, ondate, ents, role)
    contribution = target.contribution
    agent = contribution.get_agent()
    assert agent.get_id() == IRI_BASE + '/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/'
    assert agent.get_name() == 'Jed Clampett'
    assert contribution.get_ondate() == ondate
    assert contribution.get_role() == role   

def test_add_contrib_IRI_only():
    tlist = []
    ondate = 'DATESTRING'
    role = 'ROLESTRING'
    source = {}
    sourcetxt = '''
    {
    "submitted_by":  "/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/"
    }
    '''
    source=json.loads(sourcetxt)
    ents = EntityMap(source)
    subby=source['submitted_by']
    target = MockNode()
    add_contribution( subby, target, ondate, ents, role)
    contribution = target.contribution
    agent = contribution.get_agent()
    assert agent.get_id() == IRI_BASE + '/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/'
    assert DMWG_AGENT_NAME_KEY not in agent.data
    assert contribution.get_ondate() == ondate
    assert contribution.get_role() == role   


def test_canonical_have_car():
    carid='http://reg.genome.network/allele/CA229854'
    sourcetxt='''
    {
    "variant": {
        "date_created": "2017-03-22T23:03:46.771883+00:00",
        "clinVarRCVs": [],
        "hgvsNames": {
            "GRCh38": "NC_000012.12:g.102846935G>T",
            "others": [
                "NG_008690.1:g.75668C>A",
                "NM_000277.1:c.929C>A",
                "NP_000268.1:p.Ser310Tyr",
                "P00439:p.Ser310Tyr"
            ],
            "GRCh37": "NC_000012.11:g.103240713G>T"
        },
        "last_modified": "2017-03-22T23:03:46.771896+00:00",
        "clinvarVariantTitle": "NM_000277.1(PAH):c.929C>A (p.Ser310Tyr)",
        "source": "ClinVar",
        "clinvarVariantId": "102898",
        "carId": "%s",
        "variationType": "single nucleotide variant",
        "molecularConsequenceList": [
            {
                "term": "missense variant",
                "soId": "SO:0001583",
                "hgvsName": "NM_000277.1:c.929C>A"
            }
        ],
        "@id": "/variants/4cfa029b-3865-47d3-95eb-39e9dffa647b/",
        "variant_identifier": "102898",
        "submitted_by": "/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/",
        "otherDescription": "",
        "@type": [
            "variant",
            "item"
        ],
        "status": "in progress",
        "variation_type": "single nucleotide variant",
        "clinVarSCVs": [],
        "dbSNPIds": [
            "62642913"
        ],
        "schema_version": "3",
        "associatedPathogenicities": [],
        "molecular_consequence": "missense variant",
        "uuid": "4cfa029b-3865-47d3-95eb-39e9dffa647b"
    }
    }''' % carid
    source=json.loads(sourcetxt)
    ents = EntityMap(source)
    rep = canonicalizeVariant(source['variant'])
    canonical_id = rep['@id']
    assert canonical_id == carid

def test_canonical_no_car():
    carid=''
    sourcetxt='''
    {
    "variant": {
        "hgvsNames": {
            "GRCh38": "NC_000012.12:g.102846935G>T",
            "others": [
                "NG_008690.1:g.75668C>A",
                "NM_000277.1:c.929C>A",
                "NP_000268.1:p.Ser310Tyr",
                "P00439:p.Ser310Tyr"
            ],
            "GRCh37": "NC_000012.11:g.103240713G>T"
        },
        "clinvarVariantTitle": "NM_000277.1(PAH):c.929C>A (p.Ser310Tyr)",
        "source": "ClinVar",
        "clinvarVariantId": "102898",
        "carId": "%s",
        "@id": "/variants/4cfa029b-3865-47d3-95eb-39e9dffa647b/",
        "variant_identifier": "102898",
        "@type": [
            "variant",
            "item"
        ]
        }
    }''' % carid
    source=json.loads(sourcetxt)
    ents = EntityMap(source)
    rep = canonicalizeVariant(source['variant'])
    canonical_id = rep['@id']
    #This is the id that the Baylor Allele Registry returns for this allele.
    assert canonical_id == 'http://reg.genome.network/allele/CA229854'


def test_freq_esp():
    sourcetxt=''' {
    "population": {
     "variant": {"@id": "fakeid", "@type": [ "variant", "item" ], "carId": "", "hgvsNames": { "others": [ "NG_009830.1:g.108338C>T", "NM_000051.3:c.6919C>T", "NP_000042.3:p.Leu2307Phe", "Q13315:p.Leu2307Phe", "LRG_135p1:p.Leu2307Phe", "LRG_135t1:c.6919C>T", "LRG_135:g.108338C>T" ], "GRCh37": "NC_000011.9:g.108196896C>T", "GRCh38": "NC_000011.10:g.108326169C>T" } },

    "submitted_by": "/users/123",
    "last_modified":  "2017-03-22T23:05:18.070862+00:00",
    "schema_version": "1",
    "populationData": {
    "desiredCI": 95,
    "esp": {
        "_tot": { "ac": { "C": 12974, "T": 24 },
                  "gc": { "CC": 6475, "TC": 24, "TT": 0 } },
        "_extra": {
            "avg_sample_read": 106,
            "rsid": "rs56009889",
            "chrom": "11",
            "alt": "T",
            "hg19_start": 108196896,
            "ref": "C" } 
    }}}}'''
    source =json.loads(sourcetxt)
    ents = EntityMap(source)
    rep = canonicalizeVariant(source['population']['variant'])
    canonical_id = rep['@id']
    frequencies = transform_frequency( source['population'],  ents )
    assert len(frequencies) == 1
    freqnode = frequencies[0]
    #Add getters instead?
    assert freqnode.data[DMWG_TYPE_KEY] == DMWG_ALLELE_FREQUENCY_TYPE
    assert freqnode.data[DMWG_ASCERTAINMENT_KEY] == DMWG_ESP
    assert freqnode.data[DMWG_POPULATION_KEY] == DMWG_COMBINED_POP
    assert freqnode.data[DMWG_VARIANT_KEY].get_id() == canonical_id
    assert freqnode.data[DMWG_ALLELE_COUNT_KEY] == 24
    assert freqnode.data[DMWG_ALLELE_NUMBER_KEY] == 24 + 12974
    assert abs(freqnode.data[DMWG_ALLELE_FREQUENCY_KEY] - 24./(24+12974)) < 0.0001
    assert freqnode.data[DMWG_HOMOZYGOUS_INDIVIDUAL_COUNT_KEY] == 0
    assert freqnode.data[DMWG_HETEROZYGOUS_INDIVIDUAL_COUNT_KEY] == 24


def test_freq_esp_empty():
    sourcetxt=''' {
     "population": {
"variant": {"@id": "fakeid", "@type": [ "item", "variant" ], "carId": "", "hgvsNames": { "others": [ "NG_009830.1:g.108338C>T", "NM_000051.3:c.6919C>T", "NP_000042.3:p.Leu2307Phe", "Q13315:p.Leu2307Phe", "LRG_135p1:p.Leu2307Phe", "LRG_135t1:c.6919C>T", "LRG_135:g.108338C>T" ], "GRCh37": "NC_000011.9:g.108196896C>T", "GRCh38": "NC_000011.10:g.108326169C>T" }},
    "schema_version": "1",
    "submitted_by": "/users/123",
    "last_modified":  "2017-03-22T23:05:18.070862+00:00",
    "populationData": {
    "desiredCI": 95,
    "esp": {
        "_tot": { "ac": {}, "gc":{} },       
        "_extra": {}
         }}}}'''
    source =json.loads(sourcetxt)
    ents = EntityMap(source)
    rep = canonicalizeVariant(source['population']['variant'])
    canonical_id = rep['@id']
    freqs = transform_frequency( source['population'], ents )
    assert len(freqs) == 1
    freq = freqs[0]
    assert isinstance(freq, AlleleFrequency)
    assert freq.data[DMWG_VARIANT_KEY].get_id() == canonical_id
    assert freq.data[DMWG_ALLELE_COUNT_KEY] == 0 
    assert freq.data[DMWG_ALLELE_NUMBER_KEY] == 0

def test_freq_exac():
    sourcetxt=''' {
        "population": {
        "variant": {"@id": "fakeid", "@type": [ "item", "variant" ] , "carId": "", "hgvsNames": { "others": [ "NG_009830.1:g.108338C>T", "NM_000051.3:c.6919C>T", "NP_000042.3:p.Leu2307Phe", "Q13315:p.Leu2307Phe", "LRG_135p1:p.Leu2307Phe", "LRG_135t1:c.6919C>T", "LRG_135:g.108338C>T" ], "GRCh37": "NC_000011.9:g.108196896C>T", "GRCh38": "NC_000011.10:g.108326169C>T" }},
        "submitted_by": "/users/123",
        "last_modified":  "2017-03-22T23:05:18.070862+00:00",
        "schema_version": "1",
        "populationData": {
        "desiredCI": 95,
        "exac": {
            "amr": {
            "an": 11574,
            "af": 0.00017280110592707794,
            "ac": 2,
            "hom": 0
            },
            "eas": {
            "an": 8642,
            "af": 0,
            "ac": 0,
            "hom": 0
            },
            "afr": {
            "an": 10404,
            "af": 0.00009611687812379854,
            "ac": 1,
            "hom": 0
            },
            "_tot": {
            "an": 121374,
            "af": 0.0011040255738461284,
            "ac": 134,
            "hom": 1
            },
            "_extra": {
            "pos": 108196896,
            "alt": "T",
            "ref": "C",
            "chrom": "11"
            }
        }}}}'''
    source =json.loads(sourcetxt)
    ents = EntityMap(source)
    rep = canonicalizeVariant(source['population']['variant'])
    canonical_id = rep['@id']
    freqs = transform_frequency( source['population'], ents )
    assert len(freqs) == 4
    chex = { DMWG_LATINO_POP: (11574, 0.000172, 2, 0),\
             DMWG_EAST_ASIAN_POP: (8642, 0., 0, 0), \
             DMWG_AFRICAN_POP: (10404, 0.000096, 1, 0), \
             DMWG_COMBINED_POP: (121374, 0.001104, 134, 1) }
    for freq in freqs:
        assert freq.data[DMWG_TYPE_KEY] == DMWG_ALLELE_FREQUENCY_TYPE
        assert freq.data[DMWG_ASCERTAINMENT_KEY] == DMWG_EXAC
        values = chex[freq.data[DMWG_POPULATION_KEY]]
        assert freq.data[DMWG_VARIANT_KEY].get_id() == canonical_id
        assert freq.data[DMWG_ALLELE_COUNT_KEY] == values[2]
        assert freq.data[DMWG_ALLELE_NUMBER_KEY] == values[0]
        assert abs(freq.data[DMWG_ALLELE_FREQUENCY_KEY] - values[1]) < 0.0001
        assert freq.data[DMWG_HOMOZYGOUS_INDIVIDUAL_COUNT_KEY] == values[3]

def test_thousandGenomes():
    src='''
    {
    "population": {
        "variant": {"@id": "fakeid", "@type": [ "item", "variant" ] , "carId": "", "hgvsNames": {  "GRCh37": "NC_000012.11:g.103310879G>C", "GRCh38": "NC_000012.12:g.102917101G>C" }},
        "schema_version": "1",
        "submitted_by": "/users/123",
        "last_modified":  "2017-03-22T23:05:18.070862+00:00",
        "populationData": {
            "desiredCI": 95,
            "tGenomes": {
                "amr": { "af": { "G": 1 }, "gc": { "G|G": 347 }, "gf": { "G|G": 1 }, "ac": { "G": 694 } },
                "eas": { "af": { "G": 1 }, "gc": { "G|G": 504 }, "gf": { "G|G": 1 }, "ac": { "G": 1008 } },
                "_extra": { "var_class": "SNP", "name": "rs1801145", "alt": "C", "ref": "G" }, 
                "afr": { "af": { "G": 1 }, "gc": { "G|G": 661 }, "gf": { "G|G": 1 }, "ac": { "G": 1322 } }, 
                "espaa": { "af": { "G": 1 }, "gc": {}, "gf": {}, "ac": { "G": 4406 } }, 
                "eur": { "af": { "G": 0.999005964214712, "C": 0.00099403578528827 }, "gc": { "G|G": 502, "C|G": 1 }, "gf": { "G|G": 0.998011928429423, "C|G": 0.00198807157057654 }, "ac": { "G": 1005, "C": 1 } }, 
                "espea": { "af": { "G": 0.999302, "C": 0.000697674 }, "gc": {}, "gf": {}, "ac": { "G": 8594, "C": 6 } }, 
                "_tot": { "af": { "G": 0.999800319488818, "C": 0.000199680511182109 }, "gc": { "G|G": 2503, "C|G": 1 }, "gf": { "G|G": 0.999600638977636, "C|G": 0.000399361022364217 }, "ac": { "G": 5007, "C": 1 } }, 
                "sas": { "af": { "G": 1 }, "gc": { "G|G": 489 }, "gf": { "G|G": 1 }, "ac": { "G": 978 } } 
            }
        }
    } }'''
    source = json.loads(src)
    ents = EntityMap(source)
    rep= canonicalizeVariant(source['population']['variant'])
    canonical_id = rep['@id']
    freqs = transform_frequency( source['population'],  ents)
    assert len(freqs) == 8
    ex = { DMWG_LATINO_POP: (DMWG_1000_GENOMES,694, 0., 0, 0),\
           DMWG_EAST_ASIAN_POP: (DMWG_1000_GENOMES,1008, 0., 0, 0), \
           DMWG_AFRICAN_POP: (DMWG_1000_GENOMES,1322, 0., 0, 0), \
           DMWG_AFRICAN_AMERICAN_POP: (DMWG_1000_GENOMES_ESP,4406, 0., 0, 0), \
           DMWG_EUROPEAN_POP: (DMWG_1000_GENOMES,1006, 0.000994, 1, 0), \
           DMWG_EUROPEAN_AMERICAN_POP: (DMWG_1000_GENOMES_ESP,8600, 0.000697, 6, 0), \
           DMWG_COMBINED_POP: (DMWG_1000_GENOMES,5008, 0.0001996, 1, 0), 
           DMWG_SOUTH_ASIAN_POP: (DMWG_1000_GENOMES,978, 0., 0, 0) }
    for freq in freqs:
        assert freq.data[DMWG_VARIANT_KEY].get_id() == canonical_id
        values = ex[freq.data[DMWG_POPULATION_KEY]]
        assert freq.data[DMWG_ASCERTAINMENT_KEY] == values[0]
        assert freq.data[DMWG_TYPE_KEY] == DMWG_ALLELE_FREQUENCY_TYPE
        assert freq.data[DMWG_ALLELE_NUMBER_KEY] == values[1]
        assert abs(freq.data[DMWG_ALLELE_FREQUENCY_KEY] - values[2]) < 0.0001
        assert freq.data[DMWG_ALLELE_COUNT_KEY] == values[3]
        assert freq.data[DMWG_HOMOZYGOUS_INDIVIDUAL_COUNT_KEY] == values[4]

def test_conservation():
    sourcetext='''
    {
    "variant": { 
        "hgvsNames": { "GRCh38": "NC_000012.12:g.102846935G>T" }, 
        "@id": "/variants/4cfa029b-3865-47d3-95eb-39e9dffa647b/",
        "@type": [ "variant", "item" ],
        "carId":""
    },
    "computational": {
        "@id": "/computational/d9d3e60d-50b1-424b-8bcc-c892bba0da22/",
        "@type": ["item", "data"],
        "last_modified": "2017-03-22T23:08:28.110287+00:00",
        "variant":  "/variants/4cfa029b-3865-47d3-95eb-39e9dffa647b/",
        "computationalData": {
            "conservation": {
                "phylop20way": "0.948",
                "phylop7way": "9.598",
                "phastconsp20way": "0.953",
                "siphy": "19.1358",
                "gerp": "5.38",
                "phastconsp7way": "1"
            }
        },
        "submitted_by": {
            "first_name": "Amanda",
            "last_name": "Thomas",
            "title": "Amanda Thomas",
            "@id": "/users/1db40a7e-e6dc-46dd-9e17-44efbf53d508/",
            "@type": ["user", "item"]
        }
    }}'''
    source = json.loads(sourcetext)
    ents = EntityMap(source)
    print source['variant'], type(source['variant'])
    rep = canonicalizeVariant(source['variant'])
    canonical_id = rep['@id']
    comps = transform_computational( source['computational'],  ents)
    assert len(comps) == 6
    chex={"phylop20way": 0.948, \
          "phylop7way": 9.598, \
          "phastconsp20way": 0.953, \
          "siphy": 19.1358, \
          "gerp": 5.38, \
          "phastconsp7way": 1 }
    for comp in comps:
        assert comp.data[DMWG_VARIANT_KEY].get_id() == canonical_id
        score = chex[comp.data[DMWG_CONSERVATION_METHOD_KEY]]
        assert comp.data[DMWG_CONSERVATION_SCORE_KEY] == score


    


