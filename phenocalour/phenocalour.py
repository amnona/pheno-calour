from logging import getLogger
from collections import defaultdict
from csv import DictReader
import webbrowser
from operator import itemgetter
from pkg_resources import resource_filename

from calour.database import Database

logger = getLogger(__name__)


class Phenotype(Database):
    '''Phenotype database module for calour, using data from:
    Hiding in Plain Sight: Mining Bacterial Species Records for Phenotypic Trait Information
    Albert Barberán, Hildamarie Caceres Velazquez, Stuart Jones, Noah Fierer
    DOI: 10.1128/mSphere.00237-17
    '''
    def __init__(self, exp=None):
        super().__init__(exp=exp, database_name=' phenotype', methods=['get'])
        if '_calour_ijsem_pheno_db' not in exp.exp_metadata:
            seqdict = {}
            datfile = resource_filename(__package__, 'data/pheno-calour-data.txt')
            with open(datfile) as dbf:
                for crow in DictReader(dbf, delimiter='\t'):
                    cseq = crow['Sequence']
                    seqdict[cseq]={}
                    for ck,cv in crow.items():
                        if ck == 'Sequence':
                            continue
                        seqdict[cseq][ck]=cv
                    if seqdict[cseq]['Habitat'] == 'other':
                        seqdict[cseq]['Habitat'] = seqdict[cseq]["If 'other' was chosen above, please enter a habitat below"]
                exp.exp_metadata['_calour_ijsem_pheno_db'] = seqdict
        self.phenotype_data = exp.exp_metadata['_calour_ijsem_pheno_db']
        self.show_fields = ['Genus name','species name','gram status','Habitat','motility','cell shape','mean width (microns)','mean length (microns)','pigment production','cell aggregation','oxygen preference','Sole carbon substrate use','Metabolism assays','spore production','optimal NaCl concentration for growth (%)','NaCl concentration range at which growth occurred (%)','temperature optimum for growth (degC)','temperature range at which growth occurred (degC)','pH optimum for growth','pH range at which growth occurred','DNA GC content (mol%)','16S rDNA accession number']

    def get_seq_annotation_strings(self, feature):
        '''Get nice string summaries of annotations for a given sequence

        Parameters
        ----------
        sequence : str
            the DNA sequence to query the annotation strings about

        Returns
        -------
        shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
            a list of:
                annotationdetails : dict
                    'seqid' : str, the sequence annotated
                    'annotationtype : str
                    ...
                annotationsummary : str
                    a short summary of the annotation
        '''
        shortdesc = []
        if self.phenotype_data is None:
            return []
        if feature not in self.phenotype_data:
            return []
        if self.show_fields is None:
            for k,v in self.phenotype_data[feature].items():
                if v == '':
                    continue
                shortdesc.append(({'annotationtype': 'other', 'feature': feature}, '%s: %s' % (k,v)))
            shortdesc = sorted(shortdesc, key=itemgetter(1))
        else:
            for cfield in self.show_fields:
                if cfield in self.phenotype_data[feature]:
                    shortdesc.append(({'annotationtype': 'other', 'feature': feature}, '%s: %s' % (cfield,self.phenotype_data[feature][cfield])))

        return shortdesc

    def show_annotation_info(self, annotation):
        '''Show the website for the sequence

        Parameters
        ----------
        annotation : dict
            should contain 'feature'
        '''
        # open in a new tab, if possible
        journal_base = 'http://ijs.microbiologyresearch.org/content/journal/ijsem/'
        new = 2
        seq = annotation['feature']
        data = self.phenotype_data[seq]

        if 'article doi' in data:
            address = journal_base + data['article doi']
        else:
            logger.warn('no doi for phenotype database - cannot open link')
            return
        webbrowser.open(address, new=new)
