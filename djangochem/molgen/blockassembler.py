import collections
import datetime
import itertools
import time
import numpy
import logging

from rdkit.Chem.AllChem import CalcExactMolWt, MolFromSmiles, MolToSmiles

from blocks.blockset import BlockSet
from .blockset_storage_json import BlockSetStorageJSON

from .blockreactor import rxn_envs, cleaner_env, old_style_rxn_envs, old_cleaner_env
from .blockreactor import make_cleaner_func
from .blockreactor import mass_under, contains_not_smarts

from .utils import put_block, put_mol, put_smiles
logger = logging.getLogger()

from .rdkit_grafting_reaction import run_grafting_reactions

SMARTS_TO_FILTER = ["[#7]-[#7]",
                    "[#8]~[#8]",
                    "[#7]-[#8]",
                    "[#16][#16]",
                    "[SX4](=O)(=O)[#8]",
                    "[SX4](=O)(=O)[#7]",
                    "[CX3H1](=O)",
                    "[#7][SX2]",
                    "[CX3](=O)[S]",
                    "[CX3](=O)O[CX3](=O)",
                    "[N]#[C][#7]",
                    "[P][#7]",
                    "[P][S]",
                    "[P][P]",
                    "[P][Si]",
                    "[#7][Si]",
                    "[!#6][#6]#[#6]",
                    "[!#6R][#6;R0]=[#6;R0]"
                    ]

LASER_SMARTS_TO_FILTER = ["[#8]~[#8]",
                          "[#16][#16]",
                          "[SX4](=O)(=O)[#8]",
                          "[SX4](=O)(=O)[#7]",
                          "[CX3H1](=O)",
                          "[#7][SX2]",
                          "[CX3](=O)[S]",
                          "[CX3](=O)O[CX3](=O)",
                          "[N]#[C][#7]",
                          "[P][#7]",
                          "[P][S]",
                          "[P][P]",
                          "[P][Si]",
                          "[#7][Si]",
                          "[!#6][#6]#[#6]"
                          ]

CEP_SMARTS_TO_FILTER = [
                    "[#7]-[#7]",
                    "[#8]~[#8]",
                    "[#7]-[#8]",
                    "[#16][#16]",
                    "[SX4](=O)(=O)[#8]",
                    "[CX3H1](=O)",
                    "[CX3](=O)[S]",
                    "[CX3](=O)O[CX3](=O)",
                    "[N]#[C][#7]"]


########################################################################
def compare_weight(a, b):
    """
    sort by REVERSE weight
    """
    a = CalcExactMolWt(a.mol)
    b = CalcExactMolWt(b.mol)
    return (a < b) - (a > b)


def single_site_right_reactant_filter(b):
    if len(b.sites(1)) == 1:
        return True

class BlockAssembler(object):
    def __init__(self, storage, project, stdout, tags=[], style='sym'):
        self.storage = storage
        self.project = project
        self.stdout = stdout
        self.tags = tags
        self.style = style

    def log_report(self, name, count, tally, duration, estimated_total):
        """ write a status line to stderr"""
        dur = datetime.timedelta(seconds=duration)
        report = {"name": name,
                  "dur": str(dur).split(".")[0],
                  "c": tally,
                  "mps": int(count/duration),
                  "est": estimated_total}
        self.stdout.write("{name} - {c} mols of ~{est:,} total in {dur} - {mps:d} mps\n".format(**report))


    def mass_binning(self, block, lower_limit=100, upper_limit=700, bin_size=100):
        mass = [CalcExactMolWt(block.cleaned().mol)]
        bins = numpy.linspace(lower_limit, upper_limit, upper_limit/bin_size)
        pos = numpy.digitize(mass, bins)
        limits = numpy.arange(lower_limit, upper_limit + bin_size, bin_size)
        under = ['under' + str(lower_limit)]
        between = [(str(limits[i]) + '_to_' + str(limits[i+1])) for i in range(0, len(limits)-1)]
        over = ['over' + str(upper_limit)]
        tags = under + between + over
        return tags[pos]

    def clean_prod(self, prod, cleaners):
        cleaned_prod = prod
        for cleaner in cleaners:
            cleaned_prod = cleaner(cleaned_prod)
        try:
            cleaned_prod = cleaned_prod.cleaned()
        except ValueError as e:
            try:
                smiles = prod.smiles()
            except ValueError:
                smiles = '<SMILES could not be determined>'
            raise ValueError("Could not clean prod '{}'".format(smiles))
        return cleaned_prod

    def check_prod_against_filters(self, prod, filters):
        if not all([f(prod) for f in filters]):
            raise ValueError("'{}' did not pass filters.".format(prod))
        return True

    def clean_and_filter_prod(self, prod, cleaners, filters):
        cleaned_prod = self.clean_prod(prod, cleaners)
        self.check_prod_against_filters(prod, filters)
        return cleaned_prod

    def clean_filter_add_to_db(self, cleaners, prod, filters, extra_tags=[],
                               add_binning_tag=False):
        """ clean product, run filters, add to database """
        try:
            cleaned_prod = self.clean_and_filter_prod(prod, cleaners, filters)
        except ValueError as e:
            logger.error(e)
            return 'bad'
        tags = ['product'] + extra_tags
        if add_binning_tag:
            tags.append(self.mass_binning(prod))
        if put_block(cleaned_prod, self.project, tags=tags):
            return "added"
        else:
            return "dupe"

    def react_and_filter(self,
                         left_set,
                         right_set,
                         env,
                         name,
                         cleaners,
                         filters,
                         return_prods=False,
                         filter_right_side=True,
                         add_to_db=True):
        """
        given two blocksets:
            - react_sym them using env (see below about all)
            - send updates to stderr periodically
            - write successful products to database (can disable via option)

            if style == 'all' or 'fill' run react_all or react_fill instead of default react_sym.
        """
        if return_prods:
            reaction_products = BlockSet()
        with env:
            if filter_right_side is True:
                single_site_right = BlockSet(right_set.filtered(single_site_right_reactant_filter))
                right_reactant = single_site_right
            else:
                right_reactant = right_set
            tallies = {"added": 0, "dupe": 0, "bad": 0, "unsaved": 0}
            count = 0
            t = time.time()

            if self.style == 'sym':
                estimated_total = left_set.estimate_sym(right_reactant)
                self.stdout.write(name + " react_sym estimate " + str(estimated_total) + "\n")
                prod_generator = left_set.react_sym(right_reactant)
            elif self.style == 'all':
                estimated_total = left_set.estimate(right_set)
                prod_generator = left_set.react_all(right_set)
            elif self.style == 'fill':
                estimated_total = left_set.estimate(right_set)
                prod_generator = left_set.react_fill(right_set)
            else:
                raise Exception("style={} must be one of 'sym', 'fill', 'all'".format(self.style))
            for prod in prod_generator:
                count += 1
                try:
                    cleaned_prod = self.clean_and_filter_prod(prod, cleaners,
                                                              filters)
                    if add_to_db:
                        if put_block(cleaned_prod, self.project,
                                     tags=(['product'] + self.tags)):
                            status = 'added'
                        else:
                            status = 'dupe'
                    else:
                        status = 'unsaved'
                except ValueError:
                    status = 'bad'
                tallies[status] += 1
                if return_prods and status != 'bad':
                    reaction_products.add(prod)
            duration = time.time() - t + 0.00001
            self.log_report( name, count, tallies, duration, estimated_total)
            if return_prods:
                return reaction_products

    def bridgify(self, seed_set, bridges, groupname, filters):
        """
        Expand pool of acceptor or donor blockset...
        given one of these sets (seed_set) and a reaction environment,
        react with the bridges to expand the set

        Also, add to the database
        """

        with cleaner_env:
            for b in seed_set:
                put_block(b.cleaned(), self.project, tags=[groupname, 'fragment'] + self.tags)

        seed_and_bridge = BlockSet()

        # a = alpha, b = beta, c = gamma, d = delta
        for pair in [("a", "c"), ("b", "c"), ("a", "d"), ("b", "d")]:
            with rxn_envs[pair]:
                one_site_bridges = BlockSet(bridges.filtered(single_site_right_reactant_filter))
                seed_and_bridge.update(BlockSet(seed_set.react_sym(one_site_bridges)))

        for pair in [("c", "a"), ("c", "b"), ("d", "a"), ("d", "b")]:
            with rxn_envs[pair]:
                one_site_seeds = BlockSet(seed_set.filtered(single_site_right_reactant_filter))
                seed_and_bridge.update(BlockSet(bridges.react_sym(one_site_seeds)))

        with cleaner_env:
            for f in filters:
                seed_and_bridge.apply_filter(f)
            for b in bridges:
                put_block(b.cleaned(), self.project, tags=["bridge"] + self.tags)
            for b in seed_and_bridge:
                put_block(b.cleaned(), self.project, tags=[groupname, "generated"] + self.tags)
        return seed_and_bridge

    def test(self):
        donors = BlockSet(self.storage.filter(groups=["donor"]))
        acceptors = BlockSet(self.storage.filter(groups=["acceptor"]))
        with old_cleaner_env:
            for b in donors:
                put_block(b.cleaned(), self.project, tags=['donor', 'fragment'] + self.tags)
            for b in acceptors:
                put_block(b.cleaned(), self.project, tags=['acceptors', 'fragment'] + self.tags)

        for pair in [("a", "b"), ("b", "a"), ("a", "a"), ("b", "b")]:
            self.react_and_filter(donors, acceptors, old_style_rxn_envs[pair], str(pair), [], [], filter_right_side=False)
            self.react_and_filter(acceptors, donors, old_style_rxn_envs[pair], str(pair), [], [], filter_right_side=False)

    def oled(self):

        # load cleaners:
        cleaner_funcs = []
        for marker_number, smarts_list in self.storage.cleaners.items():
            cleaner_funcs.append(make_cleaner_func(smarts_list, marker_number))

        donors = BlockSet(self.storage.filter(groups=["donor"]))
        acceptors = BlockSet(self.storage.filter(groups=["acceptor"]))
        steric_donors = BlockSet(self.storage.filter(groups=["steric_donor"]))
        steric_acceptors = BlockSet(self.storage.filter(groups=["steric_acceptor"]))
        bridges = BlockSet(self.storage.filter(groups=["bridge"]))
        phenyl_bridge = BlockSet(self.storage.filter(groups=["phenyl_bridge"]))

        weight_filter = mass_under(1100)
        bond_filter = contains_not_smarts(SMARTS_TO_FILTER)
        filters = [weight_filter, bond_filter]

        bridged_donors = self.bridgify(donors, bridges, "donor", filters)
        bridged_donors |= self.bridgify(steric_donors, phenyl_bridge, "donor", filters)

        bridged_acceptors = self.bridgify(acceptors, bridges, "acceptor", filters)
        bridged_acceptors |= self.bridgify(steric_acceptors, phenyl_bridge, "acceptor", filters)

        phenyl_donors = self.bridgify(donors, phenyl_bridge, "donor", filters)
        phenyl_acceptors = self.bridgify(acceptors, phenyl_bridge, "acceptor", filters)

        self.stdout.write("DONOR COUNTS {}\n".format([len(donors), len(bridged_donors)]))
        self.stdout.write("ACCEPTOR COUNTS {}\n".format([len(acceptors), len(bridged_acceptors)]))

        args = (cleaner_funcs, filters)
        for pair in [("a", "b"), ("b", "a"), ("a", "a"), ("b", "b")]:
            self.react_and_filter(donors, acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(acceptors, donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_donors, acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(acceptors, steric_donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(donors, steric_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_acceptors, donors, rxn_envs[pair], str(pair), *args)

        for pair in [("a", "c"), ("a", "d"), ("b", "c"), ("b", "d")]:
            self.react_and_filter(donors, bridged_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_donors, phenyl_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(acceptors, bridged_donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_acceptors, phenyl_donors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "a"), ("c", "b"), ("d", "a"), ("d", "b")]:
            self.react_and_filter(bridged_donors, acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(phenyl_donors, steric_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(bridged_acceptors, donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(phenyl_acceptors, steric_donors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(phenyl_donors, bridged_acceptors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(bridged_acceptors, phenyl_donors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(bridged_donors, phenyl_acceptors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(phenyl_acceptors, bridged_donors, rxn_envs[pair], str(pair), *args)

    def cbfg(self):
        #Core-bridge-functional group

        # load cleaners:
        cleaner_funcs = []
        for marker_number, smarts_list in self.storage.cleaners.items():
            cleaner_funcs.append(make_cleaner_func(smarts_list, marker_number))

        donors = BlockSet(self.storage.filter(groups=["donor"]))
        acceptors = BlockSet(self.storage.filter(groups=["acceptor"]))
        steric_donors = BlockSet(self.storage.filter(groups=["steric_donor"]))
        steric_acceptors = BlockSet(self.storage.filter(groups=["steric_acceptor"]))
        bridges = BlockSet(self.storage.filter(groups=["bridge"]))
        phenyl_bridge = BlockSet(self.storage.filter(groups=["phenyl_bridge"]))

        weight_filter = mass_under(1100)
        #bond_filter = contains_not_smarts(SMARTS_TO_FILTER)
        filters = [weight_filter]

        bridged_donors = self.bridgify(donors, bridges, "donor", filters)
        bridged_donors |= self.bridgify(steric_donors, phenyl_bridge, "donor", filters)

        bridged_acceptors = self.bridgify(acceptors, bridges, "acceptor", filters)
        bridged_acceptors |= self.bridgify(steric_acceptors, phenyl_bridge, "acceptor", filters)

        phenyl_donors = self.bridgify(donors, phenyl_bridge, "donor", filters)
        phenyl_acceptors = self.bridgify(acceptors, phenyl_bridge, "acceptor", filters)

        self.stdout.write("DONOR COUNTS {}\n".format([len(donors), len(bridged_donors)]))
        self.stdout.write("ACCEPTOR COUNTS {}\n".format([len(acceptors), len(bridged_acceptors)]))

        args = (cleaner_funcs, filters)
        for pair in [("a", "b"), ("b", "a"), ("a", "a"), ("b", "b")]:
            self.react_and_filter(donors, acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(acceptors, donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_donors, acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(acceptors, steric_donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(donors, steric_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_acceptors, donors, rxn_envs[pair], str(pair), *args)

        for pair in [("a", "c"), ("a", "d"), ("b", "c"), ("b", "d")]:
            self.react_and_filter(donors, bridged_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_donors, phenyl_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(acceptors, bridged_donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(steric_acceptors, phenyl_donors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "a"), ("c", "b"), ("d", "a"), ("d", "b")]:
            self.react_and_filter(bridged_donors, acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(phenyl_donors, steric_acceptors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(bridged_acceptors, donors, rxn_envs[pair], str(pair), *args)
            self.react_and_filter(phenyl_acceptors, steric_donors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(phenyl_donors, bridged_acceptors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(bridged_acceptors, phenyl_donors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(bridged_donors, phenyl_acceptors, rxn_envs[pair], str(pair), *args)

        for pair in [("c", "d"), ("c", "c"), ("d", "c"), ("d", "d")]:
            self.react_and_filter(phenyl_acceptors, bridged_donors, rxn_envs[pair], str(pair), *args)

    def laser(self):
        # load cleaners:
        cleaner_funcs = []
        for marker_number, smarts_list in self.storage.cleaners.items():
            cleaner_funcs.append(make_cleaner_func(smarts_list, marker_number))

        weight_filter = mass_under(2000)
        bond_filter = contains_not_smarts(LASER_SMARTS_TO_FILTER)
        filters = [weight_filter, bond_filter]

        cores = BlockSet(self.storage.filter(groups=["core"]))
        func_gs = BlockSet(self.storage.filter(groups=["functional_group"]))

        with cleaner_env:
            for b in (cores):
                put_block(b.cleaned(), self.project, tags=['core', 'fragment'] + self.tags)
            for b in (func_gs):
                put_block(b.cleaned(), self.project, tags=['functional_group', 'fragment'] + self.tags)

        single_pass = BlockSet()
        for pair in [("a", "a"), ("b", "a")]:
            with rxn_envs[pair]:
                single_pass.update(BlockSet(cores.react_sym(func_gs)))

        with cleaner_env:
            for f in filters:
                single_pass.apply_filter(f)
            for b in single_pass:
                put_block(b.cleaned(), self.project, tags=self.tags)

        for pair in [("a", "a"), ("b", "a")]:
            with rxn_envs[pair]:
                self.react_and_filter(single_pass,
                                      func_gs,
                                      rxn_envs[pair],
                                      str(pair),
                                      cleaner_funcs,
                                      filters)

    def cep(self, maximum_generation=2):
        # load cleaners:
        cleaner_funcs = []
        for marker_number, smarts_list in self.storage.cleaners.items():
            cleaner_funcs.append(make_cleaner_func(smarts_list, marker_number))

        fragments = BlockSet(self.storage.filter(groups=["core"]))

        weight_filter = mass_under(1000)
        bond_filter = contains_not_smarts(SMARTS_TO_FILTER)
        filters = [weight_filter, bond_filter]

        self.stdout.write("FRAGMENT COUNT {}\n".format([len(fragments)]))

        with cleaner_env:
            for b in fragments:
                put_block(b.cleaned(), self.project, tags=['fragment', self.mass_binning(b.cleaned())])

        gen_1 = BlockSet()

        # a = alpha, b = beta, c = gamma, d = delta

        pairs = [("a", "a"), ("a", "b"), ("b", "a"), ("a", "c"), ("c", "a"), ("a", "d"), ("d", "a"),
                 ("b", "b"), ("b", "c"), ("c", "b"), ("b", "d"), ("d", "b"),
                 ("c", "c"), ("c", "d"), ("d", "c"),
                 ("d", "d")]

        counter = 0
        for fragmentb in fragments:
            counter += 1
            gen_1 = BlockSet()
            for pair in pairs:
                gen_1.update(self.react_and_filter(BlockSet([fragmentb]),
                                                   fragments,
                                                   rxn_envs[pair],
                                                   str(pair),
                                                   cleaner_funcs,
                                                   filters,
                                                   return_prods=True,
                                                   filter_right_side=False))

            if maximum_generation > 1:
                for pair in pairs:
                    self.react_and_filter(fragments,
                                          gen_1,
                                          rxn_envs[pair],
                                          str(pair),
                                          cleaner_funcs,
                                          filters)
                    self.react_and_filter(gen_1,
                                          fragments,
                                          rxn_envs[pair],
                                          str(pair),
                                          cleaner_funcs,
                                          filters)

    def ferrocene(self):
        donors = BlockSet(self.storage.filter(groups=["donor"]))
        acceptors = BlockSet(self.storage.filter(groups=["acceptor"]))
        cp_fragment_blocksets = []
        react_and_filter_kwargs = {
            'filter_right_side': False,
            'add_to_db': False,
            'return_prods': True,
        }
        for pair in [("a", "b"), ("b", "a"), ("a", "a"), ("b", "b")]:
            cp_fragment_blocksets.append(self.react_and_filter(
                donors, acceptors, old_style_rxn_envs[pair],
                str(pair), [], [],
                **react_and_filter_kwargs))
            cp_fragment_blocksets.append(self.react_and_filter(
                acceptors, donors, old_style_rxn_envs[pair],
                str(pair), [], [],
                **react_and_filter_kwargs))
        cp_fragment_smiles_set = set()
        for cp_fragment_blockset in cp_fragment_blocksets:
            for block in cp_fragment_blockset:
                cp_fragment_smiles_set.add(block.smiles())
        tallies = {'ferrocenes': 0, 'fragments': 0}
        # Currently we just generate symmetric ferrocenes from whatever
        # reactants we provide.
        # If we want to take arbitrary fragment pairings in the future, we could
        # make the 'molgen' step only generate the fragments, and then create a
        # second step to combine pairings into ferrocenes.
        for frag_smiles in cp_fragment_smiles_set:
            common_tags = self.tags + ['product']
            # Put frag in db.
            frag_rdkit_mol = MolFromSmiles(frag_smiles)
            frag_put_result = put_mol(
                rdkit_mol=frag_rdkit_mol,
                project=self.project,
                tags=(['ferrocene_fragment'] + common_tags),
                return_mol=True)
            [frag_created, frag_pgmol] = [frag_put_result[k]
                                          for k in ['created', 'mol']]
            if frag_created:
                tallies['fragments'] += 1

            # Assmeble ferrocene from frag.
            ferrocene_smiles = "{frag_smiles}.{frag_smiles}.{iron}".format(
                frag_smiles=frag_smiles,
                iron='[Fe+2]')
            # (adorsk) We calculate ferrocene molecular charges by combining
            # charges of components, because rdkit can't add hydrogens to an
            # entire ferrocene.
            fe_atomic_num = 26
            fe_charge = -2
            ferrocene_charge_props = {
                'electron_count': (
                    2 * frag_pgmol.details['electron_count']
                    + (fe_atomic_num + fe_charge)),
                'molecular_charge': (
                    2 * frag_pgmol.details['molecular_charge']
                    + fe_charge)
            }
            ferrocene_put_result = put_smiles(
                smiles=ferrocene_smiles,
                project=self.project,
                tags=(['product', 'ferrocene'] + common_tags),
                details={
                    'fragments': (2 * [{'inchikey': frag_pgmol.inchikey}]),
                    **ferrocene_charge_props},
                calc_charge=False,
                return_mol=True)
            [ferrocene_mol, ferrocene_created] = [ferrocene_put_result[k]
                                                  for k in ['mol', 'created']]
            if ferrocene_created:
                tallies['ferrocenes'] += 1

        self.stdout.write("Created {} fragments. Fragments are tagged with"
                          " 'ferrocene_fragment'.".format(tallies['fragments']))
        self.stdout.write("Created {} molecules. Molecules are tagged with"
                          " 'ferrocene'.".format(tallies['ferrocenes']))

    def symmetric_terminal_spacer_core(self):
        cleaner_funcs = []
        for marker_number, smarts_list in self.storage.cleaners.items():
            cleaner_funcs.append(make_cleaner_func(smarts_list, marker_number))

        cores = BlockSet(self.storage.filter(groups=["core"]))
        terminals = BlockSet(self.storage.filter(groups=["terminal"]))
        spacers = BlockSet(self.storage.filter(groups=["spacer"]))

        weight_filter = mass_under(1100)
        bond_filter = contains_not_smarts(SMARTS_TO_FILTER)
        filters = [weight_filter, bond_filter]

        cores_w_spacers = BlockSet()
        for pair in [("a", "a")]:
            with rxn_envs[pair]:
                cores_w_spacers.update(BlockSet(cores.react_sym(spacers)))

        for pair in [("a", "a")]:
            self.react_and_filter(cores_w_spacers, terminals, rxn_envs[pair],
                                  str(pair), cleaner_funcs, filters)

    def spiro(self):
        # Define helpers to load fragments from storage.
        def load_fragments_as_mols_with_atom_maps():
            def replace_handles_with_atom_maps(smiles):
                handles = ['[He]', '[Ne]', '[Ar]', '[Kr]']
                for i, handle in enumerate(handles):
                    smiles = smiles.replace(handle, '[*:%s]' % (i + 1))
                return smiles

            def storage_group_to_mols_with_atom_maps(group_name=None):
                smiles_set = self.storage.filter(groups=[group_name],
                                                 output_format='raw_smiles')
                mapped_smiles_set = [replace_handles_with_atom_maps(smiles)
                                     for smiles in smiles_set]
                mols_with_atom_maps = [MolFromSmiles(mapped_smiles)
                                       for mapped_smiles in mapped_smiles_set]
                return mols_with_atom_maps

            cores = storage_group_to_mols_with_atom_maps(group_name='core')
            rgroups = {}
            for i in [1, 2, 3, 4]:
                rgroup_key = 'rgroup_type_%s' % i
                rgroups[rgroup_key] = storage_group_to_mols_with_atom_maps(
                    group_name=rgroup_key)
            return {'cores': cores, 'rgroups': rgroups}

        def squash_duplicate_mols(mols):
            return [MolFromSmiles(smiles)
                    for smiles in set([MolToSmiles(mol) for mol in mols])]

        # Load the fragments.
        fragments_with_atom_maps = load_fragments_as_mols_with_atom_maps()
        cores = fragments_with_atom_maps['cores']
        rgroups = fragments_with_atom_maps['rgroups']

        # Generate spiros via grafting reactions.
        scion_combos = itertools.product(*rgroups.values())
        core_scion_combo_combos = itertools.product(cores, scion_combos)
        products = []
        for core, scion_combo in core_scion_combo_combos:
            working_stocks = [core]
            for scion in scion_combo:
                next_working_stocks = []
                for working_stock in working_stocks:
                    intermediate_products = run_grafting_reactions(
                        stock=working_stock,
                        scion=scion)
                    next_working_stocks.extend(list(intermediate_products))
                working_stocks = squash_duplicate_mols(next_working_stocks)
            products.extend(working_stocks)
        products = squash_duplicate_mols(products)

        # Ingest into DB.
        tags = self.tags + ['spiro']
        tallies = collections.defaultdict(int)
        for product in products:
            put_result = put_smiles(
                smiles=MolToSmiles(product),
                project=self.project,
                tags=(tags),
                return_mol=True)
            if put_result['created']:
                tallies['spiros'] += 1
        self.stdout.write(
            "Created {tally} new molecules. Molecules are tagged with"
            " '{tags}'.".format(
                tally=tallies['spiros'],
                tags=tags
            ))

    def noble_graft(self):
        # Define helpers to load fragments from storage.
        def load_fragments_as_mols_with_atom_maps():
            def replace_handles_with_atom_maps(smiles):
                handles = ['[He]', '[Ne]', '[Ar]', '[Kr]']
                for i, handle in enumerate(handles):
                    smiles = smiles.replace(handle, '[*:%s]' % (i + 1))
                return smiles

            def storage_group_to_mols_with_atom_maps(group_name=None):
                smiles_set = self.storage.filter(groups=[group_name],
                                                 output_format='raw_smiles')
                mapped_smiles_set = [replace_handles_with_atom_maps(smiles)
                                     for smiles in smiles_set]
                mols_with_atom_maps = [MolFromSmiles(mapped_smiles)
                                       for mapped_smiles in mapped_smiles_set]
                return mols_with_atom_maps

            cores = storage_group_to_mols_with_atom_maps(group_name='core')
            rgroups = {}
            for blockgroup in self.storage.blocks:
                group = blockgroup['group']
                if group.startswith('rgroup'):  
                    rgroups[group] = storage_group_to_mols_with_atom_maps(
                            group_name=group)
            return {'cores': cores, 'rgroups': rgroups}

        def squash_duplicate_mols(mols):
            return [MolFromSmiles(smiles)
                    for smiles in set([MolToSmiles(mol) for mol in mols])]

        # Load the fragments.
        fragments_with_atom_maps = load_fragments_as_mols_with_atom_maps()
        cores = fragments_with_atom_maps['cores']
        rgroups = fragments_with_atom_maps['rgroups']

        # Generate via grafting reactions.
        scion_combos = itertools.product(*rgroups.values())
        core_scion_combo_combos = itertools.product(cores, scion_combos)
        products = []
        for core, scion_combo in core_scion_combo_combos:
            working_stocks = [core]
            for scion in scion_combo:
                next_working_stocks = []
                for working_stock in working_stocks:
                    intermediate_products = run_grafting_reactions(
                        stock=working_stock,
                        scion=scion)
                    next_working_stocks.extend(list(intermediate_products))
                working_stocks = squash_duplicate_mols(next_working_stocks)
            products.extend(working_stocks)
        products = squash_duplicate_mols(products)

        # Ingest into DB.
        tags = self.tags 
        tallies = collections.defaultdict(int)
        for product in products:
            put_result = put_smiles(
                smiles=MolToSmiles(product),
                project=self.project,
                tags=(tags),
                return_mol=True)
            if put_result['created']:
                tallies['mols'] += 1
        self.stdout.write(
            "Created {tally} new molecules. Molecules are tagged with"
            " '{tags}'.".format(
                tally=tallies['mols'],
                tags=tags
                ))

def run(inputfilename, project, recipe_name, stdout, tags, style, **kwargs):

    storage = BlockSetStorageJSON(inputfilename)
    blockassembler = BlockAssembler(storage, project, stdout, tags, style)

    getattr(blockassembler, recipe_name)(**kwargs)

