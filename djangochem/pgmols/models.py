from __future__ import unicode_literals

from django.contrib.contenttypes.fields import GenericRelation
from django.utils import timezone
from django.contrib.auth.models import User, Group
from django.db import models
from django.contrib.postgres.fields import JSONField, ArrayField
from django.db.models.signals import pre_save

from guardian.models import UserObjectPermissionBase
from guardian.models import GroupObjectPermissionBase
# from django_rdkit import models as rdmodels
from rdkit.Chem import AllChem as Chem

PERIODICTABLE = Chem.GetPeriodicTable()
INCHI_OPTION_STRING = " -RecMet  -FixedH "


# class Project(models.Model):
#     name = models.CharField(max_length=256)


class MolManager(models.Manager):
    def get_by_natural_key(self, project, inchikey):
        return self.get(project=project, inchikey=inchikey)


class Mol(models.Model):
    objects = MolManager()

    """ analog to molecule, but we don't want to reuse that class name.  too confusing!"""
    # rdmol = rdmodels.MolField()
    inchikey = models.CharField(max_length=27, null=True, db_index=True)
    smiles = models.CharField(max_length=1000, db_index=True)
    nicknames = ArrayField(models.CharField(max_length=200), null=True)
    mass = models.FloatField(null=True)
    batches = models.ManyToManyField('Batch')
    tags = ArrayField(models.CharField(max_length=200), null=True)
    createtime = models.DateTimeField(null=True, default=timezone.now, db_index=True)
    details = JSONField(default={})
    parents = models.ManyToManyField('self',
                                     symmetrical=False,
                                     related_name='children')
    group = models.ForeignKey(Group,
                              null=True,
                              on_delete=models.SET_NULL,
                              db_index=True)
    parentjob = models.ForeignKey('jobs.Job',
                                  null=True,
                                  related_name='childmols',
                                  on_delete=models.SET_NULL)
    childjobs = GenericRelation('jobs.Job',
                                content_type_field="parentct",
                                object_id_field="parentid",
                                )

    class Meta:
        permissions = (
            ('view_mol', 'View Mol'),
        )
        unique_together = (('group', 'inchikey'),)
        index_together = (('group', 'inchikey'),)

    def __str__(self):
        return self.inchikey

    def __repr__(self):
        return '<Mol proj={} inchi={}>'.format(self.group, self.inchikey)

    # def auto_fill(self):
    #     if not self.rdmol and self.smiles is not None:
    #         self.rdmol = Chem.MolFromSmiles(str(self.smiles))
    #     if isinstance(self.rdmol, six.string_types):
    #         self.rdmol = Chem.MolFromSmiles(str(self.rdmol))
    #     if not self.smiles:
    #         self.smiles = Chem.MolToSmiles(self.rdmol)
    #     if not self.inchikey:
    #         non_std_inchi = Chem.MolToInchi(self.rdmol,
    #                                         options=str(INCHI_OPTION_STRING))
    #         self.inchikey = Chem.InchiToInchiKey(non_std_inchi)

    def auto_fill(self):
        if self.smiles is not None:
            rdmol = None
            rdmol = Chem.MolFromSmiles(str(self.smiles))
            if rdmol is None:
                raise Exception("Invalid smiles {}".format(self.smiles))
            if not self.inchikey:
                non_std_inchi = Chem.MolToInchi(rdmol,
                                                options=str(INCHI_OPTION_STRING))
                self.inchikey = Chem.InchiToInchiKey(non_std_inchi)

    @property
    def molecular_charge(self):
        """I'm the 'x' property."""
        return Chem.GetFormalCharge(Chem.MolFromSmiles(str(self.smiles)))

    def to_dict(self):
        return {**self.details, **{'mol': Chem.MolFromSmiles(self.smiles), 'smiles': self.smiles,**{'tags':self.tags}}}


def fill_out_mol(sender, **kwargs):
    kwargs["instance"].auto_fill()


pre_save.connect(fill_out_mol, sender=Mol)


class MolUserObjectPermission(UserObjectPermissionBase):
    content_object = models.ForeignKey(Mol, on_delete=models.CASCADE) #TODO:check on delete


class MolGroupObjectPermission(GroupObjectPermissionBase):
    content_object = models.ForeignKey(Mol, on_delete=models.CASCADE) #TODO:check on delete


class Reaction(models.Model):
    name = models.CharField(max_length=128)
    reactants = models.ManyToManyField(Mol, related_name='reactantto')
    products = models.ManyToManyField(Mol, related_name='productof')
    details = JSONField(null=True)

    def __str__(self):
        return self.name


class Method(models.Model):
    name = models.CharField(max_length=255)
    description = models.CharField(max_length=256, null=True)
    details = JSONField(null=True)

    def __str__(self):
        return self.name


class Geom(models.Model):
    xyz = ArrayField(ArrayField(models.FloatField(), size=4))
    mol = models.ForeignKey(Mol, on_delete=models.CASCADE, db_index=True)
    method = models.ForeignKey(Method,
                               on_delete=models.PROTECT,
                               related_name="+",
                               db_index=True)
    parents = models.ManyToManyField('self',
                                     symmetrical=False,
                                     related_name='children')
    parentjob = models.ForeignKey('jobs.Job',
                                  null=True,
                                  related_name='childgeoms',
                                  on_delete=models.SET_NULL)
    childjobs = GenericRelation('jobs.Job',
                                content_type_field="parentct",
                                object_id_field="parentid",
                                )
    details = JSONField(null=True)

    def get_coords(self):
        return [dict(element=PERIODICTABLE.GetElementSymbol(int(l[0])),
                     x=l[1],
                     y=l[2],
                     z=l[3]
                     ) for l in self.xyz]

    def set_coords(self, coord_list):
        self.xyz = []
        for atom_dict in coord_list:
            row = [PERIODICTABLE.GetAtomicNumber(str(atom_dict['element'])),
                   atom_dict['x'],
                   atom_dict['y'],
                   atom_dict['z']]
            self.xyz.append(row)

    def as_xyz(self):
        coords = self.get_coords()
        output = str(len(coords)) + "\n\n"
        for c in coords:
            output += " ".join([c["element"], str(c["x"]), str(c["y"]), str(c["z"])]) + "\n"
        return output

    def __str__(self):
        return str(self.id)


class Calc(models.Model):
    method = models.ForeignKey(Method, on_delete=models.CASCADE, db_index=True) #TODO:check on delete
    props = JSONField(null=True)
    parents = models.ManyToManyField('self',
                                     symmetrical=False,
                                     related_name='children')
    geoms = models.ManyToManyField(Geom,
                                   related_name='calcs')
    mol = models.ForeignKey(Mol,
                            on_delete=models.CASCADE,
                            related_name='calcs',
                            db_index=True,
                            null=True)
    reaction = models.ForeignKey(Reaction,
                                 on_delete=models.CASCADE,
                                 related_name='calcs',
                                 db_index=True,
                                 null=True)
    parentjob = models.ForeignKey('jobs.Job',
                                  null=True,
                                  related_name='childcalcs',
                                  on_delete=models.SET_NULL)
    childjobs = GenericRelation('jobs.Job',
                                content_type_field="parentct",
                                object_id_field="parentid",
                                )

    def __str__(self):
        return str(self.id) + ' : ' + self.mol.inchikey + ' : ' +self.method.name
        #return str(self.id)


class Batch(models.Model):
    name = models.CharField(max_length=255)
    createtime = models.DateTimeField(null=True, default=timezone.now)
    announced = models.BooleanField(default=False)
    creator = models.ForeignKey(User,
                                null=True,
                                on_delete=models.SET_NULL,
                                related_name='batches')

    def __str__(self):
        return self.name

    def size(self):
        return Mol.objects.filter(batch=self).count()

    def announce(self):
        self.announced = True
        self.release()

    def release(self):
        self.released = True
        self.save()
        for c in Mol.objects.filter(batch=self):
            c.released = True
            c.save()

# def add_view_permissions(sender, **kwargs):
#     """
#     This syncdb hooks takes care of adding a view permission too all our
#     content types.
#     """
#     # for each of our content types
#     for content_type in ContentType.objects.filter(app_label='pgmols'):
#         # build our permission slug
#         codename = "view_" + content_type.model

#         # if it doesn't exist..
#         if not Permission.objects.filter(content_type=content_type, codename=codename):
#             # add it
#             Permission.objects.create(content_type=content_type,
#                                       codename=codename,
#                                       name="Can view %s" % content_type.name)
#             print "Added view permission for %s" % content_type.name

# # check for all our view permissions after a syncdb
# post_migrate.connect(add_view_permissions)
