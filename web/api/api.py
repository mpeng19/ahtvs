import socket

from tastypie import authorization
from tastypie.serializers import Serializer
from tastypie_mongoengine import resources, fields
from tastypie.resources import ALL
from jinja2 import Template
from django.conf import settings

from mongoengine import Q

from blocks import block
from aag_python.molecular_storage import molecular_data_models as mdm

from blockassembler import blockdraw


from web.api.security import ApiKeyPlusWebAuthentication, MongoAuthorization

IMG_CACHE = settings.IMG_CACHE_DIR

HTML_KEYS = []


def nav_header(data):
    lower = str(data["meta"]["offset"])
    upper = str(data["meta"]["offset"] + min(len(data["objects"]), data["meta"]["limit"]))
    out = "<p>"
    if data["meta"].get("previous") is not None:
        out += "<a href='{0}'>Previous</a> --- ".format(data["meta"]["previous"])
    out += "(" + lower + " to " + upper + ")"
    if data["meta"].get("next") is not None:
        out += " --- <a href='{0}'>next</a>".format(data["meta"]["next"])
    out += "</p>"
    return out

class MoleculeSerializer(Serializer):
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'plist', 'smi', 'svg']
    content_types = {
        'json': 'application/json',
        'jsonp': 'text/javascript',
        'xml': 'application/xml',
        'yaml': 'text/yaml',
        'html': 'text/html',
        'plist': 'application/x-plist',
        'smi': 'text/plain',
        'svg': 'image/svg+xml'
    }

    def to_smi(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data:
            return "\n".join([i["meta_data"]["smiles"] for i in data["objects"]]) + "\n"
        else:
            return data["meta_data"]["smiles"] + "\n"

    def to_html(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)

        if "objects" in data and len(data["objects"]) > 0:
            out = nav_header(data)
            out += "\n".join(["<a href='{0[resource_uri]}'><img width=440px height=440px src='{0[resource_uri]}?format=svg'></a>".format(i) for i in data["objects"]])
            return out

        out = "<h1>Molecule</h1>"  #"<h3>{0}</h3>".format(data["smiles"])
        out += "<img width=440px height=440px src='{0}?format=svg'>\n".format(data["resource_uri"])
        out += "<ul>"
        out += "\n".join(["<li><strong>{0}</strong>: {1}</li>".format(k, data[k]) for k in HTML_KEYS])
        out += "</ul>"
        out += "<h3>Reactive Molecules</h3>"
        out += "\n".join(["<a href='{0}?format=html'><img width=300 height=300px src='{0}?format=svg'></a>".format(k) for k in data.get("reactive_mol_list", [])])
        out += "<h3>Parent Sets</h3><ol>"
        out += "\n".join(["<li><a href='{0}?format=html'>parent_group</a></li>".format(k) for k in data.get("parent_list", [])])
        out += "</ol>"
        # out += "<h3>Children</h3>"
        # out += "\n".join(["<a href='{0}?format=html'><img width=300 height=300px src='{0}?format=svg'></a>".format(k) for k in data.get("child_list", [])])

        out += "<h3>Calculations</h3><ul>"
        out += "\n".join(["<li><a href='{0}'>{0}</a></li>".format(k) for k in data.get("calculation_list", [])])
        out += "</ul>"

        out += "<table>"
        for k,v in data.iteritems():
            if k in ["theory", "meta_data"]:
                out += "<tr><td ><b>" + k + "</b></td><td><table>"
                for mk,mv in data[k].iteritems():
                    out += "<tr><td><b>" + mk + "</b></td><td>" + str(mv) + "</td></tr>"
                out += "</table></td></tr>"
            else:
                out += "<tr><td><b>" + k + "</b></td><td>" + str(v) + "</td></tr>"
        out += "</table>"

        return out

    def to_svg(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data:
            data = data["objects"][0]
        if "meta_data" not in data:
            raise StandardError(str(data))
        filename = blockdraw.write_svg(block.Block(data["meta_data"]["smiles"]), path=IMG_CACHE, create_groups=data["meta_data"]["project_name"])

        return open(filename, 'r').read()

class CalculationSerializer(Serializer):
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'plist', 'xyz']
    content_types = {
        'json': 'application/json',
        'jsonp': 'text/javascript',
        'xml': 'application/xml',
        'yaml': 'text/yaml',
        'html': 'text/html',
        'plist': 'application/x-plist',
        'xyz': 'text/plain',
    }

    def to_html(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data and len(data["objects"]) >= 1:
            data = data["objects"][0]
        html = """<meta http-equiv="X-UA-Compatible" content="chrome=1">
            <link rel="stylesheet" href="/static/css/ChemDoodleWeb.css" type="text/css">
            <script src="//ajax.googleapis.com/ajax/libs/jquery/1.10.2/jquery.min.js" ></script>
            <script type="text/javascript" src="/static/js/ChemDoodleWeb.js"></script>
            <script type="text/javascript" src="/static/js/ChemDoodleWeb-uis.js"></script>

            <script type="text/javascript">
              // Create 3D SketcherCanvas variable
              var sketcher3D = new ChemDoodle.TransformCanvas3D('sketcher3D', 350, 250);
              sketcher3D.specs.set3DRepresentation('Ball and Stick');
              // Read xyz from xyz
              $.ajax({
                  url:'{{resource_uri}}?format=xyz',
                  success: function(data) {
                      var xyz = new ChemDoodle.readXYZ(data);
                      // Load molecule in sketcher3D
                      sketcher3D.loadMolecule(xyz);
              }
            });
            </script>
        """

        out = "<h1>Calculation</h1>"
        if "theory" in data:
            out += "<h2>" + data["theory"]["theory_description"] + "</h2>"

        if "coord_list" in data and "resource_uri" in data and len(data["coord_list"]) > 0:
            t = Template(html)
            out  += t.render(resource_uri=data["resource_uri"])


        out += "<ul>Parent Calculations"
        out += "\n".join(["<li><a href='{0}'>{0}</a></li>".format(k) for k in data.get("parent_calculation_list",[])])
        out += "</ul>"
        out += "<ul>Child Calculations"
        out += "\n".join(["<li><a href='{0}'>{0}</a></li>".format(k) for k in data.get("child_calculation_list",[])])
        out += "</ul>"

        out += "<table>"
        for k,v in data.iteritems():
            if k in ["theory", "meta_data"]:
                out += "<tr><td ><b>" + k + "</b></td><td><table>"
                for mk,mv in data[k].iteritems():
                    out += "<tr><td><b>" + mk + "</b></td><td>" + str(mv) + "</td></tr>"
                out += "</table></td></tr>"
            else:
                out += "<tr><td><b>" + k + "</b></td><td>" + str(v) + "</td></tr>"
        out += "</table>"

        return out


    def to_xyz(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data and len(data["objects"]) == 1:
            data = data["objects"][0]
        if "coord_list" in data:
            coords = data["coord_list"]
            output = str(len(coords))+ "\n\n"
            for c in coords:
                output += " ".join([c["element"], str(c["x"]), str(c["y"]), str(c["z"])]) + "\n"
            return output
        else:
            return data


class MolecularLinkageSerializer(Serializer):
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'plist']
    content_types = {
        'json': 'application/json',
        'jsonp': 'text/javascript',
        'xml': 'application/xml',
        'yaml': 'text/yaml',
        'html': 'text/html',
        'plist': 'application/x-plist',
    }

    def to_html(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "meta" in data:
            out = nav_header(data)

        if "objects" in data:
            for ml in data["objects"]:
                out += "<div style='padding:3px;display:inline-block;border-width:1px; border-style:dotted; border-color:#999;margin:2px;'>"
                out += "<a style='color:#ccc' href='{0}?format=html'>".format(ml["resource_uri"])
                out += "".join(["<img style='margin:5px' height=100px src='{0}?format=svg'>".format(k) for k in ml.get("parent_molecule_list", [])])
                out += '<div style="display:inline-block;height:60px;margin:20px;border-width=1px;width: 0px;border-style: dotted;border-width: 1px;"></div>'
                out += "".join(["<img style='margin:5px' height=100px src='{0}?format=svg'>".format(k) for k in ml.get("child_molecule_list",[])])
                out += "</a></div>"
        else:
            out = "<h1>Molecular Linkage</h1>"
            out += '<h3>Molecules</h3>'
            out += "<h4>Parents</h4>"
            out += "\n".join(["<a href='{0}?format=html'><img width=300 height=300px src='{0}?format=svg'></a>".format(k) for k in data.get("parent_molecule_list", [])])
            out += "<h4>Children</h4>"
            out += "\n".join(["<a href='{0}?format=html'><img width=300 height=300px src='{0}?format=svg'></a>".format(k) for k in data.get("child_molecule_list",[])])

            out += "<h3>Reactive Molecules</h3>"
            out += "<h4>Parents</h4>"
            out += "\n".join(["<a href='{0}?format=html'><img width=300 height=300px src='{0}?format=svg'></a>".format(k) for k in data.get("parent_reactive_mol_list", [])])
            out += "<h4>Children</h4>"
            out += "\n".join(["<a href='{0}?format=html'><img width=300 height=300px src='{0}?format=svg'></a>".format(k) for k in data.get("child_reactive_mol_list",[])])


            out += "<table>"
            for k,v in data.iteritems():
                if k in ["info", "meta_data"]:
                    out += "<tr><td ><b>" + k + "</b></td><td><table>"
                    for mk,mv in data[k].iteritems():
                        out += "<tr><td><b>" + mk + "</b></td><td>" + str(mv) + "</td></tr>"
                    out += "</table></td></tr>"
                else:
                    out += "<tr><td><b>" + k + "</b></td><td>" + str(v) + "</td></tr>"
            out += "</table>"
        return out

class TheoryResource(resources.MongoEngineResource):
    class Meta:
        object_class = mdm.Theory
        authorization = authorization.Authorization()
        filtering = {
            'theory_level': ALL
        }

class PropertiesResource(resources.MongoEngineResource):
    class Meta:
        object_class = mdm.Properties
        authorization = authorization.Authorization()
        ordering = ["total_energy"]



class ReactiveMoleculeSerializer(Serializer):
    formats = ['json', 'html', 'svg', 'smi']
    content_types = {
        'json': 'application/json',
        'html': 'text/html',
        'smi': 'text/plain',
        'svg': 'image/svg+xml'
    }

    def to_html(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)

        if "objects" in data and len(data["objects"]) > 0:
            out = nav_header(data)

            out += "\n".join(["<a href='{0[resource_uri]}'><img width=440px height=440px src='{0[resource_uri]}?format=svg'></a>".format(i) for i in data["objects"]])
            return out

        out = "<h1>Reactive Molecule</h1><img width=440px height=440px src='{0[resource_uri]}?format=svg'>".format(data)
        out += "<h3>Parent Sets</h3><ol>"
        out += "\n".join(["<li><a href='{0}?format=html'>parent_group</a></li>".format(k) for k in data.get("parent_list", [])])
        out += "</ol>"
        # out += "<h3>Children</h3>"
        # out += "\n".join(["<a href='{0}?format=html'><img width=300 height=300px src='{0}?format=svg'></a>".format(k) for k in data.get("child_list", [])])

        out += "<table>"
        for k,v in data.iteritems():
            if k in ["meta_data"]:
                out += "<tr><td ><b>" + k + "</b></td><td><table>"
                for mk,mv in data[k].iteritems():
                    out += "<tr><td><b>" + mk + "</b></td><td>" + str(mv) + "</td></tr>"
                out += "</table></td></tr>"
            else:
                out += "<tr><td><b>" + k + "</b></td><td>" + str(v) + "</td></tr>"
        out += "</table>"
        return out

    def to_smi(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        if "objects" in data:
            return "\n".join([i["meta_data"]["smiles"] for i in data["objects"]]) + '\n'
        else:
            return data["meta_data"]["smiles"] + '\n'

    def to_svg(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        filename = blockdraw.write_svg(block.Block(data["meta_data"]["smiles"]), path=IMG_CACHE, create_groups=data["meta_data"]["project_name"])
        return open(filename, 'r').read()


class MetaResource(resources.MongoEngineResource):
    class Meta:
        object_class = mdm.MetaData
        filtering = {
            'tag_list': ALL
        }


class ReactiveMolResource(resources.MongoEngineResource):
    class Meta:
        queryset = mdm.ReactiveMolecule.objects.all()
        authentication = ApiKeyPlusWebAuthentication()
        authorization = MongoAuthorization()
        serializer = ReactiveMoleculeSerializer()
        filtering = {
            'tag_list': ALL
        }
        excludes = ["child_list"]
    meta_data = fields.EmbeddedDocumentField(embedded='api.api.MetaResource', attribute='meta_data', null=True)
    #child_list = fields.ReferencedListField(of='api.api.ReactiveMolResource', attribute='child_list', full=False, null=True)
    parent_list = fields.ReferencedListField(of='api.api.MolecularLinkageResource', attribute='parent_list', full=False, null=True)

class CalculationResource(resources.MongoEngineResource):
    class Meta:
        queryset = mdm.Calculation.objects.all()
        authentication = ApiKeyPlusWebAuthentication()
        authorization = MongoAuthorization()
        serializer = CalculationSerializer()
        allowed_methods = ('get',)
        ordering = ("properties",)
    meta_data = fields.EmbeddedDocumentField(embedded='api.api.MetaResource', attribute='meta_data', null=True)
    child_calculation_list = fields.ReferencedListField(of='api.api.CalculationResource', attribute='child_calculation_list', full=False, null=True)
    parent_calculation_list = fields.ReferencedListField(of='api.api.CalculationResource', attribute='parent_calculation_list', full=False, null=True)
    theory = fields.EmbeddedDocumentField(embedded='api.api.TheoryResource', attribute='theory', null=True)
    properties = fields.EmbeddedDocumentField(embedded='api.api.PropertiesResource', attribute='properties', null=True)

    def build_filters(self, filters=None):
        """
        for any filters specified that are in meta_data, use them with meta_data
        """
        if filters is None:
            filters = {}

        orm_filters = super(CalculationResource, self).build_filters(filters)

        meta_data_keys = set(filters.keys()) & set(mdm.MetaData._fields)
        for k in meta_data_keys:
            if "," in filters[k]:
                orm_filters.update({"meta_data__{}__all".format(k): filters[k].split(",")})
            elif "|" in filters[k]:
                orm_filters.update({"meta_data__{}__in".format(k): filters[k].split("|")})
            else:
                orm_filters.update({"meta_data__{}".format(k): filters[k]})

        theory_keys = set(filters.keys()) & set(mdm.Theory._fields)
        orm_filters.update(dict(("theory__{}".format(k), filters[k]) for k in theory_keys))

        properties_keys = set(filters.keys()) & set(mdm.Properties._fields)
        orm_filters.update(dict(("properties__{}".format(k), filters[k]) for k in properties_keys))


        calculation_keys = set(filters.keys()) & set(mdm.Calculation._fields)
        orm_filters.update(dict((k, filters[k]) for k in calculation_keys))

        return orm_filters


class MoleculeResource(resources.MongoEngineResource):
    class Meta:
        queryset = mdm.Molecule.objects.all()
        #allowed_methods = ('get', 'post', 'put', 'delete')
        allowed_methods = ('get',)

        authentication = ApiKeyPlusWebAuthentication()
        authorization = MongoAuthorization()

        serializer = MoleculeSerializer()
        filtering = {"id": ALL}
        limit = 20
        max_limit = 500
        excludes = ["child_list"]

    meta_data = fields.EmbeddedDocumentField(embedded='api.api.MetaResource', attribute='meta_data', null=True)
    reactive_mol_list = fields.ReferencedListField(of='api.api.ReactiveMolResource', attribute='reactive_mol_list', full=False, null=True)
    calculation_list = fields.ReferencedListField(of='api.api.CalculationResource', attribute='calculation_list', full=False, null=True)

    #child_list = fields.ReferencedListField(of='api.api.MoleculeResource', attribute='child_list', full=False, null=True)
    parent_list = fields.ReferencedListField(of='api.api.MolecularLinkageResource', attribute='parent_list', full=False, null=True)



    def build_filters(self, filters=None):
        """
        for any filters specified that are in meta_data, use them with meta_data
        """
        if filters is None:
            filters = {}

        orm_filters = super(MoleculeResource, self).build_filters(filters)

        meta_data_keys = set(filters.keys()) & set(mdm.MetaData._fields)
        for k in meta_data_keys:
            if "," in filters[k]:
                orm_filters.update({"meta_data__{}__all".format(k): filters[k].split(",")})
            elif "|" in filters[k]:
                orm_filters.update({"meta_data__{}__in".format(k): filters[k].split("|")})
            else:
                orm_filters.update({"meta_data__{}".format(k): filters[k]})

        molecule_keys = set(filters.keys()) & set(mdm.Molecule._fields)
        orm_filters.update(dict((k, filters[k]) for k in molecule_keys))

        return orm_filters


    #     query = filters.get('query')
    #     if query:
    #         qset = (
    #             Q(comment__icontains=query) |
    #             Q(media_text__icontains=query)
    #         )
    #         applicable_filters['custom'] = qset

    #     return applicable_filters

    # def apply_filters(self, request, applicable_filters):
    #     custom = None
    #     if 'custom' in applicable_filters:
    #         custom = applicable_filters.pop('custom')

    #     semi_filtered = super(TaggedResource, self).apply_filters(request, applicable_filters)

    #     return semi_filtered.filter(**custom) if custom else semi_filtered


class MolecularLinkageResource(resources.MongoEngineResource):
    class Meta:
        queryset = mdm.MolecularLinkage.objects.all()
        #allowed_methods = ('get', 'post', 'put', 'delete')
        allowed_methods = ('get',)
        authentication = ApiKeyPlusWebAuthentication()
        authorization = MongoAuthorization()
        serializer = MolecularLinkageSerializer()
        limit = 100
        max_limit = 1000

    meta_data = fields.EmbeddedDocumentField(embedded='api.api.MetaResource', attribute='meta_data', null=True)

    child_reactive_mol_list = fields.ReferencedListField(of='api.api.ReactiveMolResource', attribute='child_reactive_mol_list', full=False, null=True)
    parent_reactive_mol_list = fields.ReferencedListField(of='api.api.ReactiveMolResource', attribute='parent_reactive_mol_list', full=False, null=True)

    child_molecule_list = fields.ReferencedListField(of='api.api.MoleculeResource', attribute='child_molecule_list', full=False, null=True)
    parent_molecule_list = fields.ReferencedListField(of='api.api.MoleculeResource', attribute='parent_molecule_list', full=False, null=True)

    def build_filters(self, filters=None):
        """
        for any filters specified that are in meta_data, use them with meta_data
        """
        if filters is None:
            filters = {}

        orm_filters = super(MolecularLinkageResource, self).build_filters(filters)

        meta_data_keys = set(filters.keys()) & set(mdm.MetaData._fields)
        for k in meta_data_keys:
            if "," in filters[k]:
                orm_filters.update({"meta_data__{}__all".format(k): filters[k].split(",")})
            elif "|" in filters[k]:
                orm_filters.update({"meta_data__{}__in".format(k): filters[k].split("|")})
            else:
                orm_filters.update({"meta_data__{}".format(k): filters[k]})

        print orm_filters
        molecule_keys = set(filters.keys()) & set(mdm.Molecule._fields)
        orm_filters.update(dict((k, filters[k]) for k in molecule_keys))

        return orm_filters
