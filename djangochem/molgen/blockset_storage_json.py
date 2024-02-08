import json
from blocks import block, blockset

GROUP = "group"
LABEL = "label"

class BlockSetStorageJSON:
    """
    store blocks in json file as a list of dicts of this form:
    {
    "smile": "c1ccccc1",
    "label": "benzene",
    "group": "neutral",
    "reaction_sites": [
          {
            "smile": "[Be]c1ccc([Be])cc1",
            "total_count": 2,
            "mg_count": 0,
            "be_count": 2,
            "fe_count": 0
          },
          {
            "smile": "[Be]c1cc([Be])cc([Be])c1",
            "total_count": 3,
            "mg_count": 0,
            "be_count": 3,
            "fe_count": 0
          }
        ]
    }
    """
    def __init__(self, filepath):
        self.filepath = filepath
        with open(self.filepath, 'r') as input:
            root_dict = json.load(input)
            self.blocks = root_dict["molecules"]
            self.cleaners = root_dict.get("cleaners", {})

    def filter(self, groups=[], require_symmetric=False, labels=[],
               output_format=None):
        lowered_groups = [g.lower() for g in groups]
        lowered_labels = [n.lower() for n in labels]
        for m in self.blocks:
            if m[GROUP].lower() in lowered_groups:
                for b in m["blocks"]:
                    if not require_symmetric or b["symmetric"]:
                        if labels == [] or m[LABEL].lower() in lowered_labels:
                            if output_format == 'raw_smiles':
                                out = b['smiles']
                            else:
                                out = block.Block(smiles=b["smiles"])
                                out.label = m.get("label")
                            yield out
