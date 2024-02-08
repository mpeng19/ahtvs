# myapp/mongoadmin.py

# Import the MongoAdmin base class
from mongonaut.sites import MongoAdmin

# Import your custom models
from aag_python.molecular_storage import molecular_data_models as mdm

# Instantiate the MongoAdmin class
# Then attach the mongoadmin to your model
mdm.Molecule.mongoadmin = MongoAdmin()
