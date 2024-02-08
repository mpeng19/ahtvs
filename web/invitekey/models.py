import os
import hashlib
import base64

from django.db import models

# Create your models here.

def make_key():
    return base64.urlsafe_b64encode(os.urandom(20))

class Invitation(models.Model):
    key = models.CharField(max_length=50, unique=True, default=make_key)
    group = models.CharField(max_length=50)
    email = models.CharField(max_length=100, default="")

    def __unicode__(self):
        return self.group + "-" + str(self.pk)


