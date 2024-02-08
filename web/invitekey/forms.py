import urllib

from django import forms
from django.core.validators import validate_email
from django.core.mail import send_mail
from django.template import Context
from django.template.loader import render_to_string
from django.contrib.sites.shortcuts import get_current_site
from django.contrib.auth.models import Group
from django.core.urlresolvers import reverse

from registration.forms import RegistrationForm
from .models import Invitation

class InviteKeyRegistrationForm(RegistrationForm):
    key = forms.CharField(widget=forms.HiddenInput())
    first_name = forms.CharField()
    last_name = forms.CharField()


class MultiEmailField(forms.Field):
    def to_python(self, value):
        "Normalize data to a list of strings."
        # Return an empty list if no input was given.
        if not value:
            return []
        return [v.strip() for v in value.split(',')]

    def validate(self, value):
        "Check if value consists only of valid emails."

        # Use the parent's handling of required fields, etc.
        super(MultiEmailField, self).validate(value)

        for email in value:
            validate_email(email)

class InviteForm(forms.Form):
    group = forms.ModelChoiceField(Group.objects.all())
    message = forms.CharField(widget=forms.Textarea)
    emails = MultiEmailField()

    def send_invites(self, request):
        "create a bunch of invites and send emails with the corresponding keys"
        current_site = get_current_site(request)
        domain = current_site.domain.rstrip("/") + "/"
        site_name = current_site.name
        i = 0
        for email in self.cleaned_data["emails"]:
            inv = Invitation(group=self.cleaned_data["group"], email=email)
            inv.save()
            invite_url = domain + reverse("registration_register")[1:]
            invite_url += '?' + urllib.urlencode({"key": inv.key, "email": email})

            c = Context({"invite_url": invite_url, "message": self.cleaned_data["message"]})
            content = render_to_string("invitekey/invite_mail.txt", c)
            send_mail('Invitation to Join ' + site_name, content, 'noreply@' + domain,
                      [email], fail_silently=False)
            i += 1
        return i
