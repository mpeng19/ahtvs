from django.shortcuts import render
from django.contrib.auth.models import Group, User
from django.contrib.auth import authenticate
from django.contrib.auth import login, logout
from django.contrib.sites.models import RequestSite
from django.contrib.sites.models import Site
from django.contrib import messages
from django.shortcuts import get_object_or_404
from django.views.generic.edit import FormView
from django.conf import settings

from registration import signals
from registration.models import RegistrationProfile
from registration.backends.default.views import ActivationView
from registration.backends.default.views import RegistrationView

from utils.loginhelpers import StaffOnlyMixin

from .models import Invitation
from .forms import InviteKeyRegistrationForm, InviteForm


class InvitedRegistrationView(RegistrationView):
    form_class = InviteKeyRegistrationForm
    """
    A registration backend which follows a simple workflow:

    0. User follows link to registratioin page with activation key
    1. User signs up, inactive account is created.

    2. Email is sent to user with activation link.

    3. User clicks activation link, account is now active.



    """
    def register(self, request, form):
        """
        Given a username, email address and password, register a new
        user account, which will initially be inactive.

        Along with the new ``User`` object, a new
        ``registration.models.RegistrationProfile`` will be created,
        tied to that ``User``, containing the activation key which
        will be used for this account.

        An email will be sent to the supplied email address; this
        email should contain an activation link. The email will be
        rendered using two templates. See the documentation for
        ``RegistrationProfile.send_activation_email()`` for
        information about these templates and the contexts provided to
        them.

        After the ``User`` and ``RegistrationProfile`` are created and
        the activation email is sent, the signal
        ``registration.signals.user_registered`` will be sent, with
        the new ``User`` as the keyword argument ``user`` and the
        class of this backend as the sender.

        """

        cleaned_data = form.clean()
        username = cleaned_data['username']
        password = cleaned_data['password1']
        email = cleaned_data['email']
        key = cleaned_data["key"]
        if Site._meta.installed:
            site = Site.objects.get_current()
        else:
            site = RequestSite(request)

        inv = get_object_or_404(Invitation, key=key)
        if inv.email == email:
            User.objects.create_user(username, email, password)
            new_user = authenticate(username=username, password=password)
            login(request, new_user)

        else:
            # unusual edge case of user creating second account.  log them out first.
            logout(request)
            new_user = RegistrationProfile.objects.create_inactive_user(username, email,
                                                                        password, site)

        g = Group.objects.get(name=inv.group)

        new_user.groups.add(g)
        new_user.first_name = cleaned_data["first_name"]
        new_user.last_name = cleaned_data["last_name"]
        new_user.save()

        inv.delete()  # remove invitation so it cant be re-used

        signals.user_registered.send(sender=self.__class__,
                                     user=new_user,
                                     request=request)
        return new_user

    def registration_allowed(self, request):
        if request.POST.get("key") is not None:
            key = request.POST.get("key")
        else:
            key = request.GET.get("key")

        return Invitation.objects.filter(key=key).count() > 0 and getattr(settings, 'REGISTRATION_OPEN', True)

    def get_initial(self):
        if self.request.GET:
            initial = self.request.GET.dict()
            return initial
        else:
            return super(InvitedRegistrationView, self).get_initial()

    def get_success_url(self, request, user):
        if request.user.is_authenticated():
            return ('registration_activation_complete_preapproved', (), {})
        else:
            return ('registration_complete', (), {})



class InviteView(StaffOnlyMixin, FormView):
    template_name = 'registration/invite.html'
    form_class = InviteForm
    success_url = 'invite'

    def form_valid(self, form):
        # This method is called when valid form data has been POSTed.
        # It should return an HttpResponse.
        count = form.send_invites(self.request)
        messages.success(self.request, "{} Invitations sent".format(count))
        return super(InviteView, self).form_valid(form)