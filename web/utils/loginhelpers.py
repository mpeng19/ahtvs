import json
import re
import functools

from django.utils import timezone
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.utils.decorators import method_decorator
from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404
from django.core.exceptions import PermissionDenied
from django.views.generic import View


def only_staff_allowed(fn):
    '''decorator'''
    def wrapped(request, *args, **kwargs):
        if request.user.is_authenticated():
            if request.user.is_staff:
                return fn(request, *args, **kwargs)
            else:
                return HttpResponseRedirect(reverse('login') + "?next=/")
        else:
            return HttpResponseRedirect(reverse('login') + "?next=" + request.get_full_path())
    return wrapped


class StaffOnlyMixin(object):
    @method_decorator(only_staff_allowed)
    def dispatch(self, *args, **kwargs):
        return super(StaffOnlyMixin, self).dispatch(*args, **kwargs)


class LoginRequiredMixin(object):
    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super(LoginRequiredMixin, self).dispatch(*args, **kwargs)
