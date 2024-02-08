from tastypie.authentication import ApiKeyAuthentication
from tastypie.authorization import Authorization
from tastypie.exceptions import Unauthorized


class ApiKeyPlusWebAuthentication(ApiKeyAuthentication):
    def is_authenticated(self, request, **kwargs):
        if request.user.is_authenticated():
            return True
        return super(ApiKeyPlusWebAuthentication, self).is_authenticated(request, **kwargs)

    def get_identifier(self, request):
        if request.user.is_authenticated():
            return request.user.username
        else:
            return super(ApiKeyPlusWebAuthentication, self).get_identifier(request)


class MongoAuthorization(Authorization):
    def read_list(self, object_list, bundle):
        # This assumes a ``QuerySet`` from ``ModelResource``.
        if bundle.request.user.is_staff:
            return object_list
        groups = bundle.request.user.groups.values_list('name', flat=True)
        result = object_list.filter(meta_data__project_name__in=groups, meta_data__status__released=True)
        if object_list.count() > 0 and result.count() == 0:
            raise Unauthorized("You are not a member of the project for these items")
        return result

    def read_detail(self, object_list, bundle):
        # Is the requested object owned by the user?
        if bundle.request.user.is_staff:
            return True
        groups = bundle.request.user.groups.values_list('name', flat=True)
        if not bundle.obj.meta_data.status.released:
            raise Unauthorized("No such data")
        if not (bundle.obj.meta_data.project_name in groups):
            raise Unauthorized("You are not a member of the project for this item")
        return True

    def create_list(self, object_list, bundle):
        # Assuming they're auto-assigned to ``user``.
        return object_list

    def create_detail(self, object_list, bundle):
        return bundle.obj.user == bundle.request.user

    def update_list(self, object_list, bundle):
        allowed = []

        # Since they may not all be saved, iterate over them.
        for obj in object_list:
            if obj.user == bundle.request.user:
                allowed.append(obj)

        return allowed

    def update_detail(self, object_list, bundle):
        return bundle.obj.user == bundle.request.user

    def delete_list(self, object_list, bundle):
        # Sorry user, no deletes for you!
        raise Unauthorized("Sorry, no deletes.")

    def delete_detail(self, object_list, bundle):
        raise Unauthorized("Sorry, no deletes.")


class SQLAuthorization(Authorization):
    def read_list(self, object_list, bundle):
        if bundle.request.user.is_staff:
            return object_list
        groups = bundle.request.user.groups.values_list('name', flat=True)
        result = object_list.filter(released=True, project__in=groups)

        if object_list.count() > 0 and result.count() == 0:
            raise Unauthorized("You are not a member of the project for these items")

        return result

    def read_detail(self, object_list, bundle):
        # Is the requested object owned by the user?
        groups = bundle.request.user.groups.values_list('name', flat=True)
        if bundle.request.user.is_staff:
            return True

        if not bundle.obj.project in groups:
            raise Unauthorized("You are not a member of the project for this item")
        if not bundle.obj.released:
            raise Unauthorized("No such data")
        return True

    def create_list(self, object_list, bundle):
        # Assuming they're auto-assigned to ``user``.
        return object_list

    def create_detail(self, object_list, bundle):
        return bundle.obj.user == bundle.request.user

    def update_list(self, object_list, bundle):
        allowed = []

        # Since they may not all be saved, iterate over them.
        for obj in object_list:
            if obj.user == bundle.request.user:
                allowed.append(obj)

        return allowed

    def update_detail(self, object_list, bundle):
        return bundle.obj.user == bundle.request.user

    def delete_list(self, object_list, bundle):
        # Sorry user, no deletes for you!
        raise Unauthorized("Sorry, no deletes.")

    def delete_detail(self, object_list, bundle):
        raise Unauthorized("Sorry, no deletes.")


class VoteAuthorization(Authorization):
    def read_list(self, object_list, bundle):
        if bundle.request.user.is_staff:
            return object_list
        result = object_list.filter(ballot__voter=bundle.request.user)
        return result

    def read_detail(self, object_list, bundle):
        # Is the requested object owned by the user?
        if bundle.request.user.is_staff:
            return True

        if not bundle.obj.ballot.voter == bundle.request.user:
            raise Unauthorized("This is not your Vote to see")

        return True

    def create_list(self, object_list, bundle):
        # Assuming they're auto-assigned to ``user``.
        return object_list

    def create_detail(self, object_list, bundle):
        return bundle.obj.user == bundle.request.user

    def update_list(self, object_list, bundle):
        allowed = []

        # Since they may not all be saved, iterate over them.
        for obj in object_list:
            if obj.user == bundle.request.user:
                allowed.append(obj)

        return allowed

    def update_detail(self, object_list, bundle):
        return bundle.obj.user == bundle.request.user

    def delete_list(self, object_list, bundle):
        # Sorry user, no deletes for you!
        raise Unauthorized("Sorry, no deletes.")

    def delete_detail(self, object_list, bundle):
        raise Unauthorized("Sorry, no deletes.")