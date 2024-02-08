from django.contrib import admin

from models import Candidate, CandidateSet, Ballot, Vote, Comment, RedoxPair, Property, Method
# Register your models here.


# class CandidateSetAdmin(admin.ModelAdmin):
#     # it's important to override the default behavior or the
#     # admin will handing trying to make a select menu with 100k entries
#     raw_id_fields = ('members',)


class VoteAdmin(admin.ModelAdmin):
    # it's important to override the default behavior or the
    # admin will handing trying to make a select menu with 100k entries
    raw_id_fields = ('candidate',)

class CommentAdmin(admin.ModelAdmin):
    # it's important to override the default behavior or the
    # admin will handing trying to make a select menu with 100k entries
    raw_id_fields = ('candidate',)



admin.site.register(Candidate)
admin.site.register(CandidateSet)
admin.site.register(RedoxPair)
admin.site.register(Ballot)
admin.site.register(Vote, VoteAdmin)
admin.site.register(Comment, CommentAdmin)

admin.site.register(Property)
admin.site.register(Method)
