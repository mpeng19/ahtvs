from django import forms
from .models import Ballot
from datetimewidget.widgets import DateTimeWidget

dateTimeOptions = {
'format': 'yyyy/mm/dd hh:00',
'autoclose': 'true',
'showMeridian' : 'false',
'todayHighlight': 'true',
'minuteStep': 60,
'minView': 1
}


class NewBallotForm(forms.ModelForm):
    class Meta:
        model = Ballot
        fields = ['candidateset', 'voter', 'name', 'close_time']

    close_time = forms.DateTimeField(widget=DateTimeWidget(bootstrap_version=3, options = dateTimeOptions))

    def __init__(self, *args, **kwargs):
        super(NewBallotForm, self).__init__(*args, **kwargs)
        self.fields['candidateset'].choices = sorted(self.fields['candidateset'].choices, key=lambda x: x[1])
        self.fields['voter'].choices = sorted(self.fields['voter'].choices, key=lambda x: x[1])
