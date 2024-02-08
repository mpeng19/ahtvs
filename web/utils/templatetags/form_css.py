from django import template
register = template.Library()

@register.filter(name='add_class')
def add_class(field, classname):
    attrs = {'class': classname}
    return field.as_widget(attrs=attrs)

@register.filter(name='set_id')
def set_id(field, id):
    attrs = {'id': id}
    return field.as_widget(attrs=attrs)

@register.filter(name='add_attributes')
def add_attributes(field, css):
    attrs = {}
    definition = css.split(',')

    for d in definition:
        if ':' not in d:
            attrs['class'] = d
        else:
            t, v = d.split(':')
            attrs[t] = v

    return field.as_widget(attrs=attrs)