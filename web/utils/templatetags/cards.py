from django import template
register = template.Library()

@register.filter(name='format_inchi_key')
def format_inchi_key(value):
    return value.split("-")[0]
