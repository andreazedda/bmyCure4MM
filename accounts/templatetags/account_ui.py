from django import template

from accounts.ui import bilingual

register = template.Library()
register.filter("account_label", bilingual)
register.simple_tag(bilingual, name="ui")
